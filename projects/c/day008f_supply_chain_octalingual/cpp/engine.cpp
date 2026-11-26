#include "common.hpp"

static int64_t wh_total_stock(const Sim &s, int32_t w) {
    int64_t tot = 0;
    for (int k = 0; k < NUM_SKUS; ++k) tot += s.acc[w].on_hand[k];
    return tot;
}

static void do_ship(Sim &s, int32_t src, int32_t dst, int32_t sku, int64_t qty, int64_t day, int64_t tick, int64_t t_abs, int bin) {
    if (qty <= 0 || s.acc[src].on_hand[sku] <= 0) return;
    if (qty > s.acc[src].on_hand[sku]) qty = s.acc[src].on_hand[sku];
    int64_t wh_cap = 20000 * s.cap_scale;
    int64_t c_cap = 12000 * s.cap_scale;
    if (dst >= 10 && dst < 14 && wh_total_stock(s, dst) + qty > wh_cap) { s.branches[BR_WH_CAP_REJECT]++; return; }
    int32_t cid = static_cast<int32_t>((src + dst + sku) % NUM_CARRIERS);
    if (s.carrier_day_vol[cid] + qty > c_cap) {
        cid = (cid + 1) % NUM_CARRIERS;
        if (s.carrier_day_vol[cid] + qty > c_cap) { s.branches[BR_CARRIER_CAP_REJECT]++; return; }
    }
    int64_t rate = carrier_rate(s.carrier_day_vol[cid]);
    s.carrier_tier_cnt[rate == 200 ? 0 : (rate == 160 ? 1 : 2)]++;
    s.carrier_day_vol[cid] += qty;
    s.carrier_vol_tot += qty;
    if (bin >= 0 && bin < 100) s.bin_orders[sku][bin] += qty;
    int64_t goods = qty * SKU_WHOLESALE[sku], freight = qty * rate;
    s.acc[src].on_hand[sku] -= qty;
    s.acc[dst].in_transit[sku] += qty;
    transfer_cash(s, dst, src, goods, -1, day, tick);
    transfer_cash(s, dst, 17 + cid, freight, -1, day, tick);
    if (day < 6) s.daily_cost_echelon[day][2] += freight;
    int64_t delay = rand_int(&s.rng[17 + cid], 0, 25);
    if (delay > 18) s.branches[BR_CARRIER_DELAY]++;
    int64_t lt = CARRIER_BASE_LT[cid] + delay;
    if (s.lt_len < 256) { s.leadtime_pairs[s.lt_len][0] = lt; s.leadtime_pairs[s.lt_len][1] = qty; s.lt_len++; }
    if (s.log_enabled) s.ships.push_back({day, tick, src, dst, cid, sku, static_cast<int32_t>(qty), 1, goods, freight});
    heap_push(s.heap, {t_abs + lt, s.next_seq++, EV_DELIVER, src, dst, sku, static_cast<int32_t>(qty), cid});
}

static void hash_ev(Sim &s, const SimEvent &ev) {
    s.event_hash = fnv1a_i64(s.event_hash, ev.time);
    s.event_hash = fnv1a_i64(s.event_hash, ev.etype);
    s.event_hash = fnv1a_i64(s.event_hash, ev.src);
    s.event_hash = fnv1a_i64(s.event_hash, ev.dst);
    s.event_hash = fnv1a_i64(s.event_hash, ev.sku);
    s.event_hash = fnv1a_i64(s.event_hash, ev.qty);
}

void sim_run_days(Sim &s, int64_t days, int64_t ev_per_day, int safety_override, int disrupt_override, int transship_on) {
    int64_t T = s.cfg.ticks_per_day, ev_processed = 0;
    for (int64_t d = 0; d < days; ++d) {
        int64_t day_base = d * T;
        int64_t cs = s.cap_scale;
        s.carrier_cap_tot += NUM_CARRIERS * 12000 * cs;
        for (int c = 0; c < NUM_CARRIERS; ++c) s.carrier_day_vol[c] = 0;
        for (int f = 0; f < NUM_FACTORIES; ++f) {
            if (s.acc[14 + f].plant_state == PLANT_DISRUPTED) {
                if (--s.acc[14 + f].disrupt_rem <= 0) s.acc[14 + f].plant_state = PLANT_RECOVERY;
            } else if (s.acc[14 + f].plant_state == PLANT_RECOVERY) s.acc[14 + f].plant_state = PLANT_OPEN;
        }
        if (d > 0) {
            for (int w = 0; w < 3; ++w) {
                int wh = 10 + w;
                for (int sku = 0; sku < NUM_SKUS; ++sku) {
                    int32_t sts[NUM_STORES], prio[NUM_STORES]; int64_t reqs[NUM_STORES], seqs[NUM_STORES], fills[NUM_STORES]; int n = 0;
                    for (int st = 0; st < NUM_STORES; ++st) {
                        if (STORE_WH[st] != wh) continue;
                        int64_t S_lvl = (safety_override > 0 ? safety_override : 450) * cs;
                        int64_t inv = s.acc[st].on_hand[sku] + s.acc[st].in_transit[sku];
                        if (inv < S_lvl) { sts[n] = st; reqs[n] = S_lvl - inv; seqs[n] = static_cast<int64_t>(s.next_seq++); prio[n++] = (st % 3 == 0 ? 2 : 1); }
                    }
                    if (n > 0) {
                        run_alloc_auction(s, wh, sku, sts, reqs, seqs, prio, n, fills);
                        int b0 = static_cast<int>((d * 100) / days);
                        for (int i = 0; i < n; ++i) {
                            if (s.log_enabled && s.f_alloc) fprintf(s.f_alloc, "%lld,%d,%d,%d,%lld,%lld\n", static_cast<long long>(d), wh, sku, sts[i], static_cast<long long>(reqs[i]), static_cast<long long>(fills[i]));
                            if (fills[i] > 0) do_ship(s, wh, sts[i], sku, fills[i], d, 0, day_base, b0);
                        }
                    }
                }
            }
        }
        int64_t chk_step = ev_per_day / 8; if (chk_step < 1) chk_step = 1;
        int64_t shk_step = ev_per_day / 250; if (shk_step < 1) shk_step = 1;
        for (int64_t e = 0; e < ev_per_day; ++e) {
            int64_t t_abs = day_base + (e * T) / ev_per_day;
            int32_t st = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS], 0, NUM_STORES - 1));
            int32_t sku = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS], 0, NUM_SKUS - 1));
            int32_t qty = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS], 0, 9) == 0 ? rand_int(&s.rng[NUM_ACCOUNTS], 12, 24) : rand_int(&s.rng[NUM_ACCOUNTS], 1, 5));
            heap_push(s.heap, {t_abs, s.next_seq++, EV_DEMAND, ACC_CUST, st, sku, qty, 0});
            if (e % 20 == 0) { int32_t f = static_cast<int32_t>((e / 20) % NUM_FACTORIES); heap_push(s.heap, {t_abs + 1, s.next_seq++, EV_PRODUCE, 14 + f, 10 + f, f, 130, 0}); }
            if (e % 40 == 0) { int32_t kq = static_cast<int32_t>(rand_int(&s.rng[13], 15, 35)); heap_push(s.heap, {t_abs + 2, s.next_seq++, EV_KIT, 13, 13, 3, kq, 0}); }
            if (e % 60 == 0 && transship_on) { int32_t tq = static_cast<int32_t>(rand_int(&s.rng[13], 30, 50)); heap_push(s.heap, {t_abs + 3, s.next_seq++, EV_TRANSSHIP, 10, 13, 0, tq, 0}); }
            if (e % chk_step == 0) {
                int64_t dp = disrupt_override >= 0 ? disrupt_override : s.cfg.disrupt_prob_pct;
                if (rand_int(&s.rng[NUM_ACCOUNTS + 1], 1, 100) <= dp) { int32_t f = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS + 1], 0, NUM_FACTORIES - 1)); heap_push(s.heap, {t_abs + 4, s.next_seq++, EV_DISRUPT, 14 + f, 14 + f, f, 1, 0}); }
                heap_push(s.heap, {t_abs + 6, s.next_seq++, EV_RISK_CHK, 0, 0, 0, 0, 0});
            }
            if (e % shk_step == 0) {
                if (rand_int(&s.rng[NUM_ACCOUNTS + 1], 1, 100) <= s.cfg.shock_prob_pct) { int32_t st_s = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS + 1], 0, NUM_STORES - 1)); int32_t cls = static_cast<int32_t>(rand_int(&s.rng[NUM_ACCOUNTS + 1], 0, 2)); heap_push(s.heap, {t_abs + 5, s.next_seq++, EV_SHOCK, ACC_CUST, st_s, 0, (cls + 1) * 10, cls}); }
            }
        }
        int64_t day_end = (d + 1) * T;
        while (!s.heap.empty() && s.heap[0].time < day_end) {
            SimEvent ev = heap_pop(s.heap);
            int64_t tick = ev.time - day_base;
            int bin = static_cast<int>(((d * T + tick) * 100) / (days * T));
            if (bin < 0) bin = 0; else if (bin > 99) bin = 99;
            hash_ev(s, ev);
            ev_processed++;
            if (s.log_enabled && s.f_orders && ev_processed <= 5000) fprintf(s.f_orders, "%lld,%lld,%llu,%d,%d,%d,%d,%d\n", static_cast<long long>(d), static_cast<long long>(tick), static_cast<unsigned long long>(ev.seq), ev.etype, ev.src, ev.dst, ev.sku, ev.qty);
            if (ev.etype == EV_DEMAND) {
                int32_t st = ev.dst, sku = ev.sku; int64_t qty = ev.qty;
                s.dem_events_tot++; s.dem_units_tot += qty;
                s.bin_demand[sku][bin] += qty; s.bin_store_dem[st][bin] += qty;
                if (d < 6) { s.daily_demand[d][sku] += qty; s.daily_fill_store[st][0] += qty; }
                s.acc[st].ewma_milli[sku] = static_cast<int32_t>((s.cfg.ewma_alpha_milli * qty * 1000 + (1000 - s.cfg.ewma_alpha_milli) * s.acc[st].ewma_milli[sku]) / 1000);
                if (s.acc[st].status == ST_SUSPENDED) s.branches[BR_DEMAND_LOST]++;
                else if (s.acc[st].on_hand[sku] >= qty) {
                    s.acc[st].on_hand[sku] -= qty; s.total_sold[sku] += qty; s.dem_units_filled += qty; int64_t rev = qty * SKU_RETAIL[sku];
                    transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                    if (s.log_enabled) s.sales.push_back({d, tick, st, sku, static_cast<int32_t>(qty), rev});
                    s.branches[BR_DEMAND_FILL]++; if (d < 6) s.daily_fill_store[st][1] += qty;
                } else {
                    int64_t part = s.acc[st].on_hand[sku];
                    if (part > 0) {
                        s.acc[st].on_hand[sku] = 0; s.total_sold[sku] += part; s.dem_units_filled += part; int64_t rev = part * SKU_RETAIL[sku];
                        transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                        if (s.log_enabled) s.sales.push_back({d, tick, st, sku, static_cast<int32_t>(part), rev});
                        if (d < 6) s.daily_fill_store[st][1] += part;
                    }
                    s.acc[st].backorder[sku] += (qty - part); s.dem_bo_events++; s.branches[BR_DEMAND_BO]++;
                }
                int64_t S_lvl = (safety_override > 0 ? safety_override : 450) * cs;
                int64_t s_lvl = (S_lvl * 40) / 100;
                int64_t inv_pos = s.acc[st].on_hand[sku] + s.acc[st].in_transit[sku] - s.acc[st].backorder[sku];
                if (inv_pos < s_lvl && s.acc[st].status == ST_ACTIVE) {
                    s.branches[BR_STORE_REPLEN]++;
                    int32_t src_wh = (sku == 3 || st == 9) ? 13 : (10 + sku);
                    do_ship(s, src_wh, st, sku, S_lvl - inv_pos, d, tick, ev.time, bin);
                }
            } else if (ev.etype == EV_DELIVER) {
                if (s.acc[ev.dst].in_transit[ev.sku] >= ev.qty) { s.acc[ev.dst].in_transit[ev.sku] -= ev.qty; s.acc[ev.dst].on_hand[ev.sku] += ev.qty; }
                if (s.log_enabled) s.ships.push_back({d, tick, ev.src, ev.dst, ev.aux, ev.sku, ev.qty, 2, 0, 0});
            } else if (ev.etype == EV_PRODUCE) {
                int64_t qty = s.acc[ev.src].plant_state == PLANT_DISRUPTED ? 0 : (s.acc[ev.src].plant_state == PLANT_RECOVERY ? ev.qty / 2 : ev.qty);
                if (qty > 0) {
                    int64_t cost = qty * SKU_PROD[ev.sku]; s.acc[ev.src].on_hand[ev.sku] += qty; s.total_produced[ev.sku] += qty;
                    transfer_cash(s, ev.src, ACC_SUPP, cost, -1, d, tick); if (d < 6) s.daily_cost_echelon[d][0] += cost;
                    if (s.log_enabled) s.prods.push_back({d, tick, ev.src, ev.sku, static_cast<int32_t>(qty), cost});
                    do_ship(s, ev.src, ev.dst, ev.sku, qty, d, tick, ev.time, bin);
                }
            } else if (ev.etype == EV_KIT) {
                int32_t wh = 13;
                int64_t m = std::min(s.acc[wh].on_hand[0], std::min(s.acc[wh].on_hand[1], s.acc[wh].on_hand[2]));
                int64_t mx = std::max(s.acc[wh].on_hand[0], std::max(s.acc[wh].on_hand[1], s.acc[wh].on_hand[2]));
                int64_t imb = mx - m; if (s.kit_imb_len < 64) s.kit_imb_series[s.kit_imb_len++] = imb;
                int64_t q = ev.qty;
                if (m >= q) s.branches[BR_KIT_FULL]++; else { q = m > 0 ? m : 0; s.branches[BR_KIT_PARTIAL]++; }
                if (q > 0) {
                    for (int k = 0; k < 3; ++k) s.acc[wh].on_hand[k] -= q;
                    s.acc[wh].on_hand[3] += q; s.total_kits += q; int64_t fee = q * s.cfg.kit_fee;
                    transfer_cash(s, wh, ACC_FEE, fee, -1, d, tick); if (d < 6) s.daily_cost_echelon[d][1] += fee;
                    if (s.log_enabled) s.kits.push_back({d, tick, wh, static_cast<int32_t>(q), static_cast<int32_t>(imb), fee});
                }
            } else if (ev.etype == EV_TRANSSHIP) {
                int min_k = 0;
                for (int k = 1; k < 3; ++k) if (s.acc[13].on_hand[k] < s.acc[13].on_hand[min_k]) min_k = k;
                int best_w = 10;
                for (int w = 11; w < 13; ++w) if (s.acc[w].on_hand[min_k] > s.acc[best_w].on_hand[min_k]) best_w = w;
                s.branches[BR_TRANSSHIP]++; s.n_transships++;
                do_ship(s, best_w, 13, min_k, ev.qty, d, tick, ev.time, bin);
                if (s.log_enabled && s.f_trans) fprintf(s.f_trans, "%lld,%lld,%d,13,%d,%d\n", static_cast<long long>(d), static_cast<long long>(tick), best_w, min_k, ev.qty);
            } else if (ev.etype == EV_DISRUPT) {
                s.acc[ev.src].plant_state = PLANT_DISRUPTED; s.acc[ev.src].disrupt_rem = 1; s.branches[BR_DISRUPT_HIT]++;
                if (s.log_enabled && s.f_disrupt) fprintf(s.f_disrupt, "%lld,%lld,%d,1,1\n", static_cast<long long>(d), static_cast<long long>(tick), ev.src);
            } else if (ev.etype == EV_SHOCK) {
                s.branches[BR_SHOCK_HIT]++;
                int cls = ev.aux;
                int64_t mag = (cls + 1) * 15 * cs;
                s.shock_cnt[cls]++; s.shock_resp[cls] += mag;
                if (s.n_shock_evs < 128) { s.shock_ev_bin[s.n_shock_evs] = bin; s.shock_ev_cls[s.n_shock_evs] = cls; s.shock_ev_st[s.n_shock_evs] = ev.dst; s.n_shock_evs++; }
                if (s.log_enabled && s.f_shocks) fprintf(s.f_shocks, "%lld,%lld,%d,%d,%lld\n", static_cast<long long>(d), static_cast<long long>(tick), cls, ev.dst, static_cast<long long>(mag));
                heap_push(s.heap, {ev.time + 1, s.next_seq++, EV_DEMAND, ACC_CUST, ev.dst, 0, static_cast<int32_t>(mag > 5000 ? 5000 : mag), 0});
            } else if (ev.etype == EV_RISK_CHK) sim_risk_check(s, d, tick);
            if (ev_processed % s.cfg.k_chk == 0) check_invariants(s);
        }
        sim_eod_clearing(s, d); check_invariants(s);
    }
}
