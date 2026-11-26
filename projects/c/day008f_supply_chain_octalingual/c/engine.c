#include "common.h"
#include <stdlib.h>

static void push_sale(Sim *s, SaleRec r) {
    if (!s->log_enabled) return;
    if (s->n_sales == s->cap_sales) {
        s->cap_sales = s->cap_sales ? s->cap_sales * 2 : 1024;
        s->sales = (SaleRec *)realloc(s->sales, (size_t)s->cap_sales * sizeof(SaleRec));
    }
    s->sales[s->n_sales++] = r;
}

static void push_ship(Sim *s, ShipRec r) {
    if (!s->log_enabled) return;
    if (s->n_ships == s->cap_ships) {
        s->cap_ships = s->cap_ships ? s->cap_ships * 2 : 1024;
        s->ships = (ShipRec *)realloc(s->ships, (size_t)s->cap_ships * sizeof(ShipRec));
    }
    s->ships[s->n_ships++] = r;
}

static void push_prod(Sim *s, ProdRec r) {
    if (!s->log_enabled) return;
    if (s->n_prods == s->cap_prods) {
        s->cap_prods = s->cap_prods ? s->cap_prods * 2 : 1024;
        s->prods = (ProdRec *)realloc(s->prods, (size_t)s->cap_prods * sizeof(ProdRec));
    }
    s->prods[s->n_prods++] = r;
}

static void push_kit(Sim *s, KitRec r) {
    if (!s->log_enabled) return;
    if (s->n_kits == s->cap_kits) {
        s->cap_kits = s->cap_kits ? s->cap_kits * 2 : 256;
        s->kits = (KitRec *)realloc(s->kits, (size_t)s->cap_kits * sizeof(KitRec));
    }
    s->kits[s->n_kits++] = r;
}

static int64_t wh_total_stock(const Sim *s, int32_t w) {
    int64_t tot = 0;
    for (int k = 0; k < NUM_SKUS; ++k) tot += s->acc[w].on_hand[k];
    return tot;
}

static void do_ship(Sim *s, int32_t src, int32_t dst, int32_t sku, int64_t qty, int64_t day, int64_t tick, int64_t t_abs) {
    if (qty <= 0 || s->acc[src].on_hand[sku] <= 0) return;
    if (qty > s->acc[src].on_hand[sku]) qty = s->acc[src].on_hand[sku];
    if (dst >= 10 && dst < 14 && wh_total_stock(s, dst) + qty > s->cfg.wh_cap) {
        s->branches[BR_WH_CAP_REJECT]++;
        return;
    }
    int32_t cid = (int32_t)((src + dst + sku) % NUM_CARRIERS);
    if (s->carrier_day_vol[cid] + qty > s->cfg.carrier_cap_daily) {
        cid = (cid + 1) % NUM_CARRIERS;
        if (s->carrier_day_vol[cid] + qty > s->cfg.carrier_cap_daily) {
            s->branches[BR_CARRIER_CAP_REJECT]++;
            return;
        }
    }
    int64_t rate = carrier_rate(s->carrier_day_vol[cid]);
    int tier = rate == 200 ? 0 : (rate == 160 ? 1 : 2);
    s->carrier_tier_cnt[tier]++;
    s->carrier_day_vol[cid] += qty;
    int64_t goods = qty * SKU_WHOLESALE[sku];
    int64_t freight = qty * rate;
    s->acc[src].on_hand[sku] -= qty;
    s->acc[dst].in_transit[sku] += qty;
    transfer_cash(s, dst, src, goods, -1, day, tick);
    transfer_cash(s, dst, 17 + cid, freight, -1, day, tick);
    if (day < 6) s->daily_cost_echelon[day][2] += freight;
    int64_t delay = rand_int(&s->rng[17 + cid], 0, 80);
    if (delay > 50) s->branches[BR_CARRIER_DELAY]++;
    int64_t lt = CARRIER_BASE_LT[cid] + delay;
    if (s->lt_len < 256) { s->leadtime_pairs[s->lt_len][0] = lt; s->leadtime_pairs[s->lt_len][1] = qty; s->lt_len++; }
    push_ship(s, (ShipRec){day, tick, src, dst, cid, sku, (int32_t)qty, 1, goods, freight});
    heap_push(&s->heap, (SimEvent){t_abs + lt, s->next_seq++, EV_DELIVER, src, dst, sku, (int32_t)qty, cid});
}

static void hash_ev(Sim *s, const SimEvent *ev) {
    s->event_hash = fnv1a_i64(s->event_hash, ev->time);
    s->event_hash = fnv1a_i64(s->event_hash, (int64_t)ev->etype);
    s->event_hash = fnv1a_i64(s->event_hash, (int64_t)ev->src);
    s->event_hash = fnv1a_i64(s->event_hash, (int64_t)ev->dst);
    s->event_hash = fnv1a_i64(s->event_hash, (int64_t)ev->sku);
    s->event_hash = fnv1a_i64(s->event_hash, (int64_t)ev->qty);
}

void sim_run_days(Sim *s, int64_t days, int64_t ev_per_day, int safety_override, int disrupt_override, int transship_on) {
    int64_t T = s->cfg.ticks_per_day;
    int64_t step = T / ev_per_day;
    if (step < 1) step = 1;
    int64_t ev_processed = 0;

    for (int64_t d = 0; d < days; ++d) {
        int64_t day_base = d * T;
        for (int c = 0; c < NUM_CARRIERS; ++c) s->carrier_day_vol[c] = 0;
        for (int f = 0; f < NUM_FACTORIES; ++f) {
            if (s->acc[14 + f].plant_state == PLANT_DISRUPTED) {
                if (--s->acc[14 + f].disrupt_rem <= 0) s->acc[14 + f].plant_state = PLANT_RECOVERY;
            } else if (s->acc[14 + f].plant_state == PLANT_RECOVERY) {
                s->acc[14 + f].plant_state = PLANT_OPEN;
            }
        }
        /* Daily allocation auction at start of day */
        for (int w = 0; w < 3; ++w) {
            int wh = 10 + w;
            for (int sku = 0; sku < NUM_SKUS; ++sku) {
                int32_t sts[NUM_STORES]; int64_t reqs[NUM_STORES], seqs[NUM_STORES], fills[NUM_STORES]; int32_t prio[NUM_STORES];
                int n = 0;
                for (int st = 0; st < NUM_STORES; ++st) {
                    if (STORE_WH[st] != wh) continue;
                    int64_t S_lvl = (safety_override > 0 ? safety_override : s->cfg.store_S_base);
                    int64_t inv = s->acc[st].on_hand[sku] + s->acc[st].in_transit[sku];
                    if (inv < S_lvl) {
                        sts[n] = st; reqs[n] = S_lvl - inv; seqs[n] = (int64_t)s->next_seq++; prio[n] = (st % 3 == 0 ? 2 : 1);
                        n++;
                    }
                }
                if (n > 0) {
                    run_alloc_auction(s, wh, sku, sts, reqs, seqs, prio, n, fills, d, 0);
                    for (int i = 0; i < n; ++i) {
                        if (s->log_enabled && s->f_alloc) fprintf(s->f_alloc, "%lld,%d,%d,%d,%lld,%lld\n", (long long)d, wh, sku, sts[i], (long long)reqs[i], (long long)fills[i]);
                        if (fills[i] > 0) do_ship(s, wh, sts[i], sku, fills[i], d, 0, day_base);
                    }
                }
            }
        }
        /* Schedule continuous day events */
        for (int64_t e = 0; e < ev_per_day; ++e) {
            int64_t t_abs = day_base + (e * T) / ev_per_day;
            int32_t st = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS], 0, NUM_STORES - 1);
            int32_t sku = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS], 0, NUM_SKUS - 1);
            int32_t qty = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS], 1, 14);
            heap_push(&s->heap, (SimEvent){t_abs, s->next_seq++, EV_DEMAND, ACC_CUST, st, sku, qty, 0});
            if (e % 20 == 0) {
                int32_t f = (int32_t)(e / 20) % NUM_FACTORIES;
                heap_push(&s->heap, (SimEvent){t_abs + 1, s->next_seq++, EV_PRODUCE, 14 + f, 10 + f, f, 80, 0});
            }
            if (e % 40 == 0) {
                heap_push(&s->heap, (SimEvent){t_abs + 2, s->next_seq++, EV_KIT, 13, 13, 3, 25, 0});
            }
            if (e % 60 == 0 && transship_on) {
                int32_t wdst = (int32_t)((e / 60) % 3) + 10;
                int32_t sku_ts = (int32_t)((e / 60) % NUM_SKUS);
                heap_push(&s->heap, (SimEvent){t_abs + 3, s->next_seq++, EV_TRANSSHIP, 13, wdst, sku_ts, 40, 0});
            }
            if (e % 250 == 0) {
                int64_t dp = disrupt_override >= 0 ? disrupt_override : s->cfg.disrupt_prob_pct;
                if (rand_int(&s->rng[NUM_ACCOUNTS + 1], 1, 100) <= dp) {
                    int32_t f = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS + 1], 0, NUM_FACTORIES - 1);
                    heap_push(&s->heap, (SimEvent){t_abs + 4, s->next_seq++, EV_DISRUPT, 14 + f, 14 + f, f, 1, 0});
                }
                if (rand_int(&s->rng[NUM_ACCOUNTS + 1], 1, 100) <= s->cfg.shock_prob_pct) {
                    int32_t st_s = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS + 1], 0, NUM_STORES - 1);
                    int32_t cls = (int32_t)rand_int(&s->rng[NUM_ACCOUNTS + 1], 0, 2);
                    heap_push(&s->heap, (SimEvent){t_abs + 5, s->next_seq++, EV_SHOCK, ACC_CUST, st_s, 0, (cls + 1) * 10, cls});
                }
                heap_push(&s->heap, (SimEvent){t_abs + 6, s->next_seq++, EV_RISK_CHK, 0, 0, 0, 0, 0});
            }
        }
        /* Process events up to end of day */
        int64_t day_end = (d + 1) * T;
        while (s->heap.len > 0 && s->heap.data[0].time < day_end) {
            SimEvent ev = heap_pop(&s->heap);
            int64_t tick = ev.time - day_base;
            hash_ev(s, &ev);
            ev_processed++;
            if (s->log_enabled && s->f_orders && ev_processed <= 5000) {
                fprintf(s->f_orders, "%lld,%lld,%llu,%d,%d,%d,%d,%d\n", (long long)d, (long long)tick, (unsigned long long)ev.seq, ev.etype, ev.src, ev.dst, ev.sku, ev.qty);
            }
            if (ev.etype == EV_DEMAND) {
                int32_t st = ev.dst, sku = ev.sku;
                int64_t qty = ev.qty;
                if (d < 6) { s->daily_demand[d][sku] += qty; s->daily_fill_store[st][0] += qty; }
                s->acc[st].ewma_milli[sku] = (int32_t)((s->cfg.ewma_alpha_milli * qty * 1000 + (1000 - s->cfg.ewma_alpha_milli) * s->acc[st].ewma_milli[sku]) / 1000);
                if (s->acc[st].status == ST_SUSPENDED) {
                    s->branches[BR_DEMAND_LOST]++;
                } else if (s->acc[st].on_hand[sku] >= qty) {
                    s->acc[st].on_hand[sku] -= qty;
                    s->total_sold[sku] += qty;
                    int64_t rev = qty * SKU_RETAIL[sku];
                    transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                    push_sale(s, (SaleRec){d, tick, st, sku, (int32_t)qty, rev});
                    s->branches[BR_DEMAND_FILL]++;
                    if (d < 6) s->daily_fill_store[st][1] += qty;
                } else {
                    int64_t part = s->acc[st].on_hand[sku];
                    if (part > 0) {
                        s->acc[st].on_hand[sku] = 0;
                        s->total_sold[sku] += part;
                        int64_t rev = part * SKU_RETAIL[sku];
                        transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                        push_sale(s, (SaleRec){d, tick, st, sku, (int32_t)part, rev});
                        if (d < 6) s->daily_fill_store[st][1] += part;
                    }
                    s->acc[st].backorder[sku] += (qty - part);
                    s->branches[BR_DEMAND_BO]++;
                }
                int64_t s_lvl = (safety_override > 0 ? (safety_override * s->cfg.store_s_base) / s->cfg.store_S_base : s->cfg.store_s_base);
                int64_t S_lvl = (safety_override > 0 ? safety_override : s->cfg.store_S_base);
                int64_t inv_pos = s->acc[st].on_hand[sku] + s->acc[st].in_transit[sku] - s->acc[st].backorder[sku];
                if (inv_pos < s_lvl && s->acc[st].status == ST_ACTIVE) {
                    s->branches[BR_STORE_REPLEN]++;
                    do_ship(s, STORE_WH[st], st, sku, S_lvl - inv_pos, d, tick, ev.time);
                }
            } else if (ev.etype == EV_DELIVER) {
                int32_t dst = ev.dst, sku = ev.sku;
                int64_t qty = ev.qty;
                if (s->acc[dst].in_transit[sku] >= qty) {
                    s->acc[dst].in_transit[sku] -= qty;
                    s->acc[dst].on_hand[sku] += qty;
                }
                push_ship(s, (ShipRec){d, tick, ev.src, dst, ev.aux, sku, (int32_t)qty, 2, 0, 0});
            } else if (ev.etype == EV_PRODUCE) {
                int32_t fac = ev.src, sku = ev.sku;
                int64_t qty = ev.qty;
                if (s->acc[fac].plant_state == PLANT_DISRUPTED) qty = 0;
                else if (s->acc[fac].plant_state == PLANT_RECOVERY) qty /= 2;
                if (qty > 0) {
                    int64_t cost = qty * SKU_PROD[sku];
                    s->acc[fac].on_hand[sku] += qty;
                    s->total_produced[sku] += qty;
                    transfer_cash(s, fac, ACC_SUPP, cost, -1, d, tick);
                    if (d < 6) s->daily_cost_echelon[d][0] += cost;
                    push_prod(s, (ProdRec){d, tick, fac, sku, (int32_t)qty, cost});
                    do_ship(s, fac, ev.dst, sku, qty, d, tick, ev.time);
                }
            } else if (ev.etype == EV_KIT) {
                int32_t wh = 13;
                int64_t m = s->acc[wh].on_hand[0];
                if (s->acc[wh].on_hand[1] < m) m = s->acc[wh].on_hand[1];
                if (s->acc[wh].on_hand[2] < m) m = s->acc[wh].on_hand[2];
                int64_t mx = s->acc[wh].on_hand[0];
                if (s->acc[wh].on_hand[1] > mx) mx = s->acc[wh].on_hand[1];
                if (s->acc[wh].on_hand[2] > mx) mx = s->acc[wh].on_hand[2];
                int64_t imb = mx - m;
                if (s->kit_imb_len < 64) s->kit_imb_series[s->kit_imb_len++] = imb;
                int64_t q = ev.qty;
                if (m >= q) {
                    s->branches[BR_KIT_FULL]++;
                } else if (m > 0) {
                    q = m;
                    s->branches[BR_KIT_PARTIAL]++;
                } else {
                    q = 0;
                    s->branches[BR_KIT_PARTIAL]++;
                }
                if (q > 0) {
                    for (int k = 0; k < 3; ++k) s->acc[wh].on_hand[k] -= q;
                    s->acc[wh].on_hand[3] += q;
                    s->total_kits += q;
                    int64_t fee = q * s->cfg.kit_fee;
                    transfer_cash(s, wh, ACC_FEE, fee, -1, d, tick);
                    if (d < 6) s->daily_cost_echelon[d][1] += fee;
                    push_kit(s, (KitRec){d, tick, wh, (int32_t)q, (int32_t)imb, fee});
                }
            } else if (ev.etype == EV_TRANSSHIP) {
                s->branches[BR_TRANSSHIP]++;
                do_ship(s, ev.src, ev.dst, ev.sku, ev.qty, d, tick, ev.time);
                if (s->log_enabled && s->f_trans) fprintf(s->f_trans, "%lld,%lld,%d,%d,%d,%d\n", (long long)d, (long long)tick, ev.src, ev.dst, ev.sku, ev.qty);
            } else if (ev.etype == EV_DISRUPT) {
                s->acc[ev.src].plant_state = PLANT_DISRUPTED;
                s->acc[ev.src].disrupt_rem = 1;
                s->branches[BR_DISRUPT_HIT]++;
                if (s->log_enabled && s->f_disrupt) fprintf(s->f_disrupt, "%lld,%lld,%d,1,1\n", (long long)d, (long long)tick, ev.src);
            } else if (ev.etype == EV_SHOCK) {
                s->branches[BR_SHOCK_HIT]++;
                s->acc[ev.dst].on_hand[0] += 0;
                int cls = ev.aux;
                s->shock_cnt[cls]++;
                s->shock_resp[cls] += ev.qty;
                if (s->log_enabled && s->f_shocks) fprintf(s->f_shocks, "%lld,%lld,%d,%d,%d\n", (long long)d, (long long)tick, cls, ev.dst, ev.qty);
                heap_push(&s->heap, (SimEvent){ev.time + 1, s->next_seq++, EV_DEMAND, ACC_CUST, ev.dst, 0, ev.qty, 0});
            } else if (ev.etype == EV_RISK_CHK) {
                sim_risk_check(s, d, tick);
            }
            if (ev_processed % s->cfg.k_chk == 0) check_invariants(s);
        }
        sim_eod_clearing(s, d);
        check_invariants(s);
    }
}
