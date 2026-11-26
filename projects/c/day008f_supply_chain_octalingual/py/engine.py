from allocation import run_alloc_auction
from clearing import sim_eod_clearing, sim_risk_check
from inventory import NUM_ACCOUNTS, NUM_CARRIERS, NUM_STORES, ACC_CUST, ACC_SUPP, ACC_FEE, SKU_RETAIL, SKU_WHOLESALE, SKU_PROD, STORE_WH, CARRIER_BASE_LT, carrier_rate, transfer_cash, check_invariants
from rng import fnv1a_i64

def do_ship(s, src, dst, sku, qty, day, tick, t_abs, bin_idx):
    if qty <= 0 or s.acc[src]["on_hand"][sku] <= 0: return
    if qty > s.acc[src]["on_hand"][sku]: qty = s.acc[src]["on_hand"][sku]
    wh_cap = 20_000 * s.cap_scale
    c_cap = 12_000 * s.cap_scale
    if 10 <= dst < 14 and sum(s.acc[dst]["on_hand"]) + qty > wh_cap:
        s.branches[4] += 1; return
    cid = (src + dst + sku) % NUM_CARRIERS
    if s.carrier_day_vol[cid] + qty > c_cap:
        cid = (cid + 1) % NUM_CARRIERS
        if s.carrier_day_vol[cid] + qty > c_cap:
            s.branches[5] += 1; return
    rate = carrier_rate(s.carrier_day_vol[cid])
    s.carrier_tier_cnt[0 if rate == 200 else (1 if rate == 160 else 2)] += 1
    s.carrier_day_vol[cid] += qty
    s.carrier_vol_tot += qty
    if 0 <= bin_idx < 100: s.bin_orders[sku][bin_idx] += qty
    goods = qty * SKU_WHOLESALE[sku]; freight = qty * rate
    s.acc[src]["on_hand"][sku] -= qty; s.acc[dst]["in_transit"][sku] += qty
    transfer_cash(s, dst, src, goods, -1, day, tick)
    transfer_cash(s, dst, 17 + cid, freight, -1, day, tick)
    if day < 6: s.daily_cost_echelon[day][2] += freight
    delay = s.rng[17 + cid].rand_int(0, 25)
    if delay > 18: s.branches[15] += 1
    lt = CARRIER_BASE_LT[cid] + delay
    if len(s.leadtime_pairs) < 256: s.leadtime_pairs.append((lt, qty))
    if s.log_enabled: s.ships.append((day, tick, src, dst, cid, sku, qty, goods, freight, 1))
    s.heap.push((t_abs + lt, s.next_seq, 2, src, dst, sku, qty, cid)); s.next_seq += 1

def hash_ev(s, ev):
    for i in (0, 2, 3, 4, 5, 6): s.event_hash = fnv1a_i64(s.event_hash, ev[i])

def sim_run_days(s, days, ev_per_day, safety_ovr, disrupt_ovr, transship_on):
    t_day = s.cfg["ticks_per_day"]; ev_proc = 0; cs = s.cap_scale
    chk_step = max(ev_per_day // 8, 1)
    shk_step = max(ev_per_day // 250, 1)
    for d in range(days):
        day_base = d * t_day; s.carrier_day_vol = [0, 0, 0]
        s.carrier_cap_tot += NUM_CARRIERS * 12_000 * cs
        for f in range(14, 17):
            if s.acc[f]["plant_state"] == 1:
                s.acc[f]["disrupt_rem"] -= 1
                if s.acc[f]["disrupt_rem"] <= 0: s.acc[f]["plant_state"] = 2
            elif s.acc[f]["plant_state"] == 2:
                s.acc[f]["plant_state"] = 0
        if d > 0:
            for wh in range(10, 13):
                for sku in range(4):
                    sts, reqs, seqs, prio = [], [], [], []
                    s_lvl = (safety_ovr if safety_ovr > 0 else 450) * cs
                    for st in range(NUM_STORES):
                        if STORE_WH[st] != wh: continue
                        inv = s.acc[st]["on_hand"][sku] + s.acc[st]["in_transit"][sku]
                        if inv < s_lvl:
                            sts.append(st); reqs.append(s_lvl - inv); seqs.append(s.next_seq); s.next_seq += 1; prio.append(2 if st % 3 == 0 else 1)
                    if sts:
                        fills = run_alloc_auction(s, wh, sku, sts, reqs, seqs, prio)
                        fa = s.files.get("allocations.csv")
                        b0 = (d * 100) // days
                        for i, st in enumerate(sts):
                            if fa: fa.write(f"{d},{wh},{sku},{st},{reqs[i]},{fills[i]}\n")
                            if fills[i] > 0: do_ship(s, wh, st, sku, fills[i], d, 0, day_base, b0)
        for e in range(ev_per_day):
            t_abs = day_base + (e * t_day) // ev_per_day
            st = s.rng[NUM_ACCOUNTS].rand_int(0, 9)
            sku = s.rng[NUM_ACCOUNTS].rand_int(0, 3)
            qty = s.rng[NUM_ACCOUNTS].rand_int(12, 24) if s.rng[NUM_ACCOUNTS].rand_int(0, 9) == 0 else s.rng[NUM_ACCOUNTS].rand_int(1, 5)
            s.heap.push((t_abs, s.next_seq, 0, ACC_CUST, st, sku, qty, 0)); s.next_seq += 1
            if e % 20 == 0:
                f = (e // 20) % 3
                s.heap.push((t_abs + 1, s.next_seq, 3, 14 + f, 10 + f, f, 130, 0)); s.next_seq += 1
            if e % 40 == 0:
                kq = s.rng[13].rand_int(15, 35)
                s.heap.push((t_abs + 2, s.next_seq, 4, 13, 13, 3, kq, 0)); s.next_seq += 1
            if e % 60 == 0 and transship_on != 0:
                tq = s.rng[13].rand_int(30, 50)
                s.heap.push((t_abs + 3, s.next_seq, 5, 10, 13, 0, tq, 0)); s.next_seq += 1
            if e % chk_step == 0:
                dp = disrupt_ovr if disrupt_ovr >= 0 else s.cfg["disrupt_prob_pct"]
                if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= dp:
                    f = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2)
                    s.heap.push((t_abs + 4, s.next_seq, 7, 14 + f, 14 + f, f, 1, 0)); s.next_seq += 1
                s.heap.push((t_abs + 6, s.next_seq, 8, 0, 0, 0, 0, 0)); s.next_seq += 1
            if e % shk_step == 0:
                if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= s.cfg["shock_prob_pct"]:
                    st_s = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 9); cls = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2)
                    s.heap.push((t_abs + 5, s.next_seq, 6, ACC_CUST, st_s, 0, (cls + 1) * 10, cls)); s.next_seq += 1
        day_end = (d + 1) * t_day
        fo = s.files.get("orders_log.csv")
        while not s.heap.empty() and s.heap.peek_time() < day_end:
            ev = s.heap.pop(); tick = ev[0] - day_base
            bin_idx = max(0, min(99, ((d * t_day + tick) * 100) // (days * t_day)))
            hash_ev(s, ev); ev_proc += 1
            if ev_proc <= 5000 and fo: fo.write(f"{d},{tick},{ev[1]},{ev[2]},{ev[3]},{ev[4]},{ev[5]},{ev[6]}\n")
            et = ev[2]
            if et == 0:
                st, sku, qty = ev[4], ev[5], ev[6]
                s.dem_events_tot += 1; s.dem_units_tot += qty
                s.bin_demand[sku][bin_idx] += qty; s.bin_store_dem[st][bin_idx] += qty
                if d < 6: s.daily_demand[d][sku] += qty; s.daily_fill_store[st][0] += qty
                al = s.cfg["ewma_alpha_milli"]
                s.acc[st]["ewma_milli"][sku] = (al * qty * 1000 + (1000 - al) * s.acc[st]["ewma_milli"][sku]) // 1000
                if s.acc[st]["status"] == 3:
                    s.branches[2] += 1
                elif s.acc[st]["on_hand"][sku] >= qty:
                    s.acc[st]["on_hand"][sku] -= qty; s.total_sold[sku] += qty; s.dem_units_filled += qty; rev = qty * SKU_RETAIL[sku]
                    transfer_cash(s, ACC_CUST, st, rev, -1, d, tick)
                    if s.log_enabled: s.sales.append((d, tick, st, sku, qty, rev))
                    s.branches[0] += 1
                    if d < 6: s.daily_fill_store[st][1] += qty
                else:
                    part = s.acc[st]["on_hand"][sku]
                    if part > 0:
                        s.acc[st]["on_hand"][sku] = 0; s.total_sold[sku] += part; s.dem_units_filled += part; rev = part * SKU_RETAIL[sku]
                        transfer_cash(s, ACC_CUST, st, rev, -1, d, tick)
                        if s.log_enabled: s.sales.append((d, tick, st, sku, part, rev))
                        if d < 6: s.daily_fill_store[st][1] += part
                    s.acc[st]["backorder"][sku] += (qty - part); s.dem_bo_events += 1; s.branches[1] += 1
                cap_s = (safety_ovr if safety_ovr > 0 else 450) * cs
                s_lvl = (cap_s * 40) // 100
                inv_pos = s.acc[st]["on_hand"][sku] + s.acc[st]["in_transit"][sku] - s.acc[st]["backorder"][sku]
                if inv_pos < s_lvl and s.acc[st]["status"] == 0:
                    s.branches[3] += 1
                    src_wh = 13 if (sku == 3 or st == 9) else (10 + sku)
                    do_ship(s, src_wh, st, sku, cap_s - inv_pos, d, tick, ev[0], bin_idx)
            elif et == 2:
                dst, sku, qty = ev[4], ev[5], ev[6]
                if s.acc[dst]["in_transit"][sku] >= qty:
                    s.acc[dst]["in_transit"][sku] -= qty; s.acc[dst]["on_hand"][sku] += qty
                if s.log_enabled: s.ships.append((d, tick, ev[3], dst, ev[7], sku, qty, 0, 0, 2))
            elif et == 3:
                fac, sku = ev[3], ev[5]
                qty = 0 if s.acc[fac]["plant_state"] == 1 else (ev[6] // 2 if s.acc[fac]["plant_state"] == 2 else ev[6])
                if qty > 0:
                    cost = qty * SKU_PROD[sku]; s.acc[fac]["on_hand"][sku] += qty; s.total_produced[sku] += qty
                    transfer_cash(s, fac, ACC_SUPP, cost, -1, d, tick)
                    if d < 6: s.daily_cost_echelon[d][0] += cost
                    if s.log_enabled: s.prods.append((d, tick, fac, sku, qty, cost))
                    do_ship(s, fac, ev[4], sku, qty, d, tick, ev[0], bin_idx)
            elif et == 4:
                m = min(s.acc[13]["on_hand"][:3]); mx = max(s.acc[13]["on_hand"][:3]); imb = mx - m
                if len(s.kit_imb_series) < 64: s.kit_imb_series.append(imb)
                q = ev[6]
                if m >= q: s.branches[8] += 1
                else: q = m if m > 0 else 0; s.branches[9] += 1
                if q > 0:
                    for k in range(3): s.acc[13]["on_hand"][k] -= q
                    s.acc[13]["on_hand"][3] += q; s.total_kits += q; fee = q * s.cfg["kit_fee"]
                    transfer_cash(s, 13, ACC_FEE, fee, -1, d, tick)
                    if d < 6: s.daily_cost_echelon[d][1] += fee
                    if s.log_enabled: s.kits.append((d, tick, 13, q, fee, imb))
            elif et == 5:
                min_k = min(range(3), key=lambda k: s.acc[13]["on_hand"][k])
                best_w = max(range(10, 13), key=lambda w: s.acc[w]["on_hand"][min_k])
                s.branches[10] += 1; s.n_transships += 1
                do_ship(s, best_w, 13, min_k, ev[6], d, tick, ev[0], bin_idx)
                ft = s.files.get("transshipments.csv")
                if ft: ft.write(f"{d},{tick},{best_w},13,{min_k},{ev[6]}\n")
            elif et == 7:
                s.acc[ev[3]]["plant_state"] = 1; s.acc[ev[3]]["disrupt_rem"] = 1; s.branches[11] += 1
                fd = s.files.get("disruptions.csv")
                if fd: fd.write(f"{d},{tick},{ev[3]},1,1\n")
            elif et == 6:
                cls = ev[7]; mag = (cls + 1) * 15 * cs
                s.branches[12] += 1; s.shock_cnt[cls] += 1; s.shock_resp[cls] += mag
                if len(s.shock_evs) < 128: s.shock_evs.append((bin_idx, cls, ev[4]))
                fs = s.files.get("shocks_log.csv")
                if fs: fs.write(f"{d},{tick},{cls},{ev[4]},{mag}\n")
                s.heap.push((ev[0] + 1, s.next_seq, 0, ACC_CUST, ev[4], 0, min(mag, 5000), 0)); s.next_seq += 1
            elif et == 8:
                sim_risk_check(s, d, tick)
            if ev_proc % s.cfg["k_chk"] == 0: check_invariants(s)
        sim_eod_clearing(s, d); check_invariants(s)
