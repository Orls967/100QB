import pathlib, struct, time
from allocation import run_alloc_auction
from engine import sim_run_days
from entities import Sim
from inventory import NUM_ACCOUNTS, ACC_CUST, ACC_SUPP, ACC_FEE, CARRIER_BASE_LT, STORE_WH, check_invariants
from rng import M64, splitmix64, fnv1a_bytes, fnv1a_u32pair, fnv1a_i64, Xoshiro256

SKU_L = ("A", "B", "C", "D")
def wf(d, n, txt): (pathlib.Path(d) / n).write_text(txt)

def write_s01_to_s09(c, out_dir):
    hrng = 0xCBF29CE484222325
    pairs = ((0,0),(1,2),(2,1))
    exp_hv = (0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb)
    rng_ok = True
    for i in range(3):
        x = (c["master_seed"] & M64) ^ fnv1a_u32pair(*pairs[i])
        sd, _ = splitmix64(x); st = Xoshiro256(sd)
        vals = [st.next_u64() for _ in range(4)]
        hv = fnv1a_bytes(struct.pack("<4Q", *vals))
        if hv != exp_hv[i]: rng_ok = False
        hrng = fnv1a_i64(hrng, hv)
    s = Sim(c, False, None)
    fills = run_alloc_auction(s, 10, 0, [0,1,2], [1000,1500,1000], [1,2,3], [2,1,1])
    halloc = 0xCBF29CE484222325
    for fl in fills: halloc = fnv1a_i64(halloc, fl)
    inv_ok = check_invariants(s)
    ok = 1 if (rng_ok and inv_ok == 1 and sum(fills) == 1600) else 0
    wf(out_dir, "s01_selftest.txt", f"ok,h_rng,h_alloc,inv_ok\n{ok},{hrng:016x},{halloc:016x},{inv_ok}\n")
    wf(out_dir, "s02_config.csv", f"key,value\nmaster_seed,{c['master_seed']}\nnum_days,{c['num_days']}\nparity_days,{c['parity_days']}\nparity_events_per_day,{c['parity_events_per_day']}\nwh_cap,{c['wh_cap']}\nkit_fee,{c['kit_fee']}\n")

    ent = ["id,type,name"] + [f"{i},STORE,S{i}" for i in range(10)] + [f"{10+i},WAREHOUSE,W{i}" for i in range(4)] + [f"{14+i},FACTORY,F{i}" for i in range(3)] + [f"{17+i},CARRIER,C{i}" for i in range(3)] + ["20,CUSTOMER,ExtCustomer","21,SUPPLIER,RawSupplier","22,FEEPOOL,FeePool"]
    wf(out_dir, "s03_entities.csv", "\n".join(ent) + "\n")
    wf(out_dir, "s03_accounts_init.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n" + "\n".join(f"{i},{a['cash']},{','.join(map(str,a['on_hand']))},{a['credit_limit']},{a['status']}" for i, a in enumerate(s.acc)) + "\n")
    wf(out_dir, "s03_carriers.csv", "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n" + "\n".join(f"{i},{CARRIER_BASE_LT[i]},{c['carrier_cap_daily']},200,160,120" for i in range(3)) + "\n")
    edges = ["src,dst,edge_type"] + [f"{14+i},{10+i},FAC_WH" for i in range(3)] + [f"13,{10+i},KIT_TS" for i in range(3)] + [f"{STORE_WH[i]},{i},WH_STORE" for i in range(10)]
    wf(out_dir, "s03_network_edges.csv", "\n".join(edges) + "\n")

    pol = ["case_id,inv_pos,s_lvl,S_lvl,order_qty,branch"]
    for k in range(60):
        inv = (k * 37) % 1000 - 100; oq = (800 - inv) if (inv < 250 and k % 11 != 0) else 0
        br = "SUSPENDED" if k % 11 == 0 else ("REORDER" if oq > 0 else "NO_OP")
        pol.append(f"{k},{inv},250,800,{oq},{br}")
    wf(out_dir, "s04_policy_script.csv", "\n".join(pol) + "\n")

    al = ["case_id,avail,f0,f1,f2,f3"]
    for k in range(10):
        sa = Sim(c, False, None); sa.acc[10]["on_hand"][0] = 200 + k * 150
        fl = run_alloc_auction(sa, 10, 0, [0,1,2,3], [400,300+k*20,500,250], [10+k,5,20-k,8], [2,1,2,1])
        al.append(f"{k},{sa.acc[10]['on_hand'][0]},{','.join(map(str,fl))}")
    wf(out_dir, "s05_alloc_script.csv", "\n".join(al) + "\n")

    sess = ["day,session,start_tick,end_tick"]
    for d in range(c["num_days"]):
        sess += [f"{d},PRE_DAY,0,0", f"{d},ALLOCATION_AUCTION,0,1", f"{d},CONTINUOUS_OPS,0,{c['ticks_per_day']}", f"{d},EOD_CLEARING,{c['ticks_per_day']},{c['ticks_per_day']}"]
    wf(out_dir, "s06_sessions.csv", "\n".join(sess) + "\n")
    x0 = (c["master_seed"] & M64) ^ fnv1a_u32pair(2, 23); sd0, _ = splitmix64(x0); r0 = Xoshiro256(sd0)
    sh = ["idx,time,seq,etype,dst,sku,qty"] + [f"{i},{(i*c['ticks_per_day'])//8000},{i+1},0,{r0.rand_int(0,9)},{r0.rand_int(0,3)},{r0.rand_int(1,14)}" for i in range(200)]
    wf(out_dir, "s06_schedule_head.csv", "\n".join(sh) + "\n")
    wf(out_dir, "s07_shock_schedule.csv", "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n")

    sw = Sim(c, False, None); sim_run_days(sw, 1, 2000, 0, 0, 1)
    wf(out_dir, "s08_warmup_summary.csv", f"days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,{','.join(map(str,sw.total_sold))},{int(sw.inv_violations==0)}\n")
    for k in range(4): wf(out_dir, f"s08_warmup_sku_{SKU_L[k]}.csv", f"sku,sold,produced,end_snap\n{SKU_L[k]},{sw.total_sold[k]},{sw.total_produced[k]},{sw.day_inv_snap[0][k]}\n")
    wf(out_dir, "s09_day1_summary.csv", f"sold_tot,kits,alloc_cases,inv_violations\n{sum(sw.total_sold)},{sw.total_kits},{sw.alloc_cases},{sw.inv_violations}\n")
    cal_fr = sw.dem_units_filled / sw.dem_units_tot if sw.dem_units_tot > 0 else 1.0
    cal_bo = sw.dem_bo_events / sw.dem_events_tot if sw.dem_events_tot > 0 else 0.0
    cal_cu = sw.carrier_vol_tot / sw.carrier_cap_tot if sw.carrier_cap_tot > 0 else 0.0
    cal_ok = int(0.85 <= cal_fr <= 0.99 and cal_bo < 0.20 and 0.30 <= cal_cu <= 0.90)
    wf(out_dir, "s09_calibration.csv", f"fill_rate,bo_event_frac,carrier_util,calibrated\n{cal_fr:.6f},{cal_bo:.6f},{cal_cu:.6f},{cal_ok}\n")

def write_s11_to_s38(s, c, out_dir):
    for k in range(4):
        wf(out_dir, f"sales_{SKU_L[k]}.csv", "day,tick,store_id,sku,qty,rev_cents\n" + "\n".join(",".join(map(str,r)) for r in s.sales if r[3] == k) + "\n")
        wf(out_dir, f"shipments_{SKU_L[k]}.csv", "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status\n" + "\n".join(",".join(map(str,r)) for r in s.ships if r[5] == k) + "\n")
        wf(out_dir, f"inventory_snap_{SKU_L[k]}.csv", "day,sku,on_hand_tot\n" + "\n".join(f"{d},{SKU_L[k]},{s.day_inv_snap[d][k]}" for d in range(c["num_days"])) + "\n")
    for i in range(10): wf(out_dir, f"store_daily_{i}.csv", f"store,cash,oh_tot,req_tot,fill_tot\n{i},{s.acc[i]['cash']},{sum(s.acc[i]['on_hand'])},{s.daily_fill_store[i][0]},{s.daily_fill_store[i][1]}\n")
    for i in range(4): a = s.acc[10+i]; wf(out_dir, f"warehouse_daily_{i}.csv", f"wh,cash,oh_A,oh_B,oh_C,oh_D\n{i},{a['cash']},{','.join(map(str,a['on_hand']))}\n")
    for i in range(3): wf(out_dir, f"factory_daily_{i}.csv", f"fac,cash,produced,state\n{i},{s.acc[14+i]['cash']},{s.total_produced[i]},{s.acc[14+i]['plant_state']}\n")
    wf(out_dir, "production.csv", "day,tick,factory_id,sku,qty,cost_cents\n" + "\n".join(",".join(map(str,r)) for r in s.prods) + "\n")
    wf(out_dir, "kit_assembly.csv", "day,tick,wh_id,qty,fee_cents,imbalance\n" + "\n".join(",".join(map(str,r)) for r in s.kits) + "\n")
    wf(out_dir, "fees_ledger.csv", "day,tick,from_id,to_id,amount,kind\n" + "\n".join(",".join(map(str,r)) for r in s.fees) + "\n")
    wf(out_dir, "entities_final.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n" + "\n".join(f"{i},{a['cash']},{','.join(map(str,a['on_hand']))},{','.join(map(str,a['in_transit']))},{a['status']}" for i, a in enumerate(s.acc)) + "\n")

    wf(out_dir, "s11_invariants.csv", f"violations,ok\n{s.inv_violations},{int(s.inv_violations==0)}\n")
    wf(out_dir, "s12_clearing.csv", f"days,fee_pool_cash\n{c['num_days']},{s.acc[ACC_FEE]['cash']}\n")
    wf(out_dir, "s12_credit.csv", f"credit_holds,liquidations\n{s.branches[13]},{s.branches[14]}\n")
    wf(out_dir, "s13_risk.csv", f"credit_holds,liquidations\n{s.branches[13]},{s.branches[14]}\n")
    wf(out_dir, "s14_factories.csv", f"fac,produced\n0,{s.total_produced[0]}\n1,{s.total_produced[1]}\n2,{s.total_produced[2]}\n")
    wf(out_dir, "s15_carrier_fees.csv", f"tier,count\n0,{s.carrier_tier_cnt[0]}\n1,{s.carrier_tier_cnt[1]}\n2,{s.carrier_tier_cnt[2]}\n")
    wf(out_dir, "s16_kits.csv", f"total_kits,full,partial\n{s.total_kits},{s.branches[8]},{s.branches[9]}\n")
    wf(out_dir, "s16_kit_imbalance_series.csv", "idx,imbalance\n" + "\n".join(f"{i},{v}" for i, v in enumerate(s.kit_imb_series)) + "\n")
    wf(out_dir, "s17_allocations.csv", f"cases,nonzero,total_req,total_fill\n{s.alloc_cases},{s.alloc_nonzero},{s.alloc_total_req},{s.alloc_total_fill}\n")
    wf(out_dir, "s18_cost_entities.csv", "acc,cash_delta\n" + "\n".join(f"{i},{s.acc[i]['cash']-(2_000_000_000 if i<10 else 5_000_000_000)*s.cap_scale}" for i in range(17)) + "\n")
    wf(out_dir, "s19_cost_types.csv", "echelon,cost_cents\n" + "\n".join(f"{e},{sum(s.daily_cost_echelon[d][e] for d in range(c['num_days']))}" for e in range(4)) + "\n")
    wf(out_dir, "s19_cost_daily.csv", "day,prod,kit,freight,hold_pen\n" + "\n".join(f"{d},{','.join(map(str,s.daily_cost_echelon[d]))}" for d in range(c["num_days"])) + "\n")
    for k in range(4): wf(out_dir, f"s20_demand_{SKU_L[k]}.csv", "day,demand\n" + "\n".join(f"{d},{s.daily_demand[d][k]}" for d in range(c["num_days"])) + "\n")
    bw_rows = []
    for k in range(4):
        md = sum(s.bin_demand[k]) / 100.0; mo = sum(s.bin_orders[k]) / 100.0
        m2 = sum((x - md)**2 for x in s.bin_demand[k])
        m4 = sum((x - md)**4 for x in s.bin_demand[k])
        vo = sum((x - mo)**2 for x in s.bin_orders[k])
        cov1 = sum((s.bin_demand[k][b] - md) * (s.bin_demand[k][b - 1] - md) for b in range(1, 100))
        kurt = (m4 / 100.0) / ((m2 / 100.0)**2) if m2 > 0 else 0.0
        acf1 = cov1 / m2 if m2 > 0 else 0.0
        bwr = vo / m2 if m2 > 0 else 1.0
        bw_rows.append(f"{SKU_L[k]},{kurt:.6f},{acf1:.6f},{bwr:.6f}")
    wf(out_dir, "s21_bullwhip.csv", "sku,kurtosis,acf_lag1,bullwhip_ratio\n" + "\n".join(bw_rows) + "\n")
    wf(out_dir, "s22_service.csv", "store,req,fill,fill_rate\n" + "\n".join(f"{i},{s.daily_fill_store[i][0]},{s.daily_fill_store[i][1]},{s.daily_fill_store[i][1]/s.daily_fill_store[i][0] if s.daily_fill_store[i][0]>0 else 1.0:.6f}" for i in range(10)) + "\n")
    wf(out_dir, "s23_leadtime_impact.csv", "idx,leadtime,qty\n" + "\n".join(f"{i},{p[0]},{p[1]}" for i, p in enumerate(s.leadtime_pairs)) + "\n")
    imb_m = sum(s.kit_imb_series)/len(s.kit_imb_series) if s.kit_imb_series else 0.0
    wf(out_dir, "s24_kit_coupling.csv", f"mean_imbalance,total_kits\n{imb_m:.4f},{s.total_kits}\n")
    shk_rows = []
    for cl in range(3):
        vals_r = []
        for b, cev, st in s.shock_evs:
            if cev != cl: continue
            pre = [s.bin_store_dem[st][b + db] for db in (-2, -1) if b + db >= 0]
            post = [s.bin_store_dem[st][b + db] for db in (0, 1) if b + db < 100]
            if pre and post and sum(pre) > 0:
                vals_r.append((sum(post)/len(post) - sum(pre)/len(pre)) / (sum(pre)/len(pre)))
        cnt = len(vals_r); mn = sum(vals_r)/cnt if cnt > 0 else 0.0
        se = (max(0.0, (sum(x*x for x in vals_r)/cnt - mn*mn)/(cnt - 1)))**0.5 if cnt > 1 else 0.0
        shk_rows.append(f"{cl},{cnt},{mn:.6f},{se:.6f}")
    wf(out_dir, "s25_shock_event_study.csv", "cls,count,mean_rel_dev,se_rel_dev\n" + "\n".join(shk_rows) + "\n")
    wf(out_dir, "s26_disruptions.csv", f"disrupt_hits\n{s.branches[11]}\n")
    wf(out_dir, "s27_fees.csv", f"n_fee_entries,fee_pool_cash\n{len(s.fees)},{s.acc[ACC_FEE]['cash']}\n")
    wf(out_dir, "s28_capacity_rejects.csv", f"wh_cap_reject,carrier_cap_reject\n{s.branches[4]},{s.branches[5]}\n")

    sw1 = ["S_lvl,rep,fill_rate,avg_oh"]
    for sv in (80, 220, 350):
        for rp in range(5):
            c2 = dict(c); c2["master_seed"] += rp + 1; ss = Sim(c2, False, None); sim_run_days(ss, 2, 1600, sv, -1, 1)
            oh = sum(ss.acc[i]["on_hand"][0] for i in range(10))
            fr = ss.dem_units_filled / ss.dem_units_tot if ss.dem_units_tot > 0 else 1.0
            sw1.append(f"{sv},{rp},{fr:.6f},{oh/10.0:.2f}")
    wf(out_dir, "s29_sweep_safety_stock.csv", "\n".join(sw1) + "\n")
    sw2 = ["disrupt_pct,rep,produced,sold"]
    for dp in (0,4,12):
        for rp in range(5):
            c2 = dict(c); c2["master_seed"] += rp + 10; ss = Sim(c2, False, None); sim_run_days(ss, 2, 1600, 0, dp, 1)
            sw2.append(f"{dp},{rp},{sum(ss.total_produced)},{sum(ss.total_sold)}")
    wf(out_dir, "s30_sweep_shock_disrupt.csv", "\n".join(sw2) + "\n")
    sw3 = ["transship_on,rep,mean_imbalance,n_transships"]
    for ts in (0,1):
        for rp in range(5):
            c2 = dict(c); c2["master_seed"] += rp + 20; ss = Sim(c2, False, None); sim_run_days(ss, 2, 1600, 0, -1, ts)
            im = sum(ss.kit_imb_series)/len(ss.kit_imb_series) if ss.kit_imb_series else 0.0
            sw3.append(f"{ts},{rp},{im:.4f},{ss.n_transships}")
    wf(out_dir, "s31_sweep_transship.csv", "\n".join(sw3) + "\n")

    bn = ("DEMAND_FILL","DEMAND_BO","DEMAND_LOST","STORE_REPLEN","WH_CAP_REJECT","CARRIER_CAP_REJECT","ALLOC_FULL","ALLOC_PARTIAL","KIT_FULL","KIT_PARTIAL","TRANSSHIP","DISRUPT_HIT","SHOCK_HIT","CREDIT_HOLD","LIQUIDATE","CARRIER_DELAY")
    wf(out_dir, "s32_branches.csv", "branch,count\n" + "\n".join(f"{nm},{s.branches[i]}" for i, nm in enumerate(bn)) + "\n")

    sp1 = Sim(c, False, None, c["parity_events_per_day"]); sim_run_days(sp1, c["parity_days"], c["parity_events_per_day"], 0, -1, 1)
    sp2 = Sim(c, False, None, c["parity_events_per_day"]); sim_run_days(sp2, c["parity_days"], c["parity_events_per_day"], 0, -1, 1)
    h1 = h2 = 0xCBF29CE484222325
    for i in range(NUM_ACCOUNTS):
        h1 = fnv1a_i64(h1, sp1.acc[i]["cash"]); h2 = fnv1a_i64(h2, sp2.acc[i]["cash"])
        for k in range(4): h1 = fnv1a_i64(h1, sp1.acc[i]["on_hand"][k]); h2 = fnv1a_i64(h2, sp2.acc[i]["on_hand"][k])
    wf(out_dir, "s33_replay.csv", f"ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n{sp1.event_hash:016x},{sp2.event_hash:016x},{h1:016x},{h2:016x},{int(sp1.event_hash==sp2.event_hash and h1==h2)}\n")
    wf(out_dir, "s34_parity.csv", f"parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n{c['parity_days']},{c['parity_events_per_day']},{sp1.event_hash:016x},{h1:016x},{sum(sp1.total_sold)},{sp1.total_kits},{sp1.inv_violations}\n")

    vals = []
    for _ in range(5): sb = Sim(c, False, None); t0 = time.perf_counter_ns(); sim_run_days(sb, 1, 2000, 0, -1, 1); vals.append((time.perf_counter_ns()-t0)/2000.0)
    vals.sort(); med = vals[2]; dev = sorted(abs(v - med) for v in vals)
    wf(out_dir, "s35_bench.csv", f"events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n2000,{med:.4f},{dev[2]:.4f},{1e9/med:.2f}\n")
    wf(out_dir, "s36_sizes.csv", f"n_sales,n_ships,n_prods,n_kits,n_fees\n{len(s.sales)},{len(s.ships)},{len(s.prods)},{len(s.kits)},{len(s.fees)}\n")

    raw = [r[4] for r in s.sales[:4096]]; nq = len(raw); rle = []; i = 0
    while i < nq:
        dv = raw[0] if i == 0 else raw[i]-raw[i-1]; run = 1
        while i + run < nq and raw[i+run]-raw[i+run-1] == dv: run += 1
        rle.extend((dv, run)); i += run
    dec = []
    for j in range(0, len(rle), 2):
        for _ in range(rle[j+1]): dec.append(rle[j] if not dec else dec[-1] + rle[j])
    wf(out_dir, "s37_compress.csv", f"n_items,rle_ints,ratio,roundtrip_ok\n{nq},{len(rle)},{nq/len(rle) if rle else 1.0:.4f},{int(dec==raw)}\n")

    si = Sim(c, False, None, s.cap_scale * 2000); lc = [a["cash"] for a in si.acc]; loh = [list(a["on_hand"]) for a in si.acc]; lit = [[0]*4 for _ in range(NUM_ACCOUNTS)]
    for r in s.sales: lc[ACC_CUST] -= r[5]; lc[r[2]] += r[5]; loh[r[2]][r[3]] -= r[4]
    for r in s.ships:
        if r[9] == 1: loh[r[2]][r[5]] -= r[6]; lit[r[3]][r[5]] += r[6]; lc[r[3]] -= r[7]+r[8]; lc[r[2]] += r[7]; lc[17+r[4]] += r[8]
        else: lit[r[3]][r[5]] -= r[6]; loh[r[3]][r[5]] += r[6]
    for r in s.prods: lc[r[2]] -= r[5]; lc[ACC_SUPP] += r[5]; loh[r[2]][r[3]] += r[4]
    for r in s.kits: lc[r[2]] -= r[4]; lc[ACC_FEE] += r[4]; loh[r[2]][0] -= r[3]; loh[r[2]][1] -= r[3]; loh[r[2]][2] -= r[3]; loh[r[2]][3] += r[3]
    for r in s.fees: lc[r[2]] -= r[4]; lc[r[3]] += r[4]
    mcd = max(abs(lc[j] - s.acc[j]["cash"]) for j in range(NUM_ACCOUNTS))
    mid = max(abs(loh[j][k] - s.acc[j]["on_hand"][k]) + abs(lit[j][k] - s.acc[j]["in_transit"][k]) for j in range(NUM_ACCOUNTS) for k in range(4))
    wf(out_dir, "s38_ledger_replay.csv", f"max_cash_diff,max_inv_diff,match\n{mcd},{mid},{int(mcd==0 and mid==0)}\n")
