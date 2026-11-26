use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::Instant;
use crate::allocation::run_alloc_auction;
use crate::engine::sim_run_days;
use crate::entities::{Config, Sim};
use crate::inventory::*;
use crate::rng::*;

const SKU_L: [&str; 4] = ["A", "B", "C", "D"];
fn open_bw(dir: &str, name: &str) -> BufWriter<File> { BufWriter::new(File::create(format!("{}/{}", dir, name)).unwrap()) }

pub fn write_s01_to_s09(c: &Config, out_dir: &str) {
    let mut hrng = 0xCBF29CE484222325u64;
    let pairs = [(0u32, 0u32), (1, 2), (2, 1)];
    let exp_hv = [0xbe40b757fe0b6697u64, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb];
    let mut rng_ok = true;
    for i in 0..3 {
        let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(pairs[i].0, pairs[i].1);
        let mut st = Xoshiro256::seed(splitmix64(&mut x));
        let mut buf = [0u8; 32];
        for j in 0..4 { buf[j*8..(j+1)*8].copy_from_slice(&st.next().to_le_bytes()); }
        let hv = fnv1a_bytes(&buf);
        if hv != exp_hv[i] { rng_ok = false; }
        hrng = fnv1a_i64(hrng, hv as i64);
    }
    let mut s = Sim::new(c, false, None);
    let fills = run_alloc_auction(&mut s, 10, 0, &[0, 1, 2], &[1000, 1500, 1000], &[1, 2, 3], &[2, 1, 1]);
    let mut halloc = 0xCBF29CE484222325u64;
    for &fl in &fills { halloc = fnv1a_i64(halloc, fl); }
    let inv_ok = check_invariants(&mut s);
    let ok = if rng_ok && inv_ok == 1 && fills.iter().sum::<i64>() == 1600 { 1 } else { 0 };
    let mut f = open_bw(out_dir, "s01_selftest.txt");
    let _ = writeln!(f, "ok,h_rng,h_alloc,inv_ok\n{},{:016x},{:016x},{}", ok, hrng, halloc, inv_ok);

    let mut f = open_bw(out_dir, "s02_config.csv");
    let _ = writeln!(f, "key,value\nmaster_seed,{}\nnum_days,{}\nparity_days,{}\nparity_events_per_day,{}\nwh_cap,{}\nkit_fee,{}", c.master_seed, c.num_days, c.parity_days, c.parity_events_per_day, c.wh_cap, c.kit_fee);

    let mut f = open_bw(out_dir, "s03_entities.csv"); let _ = writeln!(f, "id,type,name");
    for i in 0..10 { let _ = writeln!(f, "{},STORE,S{}", i, i); }
    for i in 10..14 { let _ = writeln!(f, "{},WAREHOUSE,W{}", i, i - 10); }
    for i in 14..17 { let _ = writeln!(f, "{},FACTORY,F{}", i, i - 14); }
    for i in 17..20 { let _ = writeln!(f, "{},CARRIER,C{}", i, i - 17); }
    let _ = writeln!(f, "20,CUSTOMER,ExtCustomer\n21,SUPPLIER,RawSupplier\n22,FEEPOOL,FeePool");

    let mut f = open_bw(out_dir, "s03_accounts_init.csv"); let _ = writeln!(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status");
    for (i, a) in s.acc.iter().enumerate() { let _ = writeln!(f, "{},{},{},{},{},{},{},{}", i, a.cash, a.on_hand[0], a.on_hand[1], a.on_hand[2], a.on_hand[3], a.credit_limit, a.status); }
    let mut f = open_bw(out_dir, "s03_carriers.csv"); let _ = writeln!(f, "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2");
    for i in 0..3 { let _ = writeln!(f, "{},{},{},200,160,120", i, CARRIER_BASE_LT[i], c.carrier_cap_daily); }
    let mut f = open_bw(out_dir, "s03_network_edges.csv"); let _ = writeln!(f, "src,dst,edge_type");
    for i in 0..3 { let _ = writeln!(f, "{},{},FAC_WH", 14 + i, 10 + i); }
    for i in 0..3 { let _ = writeln!(f, "13,{},KIT_TS", 10 + i); }
    for i in 0..10 { let _ = writeln!(f, "{},{},WH_STORE", STORE_WH[i], i); }

    let mut f = open_bw(out_dir, "s04_policy_script.csv"); let _ = writeln!(f, "case_id,inv_pos,s_lvl,S_lvl,order_qty,branch");
    for k in 0..60i64 {
        let inv = (k * 37) % 1000 - 100; let oq = if inv < 250 && k % 11 != 0 { 800 - inv } else { 0 };
        let br = if k % 11 == 0 { "SUSPENDED" } else if oq > 0 { "REORDER" } else { "NO_OP" };
        let _ = writeln!(f, "{},{},250,800,{},{}", k, inv, oq, br);
    }
    let mut f = open_bw(out_dir, "s05_alloc_script.csv"); let _ = writeln!(f, "case_id,avail,f0,f1,f2,f3");
    for k in 0..10i64 {
        let mut sa = Sim::new(c, false, None); sa.acc[10].on_hand[0] = 200 + k * 150;
        let fl = run_alloc_auction(&mut sa, 10, 0, &[0, 1, 2, 3], &[400, 300 + k * 20, 500, 250], &[10 + k, 5, 20 - k, 8], &[2, 1, 2, 1]);
        let _ = writeln!(f, "{},{},{},{},{},{}", k, sa.acc[10].on_hand[0], fl[0], fl[1], fl[2], fl[3]);
    }
    let mut f = open_bw(out_dir, "s06_sessions.csv"); let _ = writeln!(f, "day,session,start_tick,end_tick");
    for d in 0..c.num_days { let _ = writeln!(f, "{},PRE_DAY,0,0\n{},ALLOCATION_AUCTION,0,1\n{},CONTINUOUS_OPS,0,{}\n{},EOD_CLEARING,{},{}", d, d, d, c.ticks_per_day, d, c.ticks_per_day, c.ticks_per_day); }
    let mut f = open_bw(out_dir, "s06_schedule_head.csv"); let _ = writeln!(f, "idx,time,seq,etype,dst,sku,qty");
    let mut x0 = (c.master_seed as u64) ^ fnv1a_u32pair(2, 23); let mut r0 = Xoshiro256::seed(splitmix64(&mut x0));
    for i in 0..200i64 { let _ = writeln!(f, "{},{},{},0,{},{},{}", i, (i * c.ticks_per_day) / 8000, i + 1, r0.rand_int(0, 9), r0.rand_int(0, 3), r0.rand_int(1, 14)); }
    let mut f = open_bw(out_dir, "s07_shock_schedule.csv"); let _ = writeln!(f, "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30");

    let mut sw = Sim::new(c, false, None); sim_run_days(&mut sw, 1, 2000, 0, 0, 1);
    let mut f = open_bw(out_dir, "s08_warmup_summary.csv"); let _ = writeln!(f, "days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,{},{},{},{},{}", sw.total_sold[0], sw.total_sold[1], sw.total_sold[2], sw.total_sold[3], (sw.inv_violations == 0) as i32);
    for k in 0..4 { let mut f = open_bw(out_dir, &format!("s08_warmup_sku_{}.csv", SKU_L[k])); let _ = writeln!(f, "sku,sold,produced,end_snap\n{},{},{},{}", SKU_L[k], sw.total_sold[k], sw.total_produced[k], sw.day_inv_snap[0][k]); }
    let mut f = open_bw(out_dir, "s09_day1_summary.csv"); let _ = writeln!(f, "sold_tot,kits,alloc_cases,inv_violations\n{},{},{},{}", sw.total_sold.iter().sum::<i64>(), sw.total_kits, sw.alloc_cases, sw.inv_violations);
    let cal_fr = if sw.dem_units_tot > 0 { sw.dem_units_filled as f64 / sw.dem_units_tot as f64 } else { 1.0 };
    let cal_bo = if sw.dem_events_tot > 0 { sw.dem_bo_events as f64 / sw.dem_events_tot as f64 } else { 0.0 };
    let cal_cu = if sw.carrier_cap_tot > 0 { sw.carrier_vol_tot as f64 / sw.carrier_cap_tot as f64 } else { 0.0 };
    let cal_ok = ((0.85..=0.99).contains(&cal_fr) && cal_bo < 0.20 && (0.30..=0.90).contains(&cal_cu)) as i32;
    let mut f = open_bw(out_dir, "s09_calibration.csv"); let _ = writeln!(f, "fill_rate,bo_event_frac,carrier_util,calibrated\n{:.6},{:.6},{:.6},{}", cal_fr, cal_bo, cal_cu, cal_ok);
}

pub fn write_s11_to_s38(s: &Sim, c: &Config, out_dir: &str) {
    for k in 0..4 {
        let mut f = open_bw(out_dir, &format!("sales_{}.csv", SKU_L[k])); let _ = writeln!(f, "day,tick,store_id,sku,qty,rev_cents");
        for r in &s.sales { if r.3 as usize == k { let _ = writeln!(f, "{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5); } }
        let mut f = open_bw(out_dir, &format!("shipments_{}.csv", SKU_L[k])); let _ = writeln!(f, "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status");
        for r in &s.ships { if r.5 as usize == k { let _ = writeln!(f, "{},{},{},{},{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5, r.6, r.7, r.8, r.9); } }
        let mut f = open_bw(out_dir, &format!("inventory_snap_{}.csv", SKU_L[k])); let _ = writeln!(f, "day,sku,on_hand_tot");
        for d in 0..(c.num_days as usize) { let _ = writeln!(f, "{},{},{}", d, SKU_L[k], s.day_inv_snap[d][k]); }
    }
    for i in 0..10 { let mut f = open_bw(out_dir, &format!("store_daily_{}.csv", i)); let _ = writeln!(f, "store,cash,oh_tot,req_tot,fill_tot\n{},{},{},{},{}", i, s.acc[i].cash, s.acc[i].on_hand.iter().sum::<i64>(), s.daily_fill_store[i][0], s.daily_fill_store[i][1]); }
    for i in 0..4 { let mut f = open_bw(out_dir, &format!("warehouse_daily_{}.csv", i)); let a = &s.acc[10 + i]; let _ = writeln!(f, "wh,cash,oh_A,oh_B,oh_C,oh_D\n{},{},{},{},{},{}", i, a.cash, a.on_hand[0], a.on_hand[1], a.on_hand[2], a.on_hand[3]); }
    for i in 0..3 { let mut f = open_bw(out_dir, &format!("factory_daily_{}.csv", i)); let _ = writeln!(f, "fac,cash,produced,state\n{},{},{},{}", i, s.acc[14 + i].cash, s.total_produced[i], s.acc[14 + i].plant_state); }
    let mut f = open_bw(out_dir, "production.csv"); let _ = writeln!(f, "day,tick,factory_id,sku,qty,cost_cents"); for r in &s.prods { let _ = writeln!(f, "{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5); }
    let mut f = open_bw(out_dir, "kit_assembly.csv"); let _ = writeln!(f, "day,tick,wh_id,qty,fee_cents,imbalance"); for r in &s.kits { let _ = writeln!(f, "{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5); }
    let mut f = open_bw(out_dir, "fees_ledger.csv"); let _ = writeln!(f, "day,tick,from_id,to_id,amount,kind"); for r in &s.fees { let _ = writeln!(f, "{},{},{},{},{},{}", r.0, r.1, r.2, r.3, r.4, r.5); }
    let mut f = open_bw(out_dir, "entities_final.csv"); let _ = writeln!(f, "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status");
    for (i, a) in s.acc.iter().enumerate() { let _ = writeln!(f, "{},{},{},{},{},{},{},{},{},{},{}", i, a.cash, a.on_hand[0], a.on_hand[1], a.on_hand[2], a.on_hand[3], a.in_transit[0], a.in_transit[1], a.in_transit[2], a.in_transit[3], a.status); }

    let _ = writeln!(open_bw(out_dir, "s11_invariants.csv"), "violations,ok\n{},{}", s.inv_violations, (s.inv_violations == 0) as i32);
    let _ = writeln!(open_bw(out_dir, "s12_clearing.csv"), "days,fee_pool_cash\n{},{}", c.num_days, s.acc[ACC_FEE].cash);
    let _ = writeln!(open_bw(out_dir, "s12_credit.csv"), "credit_holds,liquidations\n{},{}", s.branches[13], s.branches[14]);
    let _ = writeln!(open_bw(out_dir, "s13_risk.csv"), "credit_holds,liquidations\n{},{}", s.branches[13], s.branches[14]);
    let _ = writeln!(open_bw(out_dir, "s14_factories.csv"), "fac,produced\n0,{}\n1,{}\n2,{}", s.total_produced[0], s.total_produced[1], s.total_produced[2]);
    let _ = writeln!(open_bw(out_dir, "s15_carrier_fees.csv"), "tier,count\n0,{}\n1,{}\n2,{}", s.carrier_tier_cnt[0], s.carrier_tier_cnt[1], s.carrier_tier_cnt[2]);
    let _ = writeln!(open_bw(out_dir, "s16_kits.csv"), "total_kits,full,partial\n{},{},{}", s.total_kits, s.branches[8], s.branches[9]);
    let mut f = open_bw(out_dir, "s16_kit_imbalance_series.csv"); let _ = writeln!(f, "idx,imbalance"); for (i, v) in s.kit_imb_series.iter().enumerate() { let _ = writeln!(f, "{},{}", i, v); }
    let _ = writeln!(open_bw(out_dir, "s17_allocations.csv"), "cases,nonzero,total_req,total_fill\n{},{},{},{}", s.alloc_cases, s.alloc_nonzero, s.alloc_total_req, s.alloc_total_fill);
    let mut f = open_bw(out_dir, "s18_cost_entities.csv"); let _ = writeln!(f, "acc,cash_delta"); for i in 0..17 { let _ = writeln!(f, "{},{}", i, s.acc[i].cash - (if i < 10 { 2_000_000_000 } else { 5_000_000_000 }) * s.cap_scale); }
    let mut f = open_bw(out_dir, "s19_cost_types.csv"); let _ = writeln!(f, "echelon,cost_cents"); for e in 0..4 { let tot: i64 = (0..c.num_days as usize).map(|d| s.daily_cost_echelon[d][e]).sum(); let _ = writeln!(f, "{},{}", e, tot); }
    let mut f = open_bw(out_dir, "s19_cost_daily.csv"); let _ = writeln!(f, "day,prod,kit,freight,hold_pen"); for d in 0..c.num_days as usize { let _ = writeln!(f, "{},{},{},{},{}", d, s.daily_cost_echelon[d][0], s.daily_cost_echelon[d][1], s.daily_cost_echelon[d][2], s.daily_cost_echelon[d][3]); }
    for k in 0..4 { let mut f = open_bw(out_dir, &format!("s20_demand_{}.csv", SKU_L[k])); let _ = writeln!(f, "day,demand"); for d in 0..c.num_days as usize { let _ = writeln!(f, "{},{}", d, s.daily_demand[d][k]); } }
    let mut f = open_bw(out_dir, "s21_bullwhip.csv"); let _ = writeln!(f, "sku,kurtosis,acf_lag1,bullwhip_ratio");
    for k in 0..4 {
        let md = s.bin_demand[k].iter().sum::<i64>() as f64 / 100.0;
        let mo = s.bin_orders[k].iter().sum::<i64>() as f64 / 100.0;
        let (mut m2, mut m4, mut cov1, mut vo) = (0.0, 0.0, 0.0, 0.0);
        for b in 0..100 {
            let d = s.bin_demand[k][b] as f64 - md; let o = s.bin_orders[k][b] as f64 - mo;
            m2 += d * d; m4 += d * d * d * d; vo += o * o;
            if b > 0 { cov1 += (s.bin_demand[k][b] as f64 - md) * (s.bin_demand[k][b - 1] as f64 - md); }
        }
        let kurt = if m2 > 0.0 { (m4 / 100.0) / ((m2 / 100.0) * (m2 / 100.0)) } else { 0.0 };
        let acf1 = if m2 > 0.0 { cov1 / m2 } else { 0.0 };
        let bwr = if m2 > 0.0 { vo / m2 } else { 1.0 };
        let _ = writeln!(f, "{},{:.6},{:.6},{:.6}", SKU_L[k], kurt, acf1, bwr);
    }
    let mut f = open_bw(out_dir, "s22_service.csv"); let _ = writeln!(f, "store,req,fill,fill_rate"); for i in 0..10 { let fr = if s.daily_fill_store[i][0] > 0 { s.daily_fill_store[i][1] as f64 / s.daily_fill_store[i][0] as f64 } else { 1.0 }; let _ = writeln!(f, "{},{},{},{:.6}", i, s.daily_fill_store[i][0], s.daily_fill_store[i][1], fr); }
    let mut f = open_bw(out_dir, "s23_leadtime_impact.csv"); let _ = writeln!(f, "idx,leadtime,qty"); for (i, (lt, q)) in s.leadtime_pairs.iter().enumerate() { let _ = writeln!(f, "{},{},{}", i, lt, q); }
    let imb_m = if s.kit_imb_series.is_empty() { 0.0 } else { s.kit_imb_series.iter().sum::<i64>() as f64 / s.kit_imb_series.len() as f64 };
    let _ = writeln!(open_bw(out_dir, "s24_kit_coupling.csv"), "mean_imbalance,total_kits\n{:.4},{}", imb_m, s.total_kits);
    let mut f = open_bw(out_dir, "s25_shock_event_study.csv"); let _ = writeln!(f, "cls,count,mean_rel_dev,se_rel_dev");
    for cl in 0..3 {
        let (mut cnt, mut sum, mut sum2) = (0usize, 0.0, 0.0);
        for &(b, c_ev, st) in &s.shock_evs {
            if c_ev != cl { continue; }
            let (mut pre, mut post, mut pc, mut qc) = (0.0, 0.0, 0, 0);
            for db in -2i32..=-1 { let bi = b as i32 + db; if bi >= 0 { pre += s.bin_store_dem[st][bi as usize] as f64; pc += 1; } }
            for db in 0i32..=1 { let bi = b as i32 + db; if bi < 100 { post += s.bin_store_dem[st][bi as usize] as f64; qc += 1; } }
            if pc > 0 && qc > 0 && pre > 0.0 { let rdev = (post / qc as f64 - pre / pc as f64) / (pre / pc as f64); sum += rdev; sum2 += rdev * rdev; cnt += 1; }
        }
        let mn = if cnt > 0 { sum / cnt as f64 } else { 0.0 };
        let se = if cnt > 1 { ((sum2 / cnt as f64 - mn * mn) / (cnt - 1) as f64).max(0.0).sqrt() } else { 0.0 };
        let _ = writeln!(f, "{},{},{:.6},{:.6}", cl, cnt, mn, se);
    }
    let _ = writeln!(open_bw(out_dir, "s26_disruptions.csv"), "disrupt_hits\n{}", s.branches[11]);
    let _ = writeln!(open_bw(out_dir, "s27_fees.csv"), "n_fee_entries,fee_pool_cash\n{},{}", s.fees.len(), s.acc[ACC_FEE].cash);
    let _ = writeln!(open_bw(out_dir, "s28_capacity_rejects.csv"), "wh_cap_reject,carrier_cap_reject\n{},{}", s.branches[4], s.branches[5]);

    let mut f = open_bw(out_dir, "s29_sweep_safety_stock.csv"); let _ = writeln!(f, "S_lvl,rep,fill_rate,avg_oh");
    for sv in [80, 220, 350] { for rp in 0..5 { let mut c2 = c.clone(); c2.master_seed += rp + 1; let mut ss = Sim::new(&c2, false, None); sim_run_days(&mut ss, 2, 2000, sv, -1, 1); let oh: i64 = (0..10).map(|i| ss.acc[i].on_hand[0]).sum(); let fr = if ss.dem_units_tot > 0 { ss.dem_units_filled as f64 / ss.dem_units_tot as f64 } else { 1.0 }; let _ = writeln!(f, "{},{},{:.6},{:.2}", sv, rp, fr, oh as f64 / 10.0); } }
    let mut f = open_bw(out_dir, "s30_sweep_shock_disrupt.csv"); let _ = writeln!(f, "disrupt_pct,rep,produced,sold");
    for dp in [0, 4, 12] { for rp in 0..5 { let mut c2 = c.clone(); c2.master_seed += rp + 10; let mut ss = Sim::new(&c2, false, None); sim_run_days(&mut ss, 2, 2000, 0, dp, 1); let _ = writeln!(f, "{},{},{},{}", dp, rp, ss.total_produced.iter().sum::<i64>(), ss.total_sold.iter().sum::<i64>()); } }
    let mut f = open_bw(out_dir, "s31_sweep_transship.csv"); let _ = writeln!(f, "transship_on,rep,mean_imbalance,n_transships");
    for ts in [0, 1] { for rp in 0..5 { let mut c2 = c.clone(); c2.master_seed += rp + 20; let mut ss = Sim::new(&c2, false, None); sim_run_days(&mut ss, 2, 2000, 0, -1, ts); let im = if ss.kit_imb_series.is_empty() { 0.0 } else { ss.kit_imb_series.iter().sum::<i64>() as f64 / ss.kit_imb_series.len() as f64 }; let _ = writeln!(f, "{},{},{:.4},{}", ts, rp, im, ss.n_transships); } }

    let bn = ["DEMAND_FILL","DEMAND_BO","DEMAND_LOST","STORE_REPLEN","WH_CAP_REJECT","CARRIER_CAP_REJECT","ALLOC_FULL","ALLOC_PARTIAL","KIT_FULL","KIT_PARTIAL","TRANSSHIP","DISRUPT_HIT","SHOCK_HIT","CREDIT_HOLD","LIQUIDATE","CARRIER_DELAY"];
    let mut f = open_bw(out_dir, "s32_branches.csv"); let _ = writeln!(f, "branch,count"); for (i, nm) in bn.iter().enumerate() { let _ = writeln!(f, "{},{}", nm, s.branches[i]); }

    let mut sp1 = Sim::new_scaled(c, false, None, c.parity_events_per_day); sim_run_days(&mut sp1, c.parity_days, c.parity_events_per_day, 0, -1, 1);
    let mut sp2 = Sim::new_scaled(c, false, None, c.parity_events_per_day); sim_run_days(&mut sp2, c.parity_days, c.parity_events_per_day, 0, -1, 1);
    let (mut h1, mut h2) = (0xCBF29CE484222325u64, 0xCBF29CE484222325u64);
    for i in 0..NUM_ACCOUNTS { h1 = fnv1a_i64(h1, sp1.acc[i].cash); h2 = fnv1a_i64(h2, sp2.acc[i].cash); for k in 0..4 { h1 = fnv1a_i64(h1, sp1.acc[i].on_hand[k]); h2 = fnv1a_i64(h2, sp2.acc[i].on_hand[k]); } }
    let _ = writeln!(open_bw(out_dir, "s33_replay.csv"), "ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n{:016x},{:016x},{:016x},{:016x},{}", sp1.event_hash, sp2.event_hash, h1, h2, (sp1.event_hash == sp2.event_hash && h1 == h2) as i32);
    let _ = writeln!(open_bw(out_dir, "s34_parity.csv"), "parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n{},{},{:016x},{:016x},{},{},{}", c.parity_days, c.parity_events_per_day, sp1.event_hash, h1, sp1.total_sold.iter().sum::<i64>(), sp1.total_kits, sp1.inv_violations);

    let mut vals = Vec::new();
    for _ in 0..5 { let mut sb = Sim::new(c, false, None); let t0 = Instant::now(); sim_run_days(&mut sb, 1, 8000, 0, -1, 1); vals.push(t0.elapsed().as_secs_f64() * 1e9 / 8000.0); }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap()); let med = vals[2];
    let mut dev: Vec<f64> = vals.iter().map(|v| (v - med).abs()).collect(); dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let _ = writeln!(open_bw(out_dir, "s35_bench.csv"), "events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n8000,{:.4},{:.4},{:.2}", med, dev[2], 1e9 / med);
    let _ = writeln!(open_bw(out_dir, "s36_sizes.csv"), "n_sales,n_ships,n_prods,n_kits,n_fees\n{},{},{},{},{}", s.sales.len(), s.ships.len(), s.prods.len(), s.kits.len(), s.fees.len());

    let nq = s.sales.len().min(4096); let raw: Vec<i32> = s.sales.iter().take(nq).map(|r| r.4).collect();
    let mut rle = Vec::new(); let mut i = 0;
    while i < nq { let dv = if i == 0 { raw[0] } else { raw[i] - raw[i - 1] }; let mut run = 1; while i + run < nq && (raw[i + run] - raw[i + run - 1]) == dv { run += 1; } rle.push(dv); rle.push(run as i32); i += run; }
    let mut dec = Vec::with_capacity(nq);
    for ch in rle.chunks(2) { for _ in 0..ch[1] { let v = if dec.is_empty() { ch[0] } else { *dec.last().unwrap() + ch[0] }; dec.push(v); } }
    let _ = writeln!(open_bw(out_dir, "s37_compress.csv"), "n_items,rle_ints,ratio,roundtrip_ok\n{},{},{:.4},{}", nq, rle.len(), if rle.is_empty() { 1.0 } else { nq as f64 / rle.len() as f64 }, (dec == raw) as i32);

    let si = Sim::new_scaled(c, false, None, s.cap_scale * 2000);
    let mut l_cash: Vec<i64> = si.acc.iter().map(|a| a.cash).collect();
    let mut l_oh: Vec<[i64; 4]> = si.acc.iter().map(|a| a.on_hand).collect();
    let mut l_it = vec![[0i64; 4]; NUM_ACCOUNTS];
    for r in &s.sales { l_cash[ACC_CUST] -= r.5; l_cash[r.2 as usize] += r.5; l_oh[r.2 as usize][r.3 as usize] -= r.4 as i64; }
    for r in &s.ships {
        let (src, dst, cid, sku, qty, goods, freight, st) = (r.2 as usize, r.3 as usize, r.4 as usize, r.5 as usize, r.6 as i64, r.7, r.8, r.9);
        if st == 1 { l_oh[src][sku] -= qty; l_it[dst][sku] += qty; l_cash[dst] -= goods + freight; l_cash[src] += goods; l_cash[17 + cid] += freight; }
        else { l_it[dst][sku] -= qty; l_oh[dst][sku] += qty; }
    }
    for r in &s.prods { l_cash[r.2 as usize] -= r.5; l_cash[ACC_SUPP] += r.5; l_oh[r.2 as usize][r.3 as usize] += r.4 as i64; }
    for r in &s.kits  { l_cash[r.2 as usize] -= r.4; l_cash[ACC_FEE] += r.4; for k in 0..3 { l_oh[r.2 as usize][k] -= r.3 as i64; } l_oh[r.2 as usize][3] += r.3 as i64; }
    for r in &s.fees  { l_cash[r.2 as usize] -= r.4; l_cash[r.3 as usize] += r.4; }
    let mut max_cd = 0i64; let mut max_id = 0i64;
    for i in 0..NUM_ACCOUNTS { max_cd = max_cd.max((l_cash[i] - s.acc[i].cash).abs()); for k in 0..4 { max_id = max_id.max((l_oh[i][k] - s.acc[i].on_hand[k]).abs() + (l_it[i][k] - s.acc[i].in_transit[k]).abs()); } }
    let _ = writeln!(open_bw(out_dir, "s38_ledger_replay.csv"), "max_cash_diff,max_inv_diff,match\n{},{},{}", max_cd, max_id, (max_cd == 0 && max_id == 0) as i32);
}
