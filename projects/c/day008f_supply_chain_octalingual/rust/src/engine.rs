use std::io::Write;
use crate::allocation::run_alloc_auction;
use crate::clearing::{sim_eod_clearing, sim_risk_check};
use crate::entities::Sim;
use crate::inventory::*;
use crate::rng::*;

fn do_ship(s: &mut Sim, src: usize, dst: usize, sku: usize, mut qty: i64, day: i64, tick: i64, t_abs: i64) {
    if qty <= 0 || s.acc[src].on_hand[sku] <= 0 { return; }
    if qty > s.acc[src].on_hand[sku] { qty = s.acc[src].on_hand[sku]; }
    if (10..14).contains(&dst) && s.acc[dst].on_hand.iter().sum::<i64>() + qty > s.cfg.wh_cap { s.branches[4] += 1; return; }
    let mut cid = (src + dst + sku) % NUM_CARRIERS;
    if s.carrier_day_vol[cid] + qty > s.cfg.carrier_cap_daily {
        cid = (cid + 1) % NUM_CARRIERS;
        if s.carrier_day_vol[cid] + qty > s.cfg.carrier_cap_daily { s.branches[5] += 1; return; }
    }
    let rate = carrier_rate(s.carrier_day_vol[cid]);
    s.carrier_tier_cnt[if rate == 200 { 0 } else if rate == 160 { 1 } else { 2 }] += 1;
    s.carrier_day_vol[cid] += qty;
    let goods = qty * SKU_WHOLESALE[sku]; let freight = qty * rate;
    s.acc[src].on_hand[sku] -= qty; s.acc[dst].in_transit[sku] += qty;
    transfer_cash(s, dst, src, goods, -1, day, tick);
    transfer_cash(s, dst, 17 + cid, freight, -1, day, tick);
    if (day as usize) < 6 { s.daily_cost_echelon[day as usize][2] += freight; }
    let delay = s.rng[17 + cid].rand_int(0, 80);
    if delay > 50 { s.branches[15] += 1; }
    let lt = CARRIER_BASE_LT[cid] + delay;
    if s.leadtime_pairs.len() < 256 { s.leadtime_pairs.push((lt, qty)); }
    if s.log_enabled { s.ships.push((day, tick, src as i32, dst as i32, cid as i32, sku as i32, qty as i32, goods, freight, 1)); }
    let seq = s.next_seq; s.next_seq += 1;
    s.heap.push(SimEvent { time: t_abs + lt, seq, etype: 2, src: src as i32, dst: dst as i32, sku: sku as i32, qty: qty as i32, aux: cid as i32 });
}

fn hash_ev(s: &mut Sim, ev: &SimEvent) {
    s.event_hash = fnv1a_i64(s.event_hash, ev.time); s.event_hash = fnv1a_i64(s.event_hash, ev.etype as i64);
    s.event_hash = fnv1a_i64(s.event_hash, ev.src as i64); s.event_hash = fnv1a_i64(s.event_hash, ev.dst as i64);
    s.event_hash = fnv1a_i64(s.event_hash, ev.sku as i64); s.event_hash = fnv1a_i64(s.event_hash, ev.qty as i64);
}

pub fn sim_run_days(s: &mut Sim, days: i64, ev_per_day: i64, safety_override: i64, disrupt_override: i64, transship_on: i32) {
    let t_day = s.cfg.ticks_per_day; let mut ev_processed = 0i64;
    for d in 0..days {
        let day_base = d * t_day; s.carrier_day_vol = [0; 3];
        for f in 14..17 {
            if s.acc[f].plant_state == 1 { s.acc[f].disrupt_rem -= 1; if s.acc[f].disrupt_rem <= 0 { s.acc[f].plant_state = 2; } }
            else if s.acc[f].plant_state == 2 { s.acc[f].plant_state = 0; }
        }
        for wh in 10..13 {
            for sku in 0..NUM_SKUS {
                let mut sts = Vec::new(); let mut reqs = Vec::new(); let mut seqs = Vec::new(); let mut prio = Vec::new();
                let s_lvl = if safety_override > 0 { safety_override } else { s.cfg.store_S_base };
                for st in 0..NUM_STORES {
                    if STORE_WH[st] != wh { continue; }
                    let inv = s.acc[st].on_hand[sku] + s.acc[st].in_transit[sku];
                    if inv < s_lvl { sts.push(st); reqs.push(s_lvl - inv); seqs.push(s.next_seq as i64); s.next_seq += 1; prio.push(if st % 3 == 0 { 2 } else { 1 }); }
                }
                if !sts.is_empty() {
                    let fills = run_alloc_auction(s, wh, sku, &sts, &reqs, &seqs, &prio);
                    for i in 0..sts.len() {
                        if let Some(ref mut fw) = s.f_alloc { let _ = writeln!(fw, "{},{},{},{},{},{}", d, wh, sku, sts[i], reqs[i], fills[i]); }
                        if fills[i] > 0 { do_ship(s, wh, sts[i], sku, fills[i], d, 0, day_base); }
                    }
                }
            }
        }
        for e in 0..ev_per_day {
            let t_abs = day_base + (e * t_day) / ev_per_day;
            let st = s.rng[NUM_ACCOUNTS].rand_int(0, 9) as i32;
            let sku = s.rng[NUM_ACCOUNTS].rand_int(0, 3) as i32;
            let qty = s.rng[NUM_ACCOUNTS].rand_int(1, 14) as i32;
            s.heap.push(SimEvent { time: t_abs, seq: s.next_seq, etype: 0, src: ACC_CUST as i32, dst: st, sku, qty, aux: 0 }); s.next_seq += 1;
            if e % 20 == 0 { let f = ((e / 20) % 3) as i32; s.heap.push(SimEvent { time: t_abs + 1, seq: s.next_seq, etype: 3, src: 14 + f, dst: 10 + f, sku: f, qty: 80, aux: 0 }); s.next_seq += 1; }
            if e % 40 == 0 { s.heap.push(SimEvent { time: t_abs + 2, seq: s.next_seq, etype: 4, src: 13, dst: 13, sku: 3, qty: 25, aux: 0 }); s.next_seq += 1; }
            if e % 60 == 0 && transship_on != 0 { s.heap.push(SimEvent { time: t_abs + 3, seq: s.next_seq, etype: 5, src: 13, dst: 10 + ((e / 60) % 3) as i32, sku: ((e / 60) % 4) as i32, qty: 40, aux: 0 }); s.next_seq += 1; }
            if e % 250 == 0 {
                let dp = if disrupt_override >= 0 { disrupt_override } else { s.cfg.disrupt_prob_pct };
                if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= dp { let f = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2) as i32; s.heap.push(SimEvent { time: t_abs + 4, seq: s.next_seq, etype: 7, src: 14 + f, dst: 14 + f, sku: f, qty: 1, aux: 0 }); s.next_seq += 1; }
                if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= s.cfg.shock_prob_pct { let st_s = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 9) as i32; let cls = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2) as i32; s.heap.push(SimEvent { time: t_abs + 5, seq: s.next_seq, etype: 6, src: ACC_CUST as i32, dst: st_s, sku: 0, qty: (cls + 1) * 10, aux: cls }); s.next_seq += 1; }
                s.heap.push(SimEvent { time: t_abs + 6, seq: s.next_seq, etype: 8, src: 0, dst: 0, sku: 0, qty: 0, aux: 0 }); s.next_seq += 1;
            }
        }
        let day_end = (d + 1) * t_day;
        while let Some(top) = s.heap.peek() {
            if top.time >= day_end { break; }
            let ev = s.heap.pop().unwrap(); let tick = ev.time - day_base; hash_ev(s, &ev); ev_processed += 1;
            if ev_processed <= 5000 { if let Some(ref mut fo) = s.f_orders { let _ = writeln!(fo, "{},{},{},{},{},{},{},{}", d, tick, ev.seq, ev.etype, ev.src, ev.dst, ev.sku, ev.qty); } }
            match ev.etype {
                0 => {
                    let (st, sku, qty) = (ev.dst as usize, ev.sku as usize, ev.qty as i64);
                    if (d as usize) < 6 { s.daily_demand[d as usize][sku] += qty; s.daily_fill_store[st][0] += qty; }
                    s.acc[st].ewma_milli[sku] = ((s.cfg.ewma_alpha_milli * qty * 1000 + (1000 - s.cfg.ewma_alpha_milli) * (s.acc[st].ewma_milli[sku] as i64)) / 1000) as i32;
                    if s.acc[st].status == 3 { s.branches[2] += 1; }
                    else if s.acc[st].on_hand[sku] >= qty {
                        s.acc[st].on_hand[sku] -= qty; s.total_sold[sku] += qty; let rev = qty * SKU_RETAIL[sku];
                        transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                        if s.log_enabled { s.sales.push((d, tick, st as i32, sku as i32, qty as i32, rev)); }
                        s.branches[0] += 1; if (d as usize) < 6 { s.daily_fill_store[st][1] += qty; }
                    } else {
                        let part = s.acc[st].on_hand[sku];
                        if part > 0 {
                            s.acc[st].on_hand[sku] = 0; s.total_sold[sku] += part; let rev = part * SKU_RETAIL[sku];
                            transfer_cash(s, ACC_CUST, st, rev, -1, d, tick);
                            if s.log_enabled { s.sales.push((d, tick, st as i32, sku as i32, part as i32, rev)); }
                            if (d as usize) < 6 { s.daily_fill_store[st][1] += part; }
                        }
                        s.acc[st].backorder[sku] += qty - part; s.branches[1] += 1;
                    }
                    let s_lvl = if safety_override > 0 { (safety_override * s.cfg.store_s_base) / s.cfg.store_S_base } else { s.cfg.store_s_base };
                    let cap_s = if safety_override > 0 { safety_override } else { s.cfg.store_S_base };
                    let inv_pos = s.acc[st].on_hand[sku] + s.acc[st].in_transit[sku] - s.acc[st].backorder[sku];
                    if inv_pos < s_lvl && s.acc[st].status == 0 { s.branches[3] += 1; do_ship(s, STORE_WH[st], st, sku, cap_s - inv_pos, d, tick, ev.time); }
                }
                2 => {
                    let (dst, sku, qty) = (ev.dst as usize, ev.sku as usize, ev.qty as i64);
                    if s.acc[dst].in_transit[sku] >= qty { s.acc[dst].in_transit[sku] -= qty; s.acc[dst].on_hand[sku] += qty; }
                    if s.log_enabled { s.ships.push((d, tick, ev.src, ev.dst, ev.aux, ev.sku, ev.qty, 0, 0, 2)); }
                }
                3 => {
                    let (fac, sku) = (ev.src as usize, ev.sku as usize);
                    let qty = if s.acc[fac].plant_state == 1 { 0 } else if s.acc[fac].plant_state == 2 { ev.qty as i64 / 2 } else { ev.qty as i64 };
                    if qty > 0 {
                        let cost = qty * SKU_PROD[sku]; s.acc[fac].on_hand[sku] += qty; s.total_produced[sku] += qty;
                        transfer_cash(s, fac, ACC_SUPP, cost, -1, d, tick); if (d as usize) < 6 { s.daily_cost_echelon[d as usize][0] += cost; }
                        if s.log_enabled { s.prods.push((d, tick, ev.src, ev.sku, qty as i32, cost)); }
                        do_ship(s, fac, ev.dst as usize, sku, qty, d, tick, ev.time);
                    }
                }
                4 => {
                    let m = s.acc[13].on_hand[0].min(s.acc[13].on_hand[1].min(s.acc[13].on_hand[2]));
                    let mx = s.acc[13].on_hand[0].max(s.acc[13].on_hand[1].max(s.acc[13].on_hand[2]));
                    let imb = mx - m; if s.kit_imb_series.len() < 64 { s.kit_imb_series.push(imb); }
                    let mut q = ev.qty as i64;
                    if m >= q { s.branches[8] += 1; } else { q = if m > 0 { m } else { 0 }; s.branches[9] += 1; }
                    if q > 0 {
                        for k in 0..3 { s.acc[13].on_hand[k] -= q; }
                        s.acc[13].on_hand[3] += q; s.total_kits += q; let fee = q * s.cfg.kit_fee;
                        transfer_cash(s, 13, ACC_FEE, fee, -1, d, tick); if (d as usize) < 6 { s.daily_cost_echelon[d as usize][1] += fee; }
                        if s.log_enabled { s.kits.push((d, tick, 13, q as i32, fee, imb as i32)); }
                    }
                }
                5 => { s.branches[10] += 1; do_ship(s, ev.src as usize, ev.dst as usize, ev.sku as usize, ev.qty as i64, d, tick, ev.time); if let Some(ref mut ft) = s.f_trans { let _ = writeln!(ft, "{},{},{},{},{},{}", d, tick, ev.src, ev.dst, ev.sku, ev.qty); } }
                7 => { s.acc[ev.src as usize].plant_state = 1; s.acc[ev.src as usize].disrupt_rem = 1; s.branches[11] += 1; if let Some(ref mut fd) = s.f_disrupt { let _ = writeln!(fd, "{},{},{},1,1", d, tick, ev.src); } }
                6 => {
                    s.branches[12] += 1; s.shock_cnt[ev.aux as usize] += 1; s.shock_resp[ev.aux as usize] += ev.qty as i64;
                    if let Some(ref mut fs) = s.f_shocks { let _ = writeln!(fs, "{},{},{},{},{}", d, tick, ev.aux, ev.dst, ev.qty); }
                    s.heap.push(SimEvent { time: ev.time + 1, seq: s.next_seq, etype: 0, src: ACC_CUST as i32, dst: ev.dst, sku: 0, qty: ev.qty, aux: 0 }); s.next_seq += 1;
                }
                8 => { sim_risk_check(s, d, tick); }
                _ => {}
            }
            if ev_processed % s.cfg.k_chk == 0 { check_invariants(s); }
        }
        sim_eod_clearing(s, d); check_invariants(s);
    }
}
