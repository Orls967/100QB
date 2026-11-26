use std::io::Write;
use crate::entities::Sim;
use crate::inventory::*;

pub fn sim_risk_check(s: &mut Sim, day: i64, tick: i64) {
    for i in 0..(NUM_STORES + NUM_WAREHOUSES) {
        let mut bo_val = 0i64;
        for k in 0..NUM_SKUS { bo_val += s.acc[i].backorder[k] * SKU_WHOLESALE[k]; }
        let req_margin = (bo_val * s.cfg.credit_margin_bps) / 10000;
        if s.acc[i].cash < req_margin && s.acc[i].status == 0 {
            s.acc[i].status = 1; s.branches[13] += 1;
            if let Some(ref mut fc) = s.f_credit { let _ = writeln!(fc, "{},{},{},{},1", day, i, s.acc[i].cash, req_margin); }
        } else if s.acc[i].status == 1 && s.acc[i].cash >= req_margin { s.acc[i].status = 0; }
        if s.acc[i].cash < 0 && s.acc[i].status != 2 {
            s.acc[i].status = 2; s.branches[14] += 1;
            transfer_cash(s, i, ACC_FEE, 50000, 2, day, tick);
            if let Some(ref mut fr) = s.f_risk { let _ = writeln!(fr, "{},{},{},LIQUIDATE,50000", day, tick, i); }
        }
    }
}

pub fn sim_eod_clearing(s: &mut Sim, day: i64) {
    for i in 0..(NUM_STORES + NUM_WAREHOUSES) {
        let oh_tot: i64 = s.acc[i].on_hand.iter().sum();
        let it_tot: i64 = s.acc[i].in_transit.iter().sum();
        let bo_tot: i64 = s.acc[i].backorder.iter().sum();
        if (day as usize) < 6 { for k in 0..4 { s.day_inv_snap[day as usize][k] += s.acc[i].on_hand[k]; } }
        let hold = oh_tot * s.cfg.holding_cost_per_day; let pen = bo_tot * s.cfg.backorder_penalty;
        if hold > 0 { transfer_cash(s, i, ACC_FEE, hold, 0, day, s.cfg.ticks_per_day); if (day as usize) < 6 { s.daily_cost_echelon[day as usize][3] += hold; } }
        if pen > 0  { transfer_cash(s, i, ACC_FEE, pen, 1, day, s.cfg.ticks_per_day);  if (day as usize) < 6 { s.daily_cost_echelon[day as usize][3] += pen; } }
        if let Some(ref mut fa) = s.f_acc_ser { let _ = writeln!(fa, "{},{},{},{},{},{}", day, i, s.acc[i].cash, oh_tot, it_tot, s.acc[i].status); }
    }
    sim_risk_check(s, day, s.cfg.ticks_per_day);
}
