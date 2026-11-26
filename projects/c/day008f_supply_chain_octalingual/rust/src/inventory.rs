use crate::entities::Sim;

pub const NUM_SKUS: usize = 4;
pub const NUM_FACTORIES: usize = 3;
pub const NUM_WAREHOUSES: usize = 4;
pub const NUM_STORES: usize = 10;
pub const NUM_CARRIERS: usize = 3;
pub const NUM_ACCOUNTS: usize = 23;
pub const ACC_CUST: usize = 20;
pub const ACC_SUPP: usize = 21;
pub const ACC_FEE: usize = 22;

pub const SKU_RETAIL: [i64; 4] = [10000, 10200, 9800, 31500];
pub const SKU_WHOLESALE: [i64; 4] = [8000, 8200, 7800, 26000];
pub const SKU_PROD: [i64; 4] = [6000, 6200, 5800, 0];
pub const STORE_WH: [usize; 10] = [10, 10, 10, 11, 11, 11, 12, 12, 12, 12];
pub const CARRIER_BASE_LT: [i64; 3] = [120, 200, 160];

pub fn carrier_rate(day_vol: i64) -> i64 { if day_vol < 500 { 200 } else if day_vol < 2000 { 160 } else { 120 } }

pub fn transfer_cash(s: &mut Sim, from_id: usize, to_id: usize, amt: i64, kind: i32, day: i64, tick: i64) {
    if amt <= 0 || from_id == to_id { return; }
    s.acc[from_id].cash -= amt;
    s.acc[to_id].cash += amt;
    if kind >= 0 { s.fees.push((day, tick, from_id as i32, to_id as i32, amt, kind)); }
}

pub fn check_invariants(s: &mut Sim) -> i32 {
    let tot_cash: i64 = s.acc.iter().map(|a| a.cash).sum();
    if tot_cash != s.init_cash { s.inv_violations += 1; return 0; }
    for k in 0..NUM_SKUS {
        let mut oh = 0i64; let mut it = 0i64;
        for a in &s.acc {
            if a.on_hand[k] < 0 || a.in_transit[k] < 0 || a.backorder[k] < 0 { s.inv_violations += 1; return 0; }
            oh += a.on_hand[k]; it += a.in_transit[k];
        }
        let lhs = oh + it + s.total_sold[k] + if k < 3 { s.total_kits } else { 0 };
        let rhs = s.init_stock[k] + if k < 3 { s.total_produced[k] } else { s.total_kits };
        if lhs != rhs { s.inv_violations += 1; return 0; }
    }
    1
}
