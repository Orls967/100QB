use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use crate::inventory::*;
use crate::rng::*;

#[derive(Clone, Default)]
pub struct Config {
    pub master_seed: i64, pub num_days: i64, pub parity_days: i64, pub parity_events_per_day: i64,
    pub ticks_per_day: i64, pub ewma_alpha_milli: i64, pub store_s_base: i64, pub store_S_base: i64,
    pub wh_s_base: i64, pub wh_S_base: i64, pub wh_cap: i64, pub carrier_cap_daily: i64,
    pub kit_fee: i64, pub holding_cost_per_day: i64, pub backorder_penalty: i64,
    pub disrupt_prob_pct: i64, pub shock_prob_pct: i64, pub credit_margin_bps: i64, pub k_chk: i64,
}

#[derive(Clone, Default)]
pub struct Account {
    pub cash: i64, pub on_hand: [i64; 4], pub in_transit: [i64; 4], pub backorder: [i64; 4],
    pub credit_limit: i64, pub status: i32, pub ewma_milli: [i32; 4], pub plant_state: i32, pub disrupt_rem: i32,
}

pub struct Sim {
    pub cfg: Config,
    pub acc: Vec<Account>,
    pub rng: Vec<Xoshiro256>,
    pub heap: BinaryHeap<SimEvent>,
    pub next_seq: u64,
    pub init_cash: i64,
    pub init_stock: [i64; 4],
    pub total_produced: [i64; 4],
    pub total_sold: [i64; 4],
    pub total_kits: i64,
    pub branches: [i64; 16],
    pub carrier_day_vol: [i64; 3],
    pub carrier_tier_cnt: [i64; 3],
    pub alloc_total_req: i64, pub alloc_total_fill: i64, pub alloc_cases: i64, pub alloc_nonzero: i64,
    pub inv_violations: i64,
    pub log_enabled: bool,
    pub daily_demand: [[i64; 4]; 6],
    pub daily_cost_echelon: [[i64; 4]; 6],
    pub daily_fill_store: [[i64; 2]; 10],
    pub shock_cnt: [i64; 3], pub shock_resp: [i64; 3],
    pub day_inv_snap: [[i64; 4]; 6],
    pub kit_imb_series: Vec<i64>,
    pub leadtime_pairs: Vec<(i64, i64)>,
    pub sales: Vec<(i64, i64, i32, i32, i32, i64)>,
    pub ships: Vec<(i64, i64, i32, i32, i32, i32, i32, i64, i64, i32)>,
    pub prods: Vec<(i64, i64, i32, i32, i32, i64)>,
    pub kits: Vec<(i64, i64, i32, i32, i64, i32)>,
    pub fees: Vec<(i64, i64, i32, i32, i64, i32)>,
    pub f_orders: Option<BufWriter<File>>, pub f_alloc: Option<BufWriter<File>>,
    pub f_disrupt: Option<BufWriter<File>>, pub f_credit: Option<BufWriter<File>>,
    pub f_risk: Option<BufWriter<File>>, pub f_trans: Option<BufWriter<File>>,
    pub f_shocks: Option<BufWriter<File>>, pub f_acc_ser: Option<BufWriter<File>>,
    pub event_hash: u64,
}

pub fn load_config(path: &str) -> Config {
    let mut c = Config::default();
    let txt = fs::read_to_string(path).expect("config.kv");
    for line in txt.lines() {
        if let Some((k, v)) = line.split_once('=') {
            let val: i64 = v.trim().parse().unwrap_or(0);
            match k.trim() {
                "master_seed" => c.master_seed = val, "num_days" => c.num_days = val,
                "parity_days" => c.parity_days = val, "parity_events_per_day" => c.parity_events_per_day = val,
                "ticks_per_day" => c.ticks_per_day = val, "ewma_alpha_milli" => c.ewma_alpha_milli = val,
                "store_s_base" => c.store_s_base = val, "store_S_base" => c.store_S_base = val,
                "wh_s_base" => c.wh_s_base = val, "wh_S_base" => c.wh_S_base = val,
                "wh_cap" => c.wh_cap = val, "carrier_cap_daily" => c.carrier_cap_daily = val,
                "kit_fee" => c.kit_fee = val, "holding_cost_per_day" => c.holding_cost_per_day = val,
                "backorder_penalty" => c.backorder_penalty = val, "disrupt_prob_pct" => c.disrupt_prob_pct = val,
                "shock_prob_pct" => c.shock_prob_pct = val, "credit_margin_bps" => c.credit_margin_bps = val,
                "k_chk" => c.k_chk = val, _ => {}
            }
        }
    }
    c
}

fn open_bw(dir: &str, name: &str, hdr: &str) -> BufWriter<File> {
    let mut w = BufWriter::new(File::create(format!("{}/{}", dir, name)).unwrap());
    let _ = writeln!(w, "{}", hdr);
    w
}

impl Sim {
    pub fn new(c: &Config, log_en: bool, out_dir: Option<&str>) -> Self {
        let mut acc = vec![Account::default(); NUM_ACCOUNTS];
        let mut rng = Vec::with_capacity(NUM_ACCOUNTS + 2);
        for i in 0..(NUM_ACCOUNTS + 2) {
            let mut x = (c.master_seed as u64) ^ fnv1a_u32pair((i / 10) as u32, i as u32);
            rng.push(Xoshiro256::seed(splitmix64(&mut x)));
        }
        for i in 0..NUM_STORES {
            acc[i].cash = 2_000_000_000; acc[i].credit_limit = 500_000_000;
            acc[i].on_hand = [400, 400, 400, 120]; acc[i].ewma_milli = [12000; 4];
        }
        for i in 10..14 { acc[i].cash = 5_000_000_000; acc[i].credit_limit = 1_000_000_000; acc[i].on_hand = [2500, 2500, 2500, 800]; }
        for i in 14..17 { acc[i].cash = 5_000_000_000; acc[i].credit_limit = 2_000_000_000; acc[i].on_hand[i - 14] = 4000; }
        for i in 17..20 { acc[i].cash = 500_000_000; acc[i].credit_limit = 200_000_000; }
        acc[ACC_CUST].cash = 50_000_000_000; acc[ACC_SUPP].cash = 50_000_000_000; acc[ACC_FEE].cash = 1_000_000_000;
        let mut init_cash = 0i64; let mut init_stock = [0i64; 4];
        for a in &acc { init_cash += a.cash; for k in 0..4 { init_stock[k] += a.on_hand[k]; } }
        let od = out_dir.unwrap_or("");
        Self {
            cfg: c.clone(), acc, rng, heap: BinaryHeap::new(), next_seq: 1, init_cash, init_stock,
            total_produced: [0; 4], total_sold: [0; 4], total_kits: 0, branches: [0; 16],
            carrier_day_vol: [0; 3], carrier_tier_cnt: [0; 3], alloc_total_req: 0, alloc_total_fill: 0,
            alloc_cases: 0, alloc_nonzero: 0, inv_violations: 0, log_enabled: log_en,
            daily_demand: [[0; 4]; 6], daily_cost_echelon: [[0; 4]; 6], daily_fill_store: [[0; 2]; 10],
            shock_cnt: [0; 3], shock_resp: [0; 3], day_inv_snap: [[0; 4]; 6],
            kit_imb_series: Vec::new(), leadtime_pairs: Vec::new(),
            sales: Vec::new(), ships: Vec::new(), prods: Vec::new(), kits: Vec::new(), fees: Vec::new(),
            f_orders: if log_en { Some(open_bw(od, "orders_log.csv", "day,tick,seq,etype,src,dst,sku,qty")) } else { None },
            f_alloc: if log_en { Some(open_bw(od, "allocations.csv", "day,wh,sku,store,req,fill")) } else { None },
            f_disrupt: if log_en { Some(open_bw(od, "disruptions.csv", "day,tick,fac,state,dur")) } else { None },
            f_credit: if log_en { Some(open_bw(od, "credit_calls.csv", "day,acc,cash,req_margin,status")) } else { None },
            f_risk: if log_en { Some(open_bw(od, "risk_actions.csv", "day,tick,acc,action,amount")) } else { None },
            f_trans: if log_en { Some(open_bw(od, "transshipments.csv", "day,tick,src,dst,sku,qty")) } else { None },
            f_shocks: if log_en { Some(open_bw(od, "shocks_log.csv", "day,tick,cls,store,mag")) } else { None },
            f_acc_ser: if log_en { Some(open_bw(od, "accounts_series.csv", "day,acc,cash,oh_tot,it_tot,status")) } else { None },
            event_hash: 0xCBF29CE484222325,
        }
    }
}
