use std::fs;

pub const SLOTS_PER_DAY: i64 = 8640;
pub const SLOTS_PER_HOUR: i64 = 360;
pub const N_STN: usize = 6;
pub const MAX_SRV: usize = 16;
pub const MAX_Q: usize = 4096;
pub const MAX_POOL: usize = 8192;
pub const AUDIT_DAYS: i64 = 3;

pub const STN_CHECKIN: usize = 0; pub const STN_SECURITY: usize = 1; pub const STN_IMIGRASI: usize = 2;
pub const STN_GATE: usize = 3; pub const STN_BELT: usize = 4; pub const STN_CUSTOMS: usize = 5;
pub const CLS_ECON: i32 = 0; pub const CLS_BIZ: i32 = 1; pub const CLS_TRANSIT: i32 = 2;
pub const DEST_DOM: i32 = 0; pub const DEST_INTL: i32 = 1;
pub const SRV_IDLE: i32 = 0; pub const SRV_BUSY: i32 = 1; pub const SRV_DOWN: i32 = 2;

#[derive(Clone, Default)]
pub struct Config {
    pub master_seed: i64, pub slots_per_day: i64, pub slots_per_hour: i64, pub parity_days: i64, pub arr_scale_bp: i64,
    pub srv_cnt: [i64; N_STN], pub mu_bp: [i64; N_STN], pub security_extra_prob_bp: i64, pub security_extra_mu_bp: i64,
    pub belt_cap: i64, pub break_prob_bp: i64, pub repair_prob_bp: i64, pub patience_slots: i64, pub gate_close_lead: i64,
    pub bag_miss_bp: i64, pub intl_prob_bp: i64, pub biz_prob_bp: i64, pub transit_prob_bp: i64, pub rebook_prob_bp: i64,
}

#[derive(Clone)]
pub struct RunParams {
    pub days: i64, pub arr_scale_bp: i64, pub break_prob_bp: i64, pub repair_prob_bp: i64, pub belt_cap: i64, pub gate_close_lead: i64,
    pub srv_ovr: [i64; N_STN], pub mu_ovr: [i64; N_STN],
    pub priority_on: i32, pub breakdown_on: i32, pub roster_on: i32, pub customs_on: i32, pub permute_srv: i32, pub timescale_div2: i32, pub audit_on: i32,
}

pub fn default_rp(c: &Config, days: i64) -> RunParams {
    RunParams {
        days, arr_scale_bp: c.arr_scale_bp, break_prob_bp: c.break_prob_bp, repair_prob_bp: c.repair_prob_bp,
        belt_cap: c.belt_cap, gate_close_lead: c.gate_close_lead, srv_ovr: [-1; N_STN], mu_ovr: [-1; N_STN],
        priority_on: 1, breakdown_on: 1, roster_on: 1, customs_on: 0, permute_srv: 0, timescale_div2: 0, audit_on: 0,
    }
}

pub fn load_config(path: &str) -> Config {
    let mut c = Config::default();
    for line in fs::read_to_string(path).expect("config.kv").lines() {
        if let Some((k, v)) = line.split_once('=') {
            let val: i64 = v.trim().parse().unwrap_or(0);
            match k.trim() {
                "master_seed" => c.master_seed = val, "slots_per_day" => c.slots_per_day = val,
                "slots_per_hour" => c.slots_per_hour = val, "parity_days" => c.parity_days = val, "arr_scale_bp" => c.arr_scale_bp = val,
                "checkin_servers" => c.srv_cnt[0] = val, "security_servers" => c.srv_cnt[1] = val, "imigrasi_servers" => c.srv_cnt[2] = val,
                "gate_servers" => c.srv_cnt[3] = val, "belt_servers" => c.srv_cnt[4] = val, "customs_servers" => c.srv_cnt[5] = val,
                "checkin_mu_bp" => c.mu_bp[0] = val, "security_mu_bp" => c.mu_bp[1] = val, "security_extra_prob_bp" => c.security_extra_prob_bp = val,
                "security_extra_mu_bp" => c.security_extra_mu_bp = val, "imigrasi_mu_bp" => c.mu_bp[2] = val, "gate_mu_bp" => c.mu_bp[3] = val,
                "belt_mu_bp" => c.mu_bp[4] = val, "customs_mu_bp" => c.mu_bp[5] = val, "belt_cap" => c.belt_cap = val,
                "break_prob_bp" => c.break_prob_bp = val, "repair_prob_bp" => c.repair_prob_bp = val, "patience_slots" => c.patience_slots = val,
                "gate_close_lead" => c.gate_close_lead = val, "bag_miss_bp" => c.bag_miss_bp = val, "intl_prob_bp" => c.intl_prob_bp = val,
                "biz_prob_bp" => c.biz_prob_bp = val, "transit_prob_bp" => c.transit_prob_bp = val, "rebook_prob_bp" => c.rebook_prob_bp = val, _ => {}
            }
        }
    }
    c
}
