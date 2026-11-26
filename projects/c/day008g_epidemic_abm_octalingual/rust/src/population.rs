use std::fs;
use crate::rng::*;

pub const BP: i64 = 10000;
pub const MAX_DAYS: usize = 200;
pub const ST_S: i32 = 0; pub const ST_E: i32 = 1; pub const ST_IA: i32 = 2; pub const ST_IS: i32 = 3;
pub const ST_H: i32 = 4; pub const ST_ICU: i32 = 5; pub const ST_R: i32 = 6; pub const ST_D: i32 = 7;

#[derive(Clone, Default)]
pub struct Config {
    pub master_seed: i64, pub parity_n: i64, pub parity_days: i64, pub main_days: i64,
    pub beta_bp_init: i64, pub asymp_prob_bp: i64, pub asymp_rel_bp: i64,
    pub latent_days: i64, pub inf_days_min: i64, pub inf_days_max: i64,
    pub hosp_beds_per_10k: i64, pub icu_beds_per_10k: i64,
    pub waning_min_days: i64, pub waning_bp: i64, pub variant_day: i64, pub variant_mult_bp: i64,
    pub school_close_day: i64, pub lockdown_day: i64, pub mask_day: i64, pub mask_eff_bp: i64,
    pub test_day: i64, pub test_prob_bp: i64, pub trace_prob_bp: i64,
    pub vax_day: i64, pub vax_daily_bp: i64, pub vax_eff1_bp: i64, pub vax_eff2_bp: i64, pub vax_hosp_bp: i64,
    pub init_infected: i64,
}

#[derive(Clone, Default)]
pub struct Person {
    pub age_grp: i32, pub hh_id: i32, pub sw_id: i32, pub status: i32, pub variant: i32,
    pub state_day: i32, pub dur_days: i32, pub inf_day: i32, pub vax_doses: i32, pub vax_day: i32,
    pub compliance: i32, pub quar_until: i32, pub super_flag: i32, pub offspring: i32, pub ever_inf: i32, pub hosp_rej: i32,
}

#[derive(Clone)]
pub struct RunParams {
    pub n: i64, pub days: i64, pub beta_bp: i64,
    pub compliance_ovr: i64, pub test_ovr: i64, pub vax_day_ovr: i64, pub vax_cov_ovr: i64,
    pub school_ovr: i64, pub lock_ovr: i64, pub var_day_ovr: i64, pub mask_ovr: i64,
    pub contact_scale_bp: i64, pub vax_perfect: i64, pub interv_on: i32,
}

impl Default for RunParams {
    fn default() -> Self {
        Self { n: 0, days: 0, beta_bp: 0, compliance_ovr: -1, test_ovr: -1, vax_day_ovr: -1, vax_cov_ovr: -1,
               school_ovr: -1, lock_ovr: -1, var_day_ovr: -1, mask_ovr: -1, contact_scale_bp: 10000, vax_perfect: 0, interv_on: 0 }
    }
}

pub struct Sim {
    pub cfg: Config, pub rp: RunParams, pub n: i64, pub hosp_cap: i64, pub icu_cap: i64, pub hosp_occ: i64, pub icu_occ: i64,
    pub pop: Vec<Person>, pub n_hh: usize, pub n_sw: usize,
    pub hh_off: Vec<usize>, pub hh_mem: Vec<usize>, pub sw_off: Vec<usize>, pub sw_mem: Vec<usize>, pub seen_day: Vec<i32>,
    pub rng: Vec<Xoshiro256>,
    pub cand_tgt: Vec<usize>, pub cand_src: Vec<usize>, pub cand_lay: Vec<usize>, pub cand_var: Vec<i32>,
    pub daily_cnt: [[i64; 8]; MAX_DAYS], pub daily_inc_age: [[i64; 3]; MAX_DAYS], pub daily_deaths_age: [[i64; 3]; MAX_DAYS],
    pub daily_hosp: [[i64; 4]; MAX_DAYS], pub daily_tests: [[i64; 2]; MAX_DAYS], pub daily_quar: [i64; MAX_DAYS],
    pub daily_vax: [[i64; 2]; MAX_DAYS], pub daily_var: [[i64; 2]; MAX_DAYS],
    pub layer_inf: [i64; 3], pub age_inf: [i64; 3], pub age_tot: [i64; 3], pub branches: [i64; 21],
    pub inv_violations: i64, pub total_trans: i64, pub total_reinf: i64, pub event_hash: u64, pub log_enabled: bool,
    pub trans: Vec<(i32, usize, usize, usize, i32, i32)>,
    pub progs: Vec<(i32, usize, i32, i32)>,
}

pub fn load_config(path: &str) -> Config {
    let mut c = Config::default();
    for line in fs::read_to_string(path).expect("config.kv").lines() {
        if let Some((k, v)) = line.split_once('=') {
            let val: i64 = v.trim().parse().unwrap_or(0);
            match k.trim() {
                "master_seed" => c.master_seed = val, "parity_N" => c.parity_n = val, "parity_days" => c.parity_days = val,
                "main_days" => c.main_days = val, "beta_bp_init" => c.beta_bp_init = val, "asymp_prob_bp" => c.asymp_prob_bp = val,
                "asymp_rel_bp" => c.asymp_rel_bp = val, "latent_days" => c.latent_days = val, "inf_days_min" => c.inf_days_min = val,
                "inf_days_max" => c.inf_days_max = val, "hosp_beds_per_10k" => c.hosp_beds_per_10k = val, "icu_beds_per_10k" => c.icu_beds_per_10k = val,
                "waning_min_days" => c.waning_min_days = val, "waning_bp" => c.waning_bp = val, "variant_day" => c.variant_day = val,
                "variant_mult_bp" => c.variant_mult_bp = val, "school_close_day" => c.school_close_day = val, "lockdown_day" => c.lockdown_day = val,
                "mask_day" => c.mask_day = val, "mask_eff_bp" => c.mask_eff_bp = val, "test_day" => c.test_day = val,
                "test_prob_bp" => c.test_prob_bp = val, "trace_prob_bp" => c.trace_prob_bp = val, "vax_day" => c.vax_day = val,
                "vax_daily_bp" => c.vax_daily_bp = val, "vax_eff1_bp" => c.vax_eff1_bp = val, "vax_eff2_bp" => c.vax_eff2_bp = val,
                "vax_hosp_bp" => c.vax_hosp_bp = val, "init_infected" => c.init_infected = val, _ => {}
            }
        }
    }
    c
}

impl Sim {
    pub fn new(c: &Config, rp: &RunParams, log_en: bool) -> Self {
        let mut rng = Vec::with_capacity(5);
        for i in 0..5u32 { let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(100, i); rng.push(Xoshiro256::seed(splitmix64(&mut x))); }
        let n = rp.n as usize;
        let mut pop = vec![Person { sw_id: -1, inf_day: -1, vax_day: -1, quar_until: -1, ..Default::default() }; n];
        let mut age_tot = [0i64; 3];
        let mut hh_cnt = Vec::new();
        let (mut pid, mut hh) = (0usize, 0usize);
        while pid < n {
            let sz = (rng[0].rand_int(1, 6) as usize).min(n - pid);
            hh_cnt.push(sz);
            for _ in 0..sz {
                pop[pid].hh_id = hh as i32;
                let u = rng[0].rand_int(1, 100);
                pop[pid].age_grp = if u <= 22 { 0 } else if u <= 82 { 1 } else { 2 };
                age_tot[pop[pid].age_grp as usize] += 1;
                pop[pid].compliance = rng[0].rand_int(4000, 9500) as i32;
                pop[pid].super_flag = if rng[0].rand_int(1, 10) == 1 { 1 } else { 0 };
                pid += 1;
            }
            hh += 1;
        }
        let mut hh_off = vec![0usize; hh + 1];
        for i in 0..hh { hh_off[i + 1] = hh_off[i] + hh_cnt[i]; }
        let mut hh_mem = vec![0usize; n];
        let mut cur = vec![0usize; hh];
        for i in 0..n { let h = pop[i].hh_id as usize; hh_mem[hh_off[h] + cur[h]] = i; cur[h] += 1; }

        let n_schools = (age_tot[0] / 25) as usize + 1;
        let n_works = (age_tot[1] / 15) as usize + 1;
        let n_sw = n_schools + n_works;
        let mut sw_cnt = vec![0usize; n_sw];
        for i in 0..n {
            if pop[i].age_grp == 0 { let sw = rng[0].rand_int(0, n_schools as i64 - 1) as usize; pop[i].sw_id = sw as i32; sw_cnt[sw] += 1; }
            else if pop[i].age_grp == 1 { let sw = n_schools + rng[0].rand_int(0, n_works as i64 - 1) as usize; pop[i].sw_id = sw as i32; sw_cnt[sw] += 1; }
        }
        let mut sw_off = vec![0usize; n_sw + 1];
        for i in 0..n_sw { sw_off[i + 1] = sw_off[i] + sw_cnt[i]; }
        let mut sw_mem = vec![0usize; n];
        let mut sw_cur = vec![0usize; n_sw];
        for i in 0..n { if pop[i].sw_id >= 0 { let sw = pop[i].sw_id as usize; sw_mem[sw_off[sw] + sw_cur[sw]] = i; sw_cur[sw] += 1; } }

        Self {
            cfg: c.clone(), rp: rp.clone(), n: rp.n,
            hosp_cap: ((rp.n * c.hosp_beds_per_10k) / BP).max(1), icu_cap: ((rp.n * c.icu_beds_per_10k) / BP).max(1),
            hosp_occ: 0, icu_occ: 0, pop, n_hh: hh, n_sw, hh_off, hh_mem, sw_off, sw_mem, seen_day: vec![-1; n],
            rng, cand_tgt: Vec::new(), cand_src: Vec::new(), cand_lay: Vec::new(), cand_var: Vec::new(),
            daily_cnt: [[0; 8]; MAX_DAYS], daily_inc_age: [[0; 3]; MAX_DAYS], daily_deaths_age: [[0; 3]; MAX_DAYS],
            daily_hosp: [[0; 4]; MAX_DAYS], daily_tests: [[0; 2]; MAX_DAYS], daily_quar: [0; MAX_DAYS],
            daily_vax: [[0; 2]; MAX_DAYS], daily_var: [[0; 2]; MAX_DAYS],
            layer_inf: [0; 3], age_inf: [0; 3], age_tot, branches: [0; 21],
            inv_violations: 0, total_trans: 0, total_reinf: 0, event_hash: 0xcbf29ce484222325, log_enabled: log_en,
            trans: Vec::new(), progs: Vec::new(),
        }
    }
}
