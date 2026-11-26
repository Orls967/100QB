use crate::contact::sim_contact_step;
use crate::disease::{sim_disease_step, sim_set_state};
use crate::interventions::sim_interv_step;
use crate::population::*;

pub fn analytical_r0(c: &Config, beta_bp: i64, c_hh: f64, c_sw: f64, c_comm: f64) -> f64 {
    let beta = beta_bp as f64 / BP as f64;
    let d_inf = 0.5 * (c.inf_days_min + c.inf_days_max) as f64;
    let p_ia = c.asymp_prob_bp as f64 / BP as f64;
    let r_ia = c.asymp_rel_bp as f64 / BP as f64;
    beta * (c_hh + c_sw + c_comm) * d_inf * ((1.0 - p_ia) + p_ia * r_ia)
}

pub fn solve_final_size(r0: f64) -> f64 {
    if r0 <= 1.0 { return 0.0; }
    let mut z = 0.5;
    for _ in 0..80 { z = 1.0 - (-r0 * z).exp(); }
    z
}

pub fn check_invariants(s: &mut Sim, day: usize) -> bool {
    let mut cnt = [0i64; 8];
    for p in &s.pop { cnt[p.status as usize] += 1; }
    let mut sum = 0i64;
    for k in 0..8 { sum += cnt[k]; s.daily_cnt[day][k] = cnt[k]; }
    let mut ok = sum == s.n;
    if day > 0 && s.daily_cnt[day][ST_D as usize] < s.daily_cnt[day - 1][ST_D as usize] { ok = false; }
    if s.hosp_occ > s.hosp_cap || s.icu_occ > s.icu_cap { ok = false; }
    if s.log_enabled && s.trans.len() as i64 != s.total_trans { ok = false; }
    if !ok { s.inv_violations += 1; }
    ok
}

pub fn sim_run(s: &mut Sim) {
    let n = s.n as usize;
    let mut n_seed = s.cfg.init_infected;
    if s.rp.beta_bp == 0 && s.rp.contact_scale_bp < 0 { n_seed = 0; }
    if s.rp.vax_perfect != 0 { for p in &mut s.pop { p.vax_doses = 2; p.vax_day = 0; } }
    if n_seed > 0 {
        let mut seeded = 0i64;
        for _ in 0..n {
            if seeded >= n_seed { break; }
            let pid = s.rng[2].rand_int(0, n as i64 - 1) as usize;
            if s.pop[pid].status == ST_S {
                s.pop[pid].ever_inf = 1; s.pop[pid].inf_day = 0;
                sim_set_state(s, pid, ST_IS, 0, 4); seeded += 1;
            }
        }
    }
    for d in 0..(s.rp.days as usize) {
        let di = d as i32;
        sim_interv_step(s, di);
        sim_contact_step(s, di);
        sim_disease_step(s, di);
        s.daily_hosp[d][0] = s.hosp_occ; s.daily_hosp[d][1] = s.icu_occ;
        s.daily_hosp[d][2] = s.branches[9]; s.daily_hosp[d][3] = s.branches[11];
        s.daily_quar[d] = s.pop.iter().filter(|p| p.quar_until >= di).count() as i64;
        check_invariants(s, d);
    }
}
