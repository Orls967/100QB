use crate::population::*;
use crate::rng::fnv1a_i64;

pub fn sim_set_state(s: &mut Sim, pid: usize, new_st: i32, day: i32, dur: i32) {
    let old_st = s.pop[pid].status;
    s.pop[pid].status = new_st; s.pop[pid].state_day = day; s.pop[pid].dur_days = dur;
    if s.log_enabled { s.progs.push((day, pid, old_st, new_st)); }
}

pub fn sim_disease_step(s: &mut Sim, day: i32) {
    let n = s.n as usize;
    for i in 0..n {
        let elapsed = day - s.pop[i].state_day;
        let st = s.pop[i].status;
        if st == ST_E && elapsed >= s.pop[i].dur_days {
            let dur = s.rng[2].rand_int(s.cfg.inf_days_min, s.cfg.inf_days_max) as i32;
            if s.rng[2].rand_int(1, BP) <= s.cfg.asymp_prob_bp { s.branches[6] += 1; sim_set_state(s, i, ST_IA, day, dur); }
            else { s.branches[7] += 1; sim_set_state(s, i, ST_IS, day, dur); }
        } else if (st == ST_IA || st == ST_IS) && elapsed >= s.pop[i].dur_days {
            if st == ST_IS {
                let mut p_hosp = match s.pop[i].age_grp { 0 => 120, 1 => 450, _ => 1800 };
                if s.pop[i].vax_doses > 0 { p_hosp = (p_hosp * (BP - s.cfg.vax_hosp_bp)) / BP; }
                if s.rng[4].rand_int(1, BP) <= p_hosp {
                    if s.hosp_occ < s.hosp_cap { s.hosp_occ += 1; s.pop[i].hosp_rej = 0; s.branches[8] += 1; }
                    else { s.pop[i].hosp_rej = 1; s.branches[9] += 1; }
                    let dur = s.rng[4].rand_int(4, 8) as i32;
                    sim_set_state(s, i, ST_H, day, dur);
                    continue;
                }
            }
            s.branches[13] += 1; sim_set_state(s, i, ST_R, day, s.cfg.waning_min_days as i32);
        } else if st == ST_H && elapsed >= s.pop[i].dur_days {
            if s.pop[i].hosp_rej == 0 && s.hosp_occ > 0 { s.hosp_occ -= 1; }
            let p_icu = if s.pop[i].age_grp == 2 { 3500 } else { 1800 };
            if s.rng[4].rand_int(1, BP) <= p_icu {
                if s.icu_occ < s.icu_cap { s.icu_occ += 1; s.pop[i].hosp_rej = 0; s.branches[10] += 1; }
                else { s.pop[i].hosp_rej = 1; s.branches[11] += 1; }
                let dur = s.rng[4].rand_int(4, 7) as i32;
                sim_set_state(s, i, ST_ICU, day, dur);
            } else {
                let mut p_die = match s.pop[i].age_grp { 0 => 100, 1 => 500, _ => 2000 };
                if s.pop[i].hosp_rej != 0 { p_die *= 2; }
                if s.rng[4].rand_int(1, BP) <= p_die { s.branches[12] += 1; s.daily_deaths_age[day as usize][s.pop[i].age_grp as usize] += 1; sim_set_state(s, i, ST_D, day, 0); }
                else { s.branches[13] += 1; sim_set_state(s, i, ST_R, day, s.cfg.waning_min_days as i32); }
            }
        } else if st == ST_ICU && elapsed >= s.pop[i].dur_days {
            if s.pop[i].hosp_rej == 0 && s.icu_occ > 0 { s.icu_occ -= 1; }
            let mut p_die = match s.pop[i].age_grp { 0 => 1500, 1 => 2800, _ => 4800 };
            if s.pop[i].hosp_rej != 0 { p_die *= 2; }
            p_die = p_die.min(BP);
            if s.rng[4].rand_int(1, BP) <= p_die { s.branches[12] += 1; s.daily_deaths_age[day as usize][s.pop[i].age_grp as usize] += 1; sim_set_state(s, i, ST_D, day, 0); }
            else { s.branches[13] += 1; sim_set_state(s, i, ST_R, day, s.cfg.waning_min_days as i32); }
        } else if st == ST_R && elapsed >= s.pop[i].dur_days {
            if s.rng[2].rand_int(1, BP) <= s.cfg.waning_bp { s.branches[14] += 1; sim_set_state(s, i, ST_S, day, 0); }
        }
    }
    for c in 0..s.cand_tgt.len() {
        let (tgt, src, lay, var) = (s.cand_tgt[c], s.cand_src[c], s.cand_lay[c], s.cand_var[c]);
        if s.pop[tgt].status != ST_S { continue; }
        if s.pop[tgt].ever_inf != 0 { s.total_reinf += 1; }
        s.pop[tgt].ever_inf = 1; s.pop[tgt].variant = var; s.pop[tgt].inf_day = day;
        s.pop[src].offspring += 1; s.layer_inf[lay] += 1; s.branches[lay] += 1;
        let ag = s.pop[tgt].age_grp as usize;
        s.daily_inc_age[day as usize][ag] += 1; s.age_inf[ag] += 1;
        s.daily_var[day as usize][var as usize] += 1; s.total_trans += 1;
        let gi = (day - s.pop[src].inf_day).max(1);
        if s.log_enabled { s.trans.push((day, src, tgt, lay, var, gi)); }
        s.event_hash = fnv1a_i64(s.event_hash, (day as i64) * 1_000_000 + tgt as i64);
        sim_set_state(s, tgt, ST_E, day, s.cfg.latent_days as i32);
    }
}
