use crate::population::*;

fn vax_eff_bp(s: &Sim, p: &Person, day: i32) -> i64 {
    if p.vax_doses <= 0 || p.vax_day < 0 { return 0; }
    if s.rp.vax_perfect != 0 { return BP; }
    let base = if p.vax_doses == 1 { s.cfg.vax_eff1_bp } else { s.cfg.vax_eff2_bp };
    let elapsed = (day - p.vax_day) as i64;
    if elapsed <= 30 { base } else { (base - (elapsed - 30) * 35).max(2000) }
}

fn try_infect(s: &mut Sim, src: usize, tgt: usize, lay: usize, day: i32, beta_eff: i64) {
    if s.pop[tgt].status != ST_S || s.seen_day[tgt] == day { return; }
    if s.pop[tgt].quar_until >= day { s.branches[5] += 1; return; }
    let mut p_bp = beta_eff;
    if s.rp.interv_on != 0 {
        let m_day = if s.rp.mask_ovr >= 0 { s.rp.mask_ovr } else { s.cfg.mask_day };
        if day as i64 >= m_day {
            let red = (s.cfg.mask_eff_bp * s.pop[src].compliance as i64) / BP;
            let next_bp = (p_bp * (BP - red)) / BP;
            if next_bp < p_bp { s.branches[4] += 1; }
            p_bp = next_bp;
        }
    }
    let ve = vax_eff_bp(s, &s.pop[tgt], day);
    if ve > 0 {
        p_bp = (p_bp * (BP - ve)) / BP;
        if p_bp <= 0 { s.branches[3] += 1; return; }
    }
    if s.rng[1].rand_int(1, BP) <= p_bp {
        s.seen_day[tgt] = day;
        s.cand_tgt.push(tgt); s.cand_src.push(src); s.cand_lay.push(lay); s.cand_var.push(s.pop[src].variant);
    } else if ve > 0 {
        s.branches[3] += 1;
    }
}

pub fn sim_contact_step(s: &mut Sim, day: i32) {
    s.cand_tgt.clear(); s.cand_src.clear(); s.cand_lay.clear(); s.cand_var.clear();
    if s.rp.contact_scale_bp == 0 { return; }
    let n = s.n as usize;
    let sc_day = if s.rp.school_ovr >= 0 { s.rp.school_ovr } else { s.cfg.school_close_day };
    let lk_day = if s.rp.lock_ovr >= 0 { s.rp.lock_ovr } else { s.cfg.lockdown_day };
    let school_closed = s.rp.interv_on != 0 && (day as i64 >= sc_day);
    let locked = s.rp.interv_on != 0 && (day as i64 >= lk_day);

    for i in 0..n {
        let st = s.pop[i].status;
        if st != ST_IA && st != ST_IS { continue; }
        let mut base_bp = s.rp.beta_bp;
        if st == ST_IA { base_bp = (base_bp * s.cfg.asymp_rel_bp) / BP; }
        if s.pop[i].variant == 1 { base_bp = (base_bp * s.cfg.variant_mult_bp) / BP; }
        base_bp = base_bp.min(BP);

        let h = s.pop[i].hh_id as usize;
        for idx in s.hh_off[h]..s.hh_off[h + 1] {
            let j = s.hh_mem[idx];
            if j != i { try_infect(s, i, j, 0, day, base_bp); }
        }
        if s.pop[i].quar_until >= day { s.branches[5] += 1; continue; }

        if s.pop[i].sw_id >= 0 && !(s.pop[i].age_grp == 0 && school_closed) {
            let sw = s.pop[i].sw_id as usize;
            let (lo, hi) = (s.sw_off[sw], s.sw_off[sw + 1]);
            let sz = hi - lo;
            if sz > 1 {
                let mut n_sw_c = if s.pop[i].age_grp == 0 { 6 } else { 5 };
                if locked { n_sw_c = (n_sw_c + 1) / 2; }
                for _ in 0..n_sw_c {
                    let j = s.sw_mem[lo + s.rng[1].rand_int(0, sz as i64 - 1) as usize];
                    if j != i { try_infect(s, i, j, 1, day, base_bp); }
                }
            }
        }
        let mut n_comm = if s.pop[i].super_flag != 0 { 12 } else { 3 };
        if locked { n_comm = (n_comm + 1) / 3; }
        for _ in 0..n_comm {
            let j = s.rng[1].rand_int(0, n as i64 - 1) as usize;
            if j != i { try_infect(s, i, j, 2, day, base_bp); }
        }
    }
}
