use crate::disease::sim_set_state;
use crate::population::*;

pub fn sim_interv_step(s: &mut Sim, day: i32) {
    let n = s.n as usize;
    let vd = if s.rp.var_day_ovr >= 0 { s.rp.var_day_ovr } else { s.cfg.variant_day };
    if day as i64 == vd {
        let mut imported = 0;
        for _ in 0..n {
            if imported >= 5 { break; }
            let pid = s.rng[3].rand_int(0, n as i64 - 1) as usize;
            if s.pop[pid].status == ST_S {
                s.pop[pid].variant = 1; s.pop[pid].ever_inf = 1; s.pop[pid].inf_day = day;
                sim_set_state(s, pid, ST_IS, day, 4);
                s.branches[20] += 1; imported += 1;
            }
        }
    }
    if s.rp.interv_on == 0 { return; }

    if day as i64 >= s.cfg.test_day {
        let tp = if s.rp.test_ovr >= 0 { s.rp.test_ovr } else { s.cfg.test_prob_bp };
        let cp = s.rp.compliance_ovr;
        for i in 0..n {
            if s.pop[i].status == ST_IS && s.pop[i].quar_until < day {
                let eff_comp = if cp >= 0 { cp } else { s.pop[i].compliance as i64 };
                if s.rng[3].rand_int(1, BP) <= (tp * eff_comp) / BP {
                    s.daily_tests[day as usize][0] += 1;
                    if s.rng[3].rand_int(1, BP) <= 8500 {
                        s.daily_tests[day as usize][1] += 1; s.branches[15] += 1;
                        s.pop[i].quar_until = day + 7;
                        let h = s.pop[i].hh_id as usize;
                        for idx in s.hh_off[h]..s.hh_off[h + 1] {
                            let j = s.hh_mem[idx];
                            if j != i && s.pop[j].quar_until < day && s.rng[3].rand_int(1, BP) <= s.cfg.trace_prob_bp {
                                s.pop[j].quar_until = day + 5; s.branches[17] += 1;
                            }
                        }
                    } else { s.branches[16] += 1; }
                }
            }
        }
    }

    let vday = if s.rp.vax_day_ovr >= 0 { s.rp.vax_day_ovr } else { s.cfg.vax_day };
    if day as i64 >= vday {
        let cov_cap = if s.rp.vax_cov_ovr >= 0 { (n as i64 * s.rp.vax_cov_ovr) / BP } else { (n as i64 * 8000) / BP };
        let doses_today = ((n as i64 * s.cfg.vax_daily_bp) / BP).max(1);
        let mut total_vaxed = s.pop.iter().filter(|p| p.vax_doses > 0).count() as i64;
        if total_vaxed >= cov_cap { s.branches[19] += 1; return; }
        let mut given = 0i64;
        for grp in (0..=2).rev() {
            if given >= doses_today { break; }
            for i in 0..n {
                if given >= doses_today { break; }
                if s.pop[i].age_grp != grp || s.pop[i].status == ST_D { continue; }
                if s.pop[i].vax_doses == 0 && total_vaxed < cov_cap {
                    s.pop[i].vax_doses = 1; s.pop[i].vax_day = day; given += 1; total_vaxed += 1;
                    s.daily_vax[day as usize][0] += 1; s.branches[18] += 1;
                } else if s.pop[i].vax_doses == 1 && day - s.pop[i].vax_day >= 21 {
                    s.pop[i].vax_doses = 2; s.pop[i].vax_day = day; given += 1;
                    s.daily_vax[day as usize][1] += 1; s.branches[18] += 1;
                }
            }
        }
        if given == 0 { s.branches[19] += 1; }
    }
}
