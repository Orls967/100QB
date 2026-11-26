use crate::config::{MAX_POOL, N_STN, SLOTS_PER_DAY, STN_BELT};
use crate::network::{sim_step_slot, Sim};

pub fn sim_run(s: &mut Sim) {
    let total_slots = s.rp.days * if s.rp.timescale_div2 != 0 { SLOTS_PER_DAY * 2 } else { SLOTS_PER_DAY };
    for t in 0..total_slots { sim_step_slot(s, t); }
    let mut rem_pax = 0i64;
    for i in 0..MAX_POOL {
        if s.pool[i].active != 0 && s.pool[i].is_bag == 0 {
            rem_pax += 1; s.sys_sum_sojourn += total_slots - s.pool[i].arr_slot;
        }
    }
    let mut stn_rem = 0i64;
    for st in 0..N_STN { if st != STN_BELT { stn_rem += s.stn[st].len as i64 + s.stn[st].busy_cnt as i64; } }
    if s.total_arr != s.total_done + s.total_lost + rem_pax || stn_rem != rem_pax { s.inv_violations += 1; }
    if s.sys_sum_l != s.sys_sum_sojourn { s.inv_violations += 1; }
}
