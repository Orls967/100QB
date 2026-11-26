#include "common.hpp"

void sim_run(Sim &s) {
    int64_t total_slots = s.rp.days * (s.rp.timescale_div2 ? SLOTS_PER_DAY * 2 : SLOTS_PER_DAY);
    for (int64_t t = 0; t < total_slots; ++t) sim_step_slot(s, t);
    int64_t rem_pax = 0;
    for (int i = 0; i < MAX_POOL; ++i) {
        if (s.pool[i].active && !s.pool[i].is_bag) { rem_pax++; s.sys_sum_sojourn += (total_slots - s.pool[i].arr_slot); }
    }
    int64_t stn_rem = 0;
    for (int st = 0; st < N_STN; ++st) if (st != STN_BELT) stn_rem += (s.stn[st].len + s.stn[st].busy_cnt);
    if (s.total_arr != s.total_done + s.total_lost + rem_pax || stn_rem != rem_pax) s.inv_violations++;
    if (s.sys_sum_L != s.sys_sum_sojourn) s.inv_violations++;
}
