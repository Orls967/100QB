from config import N_STN, SLOTS_PER_DAY, STN_BELT
from network import sim_step_slot


def sim_run(s) -> None:
    total_slots = s.rp.days * (SLOTS_PER_DAY * 2 if s.rp.timescale_div2 != 0 else SLOTS_PER_DAY)
    for t in range(total_slots):
        sim_step_slot(s, t)
    rem_pax = 0
    for ag in s.pool:
        if ag.active != 0 and ag.is_bag == 0:
            rem_pax += 1
            s.sys_sum_sojourn += total_slots - ag.arr_slot
    stn_rem = sum(0 if st == STN_BELT else (s.stn[st].len + s.stn[st].busy_cnt) for st in range(N_STN))
    if s.total_arr != s.total_done + s.total_lost + rem_pax or stn_rem != rem_pax:
        s.inv_violations += 1
    if s.sys_sum_l != s.sys_sum_sojourn:
        s.inv_violations += 1
