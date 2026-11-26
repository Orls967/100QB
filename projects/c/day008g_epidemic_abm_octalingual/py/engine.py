import math
from contact import sim_contact_step
from disease import sim_disease_step, sim_set_state
from interventions import sim_interv_step
from population import BP, ST_S, ST_IS, ST_D


def analytical_r0(c: dict, beta_bp: int, c_hh: float, c_sw: float, c_comm: float) -> float:
    beta = beta_bp / float(BP)
    d_inf = 0.5 * (c["inf_days_min"] + c["inf_days_max"])
    p_ia = c["asymp_prob_bp"] / float(BP)
    r_ia = c["asymp_rel_bp"] / float(BP)
    return beta * (c_hh + c_sw + c_comm) * d_inf * ((1.0 - p_ia) + p_ia * r_ia)


def solve_final_size(r0: float) -> float:
    if r0 <= 1.0:
        return 0.0
    z = 0.5
    for _ in range(80):
        z = 1.0 - math.exp(-r0 * z)
    return z


def check_invariants(s, day: int) -> bool:
    cnt = [0] * 8
    for p in s.pop:
        cnt[p.status] += 1
    s.daily_cnt[day] = cnt
    ok = sum(cnt) == s.n
    if day > 0 and s.daily_cnt[day][ST_D] < s.daily_cnt[day - 1][ST_D]:
        ok = False
    if s.hosp_occ > s.hosp_cap or s.icu_occ > s.icu_cap:
        ok = False
    if s.log_enabled and len(s.trans) != s.total_trans:
        ok = False
    if not ok:
        s.inv_violations += 1
    return ok


def sim_run(s):
    n = s.n
    n_seed = s.cfg["init_infected"]
    if s.rp.beta_bp == 0 and s.rp.contact_scale_bp < 0:
        n_seed = 0
    if s.rp.vax_perfect != 0:
        for p in s.pop:
            p.vax_doses = 2
            p.vax_day = 0
    if n_seed > 0:
        seeded = 0
        for _ in range(n):
            if seeded >= n_seed:
                break
            pid = s.rng[2].rand_int(0, n - 1)
            if s.pop[pid].status == ST_S:
                s.pop[pid].ever_inf = 1
                s.pop[pid].inf_day = 0
                sim_set_state(s, pid, ST_IS, 0, 4)
                seeded += 1
    for d in range(s.rp.days):
        sim_interv_step(s, d)
        sim_contact_step(s, d)
        sim_disease_step(s, d)
        s.daily_hosp[d][0] = s.hosp_occ
        s.daily_hosp[d][1] = s.icu_occ
        s.daily_hosp[d][2] = s.branches[9]
        s.daily_hosp[d][3] = s.branches[11]
        s.daily_quar[d] = sum(1 for p in s.pop if p.quar_until >= d)
        check_invariants(s, d)
