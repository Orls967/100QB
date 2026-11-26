from population import BP, ST_S, ST_E, ST_IA, ST_IS, ST_H, ST_ICU, ST_R, ST_D
from rng import fnv1a_i64


def sim_set_state(s, pid: int, new_st: int, day: int, dur: int):
    old_st = s.pop[pid].status
    s.pop[pid].status = new_st
    s.pop[pid].state_day = day
    s.pop[pid].dur_days = dur
    if s.log_enabled:
        s.progs.append((day, pid, old_st, new_st))


def sim_disease_step(s, day: int):
    n = s.n
    for i in range(n):
        p = s.pop[i]
        elapsed = day - p.state_day
        st = p.status
        if st == ST_E and elapsed >= p.dur_days:
            dur = s.rng[2].rand_int(s.cfg["inf_days_min"], s.cfg["inf_days_max"])
            if s.rng[2].rand_int(1, BP) <= s.cfg["asymp_prob_bp"]:
                s.branches[6] += 1
                sim_set_state(s, i, ST_IA, day, dur)
            else:
                s.branches[7] += 1
                sim_set_state(s, i, ST_IS, day, dur)
        elif (st == ST_IA or st == ST_IS) and elapsed >= p.dur_days:
            if st == ST_IS:
                p_hosp = 120 if p.age_grp == 0 else (450 if p.age_grp == 1 else 1800)
                if p.vax_doses > 0:
                    p_hosp = (p_hosp * (BP - s.cfg["vax_hosp_bp"])) // BP
                if s.rng[4].rand_int(1, BP) <= p_hosp:
                    if s.hosp_occ < s.hosp_cap:
                        s.hosp_occ += 1
                        p.hosp_rej = 0
                        s.branches[8] += 1
                    else:
                        p.hosp_rej = 1
                        s.branches[9] += 1
                    dur = s.rng[4].rand_int(4, 8)
                    sim_set_state(s, i, ST_H, day, dur)
                    continue
            s.branches[13] += 1
            sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
        elif st == ST_H and elapsed >= p.dur_days:
            if p.hosp_rej == 0 and s.hosp_occ > 0:
                s.hosp_occ -= 1
            p_icu = 3500 if p.age_grp == 2 else 1800
            if s.rng[4].rand_int(1, BP) <= p_icu:
                if s.icu_occ < s.icu_cap:
                    s.icu_occ += 1
                    p.hosp_rej = 0
                    s.branches[10] += 1
                else:
                    p.hosp_rej = 1
                    s.branches[11] += 1
                dur = s.rng[4].rand_int(4, 7)
                sim_set_state(s, i, ST_ICU, day, dur)
            else:
                p_die = 100 if p.age_grp == 0 else (500 if p.age_grp == 1 else 2000)
                if p.hosp_rej != 0:
                    p_die *= 2
                if s.rng[4].rand_int(1, BP) <= p_die:
                    s.branches[12] += 1
                    s.daily_deaths_age[day][p.age_grp] += 1
                    sim_set_state(s, i, ST_D, day, 0)
                else:
                    s.branches[13] += 1
                    sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
        elif st == ST_ICU and elapsed >= p.dur_days:
            if p.hosp_rej == 0 and s.icu_occ > 0:
                s.icu_occ -= 1
            p_die = 1500 if p.age_grp == 0 else (2800 if p.age_grp == 1 else 4800)
            if p.hosp_rej != 0:
                p_die *= 2
            if p_die > BP:
                p_die = BP
            if s.rng[4].rand_int(1, BP) <= p_die:
                s.branches[12] += 1
                s.daily_deaths_age[day][p.age_grp] += 1
                sim_set_state(s, i, ST_D, day, 0)
            else:
                s.branches[13] += 1
                sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
        elif st == ST_R and elapsed >= p.dur_days:
            if s.rng[2].rand_int(1, BP) <= s.cfg["waning_bp"]:
                s.branches[14] += 1
                sim_set_state(s, i, ST_S, day, 0)

    for c in range(len(s.cand_tgt)):
        tgt, src, lay, var = s.cand_tgt[c], s.cand_src[c], s.cand_lay[c], s.cand_var[c]
        pt = s.pop[tgt]
        if pt.status != ST_S:
            continue
        if pt.ever_inf != 0:
            s.total_reinf += 1
        pt.ever_inf = 1
        pt.variant = var
        pt.inf_day = day
        s.pop[src].offspring += 1
        s.layer_inf[lay] += 1
        s.branches[lay] += 1
        ag = pt.age_grp
        s.daily_inc_age[day][ag] += 1
        s.age_inf[ag] += 1
        s.daily_var[day][var] += 1
        s.total_trans += 1
        gi = max(day - s.pop[src].inf_day, 1)
        if s.log_enabled:
            s.trans.append((day, src, tgt, lay, var, gi))
        s.event_hash = fnv1a_i64(s.event_hash, day * 1000000 + tgt)
        sim_set_state(s, tgt, ST_E, day, s.cfg["latent_days"])
