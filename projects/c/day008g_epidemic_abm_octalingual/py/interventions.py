from disease import sim_set_state
from population import BP, ST_S, ST_IS, ST_D


def sim_interv_step(s, day: int):
    n = s.n
    vd = s.rp.var_day_ovr if s.rp.var_day_ovr >= 0 else s.cfg["variant_day"]
    if day == vd:
        imported = 0
        for _ in range(n):
            if imported >= 5:
                break
            pid = s.rng[3].rand_int(0, n - 1)
            if s.pop[pid].status == ST_S:
                s.pop[pid].variant = 1
                s.pop[pid].ever_inf = 1
                s.pop[pid].inf_day = day
                sim_set_state(s, pid, ST_IS, day, 4)
                s.branches[20] += 1
                imported += 1
    if s.rp.interv_on == 0:
        return

    if day >= s.cfg["test_day"]:
        tp = s.rp.test_ovr if s.rp.test_ovr >= 0 else s.cfg["test_prob_bp"]
        cp = s.rp.compliance_ovr
        for i in range(n):
            p = s.pop[i]
            if p.status == ST_IS and p.quar_until < day:
                eff_comp = cp if cp >= 0 else p.compliance
                if s.rng[3].rand_int(1, BP) <= (tp * eff_comp) // BP:
                    s.daily_tests[day][0] += 1
                    if s.rng[3].rand_int(1, BP) <= 8500:
                        s.daily_tests[day][1] += 1
                        s.branches[15] += 1
                        p.quar_until = day + 7
                        h = p.hh_id
                        for idx in range(s.hh_off[h], s.hh_off[h + 1]):
                            j = s.hh_mem[idx]
                            if j != i and s.pop[j].quar_until < day and s.rng[3].rand_int(1, BP) <= s.cfg["trace_prob_bp"]:
                                s.pop[j].quar_until = day + 5
                                s.branches[17] += 1
                    else:
                        s.branches[16] += 1

    vday = s.rp.vax_day_ovr if s.rp.vax_day_ovr >= 0 else s.cfg["vax_day"]
    if day >= vday:
        cov_cap = (n * s.rp.vax_cov_ovr) // BP if s.rp.vax_cov_ovr >= 0 else (n * 8000) // BP
        doses_today = max((n * s.cfg["vax_daily_bp"]) // BP, 1)
        total_vaxed = sum(1 for p in s.pop if p.vax_doses > 0)
        if total_vaxed >= cov_cap:
            s.branches[19] += 1
            return
        given = 0
        for grp in (2, 1, 0):
            if given >= doses_today:
                break
            for p in s.pop:
                if given >= doses_today:
                    break
                if p.age_grp != grp or p.status == ST_D:
                    continue
                if p.vax_doses == 0 and total_vaxed < cov_cap:
                    p.vax_doses = 1
                    p.vax_day = day
                    given += 1
                    total_vaxed += 1
                    s.daily_vax[day][0] += 1
                    s.branches[18] += 1
                elif p.vax_doses == 1 and day - p.vax_day >= 21:
                    p.vax_doses = 2
                    p.vax_day = day
                    given += 1
                    s.daily_vax[day][1] += 1
                    s.branches[18] += 1
        if given == 0:
            s.branches[19] += 1
