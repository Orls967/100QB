#include "common.hpp"

void sim_interv_step(Sim &s, int32_t day) {
    int32_t N = static_cast<int32_t>(s.N);
    Xoshiro256 *ri = &s.rng[3];

    int64_t vd = s.rp.var_day_ovr >= 0 ? s.rp.var_day_ovr : s.cfg.variant_day;
    if (day == vd) {
        int32_t imported = 0;
        for (int32_t i = 0; i < N && imported < 5; ++i) {
            int32_t pid = static_cast<int32_t>(rand_int(ri, 0, N - 1));
            if (s.pop[pid].status == ST_S) {
                s.pop[pid].variant = 1; s.pop[pid].ever_inf = 1; s.pop[pid].inf_day = day;
                sim_set_state(s, pid, ST_IS, day, 4);
                s.branches[BR_VAR_IMPORT]++; imported++;
            }
        }
    }
    if (!s.rp.interv_on) return;

    if (day >= s.cfg.test_day) {
        int64_t tp = s.rp.test_ovr >= 0 ? s.rp.test_ovr : s.cfg.test_prob_bp;
        int64_t cp = s.rp.compliance_ovr >= 0 ? s.rp.compliance_ovr : -1;
        for (int32_t i = 0; i < N; ++i) {
            Person &p = s.pop[i];
            if (p.status == ST_IS && p.quar_until < day) {
                int64_t eff_comp = cp >= 0 ? cp : p.compliance;
                if (rand_int(ri, 1, BP) <= (tp * eff_comp) / BP) {
                    s.daily_tests[day][0]++;
                    if (rand_int(ri, 1, BP) <= 8500) {
                        s.daily_tests[day][1]++; s.branches[BR_TEST_POS]++;
                        p.quar_until = day + 7;
                        int32_t h = p.hh_id;
                        for (int32_t idx = s.hh_off[h]; idx < s.hh_off[h + 1]; ++idx) {
                            int32_t j = s.hh_mem[idx];
                            if (j != i && s.pop[j].quar_until < day && rand_int(ri, 1, BP) <= s.cfg.trace_prob_bp) {
                                s.pop[j].quar_until = day + 5; s.branches[BR_QUAR_TRACE]++;
                            }
                        }
                    } else s.branches[BR_TEST_NEG]++;
                }
            }
        }
    }

    int64_t vday = s.rp.vax_day_ovr >= 0 ? s.rp.vax_day_ovr : s.cfg.vax_day;
    if (day >= vday) {
        int64_t cov_cap = s.rp.vax_cov_ovr >= 0 ? (N * s.rp.vax_cov_ovr) / BP : (N * 8000) / BP;
        int64_t doses_today = std::max<int64_t>(1, (N * s.cfg.vax_daily_bp) / BP);
        int64_t total_vaxed = 0;
        for (int32_t i = 0; i < N; ++i) if (s.pop[i].vax_doses > 0) total_vaxed++;
        if (total_vaxed >= cov_cap) { s.branches[BR_VAX_OUT]++; return; }
        int64_t given = 0;
        for (int32_t grp = 2; grp >= 0 && given < doses_today; --grp) {
            for (int32_t i = 0; i < N && given < doses_today; ++i) {
                Person &p = s.pop[i];
                if (p.age_grp != grp || p.status == ST_D) continue;
                if (p.vax_doses == 0 && total_vaxed < cov_cap) {
                    p.vax_doses = 1; p.vax_day = day; given++; total_vaxed++;
                    s.daily_vax[day][0]++; s.branches[BR_VAX_DOSE]++;
                } else if (p.vax_doses == 1 && day - p.vax_day >= 21) {
                    p.vax_doses = 2; p.vax_day = day; given++;
                    s.daily_vax[day][1]++; s.branches[BR_VAX_DOSE]++;
                }
            }
        }
        if (given == 0) s.branches[BR_VAX_OUT]++;
    }
}
