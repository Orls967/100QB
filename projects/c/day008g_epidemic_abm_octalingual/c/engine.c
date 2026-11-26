#include "common.h"
#include <math.h>

double analytical_r0(const Config *c, int64_t beta_bp, double c_hh, double c_sw, double c_comm) {
    double beta = (double)beta_bp / (double)BP;
    double d_inf = 0.5 * (double)(c->inf_days_min + c->inf_days_max);
    double p_ia = (double)c->asymp_prob_bp / (double)BP;
    double r_ia = (double)c->asymp_rel_bp / (double)BP;
    double w_inf = (1.0 - p_ia) + p_ia * r_ia;
    return beta * (c_hh + c_sw + c_comm) * d_inf * w_inf;
}

double solve_final_size(double r0) {
    if (r0 <= 1.0) return 0.0;
    double z = 0.5;
    for (int it = 0; it < 80; ++it) z = 1.0 - exp(-r0 * z);
    return z;
}

int check_invariants(Sim *s, int32_t day) {
    int64_t cnt[8] = {0};
    for (int64_t i = 0; i < s->N; ++i) cnt[s->pop[i].status]++;
    int64_t sum = 0;
    for (int k = 0; k < 8; ++k) { sum += cnt[k]; s->daily_cnt[day][k] = cnt[k]; }
    int ok = 1;
    if (sum != s->N) ok = 0;
    if (day > 0 && s->daily_cnt[day][ST_D] < s->daily_cnt[day - 1][ST_D]) ok = 0;
    if (s->hosp_occ > s->hosp_cap || s->icu_occ > s->icu_cap) ok = 0;
    if (s->log_enabled && s->n_trans != s->total_trans) ok = 0;
    if (!ok) s->inv_violations++;
    return ok;
}

void sim_run(Sim *s) {
    int32_t N = (int32_t)s->N;
    int64_t n_seed = s->cfg.init_infected;
    if (s->rp.beta_bp == 0 && s->rp.contact_scale_bp < 0) n_seed = 0;
    if (s->rp.vax_perfect) {
        for (int32_t i = 0; i < N; ++i) { s->pop[i].vax_doses = 2; s->pop[i].vax_day = 0; }
    }
    /* Seed initial infections deterministically */
    if (n_seed > 0) {
        int32_t seeded = 0;
        for (int32_t i = 0; i < N && seeded < n_seed; ++i) {
            int32_t pid = (int32_t)rand_int(&s->rng[2], 0, N - 1);
            if (s->pop[pid].status == ST_S) {
                s->pop[pid].ever_inf = 1; s->pop[pid].inf_day = 0;
                sim_set_state(s, pid, ST_IS, 0, 4);
                seeded++;
            }
        }
    }
    for (int32_t d = 0; d < s->rp.days; ++d) {
        sim_interv_step(s, d);
        sim_contact_step(s, d);
        sim_disease_step(s, d);
        s->daily_hosp[d][0] = s->hosp_occ; s->daily_hosp[d][1] = s->icu_occ;
        s->daily_hosp[d][2] = s->branches[BR_HOSP_REJ]; s->daily_hosp[d][3] = s->branches[BR_ICU_REJ];
        int64_t qcnt = 0;
        for (int32_t i = 0; i < N; ++i) if (s->pop[i].quar_until >= d) qcnt++;
        s->daily_quar[d] = qcnt;
        check_invariants(s, d);
    }
}
