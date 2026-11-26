#include "common.h"
#include <stdlib.h>

static void push_prog(Sim *s, ProgRec r) {
    if (!s->log_enabled) return;
    if (s->n_progs == s->cap_progs) {
        s->cap_progs = s->cap_progs ? s->cap_progs * 2 : 1024;
        s->progs = (ProgRec *)realloc(s->progs, (size_t)s->cap_progs * sizeof(ProgRec));
    }
    s->progs[s->n_progs++] = r;
}

static void push_trans(Sim *s, TransRec r) {
    if (!s->log_enabled) return;
    if (s->n_trans == s->cap_trans) {
        s->cap_trans = s->cap_trans ? s->cap_trans * 2 : 1024;
        s->trans = (TransRec *)realloc(s->trans, (size_t)s->cap_trans * sizeof(TransRec));
    }
    s->trans[s->n_trans++] = r;
}

void sim_set_state(Sim *s, int32_t pid, int32_t new_st, int32_t day, int32_t dur) {
    int32_t old_st = s->pop[pid].status;
    s->pop[pid].status = new_st;
    s->pop[pid].state_day = day;
    s->pop[pid].dur_days = dur;
    push_prog(s, (ProgRec){day, pid, old_st, new_st});
}

void sim_disease_step(Sim *s, int32_t day) {
    int32_t N = (int32_t)s->N;
    Xoshiro256 *rd = &s->rng[2], *rh = &s->rng[4];

    /* 1. Progress existing states */
    for (int32_t i = 0; i < N; ++i) {
        Person *p = &s->pop[i];
        int32_t elapsed = day - p->state_day;
        if (p->status == ST_E && elapsed >= p->dur_days) {
            int32_t dur = (int32_t)rand_int(rd, s->cfg.inf_days_min, s->cfg.inf_days_max);
            if (rand_int(rd, 1, BP) <= s->cfg.asymp_prob_bp) {
                s->branches[BR_PROG_IA]++;
                sim_set_state(s, i, ST_IA, day, dur);
            } else {
                s->branches[BR_PROG_IS]++;
                sim_set_state(s, i, ST_IS, day, dur);
            }
        } else if ((p->status == ST_IA || p->status == ST_IS) && elapsed >= p->dur_days) {
            if (p->status == ST_IS) {
                int64_t p_hosp = (p->age_grp == 0) ? 120 : ((p->age_grp == 1) ? 450 : 1800);
                if (p->vax_doses > 0) p_hosp = (p_hosp * (BP - s->cfg.vax_hosp_bp)) / BP;
                if (rand_int(rh, 1, BP) <= p_hosp) {
                    if (s->hosp_occ < s->hosp_cap) {
                        s->hosp_occ++; p->hosp_rej = 0; s->branches[BR_HOSP_ADM]++;
                    } else {
                        p->hosp_rej = 1; s->branches[BR_HOSP_REJ]++;
                    }
                    sim_set_state(s, i, ST_H, day, (int32_t)rand_int(rh, 4, 8));
                    continue;
                }
            }
            s->branches[BR_RECOVER]++;
            sim_set_state(s, i, ST_R, day, (int32_t)s->cfg.waning_min_days);
        } else if (p->status == ST_H && elapsed >= p->dur_days) {
            if (!p->hosp_rej && s->hosp_occ > 0) s->hosp_occ--;
            int64_t p_icu = (p->age_grp == 2) ? 3500 : 1800;
            if (rand_int(rh, 1, BP) <= p_icu) {
                if (s->icu_occ < s->icu_cap) {
                    s->icu_occ++; p->hosp_rej = 0; s->branches[BR_ICU_ADM]++;
                } else {
                    p->hosp_rej = 1; s->branches[BR_ICU_REJ]++;
                }
                sim_set_state(s, i, ST_ICU, day, (int32_t)rand_int(rh, 4, 7));
            } else {
                int64_t p_die = (p->age_grp == 0) ? 100 : ((p->age_grp == 1) ? 500 : 2000);
                if (p->hosp_rej) p_die *= 2;
                if (rand_int(rh, 1, BP) <= p_die) {
                    s->branches[BR_DEATH]++; s->daily_deaths_age[day][p->age_grp]++;
                    sim_set_state(s, i, ST_D, day, 0);
                } else {
                    s->branches[BR_RECOVER]++;
                    sim_set_state(s, i, ST_R, day, (int32_t)s->cfg.waning_min_days);
                }
            }
        } else if (p->status == ST_ICU && elapsed >= p->dur_days) {
            if (!p->hosp_rej && s->icu_occ > 0) s->icu_occ--;
            int64_t p_die = (p->age_grp == 0) ? 1500 : ((p->age_grp == 1) ? 2800 : 4800);
            if (p->hosp_rej) p_die *= 2;
            if (p_die > BP) p_die = BP;
            if (rand_int(rh, 1, BP) <= p_die) {
                s->branches[BR_DEATH]++; s->daily_deaths_age[day][p->age_grp]++;
                sim_set_state(s, i, ST_D, day, 0);
            } else {
                s->branches[BR_RECOVER]++;
                sim_set_state(s, i, ST_R, day, (int32_t)s->cfg.waning_min_days);
            }
        } else if (p->status == ST_R && elapsed >= p->dur_days) {
            if (rand_int(rd, 1, BP) <= s->cfg.waning_bp) {
                s->branches[BR_WANING]++;
                sim_set_state(s, i, ST_S, day, 0);
            }
        }
    }

    /* 2. Apply synchronous new infections from contact step */
    for (int32_t c = 0; c < s->n_cand; ++c) {
        int32_t tgt = s->cand_tgt[c], src = s->cand_src[c], lay = s->cand_lay[c], var = s->cand_var[c];
        if (s->pop[tgt].status != ST_S) continue;
        if (s->pop[tgt].ever_inf) s->total_reinf++;
        s->pop[tgt].ever_inf = 1;
        s->pop[tgt].variant = var;
        s->pop[tgt].inf_day = day;
        s->pop[src].offspring++;
        s->layer_inf[lay]++;
        s->branches[BR_INF_HH + lay]++;
        s->daily_inc_age[day][s->pop[tgt].age_grp]++;
        s->age_inf[s->pop[tgt].age_grp]++;
        s->daily_var[day][var]++;
        s->total_trans++;
        int32_t gi = day - s->pop[src].inf_day;
        if (gi < 1) gi = 1;
        push_trans(s, (TransRec){day, src, tgt, lay, var, gi});
        s->event_hash = fnv1a_i64(s->event_hash, (int64_t)day * 1000000LL + tgt);
        sim_set_state(s, tgt, ST_E, day, (int32_t)s->cfg.latent_days);
    }
}
