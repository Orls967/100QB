#include "common.h"
#include <stdlib.h>

static void add_cand(Sim *s, int32_t tgt, int32_t src, int32_t lay, int32_t var) {
    if (s->n_cand == s->cap_cand) {
        s->cap_cand = s->cap_cand ? s->cap_cand * 2 : 1024;
        s->cand_tgt = (int32_t *)realloc(s->cand_tgt, (size_t)s->cap_cand * sizeof(int32_t));
        s->cand_src = (int32_t *)realloc(s->cand_src, (size_t)s->cap_cand * sizeof(int32_t));
        s->cand_lay = (int32_t *)realloc(s->cand_lay, (size_t)s->cap_cand * sizeof(int32_t));
        s->cand_var = (int32_t *)realloc(s->cand_var, (size_t)s->cap_cand * sizeof(int32_t));
    }
    s->cand_tgt[s->n_cand] = tgt; s->cand_src[s->n_cand] = src;
    s->cand_lay[s->n_cand] = lay; s->cand_var[s->n_cand] = var;
    s->n_cand++;
}

static int64_t vax_eff_bp(const Sim *s, const Person *p, int32_t day) {
    if (p->vax_doses <= 0 || p->vax_day < 0) return 0;
    if (s->rp.vax_perfect) return BP;
    int64_t base = (p->vax_doses == 1) ? s->cfg.vax_eff1_bp : s->cfg.vax_eff2_bp;
    int64_t elapsed = day - p->vax_day;
    if (elapsed <= 30) return base;
    int64_t waned = base - (elapsed - 30) * 35;
    return waned > 2000 ? waned : 2000;
}

static void try_infect(Sim *s, int32_t src, int32_t tgt, int32_t lay, int32_t day, int64_t beta_eff) {
    Person *pt = &s->pop[tgt];
    if (pt->status != ST_S) return;
    if (s->seen_day[tgt] == day) return;
    if (pt->quar_until >= day) { s->branches[BR_FAIL_QUAR]++; return; }
    int64_t p_bp = beta_eff;
    if (s->rp.interv_on) {
        int64_t m_day = s->rp.mask_ovr >= 0 ? s->rp.mask_ovr : s->cfg.mask_day;
        if (day >= m_day) {
            int64_t red = (s->cfg.mask_eff_bp * s->pop[src].compliance) / BP;
            int64_t next_bp = (p_bp * (BP - red)) / BP;
            if (next_bp < p_bp) s->branches[BR_FAIL_MASK]++;
            p_bp = next_bp;
        }
    }
    int64_t ve = vax_eff_bp(s, pt, day);
    if (ve > 0) {
        p_bp = (p_bp * (BP - ve)) / BP;
        if (p_bp <= 0) { s->branches[BR_FAIL_VAX]++; return; }
    }
    if (rand_int(&s->rng[1], 1, BP) <= p_bp) {
        s->seen_day[tgt] = day;
        add_cand(s, tgt, src, lay, s->pop[src].variant);
    } else if (ve > 0) {
        s->branches[BR_FAIL_VAX]++;
    }
}

void sim_contact_step(Sim *s, int32_t day) {
    s->n_cand = 0;
    if (s->rp.contact_scale_bp == 0) return;
    int32_t N = (int32_t)s->N;
    int64_t sc_day = s->rp.school_ovr >= 0 ? s->rp.school_ovr : s->cfg.school_close_day;
    int64_t lk_day = s->rp.lock_ovr >= 0 ? s->rp.lock_ovr : s->cfg.lockdown_day;
    int32_t school_closed = s->rp.interv_on && (day >= sc_day);
    int32_t locked = s->rp.interv_on && (day >= lk_day);

    for (int32_t i = 0; i < N; ++i) {
        Person *pi = &s->pop[i];
        if (pi->status != ST_IA && pi->status != ST_IS) continue;
        int64_t base_bp = s->rp.beta_bp;
        if (pi->status == ST_IA) base_bp = (base_bp * s->cfg.asymp_rel_bp) / BP;
        if (pi->variant == 1) base_bp = (base_bp * s->cfg.variant_mult_bp) / BP;
        if (base_bp > BP) base_bp = BP;

        /* Layer 0: Household */
        int32_t h = pi->hh_id;
        for (int32_t idx = s->hh_off[h]; idx < s->hh_off[h + 1]; ++idx) {
            int32_t j = s->hh_mem[idx];
            if (j != i) try_infect(s, i, j, 0, day, base_bp);
        }
        if (pi->quar_until >= day) { s->branches[BR_FAIL_QUAR]++; continue; }

        /* Layer 1: School / Workplace */
        if (pi->sw_id >= 0 && !(pi->age_grp == 0 && school_closed)) {
            int32_t sw = pi->sw_id;
            int32_t lo = s->sw_off[sw], hi = s->sw_off[sw + 1];
            int32_t sz = hi - lo;
            if (sz > 1) {
                int32_t n_sw_c = (pi->age_grp == 0) ? 6 : 5;
                if (locked) n_sw_c = (n_sw_c + 1) / 2;
                for (int32_t c = 0; c < n_sw_c; ++c) {
                    int32_t j = s->sw_mem[lo + (int32_t)rand_int(&s->rng[1], 0, sz - 1)];
                    if (j != i) try_infect(s, i, j, 1, day, base_bp);
                }
            }
        }

        /* Layer 2: Community */
        int32_t n_comm = pi->super_flag ? 12 : 3;
        if (locked) n_comm = (n_comm + 1) / 3;
        for (int32_t c = 0; c < n_comm; ++c) {
            int32_t j = (int32_t)rand_int(&s->rng[1], 0, N - 1);
            if (j != i) try_infect(s, i, j, 2, day, base_bp);
        }
    }
}
