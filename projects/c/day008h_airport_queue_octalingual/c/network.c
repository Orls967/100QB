#include "common.h"

static void push_stn_ev(Sim *s, int stn, int32_t slot, int32_t delta) {
    if (!s->rp.audit_on || slot >= s->audit_slots) return;
    if (s->n_stn_ev[stn] >= s->cap_stn_ev[stn]) {
        s->cap_stn_ev[stn] = s->cap_stn_ev[stn] ? s->cap_stn_ev[stn] * 2 : 1024;
        s->stn_ev[stn] = (StnEvent *)realloc(s->stn_ev[stn], (size_t)s->cap_stn_ev[stn] * sizeof(StnEvent));
    }
    s->stn_ev[stn][s->n_stn_ev[stn]++] = (StnEvent){slot, stn, delta};
}

static void record_wait_samp(Sim *s, int32_t w) {
    if (s->n_wait_samp >= 50000) return;
    if (s->n_wait_samp >= s->cap_wait_samp) {
        s->cap_wait_samp = s->cap_wait_samp ? s->cap_wait_samp * 2 : 4096;
        s->wait_samples = (int32_t *)realloc(s->wait_samples, (size_t)s->cap_wait_samp * sizeof(int32_t));
    }
    s->wait_samples[s->n_wait_samp++] = w;
}

static void record_soj_samp(Sim *s, int32_t soj) {
    if (s->n_soj_samp >= 50000) return;
    if (s->n_soj_samp >= s->cap_soj_samp) {
        s->cap_soj_samp = s->cap_soj_samp ? s->cap_soj_samp * 2 : 4096;
        s->sojourn_samples = (int32_t *)realloc(s->sojourn_samples, (size_t)s->cap_soj_samp * sizeof(int32_t));
    }
    s->sojourn_samples[s->n_soj_samp++] = soj;
}

Sim *sim_create(const Config *c, const RunParams *rp) {
    Sim *s = (Sim *)calloc(1, sizeof(Sim));
    s->cfg = *c; s->rp = *rp;
    for (int i = 0; i < 8; ++i) {
        uint64_t x = (uint64_t)c->master_seed ^ fnv1a_u32pair(200, (uint32_t)i);
        s->rng[i] = xoshiro_seed(splitmix64(&x));
    }
    for (int i = MAX_POOL - 1; i >= 0; --i) s->free_stk[s->free_top++] = i;
    int caps[N_STN] = {1000, 1000, 1000, 1000, (int)rp->belt_cap, 1000};
    int prios[N_STN] = {1, 1, 0, 0, 0, 0};
    for (int i = 0; i < N_STN; ++i) {
        int cm = (int)(rp->srv_ovr[i] > 0 ? rp->srv_ovr[i] : c->srv_cnt[i]);
        if (cm > MAX_SRV) cm = MAX_SRV;
        stn_init(&s->stn[i], cm, caps[i], prios[i]);
    }
    static const int prof[24] = {180,140,120,150,260,420,680,860,920,820,700,660,680,740,820,880,940,900,760,620,480,360,280,210};
    for (int h = 0; h < 24; ++h) {
        s->arr_prof[h] = prof[h];
        for (int st = 0; st < N_STN; ++st) {
            int cm = s->stn[st].c_max;
            int act = (h >= 6 && h <= 20) ? cm : (cm > 1 ? cm - 1 : 1);
            s->roster[h][st] = act;
        }
    }
    s->event_hash = 0xcbf29ce484222325ULL;
    s->audit_slots = rp->audit_on ? (int32_t)((rp->days < AUDIT_DAYS ? rp->days : AUDIT_DAYS) * SLOTS_PER_DAY) : 0;
    if (s->audit_slots > 0) {
        s->q_series = (int32_t *)calloc((size_t)s->audit_slots * N_STN, sizeof(int32_t));
        s->srv_series = (int32_t *)calloc((size_t)s->audit_slots * 3, sizeof(int32_t));
    }
    return s;
}

void sim_free(Sim *s) {
    if (!s) return;
    free(s->wait_samples); free(s->sojourn_samples);
    free(s->q_series); free(s->srv_series);
    for (int i = 0; i < N_STN; ++i) free(s->stn_ev[i]);
    free(s);
}

static int64_t eff_bp(const Sim *s, int64_t base_bp) {
    return s->rp.timescale_div2 ? (base_bp / 2) : base_bp;
}

static void enqueue_agent(Sim *s, int stn_id, int32_t aidx, int64_t slot, int to_head) {
    Station *st = &s->stn[stn_id];
    int cur_in_stn = st->len + st->busy_cnt;
    if (cur_in_stn >= st->cap) {
        s->branches[BR_ARR_REJ]++;
        st->block_cnt++;
        if (stn_id == STN_BELT) s->total_belt_block++;
        if (!s->pool[aidx].is_bag) {
            s->total_lost++;
            s->sys_sum_sojourn += (slot - s->pool[aidx].arr_slot);
        }
        pool_free(s, aidx);
        return;
    }
    s->pool[aidx].cur_stn = stn_id;
    s->pool[aidx].q_enter_slot = slot;
    st->arr_cnt++;
    if (to_head) stn_push_head(st, aidx);
    else stn_push_tail(st, aidx);
    push_stn_ev(s, stn_id, (int32_t)slot, +1);
}

void sim_step_slot(Sim *s, int64_t slot) {
    int32_t hr = (int32_t)((slot % SLOTS_PER_DAY) / SLOTS_PER_HOUR);
    int32_t completed_aidx[N_STN * MAX_SRV];
    int32_t completed_stn[N_STN * MAX_SRV];
    int n_comp = 0;

    /* Phase 1: Arrivals */
    int64_t arr_p = (s->arr_prof[hr] * s->rp.arr_scale_bp) / BP;
    if (s->rp.timescale_div2) arr_p /= 2;
    if (bernoulli(&s->rng[0], arr_p)) {
        int32_t aidx = pool_alloc(s);
        if (aidx >= 0) {
            Agent *ag = &s->pool[aidx];
            int64_t u_cls = rand_int(&s->rng[0], 1, BP);
            ag->cls = (u_cls <= s->cfg.biz_prob_bp) ? CLS_BIZ : ((u_cls <= s->cfg.biz_prob_bp + s->cfg.transit_prob_bp) ? CLS_TRANSIT : CLS_ECON);
            ag->dest = bernoulli(&s->rng[0], s->cfg.intl_prob_bp) ? DEST_INTL : DEST_DOM;
            ag->n_bags = (int32_t)rand_int(&s->rng[0], 0, 2);
            ag->arr_slot = slot;
            ag->flight_hr = hr;
            ag->dep_slot = slot + rand_int(&s->rng[0], 110, 220);
            s->total_arr++;
            s->hr_arr[hr]++;
            s->branches[BR_ARR_ACC]++;
            s->event_hash = fnv1a_i64(s->event_hash, slot * 100 + ag->cls * 10 + ag->dest);
            int start_stn = (ag->cls == CLS_TRANSIT) ? STN_SECURITY : STN_CHECKIN;
            if (ag->cls == CLS_TRANSIT) s->branches[BR_ROUTE_TRANSIT]++;
            if (ag->n_bags > 0 && ag->cls != CLS_TRANSIT) {
                for (int b = 0; b < ag->n_bags; ++b) {
                    int32_t bidx = pool_alloc(s);
                    if (bidx >= 0) {
                        Agent *bg = &s->pool[bidx];
                        bg->is_bag = 1; bg->cls = ag->cls; bg->dest = ag->dest;
                        bg->arr_slot = slot; bg->dep_slot = ag->dep_slot;
                        s->total_bags++;
                        enqueue_agent(s, STN_BELT, bidx, slot, 0);
                    }
                }
            }
            enqueue_agent(s, start_stn, aidx, slot, 0);
        }
    }

    /* Phase 2: Service completions for servers busy at start of slot */
    for (int st = 0; st < N_STN; ++st) {
        Station *sn = &s->stn[st];
        int64_t base_mu = (s->rp.mu_ovr[st] > 0) ? s->rp.mu_ovr[st] : s->cfg.mu_bp[st];
        for (int idx_k = 0; idx_k < sn->c_max; ++idx_k) {
            int k = s->rp.permute_srv ? ((idx_k + 1) % sn->c_max) : idx_k;
            Server *sv = &sn->srv[k];
            if (sv->status == SRV_BUSY && sv->agent_idx >= 0) {
                int64_t mu = sv->extra_sec ? s->cfg.security_extra_mu_bp : base_mu;
                if (bernoulli(&s->rng[1], eff_bp(s, mu))) {
                    completed_aidx[n_comp] = sv->agent_idx;
                    completed_stn[n_comp] = st;
                    n_comp++;
                    sv->agent_idx = -1;
                    sv->status = SRV_IDLE;
                    sv->extra_sec = 0;
                    sn->busy_cnt--;
                    sn->dep_cnt++;
                    s->branches[BR_SVC_DONE]++;
                    push_stn_ev(s, st, (int32_t)slot, -1);
                }
            }
        }
    }

    /* Phase 3: Routing & transitions */
    for (int i = 0; i < n_comp; ++i) {
        int32_t aidx = completed_aidx[i];
        int st = completed_stn[i];
        Agent *ag = &s->pool[aidx];
        if (ag->is_bag) {
            if (bernoulli(&s->rng[4], s->cfg.bag_miss_bp)) {
                s->total_bags_miss++; s->branches[BR_BAG_MISS]++;
            }
            pool_free(s, aidx);
            continue;
        }
        int next_st = -1;
        if (st == STN_CHECKIN) next_st = STN_SECURITY;
        else if (st == STN_SECURITY) {
            if (ag->dest == DEST_INTL) { next_st = STN_IMIGRASI; s->branches[BR_ROUTE_INTL]++; }
            else next_st = STN_GATE;
        } else if (st == STN_IMIGRASI) next_st = STN_GATE;
        else if (st == STN_GATE) {
            int64_t lead = s->rp.gate_close_lead;
            if (s->rp.timescale_div2) lead *= 2;
            if (slot >= ag->dep_slot - lead) {
                s->total_missed++; s->branches[BR_MISSED_FLT]++;
                if (bernoulli(&s->rng[2], s->cfg.rebook_prob_bp)) {
                    s->total_rebook++; s->branches[BR_REBOOK]++;
                }
                s->total_lost++;
                s->sys_sum_sojourn += (slot - ag->arr_slot);
                pool_free(s, aidx);
                continue;
            }
            if (s->rp.customs_on && ag->dest == DEST_INTL) next_st = STN_CUSTOMS;
        }
        if (next_st >= 0) {
            enqueue_agent(s, next_st, aidx, slot, 0);
        } else {
            s->total_done++;
            int64_t soj = slot - ag->arr_slot;
            s->sys_sum_sojourn += soj;
            s->cls_sojourn_sum[ag->cls] += soj;
            s->cls_sojourn_cnt[ag->cls]++;
            record_soj_samp(s, (int32_t)soj);
            pool_free(s, aidx);
        }
    }

    /* Phase 4: Server assignment <- queue */
    for (int st = 0; st < N_STN; ++st) {
        Station *sn = &s->stn[st];
        for (int idx_k = 0; idx_k < sn->c_max; ++idx_k) {
            int k = s->rp.permute_srv ? ((idx_k + 1) % sn->c_max) : idx_k;
            Server *sv = &sn->srv[k];
            if (sv->enabled && sv->status == SRV_IDLE && sn->len > 0) {
                int32_t aidx = stn_pop_next(sn, s->pool, s->rp.priority_on, &s->branches[BR_BIZ_PRIO]);
                if (aidx >= 0) {
                    int64_t w = slot - s->pool[aidx].q_enter_slot;
                    sn->sum_W += w; sn->wait_cnt++;
                    s->stn_wait_sum[st] += w; s->stn_wait_cnt[st]++;
                    if (!s->pool[aidx].is_bag) {
                        s->cls_wait_sum[s->pool[aidx].cls] += w;
                        s->cls_wait_cnt[s->pool[aidx].cls]++;
                        record_wait_samp(s, (int32_t)w);
                    }
                    sv->status = SRV_BUSY;
                    sv->agent_idx = aidx;
                    sv->extra_sec = (st == STN_SECURITY && bernoulli(&s->rng[5], s->cfg.security_extra_prob_bp)) ? 1 : 0;
                    if (sv->extra_sec) s->branches[BR_SEC_EXTRA]++;
                    sn->busy_cnt++;
                }
            }
        }
    }

    /* Phase 5: Abandonment (patience expired) */
    int64_t pat = s->rp.timescale_div2 ? (s->cfg.patience_slots * 2) : s->cfg.patience_slots;
    for (int st = 0; st < N_STN; ++st) {
        if (st == STN_BELT) continue;
        Station *sn = &s->stn[st];
        for (int pos = 0; pos < sn->len; ) {
            int32_t aidx = sn->q[(sn->head + pos) % MAX_Q];
            if (slot - s->pool[aidx].q_enter_slot >= pat) {
                stn_remove_at(sn, pos);
                sn->aban_cnt++; sn->dep_cnt++;
                s->total_aban++; s->total_lost++;
                s->branches[BR_ABANDON]++;
                s->sys_sum_sojourn += (slot - s->pool[aidx].arr_slot);
                push_stn_ev(s, st, (int32_t)slot, -1);
                pool_free(s, aidx);
            } else {
                pos++;
            }
        }
    }

    /* Phase 6: Breakdown & Repair */
    if (s->rp.breakdown_on && s->rp.break_prob_bp > 0) {
        int64_t bp = eff_bp(s, s->rp.break_prob_bp);
        int64_t rp_b = eff_bp(s, s->rp.repair_prob_bp);
        for (int st = 0; st < N_STN; ++st) {
            Station *sn = &s->stn[st];
            for (int idx_k = 0; idx_k < sn->c_max; ++idx_k) {
                int k = s->rp.permute_srv ? ((idx_k + 1) % sn->c_max) : idx_k;
                Server *sv = &sn->srv[k];
                if (sv->status == SRV_DOWN) {
                    if (bernoulli(&s->rng[3], rp_b)) { sv->status = SRV_IDLE; s->branches[BR_SRV_REPAIR]++; }
                } else if (sv->enabled && bernoulli(&s->rng[3], bp)) {
                    if (sv->status == SRV_BUSY && sv->agent_idx >= 0) {
                        stn_push_head(sn, sv->agent_idx);
                        sv->agent_idx = -1;
                        sn->busy_cnt--;
                    }
                    sv->status = SRV_DOWN;
                    s->total_breakdowns++; s->branches[BR_SRV_BREAK]++;
                }
            }
        }
    }

    /* Phase 7: Hourly staff roster */
    if (s->rp.roster_on && (slot % SLOTS_PER_HOUR) == 0) {
        for (int st = 0; st < N_STN; ++st) {
            Station *sn = &s->stn[st];
            int target = s->roster[hr][st];
            if (target < sn->c_active) s->branches[BR_ROSTER_DROP]++;
            sn->c_active = target;
            for (int idx_k = 0; idx_k < sn->c_max; ++idx_k) {
                int k = s->rp.permute_srv ? ((idx_k + 1) % sn->c_max) : idx_k;
                sn->srv[k].enabled = (idx_k < target) ? 1 : 0;
            }
        }
    }

    /* Phase 8: Recording & exact invariants */
    int32_t tot_idle = 0, tot_busy = 0, tot_down = 0;
    for (int st = 0; st < N_STN; ++st) {
        Station *sn = &s->stn[st];
        int in_st = sn->len + sn->busy_cnt;
        sn->sum_L += in_st;
        sn->busy_slots += sn->busy_cnt;
        if (st != STN_BELT) s->sys_sum_L += in_st;
        s->hr_util_sum[hr][st] += sn->busy_cnt;
        if (sn->busy_cnt > sn->c_max || in_st > sn->cap) s->inv_violations++;
        if (s->rp.audit_on && slot < s->audit_slots) {
            s->q_series[slot * N_STN + st] = in_st;
        }
        for (int k = 0; k < sn->c_max; ++k) {
            if (sn->srv[k].status == SRV_IDLE && sn->srv[k].enabled) tot_idle++;
            else if (sn->srv[k].status == SRV_BUSY) tot_busy++;
            else if (sn->srv[k].status == SRV_DOWN) tot_down++;
        }
    }
    if (s->rp.audit_on && slot < s->audit_slots) {
        s->srv_series[slot * 3 + 0] = tot_idle;
        s->srv_series[slot * 3 + 1] = tot_busy;
        s->srv_series[slot * 3 + 2] = tot_down;
    }
}
