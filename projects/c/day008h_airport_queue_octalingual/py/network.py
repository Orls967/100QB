from config import (
    AUDIT_DAYS,
    CLS_BIZ,
    CLS_ECON,
    CLS_TRANSIT,
    DEST_DOM,
    DEST_INTL,
    MAX_POOL,
    MAX_Q,
    MAX_SRV,
    N_STN,
    SLOTS_PER_DAY,
    SLOTS_PER_HOUR,
    SRV_BUSY,
    SRV_DOWN,
    SRV_IDLE,
    STN_BELT,
    STN_CHECKIN,
    STN_CUSTOMS,
    STN_GATE,
    STN_IMIGRASI,
    STN_SECURITY,
)
from passenger import Agent, pool_alloc, pool_free
from rng import BP, FNV_OFFSET, MASK64, Xoshiro256, fnv1a_i64, fnv1a_u32pair, splitmix64
from station import Station

BR_ARR_ACC = 0
BR_ARR_REJ = 1
BR_SVC_DONE = 2
BR_SEC_EXTRA = 3
BR_SRV_BREAK = 4
BR_SRV_REPAIR = 5
BR_ABANDON = 6
BR_MISSED_FLT = 7
BR_REBOOK = 8
BR_BAG_MISS = 9
BR_BELT_BLOCK = 10
BR_ROSTER_DROP = 11
BR_ROUTE_INTL = 12
BR_ROUTE_TRANSIT = 13
BR_BIZ_PRIO = 14


class Sim:
    def __init__(self, c: dict, rp) -> None:
        self.cfg = c
        self.rp = rp
        self.rng = []
        for i in range(8):
            x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(200, i)
            _, sd = splitmix64(x)
            self.rng.append(Xoshiro256(sd))
        self.pool = [Agent() for _ in range(MAX_POOL)]
        self.free_stk = list(range(MAX_POOL - 1, -1, -1))
        caps = [1000, 1000, 1000, 1000, rp.belt_cap, 1000]
        prios = [1, 1, 0, 0, 0, 0]
        self.stn = []
        for i in range(N_STN):
            cm = rp.srv_ovr[i] if rp.srv_ovr[i] > 0 else c["srv_cnt"][i]
            if cm > MAX_SRV:
                cm = MAX_SRV
            self.stn.append(Station(cm, caps[i], prios[i]))
        self.arr_prof = [
            180, 140, 120, 150, 260, 420, 680, 860, 920, 820, 700, 660,
            680, 740, 820, 880, 940, 900, 760, 620, 480, 360, 280, 210,
        ]
        self.roster = []
        for h in range(24):
            row = []
            for st in range(N_STN):
                cm = self.stn[st].c_max
                row.append(cm if (6 <= h <= 20) else (cm - 1 if cm > 1 else 1))
            self.roster.append(row)
        self.branches = [0] * 15
        self.total_arr = 0
        self.total_done = 0
        self.total_lost = 0
        self.total_missed = 0
        self.total_rebook = 0
        self.total_aban = 0
        self.total_bags = 0
        self.total_bags_miss = 0
        self.total_belt_block = 0
        self.total_breakdowns = 0
        self.sys_sum_l = 0
        self.sys_sum_sojourn = 0
        self.inv_violations = 0
        self.event_hash = FNV_OFFSET
        self.cls_wait_sum = [0, 0, 0]
        self.cls_wait_cnt = [0, 0, 0]
        self.cls_sojourn_sum = [0, 0, 0]
        self.cls_sojourn_cnt = [0, 0, 0]
        self.stn_wait_sum = [0] * N_STN
        self.stn_wait_cnt = [0] * N_STN
        self.hr_arr = [0] * 24
        self.hr_util_sum = [[0] * N_STN for _ in range(24)]
        self.wait_samples: list[int] = []
        self.sojourn_samples: list[int] = []
        self.audit_slots = min(rp.days, AUDIT_DAYS) * SLOTS_PER_DAY if rp.audit_on != 0 else 0
        self.q_series = [0] * (self.audit_slots * N_STN) if self.audit_slots > 0 else []
        self.srv_series = [0] * (self.audit_slots * 3) if self.audit_slots > 0 else []
        self.stn_ev: list[list[tuple[int, int, int]]] = [[] for _ in range(N_STN)]

    def eff_bp(self, base: int) -> int:
        return (base // 2) if self.rp.timescale_div2 != 0 else base

    def enqueue_agent(self, stn_id: int, aidx: int, slot: int, to_head: bool) -> None:
        st = self.stn[stn_id]
        if st.len + st.busy_cnt >= st.cap:
            self.branches[BR_ARR_REJ] += 1
            st.block_cnt += 1
            if stn_id == STN_BELT:
                self.total_belt_block += 1
            if self.pool[aidx].is_bag == 0:
                self.total_lost += 1
                self.sys_sum_sojourn += slot - self.pool[aidx].arr_slot
            pool_free(self, aidx)
            return
        self.pool[aidx].cur_stn = stn_id
        self.pool[aidx].q_enter_slot = slot
        st.arr_cnt += 1
        if to_head:
            st.push_head(aidx)
        else:
            st.push_tail(aidx)
        if self.rp.audit_on != 0 and slot < self.audit_slots:
            self.stn_ev[stn_id].append((slot, stn_id, 1))


def sim_step_slot(s: Sim, slot: int) -> None:
    hr = (slot % SLOTS_PER_DAY) // SLOTS_PER_HOUR
    comp: list[tuple[int, int]] = []

    arr_p = (s.arr_prof[hr] * s.rp.arr_scale_bp) // BP
    if s.rp.timescale_div2 != 0:
        arr_p //= 2
    if s.rng[0].bernoulli(arr_p) != 0:
        aidx = pool_alloc(s)
        if aidx >= 0:
            ag = s.pool[aidx]
            u_cls = s.rng[0].rand_int(1, BP)
            if u_cls <= s.cfg["biz_prob_bp"]:
                ag.cls = CLS_BIZ
            elif u_cls <= s.cfg["biz_prob_bp"] + s.cfg["transit_prob_bp"]:
                ag.cls = CLS_TRANSIT
            else:
                ag.cls = CLS_ECON
            ag.dest = DEST_INTL if s.rng[0].bernoulli(s.cfg["intl_prob_bp"]) != 0 else DEST_DOM
            ag.n_bags = s.rng[0].rand_int(0, 2)
            ag.arr_slot = slot
            ag.flight_hr = hr
            ag.dep_slot = slot + s.rng[0].rand_int(110, 220)
            s.total_arr += 1
            s.hr_arr[hr] += 1
            s.branches[BR_ARR_ACC] += 1
            s.event_hash = fnv1a_i64(s.event_hash, slot * 100 + ag.cls * 10 + ag.dest)
            if ag.cls == CLS_TRANSIT:
                s.branches[BR_ROUTE_TRANSIT] += 1
                start_stn = STN_SECURITY
            else:
                start_stn = STN_CHECKIN
            if ag.n_bags > 0 and ag.cls != CLS_TRANSIT:
                for _ in range(ag.n_bags):
                    bidx = pool_alloc(s)
                    if bidx >= 0:
                        bg = s.pool[bidx]
                        bg.is_bag = 1
                        bg.cls = ag.cls
                        bg.dest = ag.dest
                        bg.arr_slot = slot
                        bg.dep_slot = ag.dep_slot
                        s.total_bags += 1
                        s.enqueue_agent(STN_BELT, bidx, slot, False)
            s.enqueue_agent(start_stn, aidx, slot, False)

    for st in range(N_STN):
        sn = s.stn[st]
        base_mu = s.rp.mu_ovr[st] if s.rp.mu_ovr[st] > 0 else s.cfg["mu_bp"][st]
        for idx_k in range(sn.c_max):
            k = (idx_k + 1) % sn.c_max if s.rp.permute_srv != 0 else idx_k
            sv = sn.srv[k]
            if sv.status != SRV_BUSY or sv.agent_idx < 0:
                continue
            mu = s.cfg["security_extra_mu_bp"] if sv.extra_sec != 0 else base_mu
            if s.rng[1].bernoulli(s.eff_bp(mu)) != 0:
                comp.append((sv.agent_idx, st))
                sv.agent_idx = -1
                sv.status = SRV_IDLE
                sv.extra_sec = 0
                sn.busy_cnt -= 1
                sn.dep_cnt += 1
                s.branches[BR_SVC_DONE] += 1
                if s.rp.audit_on != 0 and slot < s.audit_slots:
                    s.stn_ev[st].append((slot, st, -1))

    for aidx, st in comp:
        ag = s.pool[aidx]
        if ag.is_bag != 0:
            if s.rng[4].bernoulli(s.cfg["bag_miss_bp"]) != 0:
                s.total_bags_miss += 1
                s.branches[BR_BAG_MISS] += 1
            pool_free(s, aidx)
            continue
        next_st = -1
        if st == STN_CHECKIN:
            next_st = STN_SECURITY
        elif st == STN_SECURITY:
            if ag.dest == DEST_INTL:
                next_st = STN_IMIGRASI
                s.branches[BR_ROUTE_INTL] += 1
            else:
                next_st = STN_GATE
        elif st == STN_IMIGRASI:
            next_st = STN_GATE
        elif st == STN_GATE:
            lead = s.rp.gate_close_lead * 2 if s.rp.timescale_div2 != 0 else s.rp.gate_close_lead
            if slot >= ag.dep_slot - lead:
                s.total_missed += 1
                s.branches[BR_MISSED_FLT] += 1
                if s.rng[2].bernoulli(s.cfg["rebook_prob_bp"]) != 0:
                    s.total_rebook += 1
                    s.branches[BR_REBOOK] += 1
                s.total_lost += 1
                s.sys_sum_sojourn += slot - ag.arr_slot
                pool_free(s, aidx)
                continue
            if s.rp.customs_on != 0 and ag.dest == DEST_INTL:
                next_st = STN_CUSTOMS
        if next_st >= 0:
            s.enqueue_agent(next_st, aidx, slot, False)
        else:
            s.total_done += 1
            soj = slot - ag.arr_slot
            s.sys_sum_sojourn += soj
            s.cls_sojourn_sum[ag.cls] += soj
            s.cls_sojourn_cnt[ag.cls] += 1
            if len(s.sojourn_samples) < 50000:
                s.sojourn_samples.append(soj)
            pool_free(s, aidx)

    for st in range(N_STN):
        sn = s.stn[st]
        for idx_k in range(sn.c_max):
            k = (idx_k + 1) % sn.c_max if s.rp.permute_srv != 0 else idx_k
            sv = sn.srv[k]
            if sv.enabled == 0 or sv.status != SRV_IDLE or sn.len <= 0:
                continue
            aidx, s.branches[BR_BIZ_PRIO] = sn.pop_next(s.pool, s.rp.priority_on, s.branches[BR_BIZ_PRIO])
            if aidx >= 0:
                w = slot - s.pool[aidx].q_enter_slot
                sn.sum_w += w
                sn.wait_cnt += 1
                s.stn_wait_sum[st] += w
                s.stn_wait_cnt[st] += 1
                if s.pool[aidx].is_bag == 0:
                    cl = s.pool[aidx].cls
                    s.cls_wait_sum[cl] += w
                    s.cls_wait_cnt[cl] += 1
                    if len(s.wait_samples) < 50000:
                        s.wait_samples.append(w)
                extra = 1 if (st == STN_SECURITY and s.rng[5].bernoulli(s.cfg["security_extra_prob_bp"]) != 0) else 0
                if extra != 0:
                    s.branches[BR_SEC_EXTRA] += 1
                sv.status = SRV_BUSY
                sv.agent_idx = aidx
                sv.extra_sec = extra
                sn.busy_cnt += 1

    pat = s.cfg["patience_slots"] * 2 if s.rp.timescale_div2 != 0 else s.cfg["patience_slots"]
    for st in range(N_STN):
        if st == STN_BELT:
            continue
        sn = s.stn[st]
        pos = 0
        while pos < sn.len:
            aidx = sn.q[(sn.head + pos) % MAX_Q]
            if slot - s.pool[aidx].q_enter_slot >= pat:
                sn.remove_at(pos)
                sn.aban_cnt += 1
                sn.dep_cnt += 1
                s.total_aban += 1
                s.total_lost += 1
                s.branches[BR_ABANDON] += 1
                s.sys_sum_sojourn += slot - s.pool[aidx].arr_slot
                if s.rp.audit_on != 0 and slot < s.audit_slots:
                    s.stn_ev[st].append((slot, st, -1))
                pool_free(s, aidx)
            else:
                pos += 1

    if s.rp.breakdown_on != 0 and s.rp.break_prob_bp > 0:
        bp_v = s.eff_bp(s.rp.break_prob_bp)
        rp_b = s.eff_bp(s.rp.repair_prob_bp)
        for st in range(N_STN):
            sn = s.stn[st]
            for idx_k in range(sn.c_max):
                k = (idx_k + 1) % sn.c_max if s.rp.permute_srv != 0 else idx_k
                sv = sn.srv[k]
                if sv.status == SRV_DOWN:
                    if s.rng[3].bernoulli(rp_b) != 0:
                        sv.status = SRV_IDLE
                        s.branches[BR_SRV_REPAIR] += 1
                elif sv.enabled != 0 and s.rng[3].bernoulli(bp_v) != 0:
                    if sv.status == SRV_BUSY and sv.agent_idx >= 0:
                        sn.push_head(sv.agent_idx)
                        sv.agent_idx = -1
                        sn.busy_cnt -= 1
                    sv.status = SRV_DOWN
                    s.total_breakdowns += 1
                    s.branches[BR_SRV_BREAK] += 1

    if s.rp.roster_on != 0 and (slot % SLOTS_PER_HOUR) == 0:
        for st in range(N_STN):
            sn = s.stn[st]
            target = s.roster[hr][st]
            if target < sn.c_active:
                s.branches[BR_ROSTER_DROP] += 1
            sn.c_active = target
            for idx_k in range(sn.c_max):
                k = (idx_k + 1) % sn.c_max if s.rp.permute_srv != 0 else idx_k
                sn.srv[k].enabled = 1 if idx_k < target else 0

    ti = 0
    tb = 0
    td = 0
    for st in range(N_STN):
        sn = s.stn[st]
        in_st = sn.len + sn.busy_cnt
        sn.sum_l += in_st
        sn.busy_slots += sn.busy_cnt
        if st != STN_BELT:
            s.sys_sum_l += in_st
        s.hr_util_sum[hr][st] += sn.busy_cnt
        if sn.busy_cnt > sn.c_max or in_st > sn.cap:
            s.inv_violations += 1
        if s.rp.audit_on != 0 and slot < s.audit_slots:
            s.q_series[slot * N_STN + st] = in_st
        for k in range(sn.c_max):
            sv = sn.srv[k]
            if sv.status == SRV_IDLE and sv.enabled != 0:
                ti += 1
            elif sv.status == SRV_BUSY:
                tb += 1
            elif sv.status == SRV_DOWN:
                td += 1
    if s.rp.audit_on != 0 and slot < s.audit_slots:
        s.srv_series[slot * 3] = ti
        s.srv_series[slot * 3 + 1] = tb
        s.srv_series[slot * 3 + 2] = td
