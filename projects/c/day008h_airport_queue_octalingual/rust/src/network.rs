use crate::config::*;
use crate::passenger::{pool_alloc, pool_free, Agent};
use crate::rng::*;
use crate::station::Station;

pub const BR_ARR_ACC: usize = 0; pub const BR_ARR_REJ: usize = 1; pub const BR_SVC_DONE: usize = 2; pub const BR_SEC_EXTRA: usize = 3;
pub const BR_SRV_BREAK: usize = 4; pub const BR_SRV_REPAIR: usize = 5; pub const BR_ABANDON: usize = 6; pub const BR_MISSED_FLT: usize = 7;
pub const BR_REBOOK: usize = 8; pub const BR_BAG_MISS: usize = 9; pub const BR_ROSTER_DROP: usize = 11;
pub const BR_ROUTE_INTL: usize = 12; pub const BR_ROUTE_TRANSIT: usize = 13; pub const BR_BIZ_PRIO: usize = 14;

pub struct Sim {
    pub cfg: Config, pub rp: RunParams, pub rng: Vec<Xoshiro256>,
    pub pool: Vec<Agent>, pub free_stk: Vec<i32>, pub stn: Vec<Station>,
    pub arr_prof: [i32; 24], pub roster: [[usize; N_STN]; 24], pub branches: [i64; 15],
    pub total_arr: i64, pub total_done: i64, pub total_lost: i64, pub total_missed: i64, pub total_rebook: i64, pub total_aban: i64,
    pub total_bags: i64, pub total_bags_miss: i64, pub total_belt_block: i64, pub total_breakdowns: i64,
    pub sys_sum_l: i64, pub sys_sum_sojourn: i64, pub inv_violations: i64, pub event_hash: u64,
    pub cls_wait_sum: [i64; 3], pub cls_wait_cnt: [i64; 3], pub cls_sojourn_sum: [i64; 3], pub cls_sojourn_cnt: [i64; 3],
    pub stn_wait_sum: [i64; N_STN], pub stn_wait_cnt: [i64; N_STN],
    pub hr_arr: [i64; 24], pub hr_util_sum: [[i64; N_STN]; 24],
    pub wait_samples: Vec<i32>, pub sojourn_samples: Vec<i32>,
    pub q_series: Vec<i32>, pub srv_series: Vec<i32>, pub stn_ev: [Vec<(i32, i32, i32)>; N_STN],
    pub audit_slots: i32,
}

impl Sim {
    pub fn new(c: &Config, rp: &RunParams) -> Self {
        let mut rng = Vec::with_capacity(8);
        for i in 0..8u32 { let mut x = (c.master_seed as u64) ^ fnv1a_u32pair(200, i); rng.push(Xoshiro256::seed(splitmix64(&mut x))); }
        let free_stk = (0..MAX_POOL as i32).rev().collect();
        let caps = [1000, 1000, 1000, 1000, rp.belt_cap as i32, 1000];
        let prios = [1, 1, 0, 0, 0, 0];
        let mut stn = Vec::with_capacity(N_STN);
        for i in 0..N_STN {
            let cm = (if rp.srv_ovr[i] > 0 { rp.srv_ovr[i] } else { c.srv_cnt[i] }).min(MAX_SRV as i64) as usize;
            stn.push(Station::new(cm, caps[i], prios[i]));
        }
        let prof = [180,140,120,150,260,420,680,860,920,820,700,660,680,740,820,880,940,900,760,620,480,360,280,210];
        let mut roster = [[1usize; N_STN]; 24];
        for h in 0..24 {
            for st in 0..N_STN {
                let cm = stn[st].c_max;
                roster[h][st] = if (6..=20).contains(&h) { cm } else if cm > 1 { cm - 1 } else { 1 };
            }
        }
        let audit_slots = if rp.audit_on != 0 { (rp.days.min(AUDIT_DAYS) * SLOTS_PER_DAY) as i32 } else { 0 };
        let q_series = if audit_slots > 0 { vec![0; audit_slots as usize * N_STN] } else { Vec::new() };
        let srv_series = if audit_slots > 0 { vec![0; audit_slots as usize * 3] } else { Vec::new() };
        Self {
            cfg: c.clone(), rp: rp.clone(), rng, pool: vec![Agent::default(); MAX_POOL], free_stk, stn,
            arr_prof: prof, roster, branches: [0; 15],
            total_arr: 0, total_done: 0, total_lost: 0, total_missed: 0, total_rebook: 0, total_aban: 0,
            total_bags: 0, total_bags_miss: 0, total_belt_block: 0, total_breakdowns: 0,
            sys_sum_l: 0, sys_sum_sojourn: 0, inv_violations: 0, event_hash: 0xcbf29ce484222325,
            cls_wait_sum: [0; 3], cls_wait_cnt: [0; 3], cls_sojourn_sum: [0; 3], cls_sojourn_cnt: [0; 3],
            stn_wait_sum: [0; N_STN], stn_wait_cnt: [0; N_STN], hr_arr: [0; 24], hr_util_sum: [[0; N_STN]; 24],
            wait_samples: Vec::new(), sojourn_samples: Vec::new(), q_series, srv_series,
            stn_ev: Default::default(), audit_slots,
        }
    }
    fn eff_bp(&self, base_bp: i64) -> i64 { if self.rp.timescale_div2 != 0 { base_bp / 2 } else { base_bp } }
    fn enqueue(&mut self, stn_id: usize, aidx: i32, slot: i64, to_head: bool) {
        let sn = &mut self.stn[stn_id];
        if sn.len as i32 + sn.busy_cnt >= sn.cap {
            self.branches[BR_ARR_REJ] += 1; sn.block_cnt += 1;
            if stn_id == STN_BELT { self.total_belt_block += 1; }
            if self.pool[aidx as usize].is_bag == 0 {
                self.total_lost += 1; self.sys_sum_sojourn += slot - self.pool[aidx as usize].arr_slot;
            }
            pool_free(&mut self.pool, &mut self.free_stk, aidx);
            return;
        }
        self.pool[aidx as usize].cur_stn = stn_id as i32; self.pool[aidx as usize].q_enter_slot = slot; sn.arr_cnt += 1;
        if to_head { sn.push_head(aidx); } else { sn.push_tail(aidx); }
        if self.rp.audit_on != 0 && (slot as i32) < self.audit_slots { self.stn_ev[stn_id].push((slot as i32, stn_id as i32, 1)); }
    }
}

pub fn sim_step_slot(s: &mut Sim, slot: i64) {
    let hr = ((slot % SLOTS_PER_DAY) / SLOTS_PER_HOUR) as usize;
    let mut comp = [(0i32, 0usize); N_STN * MAX_SRV];
    let mut n_comp = 0usize;

    let mut arr_p = (s.arr_prof[hr] as i64 * s.rp.arr_scale_bp) / BP;
    if s.rp.timescale_div2 != 0 { arr_p /= 2; }
    if s.rng[0].bernoulli(arr_p) != 0 {
        let aidx = pool_alloc(&mut s.pool, &mut s.free_stk);
        if aidx >= 0 {
            let u_cls = s.rng[0].rand_int(1, BP);
            let cls = if u_cls <= s.cfg.biz_prob_bp { CLS_BIZ } else if u_cls <= s.cfg.biz_prob_bp + s.cfg.transit_prob_bp { CLS_TRANSIT } else { CLS_ECON };
            let dest = if s.rng[0].bernoulli(s.cfg.intl_prob_bp) != 0 { DEST_INTL } else { DEST_DOM };
            let n_bags = s.rng[0].rand_int(0, 2) as i32;
            let dep_slot = slot + s.rng[0].rand_int(110, 220);
            { let ag = &mut s.pool[aidx as usize]; ag.cls = cls; ag.dest = dest; ag.n_bags = n_bags; ag.arr_slot = slot; ag.flight_hr = hr as i32; ag.dep_slot = dep_slot; }
            s.total_arr += 1; s.hr_arr[hr] += 1; s.branches[BR_ARR_ACC] += 1;
            s.event_hash = fnv1a_i64(s.event_hash, slot * 100 + cls as i64 * 10 + dest as i64);
            let start_stn = if cls == CLS_TRANSIT { s.branches[BR_ROUTE_TRANSIT] += 1; STN_SECURITY } else { STN_CHECKIN };
            if n_bags > 0 && cls != CLS_TRANSIT {
                for _ in 0..n_bags {
                    let bidx = pool_alloc(&mut s.pool, &mut s.free_stk);
                    if bidx >= 0 {
                        { let bg = &mut s.pool[bidx as usize]; bg.is_bag = 1; bg.cls = cls; bg.dest = dest; bg.arr_slot = slot; bg.dep_slot = dep_slot; }
                        s.total_bags += 1; s.enqueue(STN_BELT, bidx, slot, false);
                    }
                }
            }
            s.enqueue(start_stn, aidx, slot, false);
        }
    }

    for st in 0..N_STN {
        let base_mu = if s.rp.mu_ovr[st] > 0 { s.rp.mu_ovr[st] } else { s.cfg.mu_bp[st] };
        let c_max = s.stn[st].c_max;
        for idx_k in 0..c_max {
            let k = if s.rp.permute_srv != 0 { (idx_k + 1) % c_max } else { idx_k };
            if s.stn[st].srv[k].status == SRV_BUSY && s.stn[st].srv[k].agent_idx >= 0 {
                let mu = if s.stn[st].srv[k].extra_sec != 0 { s.cfg.security_extra_mu_bp } else { base_mu };
                let p_mu = s.eff_bp(mu);
                if s.rng[1].bernoulli(p_mu) != 0 {
                    comp[n_comp] = (s.stn[st].srv[k].agent_idx, st); n_comp += 1;
                    s.stn[st].srv[k].agent_idx = -1; s.stn[st].srv[k].status = SRV_IDLE; s.stn[st].srv[k].extra_sec = 0;
                    s.stn[st].busy_cnt -= 1; s.stn[st].dep_cnt += 1; s.branches[BR_SVC_DONE] += 1;
                    if s.rp.audit_on != 0 && (slot as i32) < s.audit_slots { s.stn_ev[st].push((slot as i32, st as i32, -1)); }
                }
            }
        }
    }

    for i in 0..n_comp {
        let (aidx, st) = comp[i];
        let is_bag = s.pool[aidx as usize].is_bag; let dest = s.pool[aidx as usize].dest;
        let dep_slot = s.pool[aidx as usize].dep_slot; let arr_slot = s.pool[aidx as usize].arr_slot; let cls = s.pool[aidx as usize].cls;
        if is_bag != 0 {
            if s.rng[4].bernoulli(s.cfg.bag_miss_bp) != 0 { s.total_bags_miss += 1; s.branches[BR_BAG_MISS] += 1; }
            pool_free(&mut s.pool, &mut s.free_stk, aidx); continue;
        }
        let mut next_st: i32 = -1;
        if st == STN_CHECKIN { next_st = STN_SECURITY as i32; }
        else if st == STN_SECURITY { if dest == DEST_INTL { next_st = STN_IMIGRASI as i32; s.branches[BR_ROUTE_INTL] += 1; } else { next_st = STN_GATE as i32; } }
        else if st == STN_IMIGRASI { next_st = STN_GATE as i32; }
        else if st == STN_GATE {
            let lead = if s.rp.timescale_div2 != 0 { s.rp.gate_close_lead * 2 } else { s.rp.gate_close_lead };
            if slot >= dep_slot - lead {
                s.total_missed += 1; s.branches[BR_MISSED_FLT] += 1;
                if s.rng[2].bernoulli(s.cfg.rebook_prob_bp) != 0 { s.total_rebook += 1; s.branches[BR_REBOOK] += 1; }
                s.total_lost += 1; s.sys_sum_sojourn += slot - arr_slot;
                pool_free(&mut s.pool, &mut s.free_stk, aidx); continue;
            }
            if s.rp.customs_on != 0 && dest == DEST_INTL { next_st = STN_CUSTOMS as i32; }
        }
        if next_st >= 0 { s.enqueue(next_st as usize, aidx, slot, false); }
        else {
            s.total_done += 1; let soj = slot - arr_slot; s.sys_sum_sojourn += soj;
            s.cls_sojourn_sum[cls as usize] += soj; s.cls_sojourn_cnt[cls as usize] += 1;
            if s.sojourn_samples.len() < 50000 { s.sojourn_samples.push(soj as i32); }
            pool_free(&mut s.pool, &mut s.free_stk, aidx);
        }
    }

    for st in 0..N_STN {
        let c_max = s.stn[st].c_max;
        for idx_k in 0..c_max {
            let k = if s.rp.permute_srv != 0 { (idx_k + 1) % c_max } else { idx_k };
            if s.stn[st].srv[k].enabled != 0 && s.stn[st].srv[k].status == SRV_IDLE && s.stn[st].len > 0 {
                let aidx = s.stn[st].pop_next(&s.pool, s.rp.priority_on, &mut s.branches[BR_BIZ_PRIO]);
                if aidx >= 0 {
                    let w = slot - s.pool[aidx as usize].q_enter_slot;
                    s.stn[st].sum_w += w; s.stn[st].wait_cnt += 1; s.stn_wait_sum[st] += w; s.stn_wait_cnt[st] += 1;
                    if s.pool[aidx as usize].is_bag == 0 {
                        let cl = s.pool[aidx as usize].cls as usize;
                        s.cls_wait_sum[cl] += w; s.cls_wait_cnt[cl] += 1;
                        if s.wait_samples.len() < 50000 { s.wait_samples.push(w as i32); }
                    }
                    let extra = if st == STN_SECURITY && s.rng[5].bernoulli(s.cfg.security_extra_prob_bp) != 0 { 1 } else { 0 };
                    if extra != 0 { s.branches[BR_SEC_EXTRA] += 1; }
                    s.stn[st].srv[k].status = SRV_BUSY; s.stn[st].srv[k].agent_idx = aidx; s.stn[st].srv[k].extra_sec = extra;
                    s.stn[st].busy_cnt += 1;
                }
            }
        }
    }

    let pat = if s.rp.timescale_div2 != 0 { s.cfg.patience_slots * 2 } else { s.cfg.patience_slots };
    for st in 0..N_STN {
        if st == STN_BELT { continue; }
        let mut pos = 0usize;
        while pos < s.stn[st].len {
            let aidx = s.stn[st].q[(s.stn[st].head + pos) % MAX_Q];
            if slot - s.pool[aidx as usize].q_enter_slot >= pat {
                s.stn[st].remove_at(pos); s.stn[st].aban_cnt += 1; s.stn[st].dep_cnt += 1;
                s.total_aban += 1; s.total_lost += 1; s.branches[BR_ABANDON] += 1;
                s.sys_sum_sojourn += slot - s.pool[aidx as usize].arr_slot;
                if s.rp.audit_on != 0 && (slot as i32) < s.audit_slots { s.stn_ev[st].push((slot as i32, st as i32, -1)); }
                pool_free(&mut s.pool, &mut s.free_stk, aidx);
            } else { pos += 1; }
        }
    }

    if s.rp.breakdown_on != 0 && s.rp.break_prob_bp > 0 {
        let bp = s.eff_bp(s.rp.break_prob_bp); let rp_b = s.eff_bp(s.rp.repair_prob_bp);
        for st in 0..N_STN {
            let c_max = s.stn[st].c_max;
            for idx_k in 0..c_max {
                let k = if s.rp.permute_srv != 0 { (idx_k + 1) % c_max } else { idx_k };
                if s.stn[st].srv[k].status == SRV_DOWN {
                    if s.rng[3].bernoulli(rp_b) != 0 { s.stn[st].srv[k].status = SRV_IDLE; s.branches[BR_SRV_REPAIR] += 1; }
                } else if s.stn[st].srv[k].enabled != 0 && s.rng[3].bernoulli(bp) != 0 {
                    if s.stn[st].srv[k].status == SRV_BUSY && s.stn[st].srv[k].agent_idx >= 0 {
                        let aidx = s.stn[st].srv[k].agent_idx; s.stn[st].push_head(aidx);
                        s.stn[st].srv[k].agent_idx = -1; s.stn[st].busy_cnt -= 1;
                    }
                    s.stn[st].srv[k].status = SRV_DOWN; s.total_breakdowns += 1; s.branches[BR_SRV_BREAK] += 1;
                }
            }
        }
    }

    if s.rp.roster_on != 0 && (slot % SLOTS_PER_HOUR) == 0 {
        for st in 0..N_STN {
            let target = s.roster[hr][st];
            if target < s.stn[st].c_active { s.branches[BR_ROSTER_DROP] += 1; }
            s.stn[st].c_active = target;
            let c_max = s.stn[st].c_max;
            for idx_k in 0..c_max {
                let k = if s.rp.permute_srv != 0 { (idx_k + 1) % c_max } else { idx_k };
                s.stn[st].srv[k].enabled = if idx_k < target { 1 } else { 0 };
            }
        }
    }

    let (mut ti, mut tb, mut td) = (0i32, 0i32, 0i32);
    for st in 0..N_STN {
        let in_st = s.stn[st].len as i32 + s.stn[st].busy_cnt;
        s.stn[st].sum_l += in_st as i64; s.stn[st].busy_slots += s.stn[st].busy_cnt as i64;
        if st != STN_BELT { s.sys_sum_l += in_st as i64; }
        s.hr_util_sum[hr][st] += s.stn[st].busy_cnt as i64;
        if s.stn[st].busy_cnt > s.stn[st].c_max as i32 || in_st > s.stn[st].cap { s.inv_violations += 1; }
        if s.rp.audit_on != 0 && (slot as i32) < s.audit_slots { s.q_series[slot as usize * N_STN + st] = in_st; }
        for k in 0..s.stn[st].c_max {
            if s.stn[st].srv[k].status == SRV_IDLE && s.stn[st].srv[k].enabled != 0 { ti += 1; }
            else if s.stn[st].srv[k].status == SRV_BUSY { tb += 1; }
            else if s.stn[st].srv[k].status == SRV_DOWN { td += 1; }
        }
    }
    if s.rp.audit_on != 0 && (slot as i32) < s.audit_slots {
        s.srv_series[slot as usize * 3] = ti; s.srv_series[slot as usize * 3 + 1] = tb; s.srv_series[slot as usize * 3 + 2] = td;
    }
}
