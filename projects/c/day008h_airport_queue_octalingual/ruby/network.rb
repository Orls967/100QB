# frozen_string_literal: true

BR_ARR_ACC = 0; BR_ARR_REJ = 1; BR_SVC_DONE = 2; BR_SEC_EXTRA = 3; BR_SRV_BREAK = 4
BR_SRV_REPAIR = 5; BR_ABANDON = 6; BR_MISSED_FLT = 7; BR_REBOOK = 8; BR_BAG_MISS = 9
BR_BELT_BLOCK = 10; BR_ROSTER_DROP = 11; BR_ROUTE_INTL = 12; BR_ROUTE_TRANSIT = 13; BR_BIZ_PRIO = 14

class Sim
  attr_accessor :cfg, :rp, :rng, :pool, :free_stk, :stn, :arr_prof, :roster, :branches,
                :total_arr, :total_done, :total_lost, :total_missed, :total_rebook, :total_aban,
                :total_bags, :total_bags_miss, :total_belt_block, :total_breakdowns,
                :sys_sum_l, :sys_sum_sojourn, :inv_violations, :event_hash,
                :cls_wait_sum, :cls_wait_cnt, :cls_sojourn_sum, :cls_sojourn_cnt,
                :stn_wait_sum, :stn_wait_cnt, :hr_arr, :hr_util_sum,
                :wait_samples, :sojourn_samples, :q_series, :srv_series, :stn_ev, :audit_slots

  def initialize(c, rp)
    @cfg = c
    @rp = rp
    @rng = Array.new(8) do |i|
      x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(200, i)
      _, sd = splitmix64(x)
      Xoshiro256.new(sd)
    end
    @pool = Array.new(MAX_POOL) { Agent.new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) }
    @free_stk = (0...MAX_POOL).to_a.reverse
    caps = [1000, 1000, 1000, 1000, rp.belt_cap, 1000]
    prios = [1, 1, 0, 0, 0, 0]
    @stn = Array.new(N_STN) do |i|
      cm = rp.srv_ovr[i] > 0 ? rp.srv_ovr[i] : c["srv_cnt"][i]
      cm = MAX_SRV if cm > MAX_SRV
      Station.new(cm, caps[i], prios[i])
    end
    @arr_prof = [180, 140, 120, 150, 260, 420, 680, 860, 920, 820, 700, 660, 680, 740, 820, 880, 940, 900, 760, 620, 480, 360, 280, 210]
    @roster = Array.new(24) do |h|
      Array.new(N_STN) do |st|
        cm = @stn[st].c_max
        (h >= 6 && h <= 20) ? cm : (cm > 1 ? cm - 1 : 1)
      end
    end
    @branches = [0] * 15
    @total_arr = 0; @total_done = 0; @total_lost = 0; @total_missed = 0; @total_rebook = 0; @total_aban = 0
    @total_bags = 0; @total_bags_miss = 0; @total_belt_block = 0; @total_breakdowns = 0
    @sys_sum_l = 0; @sys_sum_sojourn = 0; @inv_violations = 0; @event_hash = FNV_OFFSET
    @cls_wait_sum = [0, 0, 0]; @cls_wait_cnt = [0, 0, 0]; @cls_sojourn_sum = [0, 0, 0]; @cls_sojourn_cnt = [0, 0, 0]
    @stn_wait_sum = [0] * N_STN; @stn_wait_cnt = [0] * N_STN
    @hr_arr = [0] * 24; @hr_util_sum = Array.new(24) { [0] * N_STN }
    @wait_samples = []; @sojourn_samples = []
    @audit_slots = rp.audit_on != 0 ? ([rp.days, AUDIT_DAYS].min * SLOTS_PER_DAY) : 0
    @q_series = @audit_slots > 0 ? [0] * (@audit_slots * N_STN) : []
    @srv_series = @audit_slots > 0 ? [0] * (@audit_slots * 3) : []
    @stn_ev = Array.new(N_STN) { [] }
  end

  def eff_bp(base)
    @rp.timescale_div2 != 0 ? (base / 2) : base
  end

  def enqueue_agent(stn_id, aidx, slot, to_head)
    st = @stn[stn_id]
    if st.len + st.busy_cnt >= st.cap
      @branches[BR_ARR_REJ] += 1; st.block_cnt += 1
      @total_belt_block += 1 if stn_id == STN_BELT
      if @pool[aidx].is_bag == 0
        @total_lost += 1; @sys_sum_sojourn += slot - @pool[aidx].arr_slot
      end
      pool_free(self, aidx)
      return
    end
    @pool[aidx].cur_stn = stn_id; @pool[aidx].q_enter_slot = slot; st.arr_cnt += 1
    to_head ? st.push_head(aidx) : st.push_tail(aidx)
    @stn_ev[stn_id] << [slot, stn_id, 1] if @rp.audit_on != 0 && slot < @audit_slots
  end
end

def sim_step_slot(s, slot)
  hr = (slot % SLOTS_PER_DAY) / SLOTS_PER_HOUR
  comp = []

  arr_p = (s.arr_prof[hr] * s.rp.arr_scale_bp) / BP
  arr_p /= 2 if s.rp.timescale_div2 != 0
  if s.rng[0].bernoulli(arr_p) != 0
    aidx = pool_alloc(s)
    if aidx >= 0
      ag = s.pool[aidx]
      u_cls = s.rng[0].rand_int(1, BP)
      ag.cls = u_cls <= s.cfg["biz_prob_bp"] ? CLS_BIZ : (u_cls <= s.cfg["biz_prob_bp"] + s.cfg["transit_prob_bp"] ? CLS_TRANSIT : CLS_ECON)
      ag.dest = s.rng[0].bernoulli(s.cfg["intl_prob_bp"]) != 0 ? DEST_INTL : DEST_DOM
      ag.n_bags = s.rng[0].rand_int(0, 2)
      ag.arr_slot = slot; ag.flight_hr = hr
      ag.dep_slot = slot + s.rng[0].rand_int(110, 220)
      s.total_arr += 1; s.hr_arr[hr] += 1; s.branches[BR_ARR_ACC] += 1
      s.event_hash = fnv1a_i64(s.event_hash, slot * 100 + ag.cls * 10 + ag.dest)
      start_stn = ag.cls == CLS_TRANSIT ? (s.branches[BR_ROUTE_TRANSIT] += 1; STN_SECURITY) : STN_CHECKIN
      if ag.n_bags > 0 && ag.cls != CLS_TRANSIT
        ag.n_bags.times do
          bidx = pool_alloc(s)
          if bidx >= 0
            bg = s.pool[bidx]
            bg.is_bag = 1; bg.cls = ag.cls; bg.dest = ag.dest; bg.arr_slot = slot; bg.dep_slot = ag.dep_slot
            s.total_bags += 1; s.enqueue_agent(STN_BELT, bidx, slot, false)
          end
        end
      end
      s.enqueue_agent(start_stn, aidx, slot, false)
    end
  end

  N_STN.times do |st|
    sn = s.stn[st]
    base_mu = s.rp.mu_ovr[st] > 0 ? s.rp.mu_ovr[st] : s.cfg["mu_bp"][st]
    sn.c_max.times do |idx_k|
      k = s.rp.permute_srv != 0 ? ((idx_k + 1) % sn.c_max) : idx_k
      sv = sn.srv[k]
      next unless sv.status == SRV_BUSY && sv.agent_idx >= 0

      mu = sv.extra_sec != 0 ? s.cfg["security_extra_mu_bp"] : base_mu
      if s.rng[1].bernoulli(s.eff_bp(mu)) != 0
        comp << [sv.agent_idx, st]
        sv.agent_idx = -1; sv.status = SRV_IDLE; sv.extra_sec = 0
        sn.busy_cnt -= 1; sn.dep_cnt += 1; s.branches[BR_SVC_DONE] += 1
        s.stn_ev[st] << [slot, st, -1] if s.rp.audit_on != 0 && slot < s.audit_slots
      end
    end
  end

  comp.each do |aidx, st|
    ag = s.pool[aidx]
    if ag.is_bag != 0
      if s.rng[4].bernoulli(s.cfg["bag_miss_bp"]) != 0
        s.total_bags_miss += 1; s.branches[BR_BAG_MISS] += 1
      end
      pool_free(s, aidx); next
    end
    next_st = -1
    if st == STN_CHECKIN
      next_st = STN_SECURITY
    elsif st == STN_SECURITY
      if ag.dest == DEST_INTL; next_st = STN_IMIGRASI; s.branches[BR_ROUTE_INTL] += 1; else next_st = STN_GATE; end
    elsif st == STN_IMIGRASI
      next_st = STN_GATE
    elsif st == STN_GATE
      lead = s.rp.timescale_div2 != 0 ? s.rp.gate_close_lead * 2 : s.rp.gate_close_lead
      if slot >= ag.dep_slot - lead
        s.total_missed += 1; s.branches[BR_MISSED_FLT] += 1
        if s.rng[2].bernoulli(s.cfg["rebook_prob_bp"]) != 0; s.total_rebook += 1; s.branches[BR_REBOOK] += 1; end
        s.total_lost += 1; s.sys_sum_sojourn += slot - ag.arr_slot
        pool_free(s, aidx); next
      end
      next_st = STN_CUSTOMS if s.rp.customs_on != 0 && ag.dest == DEST_INTL
    end
    if next_st >= 0
      s.enqueue_agent(next_st, aidx, slot, false)
    else
      s.total_done += 1; soj = slot - ag.arr_slot; s.sys_sum_sojourn += soj
      s.cls_sojourn_sum[ag.cls] += soj; s.cls_sojourn_cnt[ag.cls] += 1
      s.sojourn_samples << soj if s.sojourn_samples.length < 50_000
      pool_free(s, aidx)
    end
  end

  N_STN.times do |st|
    sn = s.stn[st]
    sn.c_max.times do |idx_k|
      k = s.rp.permute_srv != 0 ? ((idx_k + 1) % sn.c_max) : idx_k
      sv = sn.srv[k]
      next unless sv.enabled != 0 && sv.status == SRV_IDLE && sn.len > 0

      aidx, s.branches[BR_BIZ_PRIO] = sn.pop_next(s.pool, s.rp.priority_on, s.branches[BR_BIZ_PRIO])
      if aidx >= 0
        w = slot - s.pool[aidx].q_enter_slot
        sn.sum_w += w; sn.wait_cnt += 1; s.stn_wait_sum[st] += w; s.stn_wait_cnt[st] += 1
        if s.pool[aidx].is_bag == 0
          cl = s.pool[aidx].cls; s.cls_wait_sum[cl] += w; s.cls_wait_cnt[cl] += 1
          s.wait_samples << w if s.wait_samples.length < 50_000
        end
        extra = (st == STN_SECURITY && s.rng[5].bernoulli(s.cfg["security_extra_prob_bp"]) != 0) ? 1 : 0
        s.branches[BR_SEC_EXTRA] += 1 if extra != 0
        sv.status = SRV_BUSY; sv.agent_idx = aidx; sv.extra_sec = extra; sn.busy_cnt += 1
      end
    end
  end

  pat = s.rp.timescale_div2 != 0 ? s.cfg["patience_slots"] * 2 : s.cfg["patience_slots"]
  N_STN.times do |st|
    next if st == STN_BELT

    sn = s.stn[st]
    pos = 0
    while pos < sn.len
      aidx = sn.q[(sn.head + pos) % MAX_Q]
      if slot - s.pool[aidx].q_enter_slot >= pat
        sn.remove_at(pos); sn.aban_cnt += 1; sn.dep_cnt += 1
        s.total_aban += 1; s.total_lost += 1; s.branches[BR_ABANDON] += 1
        s.sys_sum_sojourn += slot - s.pool[aidx].arr_slot
        s.stn_ev[st] << [slot, st, -1] if s.rp.audit_on != 0 && slot < s.audit_slots
        pool_free(s, aidx)
      else
        pos += 1
      end
    end
  end

  if s.rp.breakdown_on != 0 && s.rp.break_prob_bp > 0
    bp = s.eff_bp(s.rp.break_prob_bp); rp_b = s.eff_bp(s.rp.repair_prob_bp)
    N_STN.times do |st|
      sn = s.stn[st]
      sn.c_max.times do |idx_k|
        k = s.rp.permute_srv != 0 ? ((idx_k + 1) % sn.c_max) : idx_k
        sv = sn.srv[k]
        if sv.status == SRV_DOWN
          if s.rng[3].bernoulli(rp_b) != 0; sv.status = SRV_IDLE; s.branches[BR_SRV_REPAIR] += 1; end
        elsif sv.enabled != 0 && s.rng[3].bernoulli(bp) != 0
          if sv.status == SRV_BUSY && sv.agent_idx >= 0
            sn.push_head(sv.agent_idx); sv.agent_idx = -1; sn.busy_cnt -= 1
          end
          sv.status = SRV_DOWN; s.total_breakdowns += 1; s.branches[BR_SRV_BREAK] += 1
        end
      end
    end
  end

  if s.rp.roster_on != 0 && (slot % SLOTS_PER_HOUR) == 0
    N_STN.times do |st|
      sn = s.stn[st]; target = s.roster[hr][st]
      s.branches[BR_ROSTER_DROP] += 1 if target < sn.c_active
      sn.c_active = target
      sn.c_max.times do |idx_k|
        k = s.rp.permute_srv != 0 ? ((idx_k + 1) % sn.c_max) : idx_k
        sn.srv[k].enabled = idx_k < target ? 1 : 0
      end
    end
  end

  ti = 0; tb = 0; td = 0
  N_STN.times do |st|
    sn = s.stn[st]; in_st = sn.len + sn.busy_cnt
    sn.sum_l += in_st; sn.busy_slots += sn.busy_cnt
    s.sys_sum_l += in_st if st != STN_BELT
    s.hr_util_sum[hr][st] += sn.busy_cnt
    s.inv_violations += 1 if sn.busy_cnt > sn.c_max || in_st > sn.cap
    s.q_series[slot * N_STN + st] = in_st if s.rp.audit_on != 0 && slot < s.audit_slots
    sn.c_max.times do |k|
      if sn.srv[k].status == SRV_IDLE && sn.srv[k].enabled != 0; ti += 1
      elsif sn.srv[k].status == SRV_BUSY; tb += 1
      elsif sn.srv[k].status == SRV_DOWN; td += 1
      end
    end
  end
  if s.rp.audit_on != 0 && slot < s.audit_slots
    s.srv_series[slot * 3] = ti; s.srv_series[slot * 3 + 1] = tb; s.srv_series[slot * 3 + 2] = td
  end
end
