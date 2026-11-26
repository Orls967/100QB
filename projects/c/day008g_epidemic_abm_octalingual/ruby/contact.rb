# frozen_string_literal: true

def vax_eff_bp(s, p, day)
  return 0 if p.vax_doses <= 0 || p.vax_day < 0
  return BP if s.rp.vax_perfect != 0

  base = p.vax_doses == 1 ? s.cfg["vax_eff1_bp"] : s.cfg["vax_eff2_bp"]
  elapsed = day - p.vax_day
  elapsed <= 30 ? base : [base - (elapsed - 30) * 35, 2000].max
end

def try_infect(s, src, tgt, lay, day, beta_eff)
  pt = s.pop[tgt]
  return if pt.status != ST_S || s.seen_day[tgt] == day
  if pt.quar_until >= day
    s.branches[5] += 1
    return
  end

  p_bp = beta_eff
  if s.rp.interv_on != 0
    m_day = s.rp.mask_ovr >= 0 ? s.rp.mask_ovr : s.cfg["mask_day"]
    if day >= m_day
      red = (s.cfg["mask_eff_bp"] * s.pop[src].compliance) / BP
      next_bp = (p_bp * (BP - red)) / BP
      s.branches[4] += 1 if next_bp < p_bp
      p_bp = next_bp
    end
  end

  ve = vax_eff_bp(s, pt, day)
  if ve > 0
    p_bp = (p_bp * (BP - ve)) / BP
    if p_bp <= 0
      s.branches[3] += 1
      return
    end
  end

  if s.rng[1].rand_int(1, BP) <= p_bp
    s.seen_day[tgt] = day
    s.cand_tgt << tgt
    s.cand_src << src
    s.cand_lay << lay
    s.cand_var << s.pop[src].variant
  elsif ve > 0
    s.branches[3] += 1
  end
end

def sim_contact_step(s, day)
  s.cand_tgt.clear
  s.cand_src.clear
  s.cand_lay.clear
  s.cand_var.clear
  return if s.rp.contact_scale_bp == 0

  n = s.n
  sc_day = s.rp.school_ovr >= 0 ? s.rp.school_ovr : s.cfg["school_close_day"]
  lk_day = s.rp.lock_ovr >= 0 ? s.rp.lock_ovr : s.cfg["lockdown_day"]
  school_closed = s.rp.interv_on != 0 && day >= sc_day
  locked = s.rp.interv_on != 0 && day >= lk_day

  n.times do |i|
    pi = s.pop[i]
    st = pi.status
    next if st != ST_IA && st != ST_IS

    base_bp = s.rp.beta_bp
    base_bp = (base_bp * s.cfg["asymp_rel_bp"]) / BP if st == ST_IA
    base_bp = (base_bp * s.cfg["variant_mult_bp"]) / BP if pi.variant == 1
    base_bp = BP if base_bp > BP

    h = pi.hh_id
    (s.hh_off[h]...s.hh_off[h + 1]).each do |idx|
      j = s.hh_mem[idx]
      try_infect(s, i, j, 0, day, base_bp) if j != i
    end
    if pi.quar_until >= day
      s.branches[5] += 1
      next
    end

    if pi.sw_id >= 0 && !(pi.age_grp == 0 && school_closed)
      sw = pi.sw_id
      lo = s.sw_off[sw]
      hi = s.sw_off[sw + 1]
      sz = hi - lo
      if sz > 1
        n_sw_c = pi.age_grp == 0 ? 6 : 5
        n_sw_c = (n_sw_c + 1) / 2 if locked
        n_sw_c.times do
          j = s.sw_mem[lo + s.rng[1].rand_int(0, sz - 1)]
          try_infect(s, i, j, 1, day, base_bp) if j != i
        end
      end
    end

    n_comm = pi.super_flag != 0 ? 12 : 3
    n_comm = (n_comm + 1) / 3 if locked
    n_comm.times do
      j = s.rng[1].rand_int(0, n - 1)
      try_infect(s, i, j, 2, day, base_bp) if j != i
    end
  end
end
