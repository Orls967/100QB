# frozen_string_literal: true

def analytical_r0(c, beta_bp, c_hh, c_sw, c_comm)
  beta = beta_bp.to_f / BP.to_f
  d_inf = 0.5 * (c["inf_days_min"] + c["inf_days_max"])
  p_ia = c["asymp_prob_bp"].to_f / BP.to_f
  r_ia = c["asymp_rel_bp"].to_f / BP.to_f
  beta * (c_hh + c_sw + c_comm) * d_inf * ((1.0 - p_ia) + p_ia * r_ia)
end

def solve_final_size(r0)
  return 0.0 if r0 <= 1.0

  z = 0.5
  80.times { z = 1.0 - Math.exp(-r0 * z) }
  z
end

def check_invariants(s, day)
  cnt = [0] * 8
  s.pop.each { |p| cnt[p.status] += 1 }
  sum = 0
  8.times do |k|
    sum += cnt[k]
    s.daily_cnt[day][k] = cnt[k]
  end
  ok = (sum == s.n)
  ok = false if day > 0 && s.daily_cnt[day][ST_D] < s.daily_cnt[day - 1][ST_D]
  ok = false if s.hosp_occ > s.hosp_cap || s.icu_occ > s.icu_cap
  ok = false if s.log_enabled && s.trans.length != s.total_trans
  s.inv_violations += 1 unless ok
  ok
end

def sim_run(s)
  n = s.n
  n_seed = s.cfg["init_infected"]
  n_seed = 0 if s.rp.beta_bp == 0 && s.rp.contact_scale_bp < 0
  if s.rp.vax_perfect != 0
    s.pop.each do |p|
      p.vax_doses = 2
      p.vax_day = 0
    end
  end
  if n_seed > 0
    seeded = 0
    n.times do
      break if seeded >= n_seed

      pid = s.rng[2].rand_int(0, n - 1)
      if s.pop[pid].status == ST_S
        s.pop[pid].ever_inf = 1
        s.pop[pid].inf_day = 0
        sim_set_state(s, pid, ST_IS, 0, 4)
        seeded += 1
      end
    end
  end
  s.rp.days.times do |d|
    sim_interv_step(s, d)
    sim_contact_step(s, d)
    sim_disease_step(s, d)
    s.daily_hosp[d][0] = s.hosp_occ
    s.daily_hosp[d][1] = s.icu_occ
    s.daily_hosp[d][2] = s.branches[9]
    s.daily_hosp[d][3] = s.branches[11]
    s.daily_quar[d] = s.pop.count { |p| p.quar_until >= d }
    check_invariants(s, d)
  end
end
