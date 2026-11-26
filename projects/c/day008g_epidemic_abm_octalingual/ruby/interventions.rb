# frozen_string_literal: true

def sim_interv_step(s, day)
  n = s.n
  vd = s.rp.var_day_ovr >= 0 ? s.rp.var_day_ovr : s.cfg["variant_day"]
  if day == vd
    imported = 0
    n.times do
      break if imported >= 5

      pid = s.rng[3].rand_int(0, n - 1)
      if s.pop[pid].status == ST_S
        s.pop[pid].variant = 1
        s.pop[pid].ever_inf = 1
        s.pop[pid].inf_day = day
        sim_set_state(s, pid, ST_IS, day, 4)
        s.branches[20] += 1
        imported += 1
      end
    end
  end
  return if s.rp.interv_on == 0

  if day >= s.cfg["test_day"]
    tp = s.rp.test_ovr >= 0 ? s.rp.test_ovr : s.cfg["test_prob_bp"]
    cp = s.rp.compliance_ovr
    n.times do |i|
      p = s.pop[i]
      next unless p.status == ST_IS && p.quar_until < day

      eff_comp = cp >= 0 ? cp : p.compliance
      next unless s.rng[3].rand_int(1, BP) <= (tp * eff_comp) / BP

      s.daily_tests[day][0] += 1
      if s.rng[3].rand_int(1, BP) <= 8500
        s.daily_tests[day][1] += 1
        s.branches[15] += 1
        p.quar_until = day + 7
        h = p.hh_id
        (s.hh_off[h]...s.hh_off[h + 1]).each do |idx|
          j = s.hh_mem[idx]
          if j != i && s.pop[j].quar_until < day && s.rng[3].rand_int(1, BP) <= s.cfg["trace_prob_bp"]
            s.pop[j].quar_until = day + 5
            s.branches[17] += 1
          end
        end
      else
        s.branches[16] += 1
      end
    end
  end

  vday = s.rp.vax_day_ovr >= 0 ? s.rp.vax_day_ovr : s.cfg["vax_day"]
  if day >= vday
    cov_cap = s.rp.vax_cov_ovr >= 0 ? (n * s.rp.vax_cov_ovr) / BP : (n * 8000) / BP
    doses_today = [(n * s.cfg["vax_daily_bp"]) / BP, 1].max
    total_vaxed = s.pop.count { |p| p.vax_doses > 0 }
    if total_vaxed >= cov_cap
      s.branches[19] += 1
      return
    end
    given = 0
    [2, 1, 0].each do |grp|
      break if given >= doses_today

      n.times do |i|
        break if given >= doses_today

        p = s.pop[i]
        next if p.age_grp != grp || p.status == ST_D

        if p.vax_doses == 0 && total_vaxed < cov_cap
          p.vax_doses = 1
          p.vax_day = day
          given += 1
          total_vaxed += 1
          s.daily_vax[day][0] += 1
          s.branches[18] += 1
        elsif p.vax_doses == 1 && day - p.vax_day >= 21
          p.vax_doses = 2
          p.vax_day = day
          given += 1
          s.daily_vax[day][1] += 1
          s.branches[18] += 1
        end
      end
    end
    s.branches[19] += 1 if given == 0
  end
end
