# frozen_string_literal: true

def sim_set_state(s, pid, new_st, day, dur)
  old_st = s.pop[pid].status
  s.pop[pid].status = new_st
  s.pop[pid].state_day = day
  s.pop[pid].dur_days = dur
  s.progs << [day, pid, old_st, new_st] if s.log_enabled
end

def sim_disease_step(s, day)
  n = s.n
  n.times do |i|
    p = s.pop[i]
    elapsed = day - p.state_day
    st = p.status
    if st == ST_E && elapsed >= p.dur_days
      dur = s.rng[2].rand_int(s.cfg["inf_days_min"], s.cfg["inf_days_max"])
      if s.rng[2].rand_int(1, BP) <= s.cfg["asymp_prob_bp"]
        s.branches[6] += 1
        sim_set_state(s, i, ST_IA, day, dur)
      else
        s.branches[7] += 1
        sim_set_state(s, i, ST_IS, day, dur)
      end
    elsif (st == ST_IA || st == ST_IS) && elapsed >= p.dur_days
      if st == ST_IS
        p_hosp = p.age_grp == 0 ? 120 : (p.age_grp == 1 ? 450 : 1800)
        p_hosp = (p_hosp * (BP - s.cfg["vax_hosp_bp"])) / BP if p.vax_doses > 0
        if s.rng[4].rand_int(1, BP) <= p_hosp
          if s.hosp_occ < s.hosp_cap
            s.hosp_occ += 1
            p.hosp_rej = 0
            s.branches[8] += 1
          else
            p.hosp_rej = 1
            s.branches[9] += 1
          end
          dur = s.rng[4].rand_int(4, 8)
          sim_set_state(s, i, ST_H, day, dur)
          next
        end
      end
      s.branches[13] += 1
      sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
    elsif st == ST_H && elapsed >= p.dur_days
      s.hosp_occ -= 1 if p.hosp_rej == 0 && s.hosp_occ > 0
      p_icu = p.age_grp == 2 ? 3500 : 1800
      if s.rng[4].rand_int(1, BP) <= p_icu
        if s.icu_occ < s.icu_cap
          s.icu_occ += 1
          p.hosp_rej = 0
          s.branches[10] += 1
        else
          p.hosp_rej = 1
          s.branches[11] += 1
        end
        dur = s.rng[4].rand_int(4, 7)
        sim_set_state(s, i, ST_ICU, day, dur)
      else
        p_die = p.age_grp == 0 ? 100 : (p.age_grp == 1 ? 500 : 2000)
        p_die *= 2 if p.hosp_rej != 0
        if s.rng[4].rand_int(1, BP) <= p_die
          s.branches[12] += 1
          s.daily_deaths_age[day][p.age_grp] += 1
          sim_set_state(s, i, ST_D, day, 0)
        else
          s.branches[13] += 1
          sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
        end
      end
    elsif st == ST_ICU && elapsed >= p.dur_days
      s.icu_occ -= 1 if p.hosp_rej == 0 && s.icu_occ > 0
      p_die = p.age_grp == 0 ? 1500 : (p.age_grp == 1 ? 2800 : 4800)
      p_die *= 2 if p.hosp_rej != 0
      p_die = BP if p_die > BP
      if s.rng[4].rand_int(1, BP) <= p_die
        s.branches[12] += 1
        s.daily_deaths_age[day][p.age_grp] += 1
        sim_set_state(s, i, ST_D, day, 0)
      else
        s.branches[13] += 1
        sim_set_state(s, i, ST_R, day, s.cfg["waning_min_days"])
      end
    elsif st == ST_R && elapsed >= p.dur_days
      if s.rng[2].rand_int(1, BP) <= s.cfg["waning_bp"]
        s.branches[14] += 1
        sim_set_state(s, i, ST_S, day, 0)
      end
    end
  end

  s.cand_tgt.each_index do |c|
    tgt = s.cand_tgt[c]
    src = s.cand_src[c]
    lay = s.cand_lay[c]
    var = s.cand_var[c]
    pt = s.pop[tgt]
    next if pt.status != ST_S

    s.total_reinf += 1 if pt.ever_inf != 0
    pt.ever_inf = 1
    pt.variant = var
    pt.inf_day = day
    s.pop[src].offspring += 1
    s.layer_inf[lay] += 1
    s.branches[lay] += 1
    ag = pt.age_grp
    s.daily_inc_age[day][ag] += 1
    s.age_inf[ag] += 1
    s.daily_var[day][var] += 1
    s.total_trans += 1
    gi = day - s.pop[src].inf_day
    gi = 1 if gi < 1
    s.trans << [day, src, tgt, lay, var, gi] if s.log_enabled
    s.event_hash = fnv1a_i64(s.event_hash, day * 1_000_000 + tgt)
    sim_set_state(s, tgt, ST_E, day, s.cfg["latent_days"])
  end
end
