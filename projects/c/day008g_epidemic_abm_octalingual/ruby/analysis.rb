# frozen_string_literal: true

def net_degrees(s)
  shh = 0.0
  ssw = 0.0
  scom = 0.0
  s.pop.each do |p|
    h = p.hh_id
    shh += (s.hh_off[h + 1] - s.hh_off[h] - 1)
    ssw += (p.age_grp == 0 ? 6.0 : (p.age_grp == 1 ? 5.0 : 0.0))
    scom += (p.super_flag != 0 ? 12.0 : 3.0)
  end
  nf = s.n.to_f
  [shh / nf, ssw / nf, scom / nf]
end

def write_s01_to_s09(c, out_dir)
  hrng = FNV_OFFSET
  pairs = [[0, 0], [1, 2], [2, 1]]
  exp_hv = [0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb]
  rng_ok = true
  3.times do |i|
    x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1])
    _, sd = splitmix64(x)
    st = Xoshiro256.new(sd)
    vals = Array.new(4) { st.next_u64 }
    hv = fnv1a_bytes(vals.pack("Q<4"))
    rng_ok = false if hv != exp_hv[i]
    hrng = fnv1a_i64(hrng, hv)
  end
  st0 = Sim.new(c, default_rp(200, 10, 450, 0), false)
  sim_run(st0)
  inv_ok = st0.inv_violations == 0

  File.write("#{out_dir}/s01_selftest.txt", format("ok,h_rng,inv_ok,tot_trans\n%d,%016x,%d,%d\n", (rng_ok && inv_ok) ? 1 : 0, hrng, inv_ok ? 1 : 0, st0.total_trans))
  File.write("#{out_dir}/s02_config.csv", "key,value\nmaster_seed,#{c['master_seed']}\nparity_N,#{c['parity_N']}\nparity_days,#{c['parity_days']}\nmain_days,#{c['main_days']}\nbeta_bp_init,#{c['beta_bp_init']}\n")

  sp = Sim.new(c, default_rp(1000, 1, 450, 0), false)
  File.write("#{out_dir}/s03_population.csv", "age_grp,count\n0,#{sp.age_tot[0]}\n1,#{sp.age_tot[1]}\n2,#{sp.age_tot[2]}\n")
  File.open("#{out_dir}/s03_households.csv", "w") { |f| f.puts("hh_id,size"); [sp.n_hh, 50].min.times { |h| f.puts("#{h},#{sp.hh_off[h + 1] - sp.hh_off[h]}") } }
  File.open("#{out_dir}/s03_schools_work.csv", "w") { |f| f.puts("sw_id,size"); [sp.n_sw, 50].min.times { |w| f.puts("#{w},#{sp.sw_off[w + 1] - sp.sw_off[w]}") } }
  chh, csw, ccom = net_degrees(sp)
  File.write("#{out_dir}/s04_network.csv", format("n_hh,n_sw,deg_hh,deg_sw,deg_comm\n%d,%d,%.4f,%.4f,%.4f\n", sp.n_hh, sp.n_sw, chh, csw, ccom))

  File.open("#{out_dir}/s05_contact_script.csv", "w") do |f|
    f.puts("case_id,vax,mask,quar,eff_bp")
    12.times do |k|
      vax = k % 3
      mask = (k / 3) % 2
      quar = (k / 6) % 2
      p = quar != 0 ? 0 : 450
      p = (p * (BP - c["mask_eff_bp"])) / BP if mask != 0
      p = (p * (BP - c["vax_eff1_bp"])) / BP if vax == 1
      p = (p * (BP - c["vax_eff2_bp"])) / BP if vax == 2
      f.puts("#{k},#{vax},#{mask},#{quar},#{p}")
    end
  end
  File.open("#{out_dir}/s06_progression_script.csv", "w") do |f|
    f.puts("case_id,age_grp,hosp_rej,p_hosp_bp,p_die_bp")
    12.times do |k|
      ag = k % 3
      rej = (k / 3) % 2
      ph = ag == 0 ? 120 : (ag == 1 ? 450 : 1800)
      pd = ag == 0 ? 100 : (ag == 1 ? 500 : 2000)
      pd *= 2 if rej != 0
      f.puts("#{k},#{ag},#{rej},#{ph},#{pd}")
    end
  end
  File.write("#{out_dir}/s07_interventions.csv", "name,start_day,param_bp\nschool_close,#{c['school_close_day']},10000\nlockdown,#{c['lockdown_day']},5000\nmask,#{c['mask_day']},#{c['mask_eff_bp']}\ntest,#{c['test_day']},#{c['test_prob_bp']}\nvax,#{c['vax_day']},#{c['vax_daily_bp']}\nvariant,#{c['variant_day']},#{c['variant_mult_bp']}\n")

  File.open("#{out_dir}/s08_calibration.csv", "w") do |f|
    f.puts("beta_bp,r0_analytic,attack_rate,calibrated")
    [300, 380, 450, 520, 600].each do |bv|
      sc = Sim.new(c, default_rp(1000, 120, bv, 0), false)
      sim_run(sc)
      ch, cs2, cc2 = net_degrees(sc)
      r0 = analytical_r0(c, bv, ch, cs2, cc2)
      ar = sc.total_trans.to_f / sc.n.to_f
      cal = (r0 >= 1.5 && r0 <= 3.5 && ar >= 0.20 && ar <= 0.80) ? 1 : 0
      f.puts(format("%d,%.4f,%.4f,%d", bv, r0, ar, cal))
    end
  end
  s9 = Sim.new(c, default_rp(500, 30, 450, 0), false)
  sim_run(s9)
  File.write("#{out_dir}/s09_baseline_small.csv", "N,days,total_trans,deaths,inv_ok\n500,30,#{s9.total_trans},#{s9.daily_cnt[29][ST_D]},#{s9.inv_violations == 0 ? 1 : 0}\n")
end

def write_daily_csv(s, out_dir, fname)
  File.open("#{out_dir}/#{fname}", "w") do |f|
    f.puts("day,S,E,Ia,Is,H,ICU,R,D")
    s.rp.days.times { |d| f.puts("#{d},#{s.daily_cnt[d].join(',')}") }
  end
end

def write_trans_csv(s, out_dir, fname)
  File.open("#{out_dir}/#{fname}", "w") do |f|
    f.puts("day,src_id,dst_id,layer,variant,gen_int")
    s.trans.each { |r| f.puts(r.join(",")) }
  end
end

def write_s10_to_s37(s_no, s_main, c, out_dir, beta_cal)
  d_cnt = s_main.rp.days
  write_daily_csv(s_no, out_dir, "daily_counts_nointerv.csv")
  write_trans_csv(s_no, out_dir, "transmissions_log_nointerv.csv")
  write_daily_csv(s_main, out_dir, "daily_counts.csv")
  write_trans_csv(s_main, out_dir, "transmissions_log.csv")

  File.open("#{out_dir}/incidence_by_age.csv", "w") { |f| f.puts("day,age_0_17,age_18_59,age_60p"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_inc_age[d].join(',')}") } }
  File.open("#{out_dir}/hospital_daily.csv", "w") { |f| f.puts("day,hosp_occ,icu_occ,hosp_rej_cum,icu_rej_cum"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_hosp[d].join(',')}") } }
  File.open("#{out_dir}/deaths_by_age.csv", "w") { |f| f.puts("day,d_0_17,d_18_59,d_60p"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_deaths_age[d].join(',')}") } }
  File.open("#{out_dir}/tests.csv", "w") { |f| f.puts("day,tests,positives"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_tests[d][0]},#{s_main.daily_tests[d][1]}") } }
  File.open("#{out_dir}/quarantine.csv", "w") { |f| f.puts("day,quarantined"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_quar[d]}") } }
  File.open("#{out_dir}/vaccination.csv", "w") { |f| f.puts("day,dose1,dose2"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_vax[d][0]},#{s_main.daily_vax[d][1]}") } }
  File.open("#{out_dir}/variant_share.csv", "w") { |f| f.puts("day,v0_inf,v1_inf"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_var[d][0]},#{s_main.daily_var[d][1]}") } }
  File.open("#{out_dir}/generation_interval.csv", "w") { |f| f.puts("idx,gen_int"); s_main.trans.first(2000).each_with_index { |r, i| f.puts("#{i},#{r[5]}") } }
  File.open("#{out_dir}/agents_final.csv", "w") { |f| f.puts("id,age_grp,status,vax_doses,offspring,ever_inf"); s_main.pop.first(2000).each_with_index { |p, i| f.puts("#{i},#{p.age_grp},#{p.status},#{p.vax_doses},#{p.offspring},#{p.ever_inf}") } }
  File.open("#{out_dir}/progression_log.csv", "w") { |f| f.puts("day,pid,from_st,to_st"); s_main.progs.each { |r| f.puts(r.join(",")) } }

  3.times { |a| File.open("#{out_dir}/daily_age_#{a}.csv", "w") { |f| f.puts("day,inc,deaths"); d_cnt.times { |d| f.puts("#{d},#{s_main.daily_inc_age[d][a]},#{s_main.daily_deaths_age[d][a]}") } } }
  %w[layer_hh.csv layer_sw.csv layer_comm.csv].each_with_index { |nm, l| File.write("#{out_dir}/#{nm}", "layer,infections\n#{l},#{s_main.layer_inf[l]}\n") }
  File.open("#{out_dir}/hospital_cap_nointerv.csv", "w") { |f| f.puts("day,hosp_occ,icu_occ"); d_cnt.times { |d| f.puts("#{d},#{s_no.daily_hosp[d][0]},#{s_no.daily_hosp[d][1]}") } }
  File.open("#{out_dir}/deaths_nointerv.csv", "w") { |f| f.puts("day,deaths"); d_cnt.times { |d| f.puts("#{d},#{s_no.daily_cnt[d][ST_D]}") } }
  oh_bins = [0] * 6
  s_no.pop.each { |p| oh_bins[[p.offspring, 5].min] += 1 if p.ever_inf != 0 }
  File.open("#{out_dir}/offspring_hist.csv", "w") { |f| f.puts("bin,count"); 6.times { |b| f.puts("#{b},#{oh_bins[b]}") } }
  File.open("#{out_dir}/rt_nointerv.csv", "w") { |f| f.puts("day,inc"); d_cnt.times { |d| f.puts("#{d},#{s_no.daily_inc_age[d].sum}") } }

  tot_viol = s_no.inv_violations + s_main.inv_violations
  File.write("#{out_dir}/s12_invariants.csv", "violations,ok\n#{tot_viol},#{tot_viol == 0 ? 1 : 0}\n")

  File.open("#{out_dir}/s13_rt.csv", "w") do |f|
    f.puts("day,rt")
    d_cnt.times do |d|
      num = (0..4).select { |k| d >= k }.sum { |k| s_main.daily_inc_age[d - k].sum }
      den = (4..8).select { |k| d >= k }.sum { |k| s_main.daily_inc_age[d - k].sum }
      f.puts(format("%d,%.4f", d, den > 0 ? num.to_f / den.to_f : 0.0))
    end
  end

  chh, csw, ccom = net_degrees(s_no)
  r0_ana = analytical_r0(c, beta_cal, chh, csw, ccom)
  i_early = (s_no.daily_inc_age[3].sum + s_no.daily_inc_age[4].sum + 1) / 2.0
  i_mid = (10...[16, d_cnt].min).sum { |d| s_no.daily_inc_age[d].sum.to_f } / 6.0
  r_grow = i_mid > i_early ? Math.log(i_mid / i_early) / 9.0 : 0.0
  File.write("#{out_dir}/s14_growth.csv", format("growth_r,doubling_days\n%.6f,%.4f\n", r_grow, r_grow > 0.0 ? Math.log(2.0) / r_grow : 0.0))
  r0_emp = 1.0 + r_grow * 4.5
  File.write("#{out_dir}/s15_r0.csv", format("r0_analytic,r0_empirical,ratio\n%.4f,%.4f,%.4f\n", r0_ana, r0_emp, r0_ana > 0.0 ? r0_emp / r0_ana : 0.0))
  z_ana = solve_final_size(r0_ana)
  herd = r0_ana > 1.0 ? 1.0 - 1.0 / r0_ana : 0.0
  File.write("#{out_dir}/s16_final_size.csv", format("z_analytic,herd_thresh,attack_nointerv,attack_interv\n%.4f,%.4f,%.4f,%.4f\n", z_ana, herd, s_no.total_trans.to_f / s_no.n, s_main.total_trans.to_f / s_main.n))

  File.write("#{out_dir}/s17_attack.csv", format(
    "category,count,frac\nHH,%d,%.4f\nSW,%d,%.4f\nCOMM,%d,%.4f\nAGE0,%d,%.4f\nAGE1,%d,%.4f\nAGE2,%d,%.4f\n",
    s_main.layer_inf[0], s_main.layer_inf[0].to_f / s_main.n,
    s_main.layer_inf[1], s_main.layer_inf[1].to_f / s_main.n,
    s_main.layer_inf[2], s_main.layer_inf[2].to_f / s_main.n,
    s_main.age_inf[0], s_main.age_inf[0].to_f / [s_main.age_tot[0], 1].max,
    s_main.age_inf[1], s_main.age_inf[1].to_f / [s_main.age_tot[1], 1].max,
    s_main.age_inf[2], s_main.age_inf[2].to_f / [s_main.age_tot[2], 1].max
  ))

  m1 = 0.0
  m2 = 0.0
  n_inf = 0
  s_no.pop.each do |p|
    next if p.ever_inf == 0

    o = p.offspring.to_f
    m1 += o
    m2 += o * o
    n_inf += 1
  end
  if n_inf > 0
    m1 /= n_inf
    m2 /= n_inf
  end
  var_o = m2 - m1 * m1
  k_disp = (var_o > m1 && m1 > 0.0) ? (m1 * m1) / (var_o - m1) : 10.0
  File.write("#{out_dir}/s18_offspring.csv", format("n_inf,mean_offspring,var_offspring,k_disp\n%d,%.4f,%.4f,%.4f\n", n_inf, m1, var_o, k_disp))

  pk_h_no = (0...d_cnt).map { |d| s_no.daily_hosp[d][0] }.max || 0
  pk_i_no = (0...d_cnt).map { |d| s_no.daily_hosp[d][1] }.max || 0
  pk_h_m = (0...d_cnt).map { |d| s_main.daily_hosp[d][0] }.max || 0
  pk_i_m = (0...d_cnt).map { |d| s_main.daily_hosp[d][1] }.max || 0
  File.write("#{out_dir}/s19_hospital.csv", "hosp_cap,icu_cap,peak_h_no,peak_i_no,rej_no,peak_h_main,peak_i_main,rej_main\n#{s_main.hosp_cap},#{s_main.icu_cap},#{pk_h_no},#{pk_i_no},#{s_no.branches[9] + s_no.branches[11]},#{pk_h_m},#{pk_i_m},#{s_main.branches[9] + s_main.branches[11]}\n")

  File.open("#{out_dir}/s20_ifr.csv", "w") do |f|
    f.puts("age_grp,infections,deaths,ifr")
    3.times do |a|
      dt = (0...d_cnt).sum { |d| s_main.daily_deaths_age[d][a] }
      f.puts(format("%d,%d,%d,%.6f", a, s_main.age_inf[a], dt, s_main.age_inf[a] > 0 ? dt.to_f / s_main.age_inf[a] : 0.0))
    end
  end
  t_tot = (0...d_cnt).sum { |d| s_main.daily_tests[d][0] }
  t_pos = (0...d_cnt).sum { |d| s_main.daily_tests[d][1] }
  File.write("#{out_dir}/s21_testing.csv", "tests,positives,quar_traces\n#{t_tot},#{t_pos},#{s_main.branches[17]}\n")

  vn = [0, 0, 0]
  vi = [0, 0, 0]
  s_main.pop.each do |p|
    vd = [p.vax_doses, 2].min
    vn[vd] += 1
    vi[vd] += 1 if p.ever_inf != 0
  end
  File.open("#{out_dir}/s22_vaccine.csv", "w") do |f|
    f.puts("doses,n_people,n_inf,attack_rate")
    3.times { |d| f.puts(format("%d,%d,%d,%.6f", d, vn[d], vi[d], vn[d] > 0 ? vi[d].to_f / vn[d] : 0.0)) }
  end

  v0 = (0...d_cnt).sum { |d| s_main.daily_var[d][0] }
  v1 = (0...d_cnt).sum { |d| s_main.daily_var[d][1] }
  File.write("#{out_dir}/s23_variant.csv", format("v0_total,v1_total,v1_share\n%d,%d,%.6f\n", v0, v1, v0 + v1 > 0 ? v1.to_f / (v0 + v1) : 0.0))
  File.write("#{out_dir}/s24_waning.csv", "waning_events,reinfections\n#{s_main.branches[14]},#{s_main.total_reinf}\n")

  r1 = default_rp(500, 30, beta_cal, 0); r1.contact_scale_bp = 0; nc1 = Sim.new(c, r1, false); sim_run(nc1)
  r2 = default_rp(500, 30, 0, 0); r2.contact_scale_bp = -1; nc2 = Sim.new(c, r2, false); sim_run(nc2)
  r3 = default_rp(500, 30, beta_cal, 1); r3.vax_perfect = 1; nc3 = Sim.new(c, r3, false); sim_run(nc3)
  File.write("#{out_dir}/s25_negative_controls.csv", "control,infections,ok\nzero_contact,#{nc1.total_trans},#{nc1.total_trans == 0 ? 1 : 0}\nzero_seed,#{nc2.total_trans},#{nc2.total_trans == 0 ? 1 : 0}\nperfect_vax,#{nc3.total_trans},#{nc3.total_trans == 0 ? 1 : 0}\n")

  nsw = 500
  dsw = 90
  File.open("#{out_dir}/s26_sweep_comp_test.csv", "w") do |f|
    f.puts("comp_bp,test_bp,rep,attack_rate")
    [2000, 6000, 9500].each { |cv| [1000, 5000, 9000].each { |tv| 3.times { |rp| c2 = c.dup; c2["master_seed"] += rp + 1; r = default_rp(nsw, dsw, beta_cal, 1); r.compliance_ovr = cv; r.test_ovr = tv; ss = Sim.new(c2, r, false); sim_run(ss); f.puts(format("%d,%d,%d,%.6f", cv, tv, rp, ss.total_trans.to_f / nsw)) } } }
  end
  File.open("#{out_dir}/s27_sweep_vax.csv", "w") do |f|
    f.puts("vax_day,vax_cov_bp,rep,attack_rate")
    [10, 45, 80].each { |vd| [2000, 5500, 9000].each { |vc| 3.times { |rp| c2 = c.dup; c2["master_seed"] += rp + 10; r = default_rp(nsw, dsw, beta_cal, 1); r.vax_day_ovr = vd; r.vax_cov_ovr = vc; ss = Sim.new(c2, r, false); sim_run(ss); f.puts(format("%d,%d,%d,%.6f", vd, vc, rp, ss.total_trans.to_f / nsw)) } } }
  end
  File.open("#{out_dir}/s28_sweep_school_lock.csv", "w") do |f|
    f.puts("school_day,lock_day,rep,attack_rate")
    [15, 200].each { |sc| [15, 45, 200].each { |lk| 3.times { |rp| c2 = c.dup; c2["master_seed"] += rp + 20; r = default_rp(nsw, dsw, beta_cal, 1); r.school_ovr = sc; r.lock_ovr = lk; ss = Sim.new(c2, r, false); sim_run(ss); f.puts(format("%d,%d,%d,%.6f", sc, lk, rp, ss.total_trans.to_f / nsw)) } } }
  end
  File.open("#{out_dir}/s29_sweep_variant_mask.csv", "w") do |f|
    f.puts("var_day,mask_day,rep,attack_rate")
    [15, 45, 200].each { |vr| [10, 40, 200].each { |mk| 3.times { |rp| c2 = c.dup; c2["master_seed"] += rp + 30; r = default_rp(nsw, dsw, beta_cal, 1); r.var_day_ovr = vr; r.mask_ovr = mk; ss = Sim.new(c2, r, false); sim_run(ss); f.puts(format("%d,%d,%d,%.6f", vr, mk, rp, ss.total_trans.to_f / nsw)) } } }
  end

  ens_s = [0.0] * 90
  ens_s2 = [0.0] * 90
  20.times do |rp|
    c2 = c.dup; c2["master_seed"] += rp + 100
    ss = Sim.new(c2, default_rp(nsw, dsw, beta_cal, 1), false)
    sim_run(ss)
    90.times { |d| inc = ss.daily_inc_age[d].sum.to_f; ens_s[d] += inc; ens_s2[d] += inc * inc }
  end
  File.open("#{out_dir}/s30_ensemble.csv", "w") do |f|
    f.puts("day,mean_inc,se_inc")
    90.times { |d| mn = ens_s[d] / 20.0; f.puts(format("%d,%.4f,%.4f", d, mn, Math.sqrt([(ens_s2[d] / 20.0 - mn * mn) / 19.0, 0.0].max))) }
  end

  bn = %w[INF_HH INF_SW INF_COMM FAIL_VAX FAIL_MASK FAIL_QUAR PROG_IA PROG_IS HOSP_ADM HOSP_REJ ICU_ADM ICU_REJ DEATH RECOVER WANING TEST_POS TEST_NEG QUAR_TRACE VAX_DOSE VAX_OUT VAR_IMPORT]
  File.open("#{out_dir}/s31_branches.csv", "w") { |f| f.puts("branch,count"); bn.each_with_index { |nm, i| f.puts("#{nm},#{s_main.branches[i]}") } }
  File.write("#{out_dir}/s32_vacuous_check.csv", "check,ok\nneg_controls,#{(nc1.total_trans == 0 && nc2.total_trans == 0 && nc3.total_trans == 0) ? 1 : 0}\nensemble_var,#{ens_s.count { |x| x > 0.0 } > 5 ? 1 : 0}\n")

  sp1 = Sim.new(c, default_rp(c["parity_N"], c["parity_days"], beta_cal, 1), false); sim_run(sp1)
  sp2 = Sim.new(c, default_rp(c["parity_N"], c["parity_days"], beta_cal, 1), false); sim_run(sp2)
  h1 = FNV_OFFSET; h2 = FNV_OFFSET
  sp1.n.times { |i| h1 = fnv1a_i64(h1, sp1.pop[i].status); h2 = fnv1a_i64(h2, sp2.pop[i].status); h1 = fnv1a_i64(h1, sp1.pop[i].offspring); h2 = fnv1a_i64(h2, sp2.pop[i].offspring) }
  File.write("#{out_dir}/s33_replay.csv", format("ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n", sp1.event_hash, sp2.event_hash, h1, h2, (sp1.event_hash == sp2.event_hash && h1 == h2) ? 1 : 0))
  File.write("#{out_dir}/s34_parity.csv", format("parity_N,parity_days,ev_hash,fin_hash,total_trans,deaths,inv_violations\n%d,%d,%016x,%016x,%d,%d,%d\n", c["parity_N"], c["parity_days"], sp1.event_hash, h1, sp1.total_trans, sp1.daily_cnt[c["parity_days"] - 1][ST_D], sp1.inv_violations))
  File.write("#{out_dir}/s36_sizes.csv", "N,days,n_trans,n_progs\n#{s_main.n},#{s_main.rp.days},#{s_main.trans.length},#{s_main.progs.length}\n")

  cnt_r = [s_main.n, 0, 0, 0, 0, 0, 0, 0]
  pi = 0
  max_diff = 0
  d_cnt.times do |d|
    while pi < s_main.progs.length && s_main.progs[pi][0] == d
      cnt_r[s_main.progs[pi][2]] -= 1
      cnt_r[s_main.progs[pi][3]] += 1
      pi += 1
    end
    8.times { |k| max_diff = [max_diff, (cnt_r[k] - s_main.daily_cnt[d][k]).abs].max }
  end
  File.write("#{out_dir}/s37_replay_counts.csv", "max_diff,match\n#{max_diff},#{max_diff == 0 ? 1 : 0}\n")
end
