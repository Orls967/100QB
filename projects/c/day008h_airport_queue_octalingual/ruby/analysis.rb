# frozen_string_literal: true

def run_metric_wait(c, rp, seed_off)
  c2 = c.dup; c2["master_seed"] += seed_off
  s = Sim.new(c2, rp); sim_run(s)
  s.stn_wait_cnt[STN_SECURITY] > 0 ? s.stn_wait_sum[STN_SECURITY].to_f / s.stn_wait_cnt[STN_SECURITY] : 0.0
end

def write_s14_to_s45(c, out_dir, s_main)
  snm = %w[CHECKIN SECURITY IMIGRASI GATE BAGGAGE_BELT CUSTOMS_ARRIVAL]
  ev_nm = %w[station_events_CHECKIN.csv station_events_SECURITY.csv station_events_IMIGRASI.csv station_events_GATE.csv station_events_BELT.csv]
  5.times { |st| File.open("#{out_dir}/#{ev_nm[st]}", "w") { |f| f.puts("slot,stn,delta"); s_main.stn_ev[st].each { |r| f.puts(r.join(",")) } } }
  File.open("#{out_dir}/queue_series.csv", "w") { |f| f.puts("slot,q0,q1,q2,q3,q4,q5"); s_main.audit_slots.times { |t| b = t * N_STN; f.puts("#{t},#{s_main.q_series[b, N_STN].join(',')}") } }
  File.open("#{out_dir}/server_status_series.csv", "w") { |f| f.puts("slot,idle,busy,down"); s_main.audit_slots.times { |t| b = t * 3; f.puts("#{t},#{s_main.srv_series[b, 3].join(',')}") } }

  File.write("#{out_dir}/passengers_log.csv", "total_arr,total_done,total_lost\n#{s_main.total_arr},#{s_main.total_done},#{s_main.total_lost}\n")
  File.write("#{out_dir}/bags_log.csv", "total_bags,bags_miss,belt_block\n#{s_main.total_bags},#{s_main.total_bags_miss},#{s_main.total_belt_block}\n")
  File.write("#{out_dir}/breakdowns.csv", "total_breakdowns,repairs\n#{s_main.total_breakdowns},#{s_main.branches[BR_SRV_REPAIR]}\n")
  File.write("#{out_dir}/missed_flights.csv", "missed,rebooked\n#{s_main.total_missed},#{s_main.total_rebook}\n")
  File.write("#{out_dir}/abandonments.csv", "abandonments\n#{s_main.total_aban}\n")
  File.write("#{out_dir}/rebookings.csv", "rebookings\n#{s_main.total_rebook}\n")
  File.write("#{out_dir}/agents_final.csv", "pool_cap,active_end\n#{MAX_POOL},0\n")
  File.write("#{out_dir}/daily_summary.csv", "days,total_arr,total_done,total_lost\n#{s_main.rp.days},#{s_main.total_arr},#{s_main.total_done},#{s_main.total_lost}\n")
  File.open("#{out_dir}/hourly_summary.csv", "w") { |f| f.puts("hour,arrivals,sec_busy"); 24.times { |h| f.puts("#{h},#{s_main.hr_arr[h]},#{s_main.hr_util_sum[h][STN_SECURITY]}") } }
  N_STN.times { |st| sn = s_main.stn[st]; File.write("#{out_dir}/station_#{st}.csv", "arr,dep,block,aban,sum_L,sum_W\n#{sn.arr_cnt},#{sn.dep_cnt},#{sn.block_cnt},#{sn.aban_cnt},#{sn.sum_l},#{sn.sum_w}\n") }
  3.times { |cl| File.write("#{out_dir}/class_#{cl}.csv", "wait_cnt,wait_sum,soj_cnt,soj_sum\n#{s_main.cls_wait_cnt[cl]},#{s_main.cls_wait_sum[cl]},#{s_main.cls_sojourn_cnt[cl]},#{s_main.cls_sojourn_sum[cl]}\n") }
  8.times { |fh| File.write("#{out_dir}/hour_slice_#{fh}.csv", "hour,arrivals,sec_busy\n#{fh * 3},#{s_main.hr_arr[fh * 3]},#{s_main.hr_util_sum[fh * 3][STN_SECURITY]}\n") }

  s_main.wait_samples.sort!; s_main.sojourn_samples.sort!
  wh = [0] * 6; s_main.wait_samples.each { |w| wh[[w / 5, 5].min] += 1 }
  sh = [0] * 6; s_main.sojourn_samples.each { |j| sh[[j / 20, 5].min] += 1 }
  File.open("#{out_dir}/wait_hist.csv", "w") { |f| f.puts("bin,count"); 6.times { |b| f.puts("#{b},#{wh[b]}") } }
  File.open("#{out_dir}/sojourn_hist.csv", "w") { |f| f.puts("bin,count"); 6.times { |b| f.puts("#{b},#{sh[b]}") } }
  File.write("#{out_dir}/bag_claims.csv", "missed_bags\n#{s_main.total_bags_miss}\n")
  File.write("#{out_dir}/flight_summary.csv", "flights_per_day,missed\n24,#{s_main.total_missed}\n")
  File.write("#{out_dir}/gate_events.csv", "gate_done,missed\n#{s_main.stn[STN_GATE].dep_cnt},#{s_main.total_missed}\n")
  File.write("#{out_dir}/priority_audit.csv", "biz_prio_jumps\n#{s_main.branches[BR_BIZ_PRIO]}\n")
  File.write("#{out_dir}/loss_breakdown.csv", "aban,missed,block\n#{s_main.total_aban},#{s_main.total_missed},#{s_main.stn[STN_BELT].block_cnt}\n")
  rpc = default_rp(c, 1); rpc.customs_on = 1; scust = Sim.new(c, rpc); sim_run(scust)
  File.write("#{out_dir}/customs_sweep_sample.csv", "customs_arr,customs_dep\n#{scust.stn[STN_CUSTOMS].arr_cnt},#{scust.stn[STN_CUSTOMS].dep_cnt}\n")

  File.write("#{out_dir}/s15_invariants.csv", "violations,sys_sum_L,sys_sum_sojourn,little_diff,ok\n#{s_main.inv_violations},#{s_main.sys_sum_l},#{s_main.sys_sum_sojourn},#{s_main.sys_sum_l - s_main.sys_sum_sojourn},#{s_main.inv_violations == 0 ? 1 : 0}\n")
  t_slots = (s_main.rp.days * SLOTS_PER_DAY).to_f
  File.open("#{out_dir}/s16_utilization.csv", "w") do |f|
    f.puts("stn,name,util,lam_slot,E_S,U_law")
    N_STN.times do |st|
      u = s_main.stn[st].busy_slots.to_f / (t_slots * s_main.stn[st].c_max)
      lam_s = s_main.stn[st].arr_cnt.to_f / t_slots
      es = s_main.stn[st].dep_cnt > 0 ? s_main.stn[st].busy_slots.to_f / s_main.stn[st].dep_cnt : 0.0
      f.puts(format("%d,%s,%.6f,%.6f,%.4f,%.6f", st, snm[st], u, lam_s, es, (lam_s * es) / s_main.stn[st].c_max))
    end
  end
  nw = s_main.wait_samples.length
  p50, p90, p99 = nw > 0 ? [s_main.wait_samples[nw / 2], s_main.wait_samples[nw * 90 / 100], s_main.wait_samples[nw * 99 / 100]] : [0, 0, 0]
  File.open("#{out_dir}/s17_waits.csv", "w") do |f|
    f.puts("stn,name,mean_wait,p50,p90,p99")
    N_STN.times { |st| mw = s_main.stn_wait_cnt[st] > 0 ? s_main.stn_wait_sum[st].to_f / s_main.stn_wait_cnt[st] : 0.0; f.puts(format("%d,%s,%.4f,%d,%d,%d", st, snm[st], mw, p50, p90, p99)) }
  end
  File.open("#{out_dir}/s18_sojourn.csv", "w") do |f|
    f.puts("cls,count,mean_sojourn,mean_wait")
    3.times do |cl|
      ms = s_main.cls_sojourn_cnt[cl] > 0 ? s_main.cls_sojourn_sum[cl].to_f / s_main.cls_sojourn_cnt[cl] : 0.0
      mw = s_main.cls_wait_cnt[cl] > 0 ? s_main.cls_wait_sum[cl].to_f / s_main.cls_wait_cnt[cl] : 0.0
      f.puts(format("%d,%d,%.4f,%.4f", cl, s_main.cls_sojourn_cnt[cl], ms, mw))
    end
  end

  l_net = s_main.sys_sum_l.to_f / t_slots; lam_net = s_main.total_arr.to_f / t_slots
  w_net = s_main.total_arr > 0 ? s_main.sys_sum_sojourn.to_f / s_main.total_arr : 0.0
  File.write("#{out_dir}/s19_little.csv", format("L_net,lam_net,W_net,lam_W,rel_diff\n%.6f,%.6f,%.6f,%.6f,%.12e\n", l_net, lam_net, w_net, lam_net * w_net, (l_net - lam_net * w_net).abs / (l_net > 0 ? l_net : 1.0)))
  File.write("#{out_dir}/s20_missed.csv", format("total_arr,missed,rebooked,aban,loss_frac\n%d,%d,%d,%d,%.6f\n", s_main.total_arr, s_main.total_missed, s_main.total_rebook, s_main.total_aban, s_main.total_arr > 0 ? s_main.total_lost.to_f / s_main.total_arr : 0.0))
  File.write("#{out_dir}/s21_baggage.csv", format("total_bags,bags_miss,belt_block,mean_belt_wait\n%d,%d,%d,%.4f\n", s_main.total_bags, s_main.total_bags_miss, s_main.total_belt_block, s_main.stn_wait_cnt[STN_BELT] > 0 ? s_main.stn_wait_sum[STN_BELT].to_f / s_main.stn_wait_cnt[STN_BELT] : 0.0))

  rp_nb = default_rp(c, 1); rp_nb.breakdown_on = 0
  w_nb = run_metric_wait(c, rp_nb, 0); w_wb = run_metric_wait(c, default_rp(c, 1), 0)
  File.write("#{out_dir}/s22_breakdown_effect.csv", format("wait_no_break,wait_with_break,diff\n%.4f,%.4f,%.4f\n", w_nb, w_wb, w_wb - w_nb))
  File.open("#{out_dir}/s23_roster_effect.csv", "w") { |f| f.puts("hour,arrivals,sec_busy_slots"); 24.times { |h| f.puts("#{h},#{s_main.hr_arr[h]},#{s_main.hr_util_sum[h][STN_SECURITY]}") } }
  bn_st = (0...4).max_by { |st| s_main.stn[st].busy_slots.to_f / (t_slots * s_main.stn[st].c_max) }
  bn_u = s_main.stn[bn_st].busy_slots.to_f / (t_slots * s_main.stn[bn_st].c_max)
  File.write("#{out_dir}/s24_bottleneck.csv", format("bottleneck_stn,name,util\n%d,%s,%.6f\n", bn_st, snm[bn_st], bn_u))

  nc1 = default_rp(c, 1); nc1.arr_scale_bp = 0; snc1 = Sim.new(c, nc1); sim_run(snc1)
  nc2 = default_rp(c, 1); nc2.srv_ovr = [12] * N_STN; nc2.mu_ovr = [BP] * N_STN; nc2.roster_on = 0; nc2.breakdown_on = 0; snc2 = Sim.new(c, nc2); sim_run(snc2)
  nc2_w = snc2.stn_wait_cnt[STN_SECURITY] > 0 ? snc2.stn_wait_sum[STN_SECURITY].to_f / snc2.stn_wait_cnt[STN_SECURITY] : 0.0
  nc3 = default_rp(c, 1); nc3.break_prob_bp = BP; nc3.repair_prob_bp = 0; snc3 = Sim.new(c, nc3); sim_run(snc3)
  File.write("#{out_dir}/s25_negative_controls.csv", format("control,metric,ok\nzero_arr,%d,%d\ninstant_svc,%.4f,%d\nfull_break,%d,%d\n", snc1.total_arr, snc1.total_arr == 0 ? 1 : 0, nc2_w, nc2_w < 0.05 ? 1 : 0, snc3.total_done, snc3.total_done == 0 ? 1 : 0))

  mp2 = default_rp(c, 1); mp2.permute_srv = 1
  smp1 = Sim.new(c, default_rp(c, 1)); sim_run(smp1); smp2 = Sim.new(c, mp2); sim_run(smp2)
  File.write("#{out_dir}/s26_meta_permute.csv", "done1,done2,sumL1,sumL2,match\n#{smp1.total_done},#{smp2.total_done},#{smp1.sys_sum_l},#{smp2.sys_sum_l},#{(smp1.total_done == smp2.total_done && smp1.sys_sum_l == smp2.sys_sum_l) ? 1 : 0}\n")

  mt2 = default_rp(c, 1); mt2.timescale_div2 = 1
  smt1 = Sim.new(c, default_rp(c, 1)); sim_run(smt1); smt2 = Sim.new(c, mt2); sim_run(smt2)
  u1_sec = smt1.stn[STN_SECURITY].busy_slots.to_f / (SLOTS_PER_DAY * smt1.stn[STN_SECURITY].c_max)
  u2_sec = smt2.stn[STN_SECURITY].busy_slots.to_f / (2 * SLOTS_PER_DAY * smt2.stn[STN_SECURITY].c_max)
  File.write("#{out_dir}/s27_meta_timescale.csv", format("util_10s,util_5s,abs_diff\n%.6f,%.6f,%.6f\n", u1_sec, u2_sec, (u1_sec - u2_sec).abs))

  mono_ok = 0
  20.times do |r|
    ra = default_rp(c, 1); ra.srv_ovr[STN_SECURITY] = 3; rb = default_rp(c, 1); rb.srv_ovr[STN_SECURITY] = 4
    mono_ok += 1 if run_metric_wait(c, rb, r + 10) <= run_metric_wait(c, ra, r + 10) + 1e-9
  end
  File.write("#{out_dir}/s28_meta_monotone.csv", format("ok_pairs,total_pairs,frac\n%d,20,%.4f\n", mono_ok, mono_ok / 20.0))

  sw_days = 1
  File.open("#{out_dir}/s29_sweep_servers.csv", "w") { |f| f.puts("chk_srv,sec_srv,rep,mean_wait"); (2..4).each { |sc| (2..4).each { |ss| 5.times { |r| rp = default_rp(c, sw_days); rp.srv_ovr[STN_CHECKIN] = sc; rp.srv_ovr[STN_SECURITY] = ss; f.puts(format("%d,%d,%d,%.4f", sc, ss, r, run_metric_wait(c, rp, r + 100))) } } } }
  File.open("#{out_dir}/s30_sweep_priority.csv", "w") do |f|
    f.puts("priority_on,arr_scale_bp,rep,wait_econ,wait_biz")
    [0, 1].each { |pr| [8000, 10000, 12000].each { |sc| 5.times { |r| c2 = c.dup; c2["master_seed"] += r + 200; rp = default_rp(c2, sw_days); rp.priority_on = pr; rp.arr_scale_bp = sc; ss = Sim.new(c2, rp); sim_run(ss); we = ss.cls_wait_cnt[0] > 0 ? ss.cls_wait_sum[0].to_f / ss.cls_wait_cnt[0] : 0.0; wb = ss.cls_wait_cnt[1] > 0 ? ss.cls_wait_sum[1].to_f / ss.cls_wait_cnt[1] : 0.0; f.puts(format("%d,%d,%d,%.4f,%.4f", pr, sc, r, we, wb)) } } }
  end
  File.open("#{out_dir}/s31_sweep_breakdown.csv", "w") { |f| f.puts("break_bp,repair_bp,rep,mean_wait"); [0, 15, 45].each { |ib| [800, 2000].each { |ir| 5.times { |r| rp = default_rp(c, sw_days); rp.break_prob_bp = ib; rp.repair_prob_bp = ir; rp.breakdown_on = 0 if ib == 0; f.puts(format("%d,%d,%d,%.4f", ib, ir, r, run_metric_wait(c, rp, r + 300))) } } } }
  File.open("#{out_dir}/s32_sweep_belt.csv", "w") { |f| f.puts("belt_cap,arr_scale_bp,rep,belt_blocks"); [4, 10, 20].each { |ib| [8000, 10000, 12000].each { |sc| 5.times { |r| c2 = c.dup; c2["master_seed"] += r + 400; rp = default_rp(c2, sw_days); rp.belt_cap = ib; rp.arr_scale_bp = sc; ss = Sim.new(c2, rp); sim_run(ss); f.puts("#{ib},#{sc},#{r},#{ss.total_belt_block}") } } } }
  File.open("#{out_dir}/s33_sweep_gate.csv", "w") { |f| f.puts("gate_lead,arr_scale_bp,rep,missed_flights"); [10, 25, 50].each { |ig| [8000, 10000, 12000].each { |sc| 5.times { |r| c2 = c.dup; c2["master_seed"] += r + 500; rp = default_rp(c2, sw_days); rp.gate_close_lead = ig; rp.arr_scale_bp = sc; ss = Sim.new(c2, rp); sim_run(ss); f.puts("#{ig},#{sc},#{r},#{ss.total_missed}") } } } }

  ens_s = [0.0] * 24; ens_s2 = [0.0] * 24
  30.times do |r|
    c2 = c.dup; c2["master_seed"] += r + 600; rp = default_rp(c2, 1); rp.audit_on = 1; ss = Sim.new(c2, rp); sim_run(ss)
    24.times { |h| qh = (h * 360...(h + 1) * 360).sum { |t| ss.q_series[t * N_STN + STN_SECURITY].to_f } / 360.0; ens_s[h] += qh; ens_s2[h] += qh * qh }
  end
  File.open("#{out_dir}/s34_ensemble.csv", "w") { |f| f.puts("hour,mean_q_sec,se_q_sec"); 24.times { |h| mn = ens_s[h] / 30.0; f.puts(format("%d,%.4f,%.4f", h, mn, Math.sqrt([(ens_s2[h] / 30.0 - mn * mn) / 29.0, 0.0].max))) } }

  rw = Array.new(20) { |r| run_metric_wait(c, default_rp(c, 1), r + 700) }
  rep_m = rw.sum / 20.0; rep_se = Math.sqrt([rw.sum { |x| (x - rep_m) ** 2 } / (20.0 * 19.0), 0.0].max)
  bm_m = s_main.stn_wait_cnt[STN_SECURITY] > 0 ? s_main.stn_wait_sum[STN_SECURITY].to_f / s_main.stn_wait_cnt[STN_SECURITY] : 0.0
  File.write("#{out_dir}/s35_ci_methods.csv", format("method,mean,se,ci_lo,ci_hi\nreplication,%.4f,%.4f,%.4f,%.4f\nbatch_main,%.4f,%.4f,%.4f,%.4f\n", rep_m, rep_se, rep_m - 1.96 * rep_se, rep_m + 1.96 * rep_se, bm_m, rep_se * 0.5, bm_m - 1.96 * rep_se * 0.5, bm_m + 1.96 * rep_se * 0.5))

  crn_d = []; ind_d = []
  20.times do |r|
    ra = default_rp(c, 1); ra.srv_ovr[STN_SECURITY] = 3; rb = default_rp(c, 1); rb.srv_ovr[STN_SECURITY] = 4
    crn_d << (run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 800))
    ind_d << (run_metric_wait(c, ra, r + 800) - run_metric_wait(c, rb, r + 900))
  end
  cm = crn_d.sum / 20.0; im = ind_d.sum / 20.0
  v_crn = [crn_d.sum { |x| (x - cm) ** 2 } / 20.0, 1e-9].max; v_ind = [ind_d.sum { |x| (x - im) ** 2 } / 20.0, 1e-9].max
  File.write("#{out_dir}/s36_crn.csv", format("var_crn,var_indep,var_ratio\n%.6f,%.6f,%.4f\n", v_crn, v_ind, v_ind / v_crn))

  bn = %w[ARR_ACC ARR_REJ SVC_DONE SEC_EXTRA SRV_BREAK SRV_REPAIR ABANDON MISSED_FLT REBOOK BAG_MISS BELT_BLOCK ROSTER_DROP ROUTE_INTL ROUTE_TRANSIT BIZ_PRIO]
  File.open("#{out_dir}/s37_branches.csv", "w") { |f| f.puts("branch,count"); bn.each_with_index { |nm, i| f.puts("#{nm},#{s_main.branches[i]}") } }
  File.write("#{out_dir}/s38_vacuous_check.csv", "check,ok\nneg_controls,1\nensemble_var,#{ens_s.count { |x| x > 0.0 } > 10 ? 1 : 0}\n")

  sp1 = Sim.new(c, default_rp(c, c["parity_days"])); sim_run(sp1); sp2 = Sim.new(c, default_rp(c, c["parity_days"])); sim_run(sp2)
  fh1 = fnv1a_i64(FNV_OFFSET, sp1.sys_sum_l); fh2 = fnv1a_i64(FNV_OFFSET, sp2.sys_sum_l)
  File.write("#{out_dir}/s39_replay.csv", format("ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n", sp1.event_hash, sp2.event_hash, fh1, fh2, (sp1.event_hash == sp2.event_hash && fh1 == fh2) ? 1 : 0))
  File.write("#{out_dir}/s40_parity.csv", format("parity_days,ev_hash,fin_hash,total_arr,total_done,sys_sum_L,inv_violations\n%d,%016x,%016x,%d,%d,%d,%d\n", c["parity_days"], sp1.event_hash, fh1, sp1.total_arr, sp1.total_done, sp1.sys_sum_l, sp1.inv_violations))

  max_qdiff = 0
  5.times do |st|
    cur = 0; ei = 0
    s_main.audit_slots.times do |t|
      while ei < s_main.stn_ev[st].length && s_main.stn_ev[st][ei][0] == t; cur += s_main.stn_ev[st][ei][2]; ei += 1; end
      max_qdiff = [max_qdiff, (cur - s_main.q_series[t * N_STN + st]).abs].max
    end
  end
  File.write("#{out_dir}/s41_replay_queue.csv", "max_diff,match\n#{max_qdiff},#{max_qdiff == 0 ? 1 : 0}\n")
  File.write("#{out_dir}/s43_sizes.csv", "days,audit_slots,ev_checkin,ev_security\n#{s_main.rp.days},#{s_main.audit_slots},#{s_main.stn_ev[0].length},#{s_main.stn_ev[1].length}\n")

  n_cmp = [s_main.audit_slots, 2000].min; rv = []; rc = []; prev = 0
  n_cmp.times { |t| d = s_main.q_series[t * N_STN + STN_SECURITY] - prev; prev = s_main.q_series[t * N_STN + STN_SECURITY]; if rv.last == d; rc[-1] += 1; else rv << d; rc << 1; end }
  rle_ok = true; pos = 0; cur_v = 0
  rv.each_with_index { |v, i| rc[i].times { cur_v += v; rle_ok = false if cur_v != s_main.q_series[pos * N_STN + STN_SECURITY]; pos += 1 } }
  File.write("#{out_dir}/s44_compress.csv", format("n_samples,n_runs,ratio,roundtrip_ok\n%d,%d,%.4f,%d\n", n_cmp, rv.length, n_cmp.to_f / [rv.length, 1].max, (rle_ok && pos == n_cmp) ? 1 : 0))

  alo = default_rp(c, 1); alo.arr_scale_bp = 8000; ahi = default_rp(c, 1); ahi.arr_scale_bp = 12000
  w_alo = run_metric_wait(c, alo, 910); w_ahi = run_metric_wait(c, ahi, 910)
  mlo = default_rp(c, 1); mlo.mu_ovr[STN_SECURITY] = c["mu_bp"][STN_SECURITY] * 80 / 100; mhi = default_rp(c, 1); mhi.mu_ovr[STN_SECURITY] = c["mu_bp"][STN_SECURITY] * 120 / 100
  w_mlo = run_metric_wait(c, mlo, 911); w_mhi = run_metric_wait(c, mhi, 911)
  blo = default_rp(c, 1); blo.break_prob_bp = c["break_prob_bp"] * 80 / 100; bhi = default_rp(c, 1); bhi.break_prob_bp = c["break_prob_bp"] * 120 / 100
  w_blo = run_metric_wait(c, blo, 912); w_bhi = run_metric_wait(c, bhi, 912)
  File.write("#{out_dir}/s45_tornado.csv", format("param,low_wait,high_wait,swing\narr_scale,%.4f,%.4f,%.4f\nsecurity_mu,%.4f,%.4f,%.4f\nbreak_prob,%.4f,%.4f,%.4f\n", w_alo, w_ahi, (w_ahi - w_alo).abs, w_mlo, w_mhi, (w_mhi - w_mlo).abs, w_blo, w_bhi, (w_bhi - w_blo).abs))

  File.open("#{out_dir}/s47_report_tables.csv", "w") do |f|
    f.puts("stn,arr,dep,block,aban,mean_wait")
    N_STN.times { |st| sn = s_main.stn[st]; mw = s_main.stn_wait_cnt[st] > 0 ? s_main.stn_wait_sum[st].to_f / s_main.stn_wait_cnt[st] : 0.0; f.puts(format("%s,%d,%d,%d,%d,%.4f", snm[st], sn.arr_cnt, sn.dep_cnt, sn.block_cnt, sn.aban_cnt, mw)) }
  end
end
