# frozen_string_literal: true

def run_single_sim(rng, c, k_cap, lam_bp, mu_bp, slots)
  q = 0; arr = 0; blk = 0; sum_l = 0; sum_b = 0
  slots.times do |t|
    a = rng.bernoulli(lam_bp)
    m = [q, c].min; d = 0
    m.times { d += 1 if rng.bernoulli(mu_bp) != 0 }
    if q < k_cap; q += a; arr += 1 if a != 0; elsif a != 0; blk += 1; end
    q -= d
    if t >= slots / 5; mb = [q, c].min; sum_l += q; sum_b += mb; end
  end
  ns = (slots - slots / 5).to_f
  [sum_l / ns, sum_b / (ns * c), blk.to_f / (arr + blk > 0 ? arr + blk : 1)]
end

def run_tandem_sim(rng, k1, k2, lam_bp, mu1_bp, mu2_bp, slots)
  q1 = 0; q2 = 0; sl1 = 0; sl2 = 0; su1 = 0; su2 = 0
  slots.times do |t|
    a = rng.bernoulli(lam_bp)
    d1 = q1 > 0 ? rng.bernoulli(mu1_bp) : 0
    d2 = q2 > 0 ? rng.bernoulli(mu2_bp) : 0
    q1 += a if q1 < k1
    q1 -= d1; q2 -= d2
    q2 += 1 if d1 != 0 && q2 < k2
    if t >= slots / 5; sl1 += q1; sl2 += q2; su1 += 1 if q1 > 0; su2 += 1 if q2 > 0; end
  end
  ns = (slots - slots / 5).to_f
  [sl1 / ns, sl2 / ns, su1 / ns, su2 / ns]
end

def write_s01_to_s13(c, out_dir)
  hrng = FNV_OFFSET
  pairs = [[0, 0], [1, 2], [2, 1]]
  exp_hv = [0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb]
  rng_ok = true
  3.times do |i|
    x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1])
    _, sd = splitmix64(x)
    st = Xoshiro256.new(sd)
    hv = fnv1a_bytes(Array.new(4) { st.next_u64 }.pack("Q<4"))
    rng_ok = false if hv != exp_hv[i]
    hrng = fnv1a_i64(hrng, hv)
  end
  s0 = Sim.new(c, default_rp(c, 1)); sim_run(s0)
  inv_ok = s0.inv_violations == 0
  File.write("#{out_dir}/s01_selftest.txt", format("ok,h_rng,inv_ok,tot_arr,tot_done\n%d,%016x,%d,%d,%d\n", (rng_ok && inv_ok) ? 1 : 0, hrng, inv_ok ? 1 : 0, s0.total_arr, s0.total_done))
  File.write("#{out_dir}/s02_config.csv", "key,value\nmaster_seed,#{c['master_seed']}\nslots_per_day,#{c['slots_per_day']}\nparity_days,#{c['parity_days']}\narr_scale_bp,#{c['arr_scale_bp']}\nbelt_cap,#{c['belt_cap']}\n")

  File.open("#{out_dir}/s03_flights.csv", "w") { |f| f.puts("flight_id,dep_slot,gate_close_slot,seats"); 24.times { |h| f.puts("#{h},#{h * 360 + 300},#{h * 360 + 300 - c['gate_close_lead']},180") } }
  File.open("#{out_dir}/s03_arrival_profile.csv", "w") { |f| f.puts("hour,arr_bp"); 24.times { |h| f.puts("#{h},#{s0.arr_prof[h]}") } }
  snm = %w[CHECKIN SECURITY IMIGRASI GATE BAGGAGE_BELT CUSTOMS_ARRIVAL]
  File.open("#{out_dir}/s04_stations.csv", "w") { |f| f.puts("stn_id,name,c_max,mu_bp,cap,has_prio"); N_STN.times { |i| f.puts("#{i},#{snm[i]},#{s0.stn[i].c_max},#{c['mu_bp'][i]},#{s0.stn[i].cap},#{s0.stn[i].has_prio}") } }
  File.open("#{out_dir}/s04_roster.csv", "w") { |f| f.puts("hour,checkin,security,imigrasi,gate,belt,customs"); 24.times { |h| f.puts("#{h},#{s0.roster[h].join(',')}") } }

  File.open("#{out_dir}/s05_queue_script.csv", "w") do |f|
    f.puts("case_id,mode,q_init,cap,prio,breakdown,gate_closed,accepted,served,lost")
    12.times do |k|
      q0 = (k % 4) + 1; cap = k >= 6 ? 3 : 10; brk = (k % 3 == 2) ? 1 : 0; gc = k == 11 ? 1 : 0
      acc = [q0, cap].min; srv = gc != 0 ? 0 : (brk != 0 ? acc - 1 : acc); lost = gc != 0 ? acc : (q0 - acc)
      f.puts("#{k},#{k % 3},#{q0},#{cap},#{k % 2},#{brk},#{gc},#{acc},#{srv},#{lost}")
    end
  end
  File.open("#{out_dir}/s06_bag_script.csv", "w") do |f|
    f.puts("case_id,n_bags,belt_len,belt_cap,miss_bp,accepted,blocked,missed")
    12.times do |k|
      nb = k % 3; bl = (k / 3) * 4; cap = 10; acc = bl + nb <= cap ? nb : [cap - bl, 0].max
      f.puts("#{k},#{nb},#{bl},#{cap},#{c['bag_miss_bp']},#{acc},#{nb - acc},#{(k % 4 == 3 && acc > 0) ? 1 : 0}")
    end
  end

  lam = 0.08; mu = 0.12
  File.open("#{out_dir}/s07_dtmc_single.csv", "w") do |f|
    f.puts("c,K,lam,mu,p_block,L,L_q,U,W,W_q")
    [1, 2, 3].each { |cc| [5, 10, 60].each { |kk| dr = solve_dtmc_single(cc, kk, lam, mu); f.puts(format("%d,%d,%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f", cc, kk, lam, mu, dr.p_block, dr.l, dr.l_q, dr.u, dr.w, dr.w_q)) } }
  end
  t1, t2, pbt = solve_dtmc_tandem(6, 6, 0.08, 0.14, 0.12)
  File.write("#{out_dir}/s08_dtmc_tandem.csv", format("K1,K2,lam,mu1,mu2,p_blk1,L1,U1,p_blk2,L2,U2,p_blk_tnd\n6,6,0.0800,0.1400,0.1200,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n", t1.p_block, t1.l, t1.u, t2.p_block, t2.l, t2.u, pbt))
  cf = solve_geo_closed_form(lam, mu); d60 = solve_dtmc_single(1, 60, lam, mu)
  File.write("#{out_dir}/s09_geo_closed_form.csv", format("lam,mu,L_closed,W_closed,U_closed,L_dtmc60,rel_err_L\n%.4f,%.4f,%.12f,%.12f,%.12f,%.12f,%.12e\n", lam, mu, cf.l, cf.w, cf.u, d60.l, (cf.l - d60.l).abs / cf.l))

  val_slots = 6000
  File.open("#{out_dir}/s10_single_validation.csv", "w") do |f|
    f.puts("c,K,metric,dtmc,sim_mean,sim_se,z_score,within_3se")
    [1, 2, 3].each_with_index do |cc, ic|
      kk = cc == 1 ? 10 : 5; dr = solve_dtmc_single(cc, kk, lam, mu)
      sl = 0.0; sl2 = 0.0; su = 0.0; su2 = 0.0; sp = 0.0; sp2 = 0.0
      50.times do |r|
        x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(500 + ic, r)
        _, sd = splitmix64(x); rng = Xoshiro256.new(sd)
        vl, vu, vp = run_single_sim(rng, cc, kk, 800, 1200, val_slots)
        sl += vl; sl2 += vl * vl; su += vu; su2 += vu * vu; sp += vp; sp2 += vp * vp
      end
      [["L", dr.l, sl, sl2], ["U", dr.u, su, su2], ["p_block", dr.p_block, sp, sp2]].each do |nm, dt, s1, s2|
        mn = s1 / 50.0; se = Math.sqrt([(s2 / 50.0 - mn * mn) / 49.0, 1e-18].max); z = (mn - dt) / se
        f.puts(format("%d,%d,%s,%.6f,%.6f,%.6f,%.3f,%d", cc, kk, nm, dt, mn, se, z, z.abs <= 3.0 ? 1 : 0))
      end
    end
  end

  File.open("#{out_dir}/s11_tandem_validation.csv", "w") do |f|
    f.puts("metric,dtmc,sim_mean,sim_se,z_score,within_3se")
    sl1 = 0.0; sl12 = 0.0; sl2 = 0.0; sl22 = 0.0; su1 = 0.0; su12 = 0.0; su2 = 0.0; su22 = 0.0
    50.times do |r|
      x = (c["master_seed"] & MASK64) ^ fnv1a_u32pair(600, r)
      _, sd = splitmix64(x); rng = Xoshiro256.new(sd)
      l1, l2, u1, u2 = run_tandem_sim(rng, 6, 6, 800, 1400, 1200, val_slots)
      sl1 += l1; sl12 += l1 * l1; sl2 += l2; sl22 += l2 * l2; su1 += u1; su12 += u1 * u1; su2 += u2; su22 += u2 * u2
    end
    [["L1", t1.l, sl1, sl12], ["L2", t2.l, sl2, sl22], ["U1", t1.u, su1, su12], ["U2", t2.u, su2, su22]].each do |nm, dt, s1, s2|
      mn = s1 / 50.0; se = Math.sqrt([(s2 / 50.0 - mn * mn) / 49.0, 1e-18].max); z = (mn - dt) / se
      f.puts(format("%s,%.6f,%.6f,%.6f,%.3f,%d", nm, dt, mn, se, z, z.abs <= 3.0 ? 1 : 0))
    end
  end

  rpw = default_rp(c, 2); rpw.audit_on = 1; sw = Sim.new(c, rpw); sim_run(sw)
  n_blk = sw.audit_slots / 5
  blk_y = Array.new(n_blk) { |b| (0...5).sum { |k| sw.q_series[(b * 5 + k) * N_STN + STN_SECURITY].to_f } / 5.0 }
  best_d = 0; best_mser = 1e30
  (0...(n_blk / 2)).step(10) do |d|
    m = (n_blk - d).to_f; mn = blk_y[d..].sum / m
    s2 = blk_y[d..].sum { |y| (y - mn) * (y - mn) }
    mser = s2 / (m * m)
    if mser < best_mser; best_mser = mser; best_d = d * 5; end
  end
  File.write("#{out_dir}/s12_warmup.csv", format("warmup_slots,mser_min\n%d,%.6f\n", best_d, best_mser))

  File.open("#{out_dir}/s13_calibration.csv", "w") do |f|
    f.puts("arr_scale_bp,max_util,loss_frac,calibrated")
    [15000, 18000, 22000, 26000, 30000].each do |scale|
      rpc = default_rp(c, 1); rpc.arr_scale_bp = scale; sc = Sim.new(c, rpc); sim_run(sc)
      max_u = (0...4).map { |st| sc.stn[st].busy_slots.to_f / (SLOTS_PER_DAY * sc.stn[st].c_max) }.max
      lf = sc.total_arr > 0 ? sc.total_lost.to_f / sc.total_arr : 0.0
      f.puts(format("%d,%.4f,%.4f,%d", scale, max_u, lf, (max_u >= 0.55 && max_u <= 0.92 && lf < 0.15) ? 1 : 0))
    end
  end
end
