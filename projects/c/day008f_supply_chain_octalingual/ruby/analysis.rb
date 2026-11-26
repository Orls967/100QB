SKU_L = %w[A B C D].freeze

def write_s01_to_s09(c, out_dir)
  hrng = 0xCBF29CE484222325
  pairs = [[0,0],[1,2],[2,1]]
  exp_hv = [0xbe40b757fe0b6697, 0x659394d6f2a43a95, 0x3abda5bac0dd49cb]
  rng_ok = true
  3.times do |i|
    x = (c['master_seed'] & M64) ^ fnv1a_u32pair(pairs[i][0], pairs[i][1])
    sd, = splitmix64(x); st = Xoshiro256.new(sd)
    vals = (0..3).map { st.next_u64 }
    hv = fnv1a_bytes(vals.pack('Q4'))
    rng_ok = false if hv != exp_hv[i]
    hrng = fnv1a_i64(hrng, hv)
  end
  s = Sim.new(c, false, nil)
  fills = run_alloc_auction(s, 10, 0, [0,1,2], [1000,1500,1000], [1,2,3], [2,1,1])
  halloc = 0xCBF29CE484222325
  fills.each { |fl| halloc = fnv1a_i64(halloc, fl) }
  inv_ok = check_invariants(s)
  ok = (rng_ok && inv_ok == 1 && fills.sum == 1600) ? 1 : 0
  File.write("#{out_dir}/s01_selftest.txt", format("ok,h_rng,h_alloc,inv_ok\n%d,%016x,%016x,%d\n", ok, hrng, halloc, inv_ok))
  File.write("#{out_dir}/s02_config.csv", "key,value\nmaster_seed,#{c['master_seed']}\nnum_days,#{c['num_days']}\nparity_days,#{c['parity_days']}\nparity_events_per_day,#{c['parity_events_per_day']}\nwh_cap,#{c['wh_cap']}\nkit_fee,#{c['kit_fee']}\n")

  ent = ["id,type,name"]
  10.times { |i| ent << "#{i},STORE,S#{i}" }
  4.times { |i| ent << "#{10+i},WAREHOUSE,W#{i}" }
  3.times { |i| ent << "#{14+i},FACTORY,F#{i}" }
  3.times { |i| ent << "#{17+i},CARRIER,C#{i}" }
  ent += ["20,CUSTOMER,ExtCustomer", "21,SUPPLIER,RawSupplier", "22,FEEPOOL,FeePool"]
  File.write("#{out_dir}/s03_entities.csv", ent.join("\n") + "\n")
  File.write("#{out_dir}/s03_accounts_init.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,credit_limit,status\n" + s.acc.each_with_index.map { |a, i| "#{i},#{a[:cash]},#{a[:on_hand].join(',')},#{a[:credit_limit]},#{a[:status]}" }.join("\n") + "\n")
  File.write("#{out_dir}/s03_carriers.csv", "carrier,base_lt,cap_daily,rate_t0,rate_t1,rate_t2\n" + (0..2).map { |i| "#{i},#{CARRIER_BASE_LT[i]},#{c['carrier_cap_daily']},200,160,120" }.join("\n") + "\n")
  edges = ["src,dst,edge_type"] + (0..2).map { |i| "#{14+i},#{10+i},FAC_WH" } + (0..2).map { |i| "13,#{10+i},KIT_TS" } + (0..9).map { |i| "#{STORE_WH[i]},#{i},WH_STORE" }
  File.write("#{out_dir}/s03_network_edges.csv", edges.join("\n") + "\n")

  pol = ["case_id,inv_pos,s_lvl,S_lvl,order_qty,branch"]
  60.times do |k|
    inv = (k * 37) % 1000 - 100; oq = (inv < 250 && k % 11 != 0) ? (800 - inv) : 0
    br = k % 11 == 0 ? 'SUSPENDED' : (oq > 0 ? 'REORDER' : 'NO_OP')
    pol << "#{k},#{inv},250,800,#{oq},#{br}"
  end
  File.write("#{out_dir}/s04_policy_script.csv", pol.join("\n") + "\n")

  al = ["case_id,avail,f0,f1,f2,f3"]
  10.times do |k|
    sa = Sim.new(c, false, nil); sa.acc[10][:on_hand][0] = 200 + k * 150
    fl = run_alloc_auction(sa, 10, 0, [0,1,2,3], [400, 300 + k * 20, 500, 250], [10 + k, 5, 20 - k, 8], [2,1,2,1])
    al << "#{k},#{sa.acc[10][:on_hand][0]},#{fl.join(',')}"
  end
  File.write("#{out_dir}/s05_alloc_script.csv", al.join("\n") + "\n")

  sess = ["day,session,start_tick,end_tick"]
  c['num_days'].times { |d| sess += ["#{d},PRE_DAY,0,0", "#{d},ALLOCATION_AUCTION,0,1", "#{d},CONTINUOUS_OPS,0,#{c['ticks_per_day']}", "#{d},EOD_CLEARING,#{c['ticks_per_day']},#{c['ticks_per_day']}"] }
  File.write("#{out_dir}/s06_sessions.csv", sess.join("\n") + "\n")

  x0 = (c['master_seed'] & M64) ^ fnv1a_u32pair(2, 23); sd0, = splitmix64(x0); r0 = Xoshiro256.new(sd0)
  sh = ["idx,time,seq,etype,dst,sku,qty"]
  200.times { |i| sh << "#{i},#{(i * c['ticks_per_day']) / 8000},#{i + 1},0,#{r0.rand_int(0, 9)},#{r0.rand_int(0, 3)},#{r0.rand_int(1, 14)}" }
  File.write("#{out_dir}/s06_schedule_head.csv", sh.join("\n") + "\n")
  File.write("#{out_dir}/s07_shock_schedule.csv", "idx,cls,prob_pct,mag\n0,SMALL,3,10\n1,MEDIUM,3,20\n2,LARGE,3,30\n")

  sw = Sim.new(c, false, nil); sim_run_days(sw, 1, 2000, 0, 0, 1)
  File.write("#{out_dir}/s08_warmup_summary.csv", "days,sold_A,sold_B,sold_C,sold_D,inv_ok\n1,#{sw.total_sold.join(',')},#{sw.inv_violations == 0 ? 1 : 0}\n")
  4.times { |k| File.write("#{out_dir}/s08_warmup_sku_#{SKU_L[k]}.csv", "sku,sold,produced,end_snap\n#{SKU_L[k]},#{sw.total_sold[k]},#{sw.total_produced[k]},#{sw.day_inv_snap[0][k]}\n") }
  File.write("#{out_dir}/s09_day1_summary.csv", "sold_tot,kits,alloc_cases,inv_violations\n#{sw.total_sold.sum},#{sw.total_kits},#{sw.alloc_cases},#{sw.inv_violations}\n")
  cal_fr = sw.dem_units_tot > 0 ? sw.dem_units_filled.to_f / sw.dem_units_tot : 1.0
  cal_bo = sw.dem_events_tot > 0 ? sw.dem_bo_events.to_f / sw.dem_events_tot : 0.0
  cal_cu = sw.carrier_cap_tot > 0 ? sw.carrier_vol_tot.to_f / sw.carrier_cap_tot : 0.0
  cal_ok = (cal_fr >= 0.85 && cal_fr <= 0.99 && cal_bo < 0.20 && cal_cu >= 0.30 && cal_cu <= 0.90) ? 1 : 0
  File.write("#{out_dir}/s09_calibration.csv", format("fill_rate,bo_event_frac,carrier_util,calibrated\n%.6f,%.6f,%.6f,%d\n", cal_fr, cal_bo, cal_cu, cal_ok))
end

def write_s11_to_s38(s, c, out_dir)
  4.times do |k|
    File.write("#{out_dir}/sales_#{SKU_L[k]}.csv", "day,tick,store_id,sku,qty,rev_cents\n" + s.sales.select { |r| r[3] == k }.map { |r| r.join(',') }.join("\n") + "\n")
    File.write("#{out_dir}/shipments_#{SKU_L[k]}.csv", "day,tick,src_id,dst_id,carrier_id,sku,qty,goods_cents,freight_cents,status\n" + s.ships.select { |r| r[5] == k }.map { |r| r.join(',') }.join("\n") + "\n")
    File.write("#{out_dir}/inventory_snap_#{SKU_L[k]}.csv", "day,sku,on_hand_tot\n" + (0...c['num_days']).map { |d| "#{d},#{SKU_L[k]},#{s.day_inv_snap[d][k]}" }.join("\n") + "\n")
  end
  10.times { |i| File.write("#{out_dir}/store_daily_#{i}.csv", "store,cash,oh_tot,req_tot,fill_tot\n#{i},#{s.acc[i][:cash]},#{s.acc[i][:on_hand].sum},#{s.daily_fill_store[i][0]},#{s.daily_fill_store[i][1]}\n") }
  4.times { |i| a = s.acc[10+i]; File.write("#{out_dir}/warehouse_daily_#{i}.csv", "wh,cash,oh_A,oh_B,oh_C,oh_D\n#{i},#{a[:cash]},#{a[:on_hand].join(',')}\n") }
  3.times { |i| File.write("#{out_dir}/factory_daily_#{i}.csv", "fac,cash,produced,state\n#{i},#{s.acc[14+i][:cash]},#{s.total_produced[i]},#{s.acc[14+i][:plant_state]}\n") }
  File.write("#{out_dir}/production.csv", "day,tick,factory_id,sku,qty,cost_cents\n" + s.prods.map { |r| r.join(',') }.join("\n") + "\n")
  File.write("#{out_dir}/kit_assembly.csv", "day,tick,wh_id,qty,fee_cents,imbalance\n" + s.kits.map { |r| r.join(',') }.join("\n") + "\n")
  File.write("#{out_dir}/fees_ledger.csv", "day,tick,from_id,to_id,amount,kind\n" + s.fees.map { |r| r.join(',') }.join("\n") + "\n")
  File.write("#{out_dir}/entities_final.csv", "acc,cash,oh_A,oh_B,oh_C,oh_D,it_A,it_B,it_C,it_D,status\n" + s.acc.each_with_index.map { |a, i| "#{i},#{a[:cash]},#{a[:on_hand].join(',')},#{a[:in_transit].join(',')},#{a[:status]}" }.join("\n") + "\n")

  File.write("#{out_dir}/s11_invariants.csv", "violations,ok\n#{s.inv_violations},#{s.inv_violations == 0 ? 1 : 0}\n")
  File.write("#{out_dir}/s12_clearing.csv", "days,fee_pool_cash\n#{c['num_days']},#{s.acc[ACC_FEE][:cash]}\n")
  File.write("#{out_dir}/s12_credit.csv", "credit_holds,liquidations\n#{s.branches[13]},#{s.branches[14]}\n")
  File.write("#{out_dir}/s13_risk.csv", "credit_holds,liquidations\n#{s.branches[13]},#{s.branches[14]}\n")
  File.write("#{out_dir}/s14_factories.csv", "fac,produced\n0,#{s.total_produced[0]}\n1,#{s.total_produced[1]}\n2,#{s.total_produced[2]}\n")
  File.write("#{out_dir}/s15_carrier_fees.csv", "tier,count\n0,#{s.carrier_tier_cnt[0]}\n1,#{s.carrier_tier_cnt[1]}\n2,#{s.carrier_tier_cnt[2]}\n")
  File.write("#{out_dir}/s16_kits.csv", "total_kits,full,partial\n#{s.total_kits},#{s.branches[8]},#{s.branches[9]}\n")
  File.write("#{out_dir}/s16_kit_imbalance_series.csv", "idx,imbalance\n" + s.kit_imb_series.each_with_index.map { |v, i| "#{i},#{v}" }.join("\n") + "\n")
  File.write("#{out_dir}/s17_allocations.csv", "cases,nonzero,total_req,total_fill\n#{s.alloc_cases},#{s.alloc_nonzero},#{s.alloc_total_req},#{s.alloc_total_fill}\n")
  File.write("#{out_dir}/s18_cost_entities.csv", "acc,cash_delta\n" + (0..16).map { |i| "#{i},#{s.acc[i][:cash] - (i < 10 ? 2_000_000_000 : 5_000_000_000) * s.cap_scale}" }.join("\n") + "\n")
  File.write("#{out_dir}/s19_cost_types.csv", "echelon,cost_cents\n" + (0..3).map { |e| "#{e},#{(0...c['num_days']).sum { |d| s.daily_cost_echelon[d][e] }}" }.join("\n") + "\n")
  File.write("#{out_dir}/s19_cost_daily.csv", "day,prod,kit,freight,hold_pen\n" + (0...c['num_days']).map { |d| "#{d},#{s.daily_cost_echelon[d].join(',')}" }.join("\n") + "\n")
  4.times { |k| File.write("#{out_dir}/s20_demand_#{SKU_L[k]}.csv", "day,demand\n" + (0...c['num_days']).map { |d| "#{d},#{s.daily_demand[d][k]}" }.join("\n") + "\n") }
  bw_lines = (0..3).map do |k|
    md = s.bin_demand[k].sum / 100.0; mo = s.bin_orders[k].sum / 100.0
    m2 = s.bin_demand[k].sum { |x| (x - md)**2 }
    m4 = s.bin_demand[k].sum { |x| (x - md)**4 }
    vo = s.bin_orders[k].sum { |x| (x - mo)**2 }
    cov1 = (1...100).sum { |b| (s.bin_demand[k][b] - md) * (s.bin_demand[k][b - 1] - md) }
    kurt = m2 > 0 ? (m4 / 100.0) / ((m2 / 100.0)**2) : 0.0
    acf1 = m2 > 0 ? cov1 / m2 : 0.0
    bwr = m2 > 0 ? vo / m2 : 1.0
    format("%s,%.6f,%.6f,%.6f", SKU_L[k], kurt, acf1, bwr)
  end
  File.write("#{out_dir}/s21_bullwhip.csv", "sku,kurtosis,acf_lag1,bullwhip_ratio\n" + bw_lines.join("\n") + "\n")
  File.write("#{out_dir}/s22_service.csv", "store,req,fill,fill_rate\n" + (0..9).map { |i| fr = s.daily_fill_store[i][0] > 0 ? s.daily_fill_store[i][1].to_f / s.daily_fill_store[i][0] : 1.0; format("%d,%d,%d,%.6f", i, s.daily_fill_store[i][0], s.daily_fill_store[i][1], fr) }.join("\n") + "\n")
  File.write("#{out_dir}/s23_leadtime_impact.csv", "idx,leadtime,qty\n" + s.leadtime_pairs.each_with_index.map { |p, i| "#{i},#{p[0]},#{p[1]}" }.join("\n") + "\n")
  imb_m = s.kit_imb_series.empty? ? 0.0 : s.kit_imb_series.sum.to_f / s.kit_imb_series.length
  File.write("#{out_dir}/s24_kit_coupling.csv", format("mean_imbalance,total_kits\n%.4f,%d\n", imb_m, s.total_kits))
  shk_lines = (0..2).map do |cl|
    vals_r = s.shock_evs.select { |ev| ev[1] == cl }.filter_map do |b, _, st|
      pre = ((-2..-1).select { |db| b + db >= 0 }).map { |db| s.bin_store_dem[st][b + db].to_f }
      post = ((0..1).select { |db| b + db < 100 }).map { |db| s.bin_store_dem[st][b + db].to_f }
      (!pre.empty? && !post.empty? && pre.sum > 0) ? ((post.sum / post.length) - (pre.sum / pre.length)) / (pre.sum / pre.length) : nil
    end
    cnt = vals_r.length
    mn = cnt > 0 ? vals_r.sum / cnt : 0.0
    se = cnt > 1 ? Math.sqrt([0.0, (vals_r.sum { |x| x * x } / cnt - mn * mn) / (cnt - 1)].max) : 0.0
    format("%d,%d,%.6f,%.6f", cl, cnt, mn, se)
  end
  File.write("#{out_dir}/s25_shock_event_study.csv", "cls,count,mean_rel_dev,se_rel_dev\n" + shk_lines.join("\n") + "\n")
  File.write("#{out_dir}/s26_disruptions.csv", "disrupt_hits\n#{s.branches[11]}\n")
  File.write("#{out_dir}/s27_fees.csv", "n_fee_entries,fee_pool_cash\n#{s.fees.length},#{s.acc[ACC_FEE][:cash]}\n")
  File.write("#{out_dir}/s28_capacity_rejects.csv", "wh_cap_reject,carrier_cap_reject\n#{s.branches[4]},#{s.branches[5]}\n")

  sw1 = ["S_lvl,rep,fill_rate,avg_oh"]
  [80, 220, 350].each { |sv| 5.times { |rp| c2 = c.dup; c2['master_seed'] += rp + 1; ss = Sim.new(c2, false, nil); sim_run_days(ss, 2, 1600, sv, -1, 1); oh = (0..9).sum { |i| ss.acc[i][:on_hand][0] }; fr = ss.dem_units_tot > 0 ? ss.dem_units_filled.to_f / ss.dem_units_tot : 1.0; sw1 << format("%d,%d,%.6f,%.2f", sv, rp, fr, oh / 10.0) } }
  File.write("#{out_dir}/s29_sweep_safety_stock.csv", sw1.join("\n") + "\n")
  sw2 = ["disrupt_pct,rep,produced,sold"]
  [0, 4, 12].each { |dp| 5.times { |rp| c2 = c.dup; c2['master_seed'] += rp + 10; ss = Sim.new(c2, false, nil); sim_run_days(ss, 2, 1600, 0, dp, 1); sw2 << "#{dp},#{rp},#{ss.total_produced.sum},#{ss.total_sold.sum}" } }
  File.write("#{out_dir}/s30_sweep_shock_disrupt.csv", sw2.join("\n") + "\n")
  sw3 = ["transship_on,rep,mean_imbalance,n_transships"]
  [0, 1].each { |ts| 5.times { |rp| c2 = c.dup; c2['master_seed'] += rp + 20; ss = Sim.new(c2, false, nil); sim_run_days(ss, 2, 1600, 0, -1, ts); im = ss.kit_imb_series.empty? ? 0.0 : ss.kit_imb_series.sum.to_f / ss.kit_imb_series.length; sw3 << format("%d,%d,%.4f,%d", ts, rp, im, ss.n_transships) } }
  File.write("#{out_dir}/s31_sweep_transship.csv", sw3.join("\n") + "\n")

  bn = %w[DEMAND_FILL DEMAND_BO DEMAND_LOST STORE_REPLEN WH_CAP_REJECT CARRIER_CAP_REJECT ALLOC_FULL ALLOC_PARTIAL KIT_FULL KIT_PARTIAL TRANSSHIP DISRUPT_HIT SHOCK_HIT CREDIT_HOLD LIQUIDATE CARRIER_DELAY]
  File.write("#{out_dir}/s32_branches.csv", "branch,count\n" + bn.each_with_index.map { |nm, i| "#{nm},#{s.branches[i]}" }.join("\n") + "\n")

  sp1 = Sim.new(c, false, nil, c['parity_events_per_day']); sim_run_days(sp1, c['parity_days'], c['parity_events_per_day'], 0, -1, 1)
  sp2 = Sim.new(c, false, nil, c['parity_events_per_day']); sim_run_days(sp2, c['parity_days'], c['parity_events_per_day'], 0, -1, 1)
  h1 = 0xCBF29CE484222325; h2 = 0xCBF29CE484222325
  NUM_ACCOUNTS.times { |i| h1 = fnv1a_i64(h1, sp1.acc[i][:cash]); h2 = fnv1a_i64(h2, sp2.acc[i][:cash]); 4.times { |k| h1 = fnv1a_i64(h1, sp1.acc[i][:on_hand][k]); h2 = fnv1a_i64(h2, sp2.acc[i][:on_hand][k]) } }
  File.write("#{out_dir}/s33_replay.csv", format("ev_hash_1,ev_hash_2,fin_hash_1,fin_hash_2,match\n%016x,%016x,%016x,%016x,%d\n", sp1.event_hash, sp2.event_hash, h1, h2, (sp1.event_hash == sp2.event_hash && h1 == h2) ? 1 : 0))
  File.write("#{out_dir}/s34_parity.csv", format("parity_days,parity_ev_per_day,ev_hash,fin_hash,sold_tot,kits,inv_violations\n%d,%d,%016x,%016x,%d,%d,%d\n", c['parity_days'], c['parity_events_per_day'], sp1.event_hash, h1, sp1.total_sold.sum, sp1.total_kits, sp1.inv_violations))

  vals = (0..4).map { sb = Sim.new(c, false, nil); t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC); sim_run_days(sb, 1, 2000, 0, -1, 1); ((Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0) * 1e9) / 2000.0 }.sort
  med = vals[2]; dev = vals.map { |v| (v - med).abs }.sort
  File.write("#{out_dir}/s35_bench.csv", format("events,ns_per_ev_med,ns_per_ev_mad,ev_per_sec_med\n2000,%.4f,%.4f,%.2f\n", med, dev[2], 1e9 / med))
  File.write("#{out_dir}/s36_sizes.csv", "n_sales,n_ships,n_prods,n_kits,n_fees\n#{s.sales.length},#{s.ships.length},#{s.prods.length},#{s.kits.length},#{s.fees.length}\n")

  raw = s.sales.first(4096).map { |r| r[4] }; nq = raw.length; rle = []; i = 0
  while i < nq; dv = i == 0 ? raw[0] : (raw[i] - raw[i - 1]); run = 1; run += 1 while i + run < nq && (raw[i + run] - raw[i + run - 1]) == dv; rle << dv << run; i += run; end
  dec = []; rle.each_slice(2) { |dv, run| run.times { dec << (dec.empty? ? dv : dec[-1] + dv) } }
  File.write("#{out_dir}/s37_compress.csv", format("n_items,rle_ints,ratio,roundtrip_ok\n%d,%d,%.4f,%d\n", nq, rle.length, rle.empty? ? 1.0 : nq.to_f / rle.length, dec == raw ? 1 : 0))

  si = Sim.new(c, false, nil, s.cap_scale * 2000); l_cash = si.acc.map { |a| a[:cash] }; l_oh = si.acc.map { |a| a[:on_hand].dup }; l_it = Array.new(NUM_ACCOUNTS) { [0,0,0,0] }
  s.sales.each { |r| l_cash[ACC_CUST] -= r[5]; l_cash[r[2]] += r[5]; l_oh[r[2]][r[3]] -= r[4] }
  s.ships.each { |r| if r[9] == 1; l_oh[r[2]][r[5]] -= r[6]; l_it[r[3]][r[5]] += r[6]; l_cash[r[3]] -= (r[7] + r[8]); l_cash[r[2]] += r[7]; l_cash[17 + r[4]] += r[8]; else l_it[r[3]][r[5]] -= r[6]; l_oh[r[3]][r[5]] += r[6]; end }
  s.prods.each { |r| l_cash[r[2]] -= r[5]; l_cash[ACC_SUPP] += r[5]; l_oh[r[2]][r[3]] += r[4] }
  s.kits.each  { |r| l_cash[r[2]] -= r[4]; l_cash[ACC_FEE] += r[4]; 3.times { |k| l_oh[r[2]][k] -= r[3] }; l_oh[r[2]][3] += r[3] }
  s.fees.each  { |r| l_cash[r[2]] -= r[4]; l_cash[r[3]] += r[4] }
  max_cd = (0...NUM_ACCOUNTS).map { |j| (l_cash[j] - s.acc[j][:cash]).abs }.max
  max_id = (0...NUM_ACCOUNTS).flat_map { |j| (0..3).map { |k| (l_oh[j][k] - s.acc[j][:on_hand][k]).abs + (l_it[j][k] - s.acc[j][:in_transit][k]).abs } }.max
  File.write("#{out_dir}/s38_ledger_replay.csv", "max_cash_diff,max_inv_diff,match\n#{max_cd},#{max_id},#{(max_cd == 0 && max_id == 0) ? 1 : 0}\n")
end
