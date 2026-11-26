def do_ship(s, src, dst, sku, qty, day, tick, t_abs)
  return if qty <= 0 || s.acc[src][:on_hand][sku] <= 0
  qty = s.acc[src][:on_hand][sku] if qty > s.acc[src][:on_hand][sku]
  if dst >= 10 && dst < 14 && s.acc[dst][:on_hand].sum + qty > s.cfg['wh_cap']
    s.branches[4] += 1; return
  end
  cid = (src + dst + sku) % NUM_CARRIERS
  if s.carrier_day_vol[cid] + qty > s.cfg['carrier_cap_daily']
    cid = (cid + 1) % NUM_CARRIERS
    if s.carrier_day_vol[cid] + qty > s.cfg['carrier_cap_daily']
      s.branches[5] += 1; return
    end
  end
  rate = carrier_rate(s.carrier_day_vol[cid])
  s.carrier_tier_cnt[rate == 200 ? 0 : (rate == 160 ? 1 : 2)] += 1
  s.carrier_day_vol[cid] += qty
  goods = qty * SKU_WHOLESALE[sku]; freight = qty * rate
  s.acc[src][:on_hand][sku] -= qty; s.acc[dst][:in_transit][sku] += qty
  transfer_cash(s, dst, src, goods, -1, day, tick)
  transfer_cash(s, dst, 17 + cid, freight, -1, day, tick)
  s.daily_cost_echelon[day][2] += freight if day < 6
  delay = s.rng[17 + cid].rand_int(0, 80)
  s.branches[15] += 1 if delay > 50
  lt = CARRIER_BASE_LT[cid] + delay
  s.leadtime_pairs << [lt, qty] if s.leadtime_pairs.length < 256
  s.ships << [day, tick, src, dst, cid, sku, qty, goods, freight, 1] if s.log_enabled
  s.heap.push([t_abs + lt, s.next_seq, 2, src, dst, sku, qty, cid]); s.next_seq += 1
end

def hash_ev(s, ev)
  [0, 2, 3, 4, 5, 6].each { |i| s.event_hash = fnv1a_i64(s.event_hash, ev[i]) }
end

def sim_run_days(s, days, ev_per_day, safety_ovr, disrupt_ovr, transship_on)
  t_day = s.cfg['ticks_per_day']; ev_proc = 0
  days.times do |d|
    day_base = d * t_day; s.carrier_day_vol = [0, 0, 0]
    (14..16).each do |f|
      if s.acc[f][:plant_state] == 1
        s.acc[f][:disrupt_rem] -= 1
        s.acc[f][:plant_state] = 2 if s.acc[f][:disrupt_rem] <= 0
      elsif s.acc[f][:plant_state] == 2
        s.acc[f][:plant_state] = 0
      end
    end
    (10..12).each do |wh|
      4.times do |sku|
        sts = []; reqs = []; seqs = []; prio = []
        s_lvl = safety_ovr > 0 ? safety_ovr : s.cfg['store_S_base']
        10.times do |st|
          next if STORE_WH[st] != wh
          inv = s.acc[st][:on_hand][sku] + s.acc[st][:in_transit][sku]
          if inv < s_lvl
            sts << st; reqs << (s_lvl - inv); seqs << s.next_seq; s.next_seq += 1; prio << (st % 3 == 0 ? 2 : 1)
          end
        end
        unless sts.empty?
          fills = run_alloc_auction(s, wh, sku, sts, reqs, seqs, prio)
          sts.each_with_index do |st, i|
            s.f_alloc.puts "#{d},#{wh},#{sku},#{st},#{reqs[i]},#{fills[i]}" if s.f_alloc
            do_ship(s, wh, st, sku, fills[i], d, 0, day_base) if fills[i] > 0
          end
        end
      end
    end
    ev_per_day.times do |e|
      t_abs = day_base + (e * t_day) / ev_per_day
      st = s.rng[NUM_ACCOUNTS].rand_int(0, 9)
      sku = s.rng[NUM_ACCOUNTS].rand_int(0, 3)
      qty = s.rng[NUM_ACCOUNTS].rand_int(1, 14)
      s.heap.push([t_abs, s.next_seq, 0, ACC_CUST, st, sku, qty, 0]); s.next_seq += 1
      if e % 20 == 0
        f = (e / 20) % 3
        s.heap.push([t_abs + 1, s.next_seq, 3, 14 + f, 10 + f, f, 80, 0]); s.next_seq += 1
      end
      if e % 40 == 0
        s.heap.push([t_abs + 2, s.next_seq, 4, 13, 13, 3, 25, 0]); s.next_seq += 1
      end
      if e % 60 == 0 && transship_on != 0
        s.heap.push([t_abs + 3, s.next_seq, 5, 13, 10 + ((e / 60) % 3), (e / 60) % 4, 40, 0]); s.next_seq += 1
      end
      if e % 250 == 0
        dp = disrupt_ovr >= 0 ? disrupt_ovr : s.cfg['disrupt_prob_pct']
        if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= dp
          f = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2)
          s.heap.push([t_abs + 4, s.next_seq, 7, 14 + f, 14 + f, f, 1, 0]); s.next_seq += 1
        end
        if s.rng[NUM_ACCOUNTS + 1].rand_int(1, 100) <= s.cfg['shock_prob_pct']
          st_s = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 9)
          cls = s.rng[NUM_ACCOUNTS + 1].rand_int(0, 2)
          s.heap.push([t_abs + 5, s.next_seq, 6, ACC_CUST, st_s, 0, (cls + 1) * 10, cls]); s.next_seq += 1
        end
        s.heap.push([t_abs + 6, s.next_seq, 8, 0, 0, 0, 0, 0]); s.next_seq += 1
      end
    end
    day_end = (d + 1) * t_day
    while !s.heap.empty? && s.heap.peek_time < day_end
      ev = s.heap.pop; tick = ev[0] - day_base; hash_ev(s, ev); ev_proc += 1
      s.f_orders.puts "#{d},#{tick},#{ev[1]},#{ev[2]},#{ev[3]},#{ev[4]},#{ev[5]},#{ev[6]}" if ev_proc <= 5000 && s.f_orders
      case ev[2]
      when 0
        st = ev[4]; sku = ev[5]; qty = ev[6]
        if d < 6; s.daily_demand[d][sku] += qty; s.daily_fill_store[st][0] += qty; end
        al = s.cfg['ewma_alpha_milli']
        s.acc[st][:ewma_milli][sku] = (al * qty * 1000 + (1000 - al) * s.acc[st][:ewma_milli][sku]) / 1000
        if s.acc[st][:status] == 3
          s.branches[2] += 1
        elsif s.acc[st][:on_hand][sku] >= qty
          s.acc[st][:on_hand][sku] -= qty; s.total_sold[sku] += qty; rev = qty * SKU_RETAIL[sku]
          transfer_cash(s, ACC_CUST, st, rev, -1, d, tick)
          s.sales << [d, tick, st, sku, qty, rev] if s.log_enabled
          s.branches[0] += 1; s.daily_fill_store[st][1] += qty if d < 6
        else
          part = s.acc[st][:on_hand][sku]
          if part > 0
            s.acc[st][:on_hand][sku] = 0; s.total_sold[sku] += part; rev = part * SKU_RETAIL[sku]
            transfer_cash(s, ACC_CUST, st, rev, -1, d, tick)
            s.sales << [d, tick, st, sku, part, rev] if s.log_enabled
            s.daily_fill_store[st][1] += part if d < 6
          end
          s.acc[st][:backorder][sku] += (qty - part); s.branches[1] += 1
        end
        s_lvl = safety_ovr > 0 ? (safety_ovr * s.cfg['store_s_base']) / s.cfg['store_S_base'] : s.cfg['store_s_base']
        cap_s = safety_ovr > 0 ? safety_ovr : s.cfg['store_S_base']
        inv_pos = s.acc[st][:on_hand][sku] + s.acc[st][:in_transit][sku] - s.acc[st][:backorder][sku]
        if inv_pos < s_lvl && s.acc[st][:status] == 0
          s.branches[3] += 1; do_ship(s, STORE_WH[st], st, sku, cap_s - inv_pos, d, tick, ev[0])
        end
      when 2
        dst = ev[4]; sku = ev[5]; qty = ev[6]
        if s.acc[dst][:in_transit][sku] >= qty
          s.acc[dst][:in_transit][sku] -= qty; s.acc[dst][:on_hand][sku] += qty
        end
        s.ships << [d, tick, ev[3], dst, ev[7], sku, qty, 0, 0, 2] if s.log_enabled
      when 3
        fac = ev[3]; sku = ev[5]
        qty = s.acc[fac][:plant_state] == 1 ? 0 : (s.acc[fac][:plant_state] == 2 ? ev[6] / 2 : ev[6])
        if qty > 0
          cost = qty * SKU_PROD[sku]; s.acc[fac][:on_hand][sku] += qty; s.total_produced[sku] += qty
          transfer_cash(s, fac, ACC_SUPP, cost, -1, d, tick)
          s.daily_cost_echelon[d][0] += cost if d < 6
          s.prods << [d, tick, fac, sku, qty, cost] if s.log_enabled
          do_ship(s, fac, ev[4], sku, qty, d, tick, ev[0])
        end
      when 4
        m = s.acc[13][:on_hand][0..2].min; mx = s.acc[13][:on_hand][0..2].max; imb = mx - m
        s.kit_imb_series << imb if s.kit_imb_series.length < 64
        q = ev[6]
        if m >= q; s.branches[8] += 1; else q = m > 0 ? m : 0; s.branches[9] += 1; end
        if q > 0
          3.times { |k| s.acc[13][:on_hand][k] -= q }
          s.acc[13][:on_hand][3] += q; s.total_kits += q; fee = q * s.cfg['kit_fee']
          transfer_cash(s, 13, ACC_FEE, fee, -1, d, tick)
          s.daily_cost_echelon[d][1] += fee if d < 6
          s.kits << [d, tick, 13, q, fee, imb] if s.log_enabled
        end
      when 5
        s.branches[10] += 1; do_ship(s, ev[3], ev[4], ev[5], ev[6], d, tick, ev[0])
        s.f_trans.puts "#{d},#{tick},#{ev[3]},#{ev[4]},#{ev[5]},#{ev[6]}" if s.f_trans
      when 7
        s.acc[ev[3]][:plant_state] = 1; s.acc[ev[3]][:disrupt_rem] = 1; s.branches[11] += 1
        s.f_disrupt.puts "#{d},#{tick},#{ev[3]},1,1" if s.f_disrupt
      when 6
        s.branches[12] += 1; s.shock_cnt[ev[7]] += 1; s.shock_resp[ev[7]] += ev[6]
        s.f_shocks.puts "#{d},#{tick},#{ev[7]},#{ev[4]},#{ev[6]}" if s.f_shocks
        s.heap.push([ev[0] + 1, s.next_seq, 0, ACC_CUST, ev[4], 0, ev[6], 0]); s.next_seq += 1
      when 8
        sim_risk_check(s, d, tick)
      end
      check_invariants(s) if ev_proc % s.cfg['k_chk'] == 0
    end
    sim_eod_clearing(s, d); check_invariants(s)
  end
end
