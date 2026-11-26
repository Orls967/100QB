def sim_risk_check(s, day, tick)
  (NUM_STORES + NUM_WAREHOUSES).times do |i|
    bo_val = (0..3).sum { |k| s.acc[i][:backorder][k] * SKU_WHOLESALE[k] }
    req_margin = (bo_val * s.cfg['credit_margin_bps']) / 10000
    if s.acc[i][:cash] < req_margin && s.acc[i][:status] == 0
      s.acc[i][:status] = 1; s.branches[13] += 1
      s.f_credit.puts "#{day},#{i},#{s.acc[i][:cash]},#{req_margin},1" if s.f_credit
    elsif s.acc[i][:status] == 1 && s.acc[i][:cash] >= req_margin
      s.acc[i][:status] = 0
    end
    if s.acc[i][:cash] < 0 && s.acc[i][:status] != 2
      s.acc[i][:status] = 2; s.branches[14] += 1
      transfer_cash(s, i, ACC_FEE, 50000, 2, day, tick)
      s.f_risk.puts "#{day},#{tick},#{i},LIQUIDATE,50000" if s.f_risk
    end
  end
end

def sim_eod_clearing(s, day)
  (NUM_STORES + NUM_WAREHOUSES).times do |i|
    oh_tot = s.acc[i][:on_hand].sum
    it_tot = s.acc[i][:in_transit].sum
    bo_tot = s.acc[i][:backorder].sum
    4.times { |k| s.day_inv_snap[day][k] += s.acc[i][:on_hand][k] } if day < 6
    hold = oh_tot * s.cfg['holding_cost_per_day']; pen = bo_tot * s.cfg['backorder_penalty']
    if hold > 0; transfer_cash(s, i, ACC_FEE, hold, 0, day, s.cfg['ticks_per_day']); s.daily_cost_echelon[day][3] += hold if day < 6; end
    if pen > 0; transfer_cash(s, i, ACC_FEE, pen, 1, day, s.cfg['ticks_per_day']); s.daily_cost_echelon[day][3] += pen if day < 6; end
    s.f_acc_ser.puts "#{day},#{i},#{s.acc[i][:cash]},#{oh_tot},#{it_tot},#{s.acc[i][:status]}" if s.f_acc_ser
  end
  sim_risk_check(s, day, s.cfg['ticks_per_day'])
end
