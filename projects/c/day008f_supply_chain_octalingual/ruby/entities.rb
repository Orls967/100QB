def load_config(path)
  c = {}
  File.readlines(path).each do |ln|
    k, v = ln.strip.split('=', 2)
    c[k] = v.to_i if v
  end
  c
end

class Sim
  attr_accessor :cfg, :acc, :rng, :heap, :next_seq, :init_cash, :init_stock,
                :total_produced, :total_sold, :total_kits, :branches, :carrier_day_vol,
                :carrier_tier_cnt, :alloc_total_req, :alloc_total_fill, :alloc_cases,
                :alloc_nonzero, :inv_violations, :log_enabled, :daily_demand,
                :daily_cost_echelon, :daily_fill_store, :shock_cnt, :shock_resp,
                :day_inv_snap, :kit_imb_series, :leadtime_pairs, :sales, :ships,
                :prods, :kits, :fees, :event_hash, :f_orders, :f_alloc, :f_disrupt,
                :f_credit, :f_risk, :f_trans, :f_shocks, :f_acc_ser,
                :bin_demand, :bin_orders, :bin_store_dem, :shock_evs, :n_transships,
                :dem_events_tot, :dem_bo_events, :dem_units_tot, :dem_units_filled,
                :carrier_cap_tot, :carrier_vol_tot, :cap_scale

  def initialize(c, log_en, out_dir, ev_per_day = 2000)
    @cfg = c; @log_enabled = log_en; @next_seq = 1; @event_hash = 0xCBF29CE484222325
    @cap_scale = [ev_per_day / 2000, 1].max
    cs = @cap_scale
    @heap = MinHeap.new
    @rng = (0...(NUM_ACCOUNTS + 2)).map do |i|
      x = (c['master_seed'] & M64) ^ fnv1a_u32pair(i / 10, i)
      sd, = splitmix64(x)
      Xoshiro256.new(sd)
    end
    @acc = Array.new(NUM_ACCOUNTS) { { cash: 0, on_hand: [0,0,0,0], in_transit: [0,0,0,0], backorder: [0,0,0,0], credit_limit: 0, status: 0, ewma_milli: [0,0,0,0], plant_state: 0, disrupt_rem: 0 } }
    10.times { |i| @acc[i][:cash] = 2_000_000_000 * cs; @acc[i][:credit_limit] = 500_000_000 * cs; @acc[i][:on_hand] = [60*cs,60*cs,60*cs,40*cs]; @acc[i][:ewma_milli] = [12000]*4 }
    (10..12).each { |i| @acc[i][:cash] = 5_000_000_000 * cs; @acc[i][:credit_limit] = 1_000_000_000 * cs; @acc[i][:on_hand] = [1600*cs,1400*cs,1200*cs,600*cs] }
    @acc[13][:cash] = 5_000_000_000 * cs; @acc[13][:credit_limit] = 1_000_000_000 * cs; @acc[13][:on_hand] = [1800*cs,1300*cs,900*cs,500*cs]
    (14..16).each { |i| @acc[i][:cash] = 5_000_000_000 * cs; @acc[i][:credit_limit] = 2_000_000_000 * cs; @acc[i][:on_hand][i - 14] = 4000 * cs }
    (17..19).each { |i| @acc[i][:cash] = 500_000_000 * cs; @acc[i][:credit_limit] = 200_000_000 * cs }
    @acc[ACC_CUST][:cash] = 50_000_000_000 * cs; @acc[ACC_SUPP][:cash] = 50_000_000_000 * cs; @acc[ACC_FEE][:cash] = 1_000_000_000 * cs
    @init_cash = @acc.sum { |a| a[:cash] }
    @init_stock = (0..3).map { |k| @acc.sum { |a| a[:on_hand][k] } }
    @total_produced = [0,0,0,0]; @total_sold = [0,0,0,0]; @total_kits = 0; @branches = Array.new(16, 0)
    @carrier_day_vol = [0,0,0]; @carrier_tier_cnt = [0,0,0]
    @alloc_total_req = 0; @alloc_total_fill = 0; @alloc_cases = 0; @alloc_nonzero = 0; @inv_violations = 0
    @daily_demand = Array.new(6) { [0,0,0,0] }; @daily_cost_echelon = Array.new(6) { [0,0,0,0] }
    @daily_fill_store = Array.new(10) { [0,0] }; @shock_cnt = [0,0,0]; @shock_resp = [0,0,0]
    @day_inv_snap = Array.new(6) { [0,0,0,0] }; @kit_imb_series = []; @leadtime_pairs = []
    @bin_demand = Array.new(4) { Array.new(100, 0) }
    @bin_orders = Array.new(4) { Array.new(100, 0) }
    @bin_store_dem = Array.new(10) { Array.new(100, 0) }
    @shock_evs = []; @n_transships = 0
    @dem_events_tot = 0; @dem_bo_events = 0; @dem_units_tot = 0; @dem_units_filled = 0
    @carrier_cap_tot = 0; @carrier_vol_tot = 0
    @sales = []; @ships = []; @prods = []; @kits = []; @fees = []
    if log_en && out_dir
      @f_orders = File.open("#{out_dir}/orders_log.csv", 'w'); @f_orders.puts 'day,tick,seq,etype,src,dst,sku,qty'
      @f_alloc = File.open("#{out_dir}/allocations.csv", 'w'); @f_alloc.puts 'day,wh,sku,store,req,fill'
      @f_disrupt = File.open("#{out_dir}/disruptions.csv", 'w'); @f_disrupt.puts 'day,tick,fac,state,dur'
      @f_credit = File.open("#{out_dir}/credit_calls.csv", 'w'); @f_credit.puts 'day,acc,cash,req_margin,status'
      @f_risk = File.open("#{out_dir}/risk_actions.csv", 'w'); @f_risk.puts 'day,tick,acc,action,amount'
      @f_trans = File.open("#{out_dir}/transshipments.csv", 'w'); @f_trans.puts 'day,tick,src,dst,sku,qty'
      @f_shocks = File.open("#{out_dir}/shocks_log.csv", 'w'); @f_shocks.puts 'day,tick,cls,store,mag'
      @f_acc_ser = File.open("#{out_dir}/accounts_series.csv", 'w'); @f_acc_ser.puts 'day,acc,cash,oh_tot,it_tot,status'
    end
  end

  def close_logs
    [@f_orders, @f_alloc, @f_disrupt, @f_credit, @f_risk, @f_trans, @f_shocks, @f_acc_ser].compact.each(&:close)
  end
end
