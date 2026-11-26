NUM_SKUS = 4
NUM_STORES = 10
NUM_WAREHOUSES = 4
NUM_CARRIERS = 3
NUM_ACCOUNTS = 23
ACC_CUST = 20
ACC_SUPP = 21
ACC_FEE = 22

SKU_RETAIL = [10000, 10200, 9800, 31500].freeze
SKU_WHOLESALE = [8000, 8200, 7800, 26000].freeze
SKU_PROD = [6000, 6200, 5800, 0].freeze
STORE_WH = [10, 10, 10, 11, 11, 11, 12, 12, 12, 12].freeze
CARRIER_BASE_LT = [150, 250, 200].freeze

def carrier_rate(v)
  v < 500 ? 200 : (v < 2000 ? 160 : 120)
end

def transfer_cash(s, from_id, to_id, amt, kind, day, tick)
  return if amt <= 0 || from_id == to_id
  s.acc[from_id][:cash] -= amt
  s.acc[to_id][:cash] += amt
  s.fees << [day, tick, from_id, to_id, amt, kind] if kind >= 0
end

def check_invariants(s)
  if s.acc.sum { |a| a[:cash] } != s.init_cash
    s.inv_violations += 1; return 0
  end
  4.times do |k|
    oh = 0; it = 0
    s.acc.each do |a|
      if a[:on_hand][k] < 0 || a[:in_transit][k] < 0 || a[:backorder][k] < 0
        s.inv_violations += 1; return 0
      end
      oh += a[:on_hand][k]; it += a[:in_transit][k]
    end
    lhs = oh + it + s.total_sold[k] + (k < 3 ? s.total_kits : 0)
    rhs = s.init_stock[k] + (k < 3 ? s.total_produced[k] : s.total_kits)
    if lhs != rhs
      s.inv_violations += 1; return 0
    end
  end
  1
end
