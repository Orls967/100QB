NUM_SKUS = 4
NUM_STORES = 10
NUM_WAREHOUSES = 4
NUM_CARRIERS = 3
NUM_ACCOUNTS = 23
ACC_CUST = 20
ACC_SUPP = 21
ACC_FEE = 22

SKU_RETAIL = (10000, 10200, 9800, 31500)
SKU_WHOLESALE = (8000, 8200, 7800, 26000)
SKU_PROD = (6000, 6200, 5800, 0)
STORE_WH = (10, 10, 10, 11, 11, 11, 12, 12, 12, 12)
CARRIER_BASE_LT = (120, 200, 160)

def carrier_rate(v):
    return 200 if v < 500 else (160 if v < 2000 else 120)

def transfer_cash(s, from_id, to_id, amt, kind, day, tick):
    if amt <= 0 or from_id == to_id: return
    s.acc[from_id]["cash"] -= amt
    s.acc[to_id]["cash"] += amt
    if kind >= 0: s.fees.append((day, tick, from_id, to_id, amt, kind))

def check_invariants(s):
    if sum(a["cash"] for a in s.acc) != s.init_cash:
        s.inv_violations += 1; return 0
    for k in range(4):
        oh = it = 0
        for a in s.acc:
            if a["on_hand"][k] < 0 or a["in_transit"][k] < 0 or a["backorder"][k] < 0:
                s.inv_violations += 1; return 0
            oh += a["on_hand"][k]; it += a["in_transit"][k]
        lhs = oh + it + s.total_sold[k] + (s.total_kits if k < 3 else 0)
        rhs = s.init_stock[k] + (s.total_produced[k] if k < 3 else s.total_kits)
        if lhs != rhs:
            s.inv_violations += 1; return 0
    return 1
