from inventory import NUM_STORES, NUM_WAREHOUSES, ACC_FEE, SKU_WHOLESALE, transfer_cash

def sim_risk_check(s, day, tick):
    for i in range(NUM_STORES + NUM_WAREHOUSES):
        bo_val = sum(s.acc[i]["backorder"][k] * SKU_WHOLESALE[k] for k in range(4))
        req_margin = (bo_val * s.cfg["credit_margin_bps"]) // 10000
        if s.acc[i]["cash"] < req_margin and s.acc[i]["status"] == 0:
            s.acc[i]["status"] = 1; s.branches[13] += 1
            fc = s.files.get("credit_calls.csv")
            if fc: fc.write(f"{day},{i},{s.acc[i]['cash']},{req_margin},1\n")
        elif s.acc[i]["status"] == 1 and s.acc[i]["cash"] >= req_margin:
            s.acc[i]["status"] = 0
        if s.acc[i]["cash"] < 0 and s.acc[i]["status"] != 2:
            s.acc[i]["status"] = 2; s.branches[14] += 1
            transfer_cash(s, i, ACC_FEE, 50000, 2, day, tick)
            fr = s.files.get("risk_actions.csv")
            if fr: fr.write(f"{day},{tick},{i},LIQUIDATE,50000\n")

def sim_eod_clearing(s, day):
    for i in range(NUM_STORES + NUM_WAREHOUSES):
        oh_tot = sum(s.acc[i]["on_hand"]); it_tot = sum(s.acc[i]["in_transit"]); bo_tot = sum(s.acc[i]["backorder"])
        if day < 6:
            for k in range(4): s.day_inv_snap[day][k] += s.acc[i]["on_hand"][k]
        hold = oh_tot * s.cfg["holding_cost_per_day"]; pen = bo_tot * s.cfg["backorder_penalty"]
        if hold > 0:
            transfer_cash(s, i, ACC_FEE, hold, 0, day, s.cfg["ticks_per_day"])
            if day < 6: s.daily_cost_echelon[day][3] += hold
        if pen > 0:
            transfer_cash(s, i, ACC_FEE, pen, 1, day, s.cfg["ticks_per_day"])
            if day < 6: s.daily_cost_echelon[day][3] += pen
        fa = s.files.get("accounts_series.csv")
        if fa: fa.write(f"{day},{i},{s.acc[i]['cash']},{oh_tot},{it_tot},{s.acc[i]['status']}\n")
    sim_risk_check(s, day, s.cfg["ticks_per_day"])
