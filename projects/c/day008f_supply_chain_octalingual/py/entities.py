from rng import M64, splitmix64, fnv1a_u32pair, Xoshiro256, MinHeap
from inventory import NUM_ACCOUNTS, NUM_STORES, ACC_CUST, ACC_SUPP, ACC_FEE

def load_config(path):
    c = {}
    with open(path) as f:
        for ln in f:
            if "=" in ln:
                k, v = ln.strip().split("=", 1)
                c[k] = int(v)
    return c

class Sim:
    def __init__(self, c, log_en, out_dir):
        self.cfg = c; self.log_enabled = log_en; self.next_seq = 1; self.event_hash = 0xCBF29CE484222325
        self.heap = MinHeap()
        self.rng = []
        for i in range(NUM_ACCOUNTS + 2):
            x = (c["master_seed"] & M64) ^ fnv1a_u32pair(i // 10, i)
            sd, _ = splitmix64(x)
            self.rng.append(Xoshiro256(sd))
        self.acc = [{"cash": 0, "on_hand": [0,0,0,0], "in_transit": [0,0,0,0], "backorder": [0,0,0,0],
                     "credit_limit": 0, "status": 0, "ewma_milli": [0,0,0,0], "plant_state": 0, "disrupt_rem": 0}
                    for _ in range(NUM_ACCOUNTS)]
        for i in range(NUM_STORES):
            self.acc[i]["cash"] = 2_000_000_000; self.acc[i]["credit_limit"] = 500_000_000
            self.acc[i]["on_hand"] = [400,400,400,120]; self.acc[i]["ewma_milli"] = [12000]*4
        for i in range(10, 14):
            self.acc[i]["cash"] = 5_000_000_000; self.acc[i]["credit_limit"] = 1_000_000_000
            self.acc[i]["on_hand"] = [2500,2500,2500,800]
        for i in range(14, 17):
            self.acc[i]["cash"] = 5_000_000_000; self.acc[i]["credit_limit"] = 2_000_000_000
            self.acc[i]["on_hand"][i - 14] = 4000
        for i in range(17, 20):
            self.acc[i]["cash"] = 500_000_000; self.acc[i]["credit_limit"] = 200_000_000
        self.acc[ACC_CUST]["cash"] = 50_000_000_000; self.acc[ACC_SUPP]["cash"] = 50_000_000_000; self.acc[ACC_FEE]["cash"] = 1_000_000_000
        self.init_cash = sum(a["cash"] for a in self.acc)
        self.init_stock = [sum(a["on_hand"][k] for a in self.acc) for k in range(4)]
        self.total_produced = [0,0,0,0]; self.total_sold = [0,0,0,0]; self.total_kits = 0; self.branches = [0]*16
        self.carrier_day_vol = [0,0,0]; self.carrier_tier_cnt = [0,0,0]
        self.alloc_total_req = self.alloc_total_fill = self.alloc_cases = self.alloc_nonzero = self.inv_violations = 0
        self.daily_demand = [[0]*4 for _ in range(6)]; self.daily_cost_echelon = [[0]*4 for _ in range(6)]
        self.daily_fill_store = [[0,0] for _ in range(10)]; self.shock_cnt = [0,0,0]; self.shock_resp = [0,0,0]
        self.day_inv_snap = [[0]*4 for _ in range(6)]; self.kit_imb_series = []; self.leadtime_pairs = []
        self.sales = []; self.ships = []; self.prods = []; self.kits = []; self.fees = []
        self.files = {}
        if log_en and out_dir:
            hdrs = {"orders_log.csv": "day,tick,seq,etype,src,dst,sku,qty", "allocations.csv": "day,wh,sku,store,req,fill",
                    "disruptions.csv": "day,tick,fac,state,dur", "credit_calls.csv": "day,acc,cash,req_margin,status",
                    "risk_actions.csv": "day,tick,acc,action,amount", "transshipments.csv": "day,tick,src,dst,sku,qty",
                    "shocks_log.csv": "day,tick,cls,store,mag", "accounts_series.csv": "day,acc,cash,oh_tot,it_tot,status"}
            for fn, h in hdrs.items():
                fp = open(f"{out_dir}/{fn}", "w"); fp.write(h + "\n"); self.files[fn] = fp

    def close_logs(self):
        for fp in self.files.values(): fp.close()
        self.files.clear()
