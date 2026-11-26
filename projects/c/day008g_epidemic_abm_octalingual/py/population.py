from dataclasses import dataclass
from rng import MASK64, FNV_OFFSET, splitmix64, fnv1a_u32pair, Xoshiro256

BP = 10000
MAX_DAYS = 200
ST_S, ST_E, ST_IA, ST_IS, ST_H, ST_ICU, ST_R, ST_D = 0, 1, 2, 3, 4, 5, 6, 7


@dataclass
class Person:
    age_grp: int = 0
    hh_id: int = 0
    sw_id: int = -1
    status: int = 0
    variant: int = 0
    state_day: int = 0
    dur_days: int = 0
    inf_day: int = -1
    vax_doses: int = 0
    vax_day: int = -1
    compliance: int = 0
    quar_until: int = -1
    super_flag: int = 0
    offspring: int = 0
    ever_inf: int = 0
    hosp_rej: int = 0


@dataclass
class RunParams:
    n: int
    days: int
    beta_bp: int
    compliance_ovr: int = -1
    test_ovr: int = -1
    vax_day_ovr: int = -1
    vax_cov_ovr: int = -1
    school_ovr: int = -1
    lock_ovr: int = -1
    var_day_ovr: int = -1
    mask_ovr: int = -1
    contact_scale_bp: int = 10000
    vax_perfect: int = 0
    interv_on: int = 0


def load_config(path: str) -> dict:
    c = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                c[k.strip()] = int(v.strip())
    return c


class Sim:
    def __init__(self, cfg: dict, rp: RunParams, log_en: bool):
        self.cfg = dict(cfg)
        self.rp = RunParams(**vars(rp))
        self.n = rp.n
        self.hosp_cap = max((self.n * cfg["hosp_beds_per_10k"]) // BP, 1)
        self.icu_cap = max((self.n * cfg["icu_beds_per_10k"]) // BP, 1)
        self.hosp_occ = 0
        self.icu_occ = 0
        self.rng = []
        for i in range(5):
            x = (cfg["master_seed"] & MASK64) ^ fnv1a_u32pair(100, i)
            _, sd = splitmix64(x)
            self.rng.append(Xoshiro256(sd))

        n = self.n
        self.pop = [Person() for _ in range(n)]
        self.seen_day = [-1] * n
        self.age_tot = [0, 0, 0]

        hh_cnt = []
        pid, hh = 0, 0
        while pid < n:
            sz = min(self.rng[0].rand_int(1, 6), n - pid)
            hh_cnt.append(sz)
            for _ in range(sz):
                p = self.pop[pid]
                p.hh_id = hh
                u = self.rng[0].rand_int(1, 100)
                ag = 0 if u <= 22 else (1 if u <= 82 else 2)
                p.age_grp = ag
                self.age_tot[ag] += 1
                p.compliance = self.rng[0].rand_int(4000, 9500)
                p.super_flag = 1 if self.rng[0].rand_int(1, 10) == 1 else 0
                pid += 1
            hh += 1
        self.n_hh = hh
        self.hh_off = [0] * (hh + 1)
        for i in range(hh):
            self.hh_off[i + 1] = self.hh_off[i] + hh_cnt[i]
        self.hh_mem = [0] * n
        cur = [0] * hh
        for i in range(n):
            h = self.pop[i].hh_id
            self.hh_mem[self.hh_off[h] + cur[h]] = i
            cur[h] += 1

        n_schools = (self.age_tot[0] // 25) + 1
        n_works = (self.age_tot[1] // 15) + 1
        self.n_sw = n_schools + n_works
        sw_cnt = [0] * self.n_sw
        for i in range(n):
            ag = self.pop[i].age_grp
            if ag == 0:
                sw = self.rng[0].rand_int(0, n_schools - 1)
                self.pop[i].sw_id = sw
                sw_cnt[sw] += 1
            elif ag == 1:
                sw = n_schools + self.rng[0].rand_int(0, n_works - 1)
                self.pop[i].sw_id = sw
                sw_cnt[sw] += 1
        self.sw_off = [0] * (self.n_sw + 1)
        for i in range(self.n_sw):
            self.sw_off[i + 1] = self.sw_off[i] + sw_cnt[i]
        self.sw_mem = [0] * n
        sw_cur = [0] * self.n_sw
        for i in range(n):
            sw = self.pop[i].sw_id
            if sw >= 0:
                self.sw_mem[self.sw_off[sw] + sw_cur[sw]] = i
                sw_cur[sw] += 1

        self.cand_tgt, self.cand_src, self.cand_lay, self.cand_var = [], [], [], []
        self.daily_cnt = [[0] * 8 for _ in range(MAX_DAYS)]
        self.daily_inc_age = [[0] * 3 for _ in range(MAX_DAYS)]
        self.daily_deaths_age = [[0] * 3 for _ in range(MAX_DAYS)]
        self.daily_hosp = [[0] * 4 for _ in range(MAX_DAYS)]
        self.daily_tests = [[0] * 2 for _ in range(MAX_DAYS)]
        self.daily_quar = [0] * MAX_DAYS
        self.daily_vax = [[0] * 2 for _ in range(MAX_DAYS)]
        self.daily_var = [[0] * 2 for _ in range(MAX_DAYS)]
        self.layer_inf = [0, 0, 0]
        self.age_inf = [0, 0, 0]
        self.branches = [0] * 21
        self.inv_violations = 0
        self.total_trans = 0
        self.total_reinf = 0
        self.event_hash = FNV_OFFSET
        self.log_enabled = log_en
        self.trans = []
        self.progs = []
