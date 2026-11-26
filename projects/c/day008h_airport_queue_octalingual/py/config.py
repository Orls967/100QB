from dataclasses import dataclass, field

SLOTS_PER_DAY = 8640
SLOTS_PER_HOUR = 360
N_STN = 6
MAX_SRV = 16
MAX_Q = 4096
MAX_POOL = 8192
AUDIT_DAYS = 3

STN_CHECKIN = 0
STN_SECURITY = 1
STN_IMIGRASI = 2
STN_GATE = 3
STN_BELT = 4
STN_CUSTOMS = 5

CLS_ECON = 0
CLS_BIZ = 1
CLS_TRANSIT = 2

DEST_DOM = 0
DEST_INTL = 1

SRV_IDLE = 0
SRV_BUSY = 1
SRV_DOWN = 2


@dataclass
class RunParams:
    days: int
    arr_scale_bp: int
    break_prob_bp: int
    repair_prob_bp: int
    belt_cap: int
    gate_close_lead: int
    srv_ovr: list[int] = field(default_factory=lambda: [-1] * N_STN)
    mu_ovr: list[int] = field(default_factory=lambda: [-1] * N_STN)
    priority_on: int = 1
    breakdown_on: int = 1
    roster_on: int = 1
    customs_on: int = 0
    permute_srv: int = 0
    timescale_div2: int = 0
    audit_on: int = 0


def default_rp(c: dict, days: int) -> RunParams:
    return RunParams(
        days=days,
        arr_scale_bp=c["arr_scale_bp"],
        break_prob_bp=c["break_prob_bp"],
        repair_prob_bp=c["repair_prob_bp"],
        belt_cap=c["belt_cap"],
        gate_close_lead=c["gate_close_lead"],
    )


def load_config(path: str) -> dict:
    c: dict = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            c[k.strip()] = int(v.strip())
    c["srv_cnt"] = [
        c["checkin_servers"],
        c["security_servers"],
        c["imigrasi_servers"],
        c["gate_servers"],
        c["belt_servers"],
        c["customs_servers"],
    ]
    c["mu_bp"] = [
        c["checkin_mu_bp"],
        c["security_mu_bp"],
        c["imigrasi_mu_bp"],
        c["gate_mu_bp"],
        c["belt_mu_bp"],
        c["customs_mu_bp"],
    ]
    return c
