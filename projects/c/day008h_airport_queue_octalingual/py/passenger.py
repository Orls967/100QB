from dataclasses import dataclass
from config import MAX_POOL


@dataclass
class Agent:
    active: int = 0
    cls: int = 0
    dest: int = 0
    flight_hr: int = 0
    n_bags: int = 0
    cur_stn: int = 0
    extra_sec: int = 0
    is_bag: int = 0
    arr_slot: int = 0
    q_enter_slot: int = 0
    dep_slot: int = 0


def pool_alloc(s) -> int:
    if not s.free_stk:
        return -1
    idx = s.free_stk.pop()
    s.pool[idx] = Agent(active=1)
    return idx


def pool_free(s, idx: int) -> None:
    if idx < 0 or idx >= MAX_POOL or s.pool[idx].active == 0:
        return
    s.pool[idx].active = 0
    s.free_stk.append(idx)
