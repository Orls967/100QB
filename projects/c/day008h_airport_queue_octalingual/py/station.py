from dataclasses import dataclass
from config import CLS_BIZ, MAX_Q, MAX_SRV


@dataclass
class Server:
    status: int = 0
    enabled: int = 0
    agent_idx: int = -1
    extra_sec: int = 0


class Station:
    def __init__(self, c_max: int, cap: int, has_prio: int) -> None:
        self.q = [0] * MAX_Q
        self.head = 0
        self.tail = 0
        self.len = 0
        self.cap = cap
        self.has_prio = has_prio
        self.c_max = c_max
        self.c_active = c_max
        self.busy_cnt = 0
        self.srv = [Server(0, 1 if i < c_max else 0, -1, 0) for i in range(MAX_SRV)]
        self.arr_cnt = 0
        self.dep_cnt = 0
        self.block_cnt = 0
        self.aban_cnt = 0
        self.sum_l = 0
        self.sum_w = 0
        self.wait_cnt = 0
        self.busy_slots = 0

    def push_tail(self, idx: int) -> None:
        if self.len >= MAX_Q:
            return
        self.q[self.tail] = idx
        self.tail = (self.tail + 1) % MAX_Q
        self.len += 1

    def push_head(self, idx: int) -> None:
        if self.len >= MAX_Q:
            return
        self.head = (self.head - 1 + MAX_Q) % MAX_Q
        self.q[self.head] = idx
        self.len += 1

    def remove_at(self, pos: int) -> int:
        idx = self.q[(self.head + pos) % MAX_Q]
        for i in range(pos, self.len - 1):
            self.q[(self.head + i) % MAX_Q] = self.q[(self.head + i + 1) % MAX_Q]
        self.tail = (self.tail - 1 + MAX_Q) % MAX_Q
        self.len -= 1
        return idx

    def pop_next(self, pool: list, priority_on: int, prio_br: int) -> tuple[int, int]:
        if self.len == 0:
            return -1, prio_br
        if priority_on != 0 and self.has_prio != 0:
            for i in range(self.len):
                if pool[self.q[(self.head + i) % MAX_Q]].cls == CLS_BIZ:
                    if i > 0:
                        return self.remove_at(i), prio_br + 1
                    break
        res = self.q[self.head]
        self.head = (self.head + 1) % MAX_Q
        self.len -= 1
        return res, prio_br
