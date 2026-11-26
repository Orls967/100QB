use crate::config::{CLS_BIZ, MAX_Q, MAX_SRV};
use crate::passenger::Agent;

#[derive(Clone)]
pub struct Server { pub status: i32, pub enabled: i32, pub agent_idx: i32, pub extra_sec: i32 }

impl Default for Server {
    fn default() -> Self { Self { status: 0, enabled: 0, agent_idx: -1, extra_sec: 0 } }
}

pub struct Station {
    pub q: Vec<i32>, pub head: usize, pub tail: usize, pub len: usize, pub cap: i32, pub has_prio: i32,
    pub c_max: usize, pub c_active: usize, pub busy_cnt: i32,
    pub srv: Vec<Server>,
    pub arr_cnt: i64, pub dep_cnt: i64, pub block_cnt: i64, pub aban_cnt: i64, pub sum_l: i64, pub sum_w: i64, pub wait_cnt: i64, pub busy_slots: i64,
}

impl Station {
    pub fn new(c_max: usize, cap: i32, has_prio: i32) -> Self {
        let mut srv = vec![Server::default(); MAX_SRV];
        for i in 0..c_max { srv[i].enabled = 1; }
        Self { q: vec![0; MAX_Q], head: 0, tail: 0, len: 0, cap, has_prio, c_max, c_active: c_max, busy_cnt: 0, srv,
               arr_cnt: 0, dep_cnt: 0, block_cnt: 0, aban_cnt: 0, sum_l: 0, sum_w: 0, wait_cnt: 0, busy_slots: 0 }
    }
    pub fn push_tail(&mut self, idx: i32) {
        if self.len >= MAX_Q { return; }
        self.q[self.tail] = idx; self.tail = (self.tail + 1) % MAX_Q; self.len += 1;
    }
    pub fn push_head(&mut self, idx: i32) {
        if self.len >= MAX_Q { return; }
        self.head = (self.head + MAX_Q - 1) % MAX_Q; self.q[self.head] = idx; self.len += 1;
    }
    pub fn remove_at(&mut self, pos: usize) -> i32 {
        let idx = self.q[(self.head + pos) % MAX_Q];
        for i in pos..(self.len - 1) { self.q[(self.head + i) % MAX_Q] = self.q[(self.head + i + 1) % MAX_Q]; }
        self.tail = (self.tail + MAX_Q - 1) % MAX_Q; self.len -= 1;
        idx
    }
    pub fn pop_next(&mut self, pool: &[Agent], priority_on: i32, prio_br: &mut i64) -> i32 {
        if self.len == 0 { return -1; }
        if priority_on != 0 && self.has_prio != 0 {
            for i in 0..self.len {
                if pool[self.q[(self.head + i) % MAX_Q] as usize].cls == CLS_BIZ {
                    if i > 0 { *prio_br += 1; return self.remove_at(i); }
                    break;
                }
            }
        }
        let res = self.q[self.head]; self.head = (self.head + 1) % MAX_Q; self.len -= 1;
        res
    }
}
