use crate::config::MAX_POOL;

#[derive(Clone, Default)]
pub struct Agent {
    pub active: i32, pub cls: i32, pub dest: i32, pub flight_hr: i32, pub n_bags: i32, pub cur_stn: i32, pub extra_sec: i32, pub is_bag: i32,
    pub arr_slot: i64, pub q_enter_slot: i64, pub dep_slot: i64,
}

pub fn pool_alloc(pool: &mut [Agent], free_stk: &mut Vec<i32>) -> i32 {
    if let Some(idx) = free_stk.pop() {
        pool[idx as usize] = Agent { active: 1, ..Default::default() };
        idx
    } else { -1 }
}

pub fn pool_free(pool: &mut [Agent], free_stk: &mut Vec<i32>, idx: i32) {
    if idx < 0 || idx as usize >= MAX_POOL || pool[idx as usize].active == 0 { return; }
    pool[idx as usize].active = 0;
    free_stk.push(idx);
}
