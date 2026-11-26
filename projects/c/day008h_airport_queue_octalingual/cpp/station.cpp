#include "common.hpp"

void stn_init(Station &st, int c_max, int cap, int has_prio) {
    st = Station{};
    st.c_max = c_max; st.c_active = c_max; st.cap = cap; st.has_prio = has_prio;
    for (int i = 0; i < c_max; ++i) { st.srv[i].enabled = 1; st.srv[i].agent_idx = -1; }
}

void stn_push_tail(Station &st, int32_t idx) {
    if (st.len >= MAX_Q) return;
    st.q[st.tail] = idx; st.tail = (st.tail + 1) % MAX_Q; st.len++;
}

void stn_push_head(Station &st, int32_t idx) {
    if (st.len >= MAX_Q) return;
    st.head = (st.head - 1 + MAX_Q) % MAX_Q; st.q[st.head] = idx; st.len++;
}

int32_t stn_remove_at(Station &st, int pos) {
    int32_t idx = st.q[(st.head + pos) % MAX_Q];
    for (int i = pos; i < st.len - 1; ++i) st.q[(st.head + i) % MAX_Q] = st.q[(st.head + i + 1) % MAX_Q];
    st.tail = (st.tail - 1 + MAX_Q) % MAX_Q; st.len--;
    return idx;
}

int32_t stn_pop_next(Station &st, const Agent *pool, int32_t priority_on, int64_t *prio_br) {
    if (st.len == 0) return -1;
    if (priority_on && st.has_prio) {
        int biz_pos = -1;
        for (int i = 0; i < st.len; ++i) {
            if (pool[st.q[(st.head + i) % MAX_Q]].cls == CLS_BIZ) { biz_pos = i; break; }
        }
        if (biz_pos > 0) { if (prio_br) (*prio_br)++; return stn_remove_at(st, biz_pos); }
    }
    int32_t res = st.q[st.head]; st.head = (st.head + 1) % MAX_Q; st.len--;
    return res;
}
