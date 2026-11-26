#include "common.h"

int32_t pool_alloc(Sim *s) {
    if (s->free_top <= 0) return -1;
    int32_t idx = s->free_stk[--s->free_top];
    memset(&s->pool[idx], 0, sizeof(Agent));
    s->pool[idx].active = 1;
    return idx;
}

void pool_free(Sim *s, int32_t idx) {
    if (idx < 0 || idx >= MAX_POOL || !s->pool[idx].active) return;
    s->pool[idx].active = 0;
    s->free_stk[s->free_top++] = idx;
}
