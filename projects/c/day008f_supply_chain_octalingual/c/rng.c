#include "common.h"
#include <stdlib.h>
#include <string.h>

uint64_t splitmix64(uint64_t *x) {
    *x += 0x9E3779B97F4A7C15ULL;
    uint64_t z = *x;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

uint64_t fnv1a_bytes(const void *buf, size_t n) {
    const uint8_t *p = (const uint8_t *)buf;
    uint64_t h = 0xCBF29CE484222325ULL;
    for (size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 0x100000001B3ULL;
    }
    return h;
}

uint64_t fnv1a_u32pair(uint32_t a, uint32_t b) {
    uint8_t buf[8];
    for (int i = 0; i < 4; ++i) { buf[i] = (uint8_t)((a >> (8 * i)) & 0xFFU); }
    for (int i = 0; i < 4; ++i) { buf[4 + i] = (uint8_t)((b >> (8 * i)) & 0xFFU); }
    return fnv1a_bytes(buf, 8);
}

uint64_t fnv1a_i64(uint64_t h, int64_t v) {
    uint64_t u = (uint64_t)v;
    for (int i = 0; i < 8; ++i) {
        h ^= (u >> (8 * i)) & 0xFFU;
        h *= 0x100000001B3ULL;
    }
    return h;
}

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

void x256_seed(Xoshiro256 *st, uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) st->s[i] = splitmix64(&x);
}

uint64_t x256_next(Xoshiro256 *st) {
    uint64_t *s = st->s;
    const uint64_t res = rotl64(s[0] + s[3], 23) + s[0];
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl64(s[3], 45);
    return res;
}

int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi) {
    if (hi <= lo) return lo;
    uint64_t span = (uint64_t)(hi - lo + 1);
    return lo + (int64_t)((x256_next(st) >> 11) % span);
}

void heap_init(MinHeap *h) {
    h->cap = 1024;
    h->len = 0;
    h->data = (SimEvent *)malloc((size_t)h->cap * sizeof(SimEvent));
}

void heap_free(MinHeap *h) {
    free(h->data);
    h->data = NULL;
    h->len = h->cap = 0;
}

static inline int ev_lt(const SimEvent *a, const SimEvent *b) {
    if (a->time != b->time) return a->time < b->time;
    return a->seq < b->seq;
}

void heap_push(MinHeap *h, SimEvent ev) {
    if (h->len == h->cap) {
        h->cap *= 2;
        h->data = (SimEvent *)realloc(h->data, (size_t)h->cap * sizeof(SimEvent));
    }
    int i = h->len++;
    while (i > 0) {
        int p = (i - 1) >> 1;
        if (!ev_lt(&ev, &h->data[p])) break;
        h->data[i] = h->data[p];
        i = p;
    }
    h->data[i] = ev;
}

SimEvent heap_pop(MinHeap *h) {
    SimEvent top = h->data[0];
    SimEvent last = h->data[--h->len];
    int i = 0;
    while (1) {
        int l = 2 * i + 1, r = l + 1, best = i;
        if (l < h->len && ev_lt(&h->data[l], &last)) best = l;
        if (r < h->len && ev_lt(&h->data[r], best == i ? &last : &h->data[best])) best = r;
        if (best == i) break;
        h->data[i] = h->data[best];
        i = best;
    }
    if (h->len > 0) h->data[i] = last;
    return top;
}
