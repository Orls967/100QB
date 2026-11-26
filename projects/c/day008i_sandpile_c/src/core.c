#include "sandpile.h"
#include <stdlib.h>
#include <string.h>

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

void rng_seed(Rng *r, uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
}

uint64_t rng_u64(Rng *r) {
    uint64_t res = rotl64(r->s[0] + r->s[3], 23) + r->s[0];
    uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0];
    r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2];
    r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl64(r->s[3], 45);
    return res;
}

int32_t rng_int(Rng *r, int32_t lo, int32_t hi) {
    if (hi <= lo) return lo;
    return lo + (int32_t)((rng_u64(r) >> 11) % (uint64_t)(hi - lo + 1));
}

uint64_t fnv1a_init(void) { return 0xcbf29ce484222325ULL; }

uint64_t fnv1a_bytes(uint64_t h, const void *buf, size_t len) {
    const uint8_t *p = (const uint8_t *)buf;
    for (size_t i = 0; i < len; i++) h = (h ^ p[i]) * 0x100000001b3ULL;
    return h;
}

uint64_t grid_state_hash(const Grid *g) {
    size_t m = (size_t)g->n * (size_t)g->n;
    uint64_t h = fnv1a_init();
    h = fnv1a_bytes(h, g->z, m * sizeof(int32_t));
    h = fnv1a_bytes(h, g->odom, m * sizeof(int64_t));
    h = fnv1a_bytes(h, &g->sink, sizeof(int64_t));
    return h;
}

void grid_init(Grid *g, int n) {
    size_t m = (size_t)n * (size_t)n;
    g->n = n;
    g->z = (int32_t *)calloc(m, sizeof(int32_t));
    g->odom = (int64_t *)calloc(m, sizeof(int64_t));
    g->q_a = (int32_t *)malloc((m + 1) * sizeof(int32_t));
    g->q_b = (int32_t *)malloc((m + 1) * sizeof(int32_t));
    g->stamp_a = (uint64_t *)calloc(m, sizeof(uint64_t));
    g->stamp_b = (uint64_t *)calloc(m, sizeof(uint64_t));
    g->sink = 0;
    g->cur_stamp = 1;
}

void grid_clear(Grid *g) {
    size_t m = (size_t)g->n * (size_t)g->n;
    memset(g->z, 0, m * sizeof(int32_t));
    memset(g->odom, 0, m * sizeof(int64_t));
    g->sink = 0;
}

void grid_free(Grid *g) {
    free(g->z);
    free(g->odom);
    free(g->q_a);
    free(g->q_b);
    free(g->stamp_a);
    free(g->stamp_b);
}

void stabilize_det(Grid *g) {
    int n = g->n, m = n * n, cap = m + 1, head = 0, tail = 0;
    uint64_t st = ++g->cur_stamp;
    for (int i = 0; i < m; i++) {
        if (g->z[i] >= 4) {
            g->q_a[tail++] = i;
            g->stamp_a[i] = st;
        }
    }
    while (head != tail) {
        int i = g->q_a[head++];
        if (head == cap) head = 0;
        g->stamp_a[i] = 0;
        int32_t k = g->z[i] >> 2;
        if (k <= 0) continue;
        g->z[i] &= 3;
        g->odom[i] += k;
        int r = i / n, c = i % n;
        int nb[4] = {r > 0 ? i - n : -1, r + 1 < n ? i + n : -1,
                     c > 0 ? i - 1 : -1, c + 1 < n ? i + 1 : -1};
        for (int d = 0; d < 4; d++) {
            int j = nb[d];
            if (j < 0) {
                g->sink += k;
            } else {
                g->z[j] += k;
                if (g->z[j] >= 4 && g->stamp_a[j] != st) {
                    g->stamp_a[j] = st;
                    g->q_a[tail++] = j;
                    if (tail == cap) tail = 0;
                }
            }
        }
    }
}

Aval stabilize_waves(Grid *g) {
    Aval av = {0, 0, 0};
    int n = g->n, m = n * n, n_cur = 0;
    uint64_t area_st = ++g->cur_stamp;
    for (int i = 0; i < m; i++) {
        if (g->z[i] >= 4) g->q_a[n_cur++] = i;
    }
    int32_t *cur = g->q_a, *nxt = g->q_b;
    while (n_cur > 0) {
        av.dur++;
        int n_nxt = 0;
        uint64_t wave_st = ++g->cur_stamp;
        for (int u = 0; u < n_cur; u++) {
            int i = cur[u];
            g->z[i] -= 4;
            g->odom[i]++;
            av.size++;
            if (g->stamp_b[i] != area_st) {
                g->stamp_b[i] = area_st;
                av.area++;
            }
        }
        for (int u = 0; u < n_cur; u++) {
            int i = cur[u];
            if (g->z[i] >= 4 && g->stamp_a[i] != wave_st) {
                g->stamp_a[i] = wave_st;
                nxt[n_nxt++] = i;
            }
            int r = i / n, c = i % n;
            int nb[4] = {r > 0 ? i - n : -1, r + 1 < n ? i + n : -1,
                         c > 0 ? i - 1 : -1, c + 1 < n ? i + 1 : -1};
            for (int d = 0; d < 4; d++) {
                int j = nb[d];
                if (j < 0) {
                    g->sink++;
                } else {
                    g->z[j]++;
                    if (g->z[j] >= 4 && g->stamp_a[j] != wave_st) {
                        g->stamp_a[j] = wave_st;
                        nxt[n_nxt++] = j;
                    }
                }
            }
        }
        int32_t *tmp = cur; cur = nxt; nxt = tmp;
        n_cur = n_nxt;
    }
    return av;
}

void stabilize_rand(Grid *g, Rng *rng) {
    int n = g->n, m = n * n, cnt = 0;
    uint64_t st = ++g->cur_stamp;
    for (int i = 0; i < m; i++) {
        if (g->z[i] >= 4) {
            g->q_a[cnt++] = i;
            g->stamp_a[i] = st;
        }
    }
    while (cnt > 0) {
        int pick = rng_int(rng, 0, cnt - 1);
        int i = g->q_a[pick];
        g->q_a[pick] = g->q_a[--cnt];
        g->stamp_a[i] = 0;
        g->z[i] -= 4;
        g->odom[i]++;
        if (g->z[i] >= 4) {
            g->stamp_a[i] = st;
            g->q_a[cnt++] = i;
        }
        int r = i / n, c = i % n;
        int nb[4] = {r > 0 ? i - n : -1, r + 1 < n ? i + n : -1,
                     c > 0 ? i - 1 : -1, c + 1 < n ? i + 1 : -1};
        for (int d = 0; d < 4; d++) {
            int j = nb[d];
            if (j < 0) {
                g->sink++;
            } else {
                g->z[j]++;
                if (g->z[j] >= 4 && g->stamp_a[j] != st) {
                    g->stamp_a[j] = st;
                    g->q_a[cnt++] = j;
                }
            }
        }
    }
}

int is_recurrent(const Grid *g, int32_t *stk, uint8_t *deg) {
    int n = g->n, m = n * n, top = 0, burnt = 0;
    for (int i = 0; i < m; i++) {
        int r = i / n, c = i % n;
        uint8_t b = (uint8_t)((r == 0) + (r == n - 1) + (c == 0) + (c == n - 1));
        if (g->z[i] + b >= 4) {
            deg[i] = 255;
            stk[top++] = i;
            burnt++;
        } else {
            deg[i] = b;
        }
    }
    while (top > 0) {
        int i = stk[--top];
        int r = i / n, c = i % n;
        int nb[4] = {r > 0 ? i - n : -1, r + 1 < n ? i + n : -1,
                     c > 0 ? i - 1 : -1, c + 1 < n ? i + 1 : -1};
        for (int d = 0; d < 4; d++) {
            int j = nb[d];
            if (j >= 0 && deg[j] != 255) {
                deg[j]++;
                if (g->z[j] + deg[j] >= 4) {
                    deg[j] = 255;
                    stk[top++] = j;
                    burnt++;
                }
            }
        }
    }
    return burnt == m;
}
