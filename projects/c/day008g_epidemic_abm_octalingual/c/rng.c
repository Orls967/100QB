#include "common.h"

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

uint64_t fnv1a_bytes(const void *buf, size_t n) {
    const uint8_t *p = (const uint8_t *)buf;
    uint64_t h = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < n; ++i) { h ^= (uint64_t)p[i]; h *= 0x100000001b3ULL; }
    return h;
}

uint64_t fnv1a_u32pair(uint32_t a, uint32_t b) {
    uint8_t buf[8];
    for (int i = 0; i < 4; ++i) { buf[i] = (uint8_t)((a >> (8 * i)) & 0xffU); buf[4 + i] = (uint8_t)((b >> (8 * i)) & 0xffU); }
    return fnv1a_bytes(buf, 8);
}

uint64_t fnv1a_i64(uint64_t h, int64_t v) {
    uint64_t uv = (uint64_t)v;
    for (int i = 0; i < 8; ++i) { h ^= (uv >> (8 * i)) & 0xffULL; h *= 0x100000001b3ULL; }
    return h;
}

void x256_seed(Xoshiro256 *st, uint64_t seed) {
    uint64_t x = seed;
    for (int i = 0; i < 4; ++i) st->s[i] = splitmix64(&x);
}

uint64_t x256_next(Xoshiro256 *st) {
    uint64_t res = rotl64(st->s[0] + st->s[3], 23) + st->s[0];
    uint64_t t = st->s[1] << 17;
    st->s[2] ^= st->s[0]; st->s[3] ^= st->s[1]; st->s[1] ^= st->s[2]; st->s[0] ^= st->s[3];
    st->s[2] ^= t; st->s[3] = rotl64(st->s[3], 45);
    return res;
}

int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi) {
    if (hi <= lo) return lo;
    return lo + (int64_t)((x256_next(st) >> 11) % (uint64_t)(hi - lo + 1));
}
