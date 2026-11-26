#include "common.hpp"

uint64_t splitmix64(uint64_t *x) {
    *x += 0x9e3779b97f4a7c15ULL;
    uint64_t z = *x;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

uint64_t fnv1a_bytes(const void *buf, size_t n) {
    const auto *p = static_cast<const uint8_t *>(buf);
    uint64_t h = 0xcbf29ce484222325ULL;
    for (size_t i = 0; i < n; ++i) h = (h ^ p[i]) * 0x100000001b3ULL;
    return h;
}

uint64_t fnv1a_u32pair(uint32_t a, uint32_t b) {
    uint8_t buf[8];
    for (int i = 0; i < 4; ++i) { buf[i] = static_cast<uint8_t>(a >> (8 * i)); buf[4 + i] = static_cast<uint8_t>(b >> (8 * i)); }
    return fnv1a_bytes(buf, 8);
}

uint64_t fnv1a_i64(uint64_t h, int64_t v) {
    auto uv = static_cast<uint64_t>(v);
    for (int i = 0; i < 8; ++i) h = (h ^ ((uv >> (8 * i)) & 0xffULL)) * 0x100000001b3ULL;
    return h;
}

static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

Xoshiro256 xoshiro_seed(uint64_t sd) {
    Xoshiro256 st{}; uint64_t x = sd;
    for (int i = 0; i < 4; ++i) st.s[i] = splitmix64(&x);
    return st;
}

uint64_t xoshiro_next(Xoshiro256 *st) {
    uint64_t res = rotl64(st->s[0] + st->s[3], 23) + st->s[0];
    uint64_t t = st->s[1] << 17;
    st->s[2] ^= st->s[0]; st->s[3] ^= st->s[1]; st->s[1] ^= st->s[2]; st->s[0] ^= st->s[3];
    st->s[2] ^= t; st->s[3] = rotl64(st->s[3], 45);
    return res;
}

int64_t rand_int(Xoshiro256 *st, int64_t lo, int64_t hi) {
    if (hi <= lo) return lo;
    return lo + static_cast<int64_t>((xoshiro_next(st) >> 11) % static_cast<uint64_t>(hi - lo + 1));
}

int32_t bernoulli(Xoshiro256 *st, int64_t p_bp) {
    if (p_bp <= 0) return 0;
    if (p_bp >= BP) return 1;
    return static_cast<int32_t>(((xoshiro_next(st) >> 11) % static_cast<uint64_t>(BP)) < static_cast<uint64_t>(p_bp));
}
