#include "md_rng.h"
#include <math.h>

static inline uint64_t rotl64(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t md_splitmix64(uint64_t *state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

void md_rng_seed(md_rng *rng, uint64_t seed) {
    uint64_t sm = seed;
    for (int i = 0; i < 4; ++i) {
        rng->s[i] = md_splitmix64(&sm);
    }
    /* Guard against degenerate all-zero state */
    if ((rng->s[0] | rng->s[1] | rng->s[2] | rng->s[3]) == 0ULL) {
        rng->s[0] = 0x9e3779b97f4a7c15ULL;
    }
    rng->has_spare_normal = false;
    rng->spare_normal = 0.0;
}

void md_rng_seed_stage(md_rng *rng, uint64_t master_seed, uint64_t stage_id, uint64_t sub_id) {
    uint64_t mix = master_seed;
    mix ^= (stage_id + 1ULL) * 0x9e3779b97f4a7c15ULL;
    mix = md_splitmix64(&mix);
    mix ^= (sub_id + 1ULL) * 0xbf58476d1ce4e5b9ULL;
    mix = md_splitmix64(&mix);
    md_rng_seed(rng, mix);
}

uint64_t md_rng_u64(md_rng *rng) {
    const uint64_t result = rotl64(rng->s[0] + rng->s[3], 23) + rng->s[0];
    const uint64_t t = rng->s[1] << 17;

    rng->s[2] ^= rng->s[0];
    rng->s[3] ^= rng->s[1];
    rng->s[1] ^= rng->s[2];
    rng->s[0] ^= rng->s[3];

    rng->s[2] ^= t;
    rng->s[3] = rotl64(rng->s[3], 45);
    return result;
}

double md_rng_uniform(md_rng *rng) {
    /* Top 53 bits mapped cleanly into [0, 1) */
    return (double)(md_rng_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}

double md_rng_range(md_rng *rng, double a, double b) {
    return a + (b - a) * md_rng_uniform(rng);
}

double md_rng_normal(md_rng *rng) {
    if (rng->has_spare_normal) {
        rng->has_spare_normal = false;
        return rng->spare_normal;
    }
    double u1, u2, s;
    do {
        u1 = 2.0 * md_rng_uniform(rng) - 1.0;
        u2 = 2.0 * md_rng_uniform(rng) - 1.0;
        s = u1 * u1 + u2 * u2;
    } while (s >= 1.0 || s <= 1e-300);
    double fac = sqrt(-2.0 * log(s) / s);
    rng->spare_normal = u2 * fac;
    rng->has_spare_normal = true;
    return u1 * fac;
}

void md_rng_jump(md_rng *rng) {
    static const uint64_t JUMP[4] = {
        0x180ec6d33cfd0abaULL,
        0xd5a61266f0c9392cULL,
        0xa9582618e03fc9aaULL,
        0x39abdc4529b1661cULL
    };
    uint64_t s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    for (int i = 0; i < 4; ++i) {
        for (int b = 0; b < 64; ++b) {
            if (JUMP[i] & (1ULL << b)) {
                s0 ^= rng->s[0];
                s1 ^= rng->s[1];
                s2 ^= rng->s[2];
                s3 ^= rng->s[3];
            }
            (void)md_rng_u64(rng);
        }
    }
    rng->s[0] = s0;
    rng->s[1] = s1;
    rng->s[2] = s2;
    rng->s[3] = s3;
    rng->has_spare_normal = false;
}
