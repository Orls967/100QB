#ifndef DAY008A_RNG_H
#define DAY008A_RNG_H

#include <math.h>
#include <stdint.h>
#include <stddef.h>

#define MASTER_SEED 20261126ULL

typedef struct {
    uint64_t s[4];
} Xoshiro256pp;

static inline uint64_t splitmix64_mix(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static inline uint64_t splitmix64_next(uint64_t *state) {
    *state += 0x9E3779B97F4A7C15ULL;
    return splitmix64_mix(*state);
}

static inline void xoshiro_seed(Xoshiro256pp *rng, uint64_t seed) {
    uint64_t st = seed;
    rng->s[0] = splitmix64_next(&st);
    rng->s[1] = splitmix64_next(&st);
    rng->s[2] = splitmix64_next(&st);
    rng->s[3] = splitmix64_next(&st);
}

static inline uint64_t rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint64_t xoshiro_next(Xoshiro256pp *rng) {
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

static inline double rng_u01(Xoshiro256pp *rng) {
    return (double)(xoshiro_next(rng) >> 11) * 0x1.0p-53;
}

static inline int rng_int(Xoshiro256pp *rng, int n) {
    double u = rng_u01(rng);
    int idx = (int)floor(u * (double)n);
    if (idx >= n) idx = n - 1;
    return idx;
}

static inline uint64_t compute_task_seed(int L, int t_index, int algo_id, int replica) {
    uint64_t h = 0x517CC1B727220A95ULL;
    h = splitmix64_mix(h ^ (uint64_t)L);
    h = splitmix64_mix(h ^ (uint64_t)t_index);
    h = splitmix64_mix(h ^ (uint64_t)algo_id);
    h = splitmix64_mix(h ^ (uint64_t)replica);
    return splitmix64_mix(MASTER_SEED ^ h);
}

static inline uint64_t fnv1a_64(const int8_t *data, size_t len) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (size_t i = 0; i < len; ++i) {
        h ^= (uint8_t)data[i];
        h *= 0x100000001B3ULL;
    }
    return h;
}

#endif
