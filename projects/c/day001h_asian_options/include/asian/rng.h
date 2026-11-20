#ifndef ASIAN_RNG_H
#define ASIAN_RNG_H

#include <math.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  uint64_t s[4];
  int has_spare;
  double spare;
} asian_rng_t;

static inline uint64_t asian_splitmix64(uint64_t *state) {
  *state += 0x9e3779b97f4a7c15ULL;
  uint64_t z = *state;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline uint64_t asian_rotl64(uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

static inline uint64_t asian_derive_chunk_seed(uint64_t seed, size_t chunk_idx) {
  return seed ^ (((uint64_t)chunk_idx + 1ULL) * 0x9e3779b97f4a7c15ULL);
}

static inline void asian_rng_seed(asian_rng_t *rng, uint64_t seed) {
  uint64_t st = seed;
  rng->s[0] = asian_splitmix64(&st);
  rng->s[1] = asian_splitmix64(&st);
  rng->s[2] = asian_splitmix64(&st);
  rng->s[3] = asian_splitmix64(&st);
  if ((rng->s[0] | rng->s[1] | rng->s[2] | rng->s[3]) == 0ULL) {
    rng->s[0] = 1ULL;
  }
  rng->has_spare = 0;
  rng->spare = 0.0;
}

static inline uint64_t asian_rng_u64(asian_rng_t *rng) {
  const uint64_t result = asian_rotl64(rng->s[1] * 5ULL, 7) * 9ULL;
  const uint64_t t = rng->s[1] << 17;

  rng->s[2] ^= rng->s[0];
  rng->s[3] ^= rng->s[1];
  rng->s[1] ^= rng->s[2];
  rng->s[0] ^= rng->s[3];
  rng->s[2] ^= t;
  rng->s[3] = asian_rotl64(rng->s[3], 45);
  return result;
}

static inline double asian_rng_uniform(asian_rng_t *rng) {
  const double inv_2_53 = 1.0 / 9007199254740992.0;
  return (double)(asian_rng_u64(rng) >> 11) * inv_2_53;
}

static inline double asian_rng_normal(asian_rng_t *rng) {
  if (rng->has_spare) {
    rng->has_spare = 0;
    return rng->spare;
  }
  double u = 0.0;
  double v = 0.0;
  double s = 0.0;
  do {
    u = 2.0 * asian_rng_uniform(rng) - 1.0;
    v = 2.0 * asian_rng_uniform(rng) - 1.0;
    s = u * u + v * v;
  } while (s >= 1.0 || s == 0.0);

  const double mul = sqrt(-2.0 * log(s) / s);
  rng->spare = v * mul;
  rng->has_spare = 1;
  return u * mul;
}

#endif /* ASIAN_RNG_H */
