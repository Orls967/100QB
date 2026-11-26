#ifndef MD_RNG_H
#define MD_RNG_H

#include "md_types.h"

/* SplitMix64 mixing function used to expand seeds cleanly */
uint64_t md_splitmix64(uint64_t *state);

/* Seed xoshiro256++ from a single 64-bit seed via SplitMix64 */
void md_rng_seed(md_rng *rng, uint64_t seed);

/* Deterministic independent stream per (master_seed, stage_id, sub_id) */
void md_rng_seed_stage(md_rng *rng, uint64_t master_seed, uint64_t stage_id, uint64_t sub_id);

/* Next 64-bit unsigned integer (xoshiro256++) */
uint64_t md_rng_u64(md_rng *rng);

/* Uniform double in [0, 1) with full 53-bit precision */
double md_rng_uniform(md_rng *rng);

/* Uniform double in (a, b) */
double md_rng_range(md_rng *rng, double a, double b);

/* Standard normal N(0, 1) via polar Box-Muller with cached spare */
double md_rng_normal(md_rng *rng);

/* Jump function equivalent to 2^128 calls to md_rng_u64 */
void md_rng_jump(md_rng *rng);

#endif /* MD_RNG_H */
