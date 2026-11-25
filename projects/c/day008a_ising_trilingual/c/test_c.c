#include "ising.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

static void test_rng_and_hash(void) {
    Xoshiro256pp rng;
    xoshiro_seed(&rng, 20261126ULL);
    /* Reference values from Python SplitMix64 + xoshiro256++ */
    uint64_t x0 = xoshiro_next(&rng);
    uint64_t x1 = xoshiro_next(&rng);
    assert(x0 == 0xA9920E0E7DAB3EF0ULL);
    assert(x1 == 0x033B9CEB4B3BA7A0ULL);

    uint64_t s0 = compute_task_seed(8, 0, 0, 0);
    uint64_t s1 = compute_task_seed(8, 0, 0, 0);
    uint64_t s2 = compute_task_seed(8, 0, 1, 0);
    assert(s0 == s1);
    assert(s0 != s2);

    int8_t cfg[4] = {1, -1, 1, 1};
    uint64_t h1 = fnv1a_64(cfg, 4);
    uint64_t h2 = fnv1a_64(cfg, 4);
    assert(h1 == h2 && h1 != 0);
}

static void test_union_find(void) {
    int parent[6];
    for (int i = 0; i < 6; ++i) parent[i] = -1;
    uf_unite(parent, 0, 1);
    uf_unite(parent, 2, 3);
    uf_unite(parent, 1, 2);
    assert(uf_find(parent, 0) == uf_find(parent, 3));
    assert(uf_find(parent, 0) != uf_find(parent, 4));
    assert(-parent[uf_find(parent, 0)] == 4);
    uf_unite(parent, 4, 5);
    uf_unite(parent, 0, 5);
    assert(-parent[uf_find(parent, 3)] == 6);
}

static void test_enum_l3(void) {
    /* L=3 has 512 states. Hand formula for 3x3 periodic Ising:
       E in {-18, -6, -2, 2, 6, 18}? Let's verify Z at beta=0 (Z=512) and beta=0.4 */
    EnumResult r0 = ising_enumerate_exact(3, 0.0);
    assert(fabs(r0.Z - 512.0) < 1e-12);
    assert(fabs(r0.e_mean - 0.0) < 1e-12);
    assert(fabs(r0.m2_mean - 1.0 / 9.0) < 1e-12);

    /* At beta=0.4, check Onsager Tc energy */
    double tc = 2.0 / log(1.0 + sqrt(2.0));
    double e_tc = ising_onsager_energy(1.0 / tc);
    assert(fabs(e_tc + sqrt(2.0)) < 1e-12);
}

static void test_determinism_and_cold_energy(void) {
    IsingLattice lat1, lat2;
    ising_init(&lat1, 8);
    ising_init(&lat2, 8);

    int64_t E0, M0;
    ising_energy_mag(&lat1, &E0, &M0);
    assert(E0 == -2 * 64 && M0 == 64);

    double beta = 0.44;
    double p_table[5] = {1.0, 1.0, 1.0, exp(-4.0 * beta), exp(-8.0 * beta)};
    double p_add = 1.0 - exp(-2.0 * beta);

    Xoshiro256pp r1, r2;
    xoshiro_seed(&r1, compute_task_seed(8, 0, 0, 0));
    xoshiro_seed(&r2, compute_task_seed(8, 0, 0, 0));
    for (int s = 0; s < 50; ++s) {
        ising_metropolis_sweep_table(&lat1, &r1, p_table);
        ising_metropolis_sweep_computed(&lat2, &r2, p_table);
    }
    assert(fnv1a_64(lat1.spins, (size_t)lat1.N) == fnv1a_64(lat2.spins, (size_t)lat2.N));

    ising_set_cold(&lat1);
    ising_set_cold(&lat2);
    xoshiro_seed(&r1, compute_task_seed(8, 0, 1, 0));
    xoshiro_seed(&r2, compute_task_seed(8, 0, 1, 0));
    for (int s = 0; s < 50; ++s) {
        ising_wolff_step(&lat1, &r1, p_add);
        ising_wolff_step(&lat2, &r2, p_add);
    }
    assert(fnv1a_64(lat1.spins, (size_t)lat1.N) == fnv1a_64(lat2.spins, (size_t)lat2.N));

    ising_set_cold(&lat1);
    ising_set_cold(&lat2);
    xoshiro_seed(&r1, compute_task_seed(8, 0, 2, 0));
    xoshiro_seed(&r2, compute_task_seed(8, 0, 2, 0));
    for (int s = 0; s < 50; ++s) {
        ising_sw_step(&lat1, &r1, p_add);
        ising_sw_step(&lat2, &r2, p_add);
    }
    assert(fnv1a_64(lat1.spins, (size_t)lat1.N) == fnv1a_64(lat2.spins, (size_t)lat2.N));

    ising_free(&lat1);
    ising_free(&lat2);
}

int main(void) {
    test_rng_and_hash();
    test_union_find();
    test_enum_l3();
    test_determinism_and_cold_energy();
    printf("C unit tests PASSED\n");
    return 0;
}
