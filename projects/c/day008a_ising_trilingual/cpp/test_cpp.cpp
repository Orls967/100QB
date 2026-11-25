#include "ising.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

static void test_rng_and_hash() {
    ising::Xoshiro256pp rng(ising::MASTER_SEED);
    uint64_t x0 = rng.next_u64();
    uint64_t x1 = rng.next_u64();
    assert(x0 == 0xA9920E0E7DAB3EF0ULL);
    assert(x1 == 0x033B9CEB4B3BA7A0ULL);

    uint64_t s0 = ising::compute_task_seed(8, 0, 0, 0);
    uint64_t s1 = ising::compute_task_seed(8, 0, 0, 0);
    uint64_t s2 = ising::compute_task_seed(8, 0, 1, 0);
    assert(s0 == s1);
    assert(s0 != s2);

    int8_t cfg[4] = {1, -1, 1, 1};
    assert(ising::fnv1a_64(cfg, 4) == ising::fnv1a_64(cfg, 4));
}

static void test_union_find() {
    std::vector<int> parent(6, -1);
    ising::uf_unite(parent, 0, 1);
    ising::uf_unite(parent, 2, 3);
    ising::uf_unite(parent, 1, 2);
    assert(ising::uf_find(parent, 0) == ising::uf_find(parent, 3));
    assert(ising::uf_find(parent, 0) != ising::uf_find(parent, 4));
    assert(-parent[static_cast<size_t>(ising::uf_find(parent, 0))] == 4);
}

static void test_enum_l3() {
    auto r0 = ising::enumerate_exact(3, 0.0);
    assert(std::fabs(r0.Z - 512.0) < 1e-12);
    assert(std::fabs(r0.e_mean - 0.0) < 1e-12);
    assert(std::fabs(r0.m2_mean - 1.0 / 9.0) < 1e-12);

    double tc = 2.0 / std::log(1.0 + std::sqrt(2.0));
    double e_tc = ising::onsager_energy(1.0 / tc);
    assert(std::fabs(e_tc + std::sqrt(2.0)) < 1e-12);
}

static void test_determinism_and_constexpr() {
    ising::Lattice lat1(8), lat2(8);
    int64_t E0 = 0, M0 = 0;
    lat1.energy_mag(E0, M0);
    assert(E0 == -128 && M0 == 64);

    double beta = 0.44;
    double p_table[5] = {1.0, 1.0, 1.0, std::exp(-4.0 * beta), std::exp(-8.0 * beta)};
    ising::Xoshiro256pp r1(ising::compute_task_seed(8, 0, 0, 0));
    ising::Xoshiro256pp r2(ising::compute_task_seed(8, 0, 0, 0));
    for (int s = 0; s < 50; ++s) {
        lat1.metropolis_sweep_table(r1, p_table);
        ising::metropolis_sweep_constexpr<8>(lat2.spins, r2, p_table);
    }
    assert(ising::fnv1a_64(lat1.spins.data(), static_cast<size_t>(lat1.N)) ==
           ising::fnv1a_64(lat2.spins.data(), static_cast<size_t>(lat2.N)));
}

int main() {
    test_rng_and_hash();
    test_union_find();
    test_enum_l3();
    test_determinism_and_constexpr();
    std::printf("C++ unit tests PASSED\n");
    return 0;
}
