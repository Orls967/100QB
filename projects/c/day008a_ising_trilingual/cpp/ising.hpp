#pragma once

#include "rng.hpp"
#include <cstdint>
#include <vector>

namespace ising {

struct EnumResult {
    double beta;
    double Z;
    double e_mean;
    double m2_mean;
};

struct Lattice {
    int L;
    int N;
    std::vector<int8_t> spins;
    std::vector<int> nbr;
    std::vector<int> uf_parent;
    std::vector<int> stack;
    std::vector<uint8_t> flip;

    explicit Lattice(int l);
    void set_cold();
    void set_hot(Xoshiro256pp &rng);
    void energy_mag(int64_t &E_out, int64_t &M_out) const;

    void metropolis_sweep_table(Xoshiro256pp &rng, const double p_table[5]);
    void metropolis_sweep_computed(Xoshiro256pp &rng, const double p_table[5]);
    void wolff_step(Xoshiro256pp &rng, double p_add);
    void sw_step(Xoshiro256pp &rng, double p_add);
};

template <int L_CONST>
void metropolis_sweep_constexpr(std::vector<int8_t> &spins, Xoshiro256pp &rng, const double p_table[5]) {
    constexpr int N = L_CONST * L_CONST;
    int8_t *s = spins.data();
    for (int y = 0; y < L_CONST; ++y) {
        const int yp = (y + 1 == L_CONST) ? 0 : y + 1;
        const int ym = (y == 0) ? L_CONST - 1 : y - 1;
        for (int x = 0; x < L_CONST; ++x) {
            const int xp = (x + 1 == L_CONST) ? 0 : x + 1;
            const int xm = (x == 0) ? L_CONST - 1 : x - 1;
            const int i = y * L_CONST + x;
            const int sum_nn = s[y * L_CONST + xp] + s[yp * L_CONST + x] + s[y * L_CONST + xm] + s[ym * L_CONST + x];
            const int dE = 2 * s[i] * sum_nn;
            const double u = rng.u01();
            const int accept = (u < p_table[(dE + 8) >> 2]);
            s[i] = accept ? static_cast<int8_t>(-s[i]) : s[i];
        }
    }
    (void)N;
}

int uf_find(std::vector<int> &parent, int x);
void uf_unite(std::vector<int> &parent, int a, int b);
EnumResult enumerate_exact(int L, double beta);
double onsager_energy(double beta);
double yang_mag(double beta);

}  // namespace ising
