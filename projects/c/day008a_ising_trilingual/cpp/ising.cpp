#include "ising.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace ising {

Lattice::Lattice(int l)
    : L(l), N(l * l), spins(static_cast<size_t>(N), 1), nbr(static_cast<size_t>(N) * 4, 0),
      uf_parent(static_cast<size_t>(N), -1), stack(static_cast<size_t>(N), 0), flip(static_cast<size_t>(N), 0) {
    for (int y = 0; y < L; ++y) {
        int yp = (y + 1 == L) ? 0 : y + 1;
        int ym = (y == 0) ? L - 1 : y - 1;
        for (int x = 0; x < L; ++x) {
            int xp = (x + 1 == L) ? 0 : x + 1;
            int xm = (x == 0) ? L - 1 : x - 1;
            int i = y * L + x;
            nbr[static_cast<size_t>(i) * 4 + 0] = y * L + xp;
            nbr[static_cast<size_t>(i) * 4 + 1] = yp * L + x;
            nbr[static_cast<size_t>(i) * 4 + 2] = y * L + xm;
            nbr[static_cast<size_t>(i) * 4 + 3] = ym * L + x;
        }
    }
}

void Lattice::set_cold() {
    std::fill(spins.begin(), spins.end(), static_cast<int8_t>(1));
}

void Lattice::set_hot(Xoshiro256pp &rng) {
    for (int i = 0; i < N; ++i) {
        spins[static_cast<size_t>(i)] = (rng.u01() < 0.5) ? static_cast<int8_t>(1) : static_cast<int8_t>(-1);
    }
}

void Lattice::energy_mag(int64_t &E_out, int64_t &M_out) const {
    int64_t E = 0;
    int64_t M = 0;
    const int8_t *s = spins.data();
    const int *nb = nbr.data();
    for (int i = 0; i < N; ++i) {
        int8_t si = s[i];
        M += si;
        E -= static_cast<int64_t>(si) * static_cast<int64_t>(s[nb[i * 4 + 0]] + s[nb[i * 4 + 1]]);
    }
    E_out = E;
    M_out = M;
}

void Lattice::metropolis_sweep_table(Xoshiro256pp &rng, const double p_table[5]) {
    int8_t *s = spins.data();
    const int *nb = nbr.data();
    for (int i = 0; i < N; ++i) {
        const int base = i * 4;
        const int sum_nn = s[nb[base]] + s[nb[base + 1]] + s[nb[base + 2]] + s[nb[base + 3]];
        const int dE = 2 * s[i] * sum_nn;
        const double u = rng.u01();
        const int accept = (u < p_table[(dE + 8) >> 2]);
        s[i] = accept ? static_cast<int8_t>(-s[i]) : s[i];
    }
}

void Lattice::metropolis_sweep_computed(Xoshiro256pp &rng, const double p_table[5]) {
    int8_t *s = spins.data();
    for (int y = 0; y < L; ++y) {
        const int yp = (y + 1 == L) ? 0 : y + 1;
        const int ym = (y == 0) ? L - 1 : y - 1;
        for (int x = 0; x < L; ++x) {
            const int xp = (x + 1 == L) ? 0 : x + 1;
            const int xm = (x == 0) ? L - 1 : x - 1;
            const int i = y * L + x;
            const int sum_nn = s[y * L + xp] + s[yp * L + x] + s[y * L + xm] + s[ym * L + x];
            const int dE = 2 * s[i] * sum_nn;
            const double u = rng.u01();
            const int accept = (u < p_table[(dE + 8) >> 2]);
            s[i] = accept ? static_cast<int8_t>(-s[i]) : s[i];
        }
    }
}

void Lattice::wolff_step(Xoshiro256pp &rng, double p_add) {
    int8_t *s = spins.data();
    const int *nb = nbr.data();
    int *stk = stack.data();

    int i0 = rng.rand_int(N);
    const int8_t old_spin = s[i0];
    const int8_t new_spin = static_cast<int8_t>(-old_spin);
    s[i0] = new_spin;
    int sp = 0;
    stk[sp++] = i0;

    while (sp > 0) {
        int curr = stk[--sp];
        const int base = curr * 4;
        for (int d = 0; d < 4; ++d) {
            int n_idx = nb[base + d];
            if (s[n_idx] == old_spin) {
                if (rng.u01() < p_add) {
                    s[n_idx] = new_spin;
                    stk[sp++] = n_idx;
                }
            }
        }
    }
}

int uf_find(std::vector<int> &parent, int x) {
    int *p = parent.data();
    while (p[x] >= 0) {
        if (p[p[x]] >= 0) {
            p[x] = p[p[x]];
        }
        x = p[x];
    }
    return x;
}

void uf_unite(std::vector<int> &parent, int a, int b) {
    int ra = uf_find(parent, a);
    int rb = uf_find(parent, b);
    if (ra == rb) return;
    if (-parent[static_cast<size_t>(ra)] < -parent[static_cast<size_t>(rb)]) {
        std::swap(ra, rb);
    }
    parent[static_cast<size_t>(ra)] += parent[static_cast<size_t>(rb)];
    parent[static_cast<size_t>(rb)] = ra;
}

void Lattice::sw_step(Xoshiro256pp &rng, double p_add) {
    std::fill(uf_parent.begin(), uf_parent.end(), -1);
    const int8_t *s = spins.data();
    const int *nb = nbr.data();

    for (int i = 0; i < N; ++i) {
        const int8_t si = s[i];
        const int r_nb = nb[i * 4 + 0];
        if (si == s[r_nb] && rng.u01() < p_add) {
            uf_unite(uf_parent, i, r_nb);
        }
        const int u_nb = nb[i * 4 + 1];
        if (si == s[u_nb] && rng.u01() < p_add) {
            uf_unite(uf_parent, i, u_nb);
        }
    }

    for (int i = 0; i < N; ++i) {
        if (uf_parent[static_cast<size_t>(i)] < 0) {
            flip[static_cast<size_t>(i)] = (rng.u01() < 0.5) ? 1 : 0;
        }
    }

    int8_t *sm = spins.data();
    for (int i = 0; i < N; ++i) {
        int r = uf_find(uf_parent, i);
        if (flip[static_cast<size_t>(r)]) {
            sm[i] = static_cast<int8_t>(-sm[i]);
        }
    }
}

EnumResult enumerate_exact(int L, double beta) {
    const int N = L * L;
    const int n_rows = 1 << L;
    const int row_mask = n_rows - 1;
    int e_row[32]{};
    int m_row[32]{};
    int e_vert[32][32]{};

    for (int r = 0; r < n_rows; ++r) {
        int er = 0, mr = 0;
        for (int x = 0; x < L; ++x) {
            int s0 = ((r >> x) & 1) ? 1 : -1;
            int s1 = ((r >> ((x + 1) % L)) & 1) ? 1 : -1;
            mr += s0;
            er -= s0 * s1;
        }
        e_row[r] = er;
        m_row[r] = mr;
    }
    for (int r1 = 0; r1 < n_rows; ++r1) {
        for (int r2 = 0; r2 < n_rows; ++r2) {
            int ev = 0;
            for (int x = 0; x < L; ++x) {
                int s1 = ((r1 >> x) & 1) ? 1 : -1;
                int s2 = ((r2 >> x) & 1) ? 1 : -1;
                ev -= s1 * s2;
            }
            e_vert[r1][r2] = ev;
        }
    }

    double w_by_E[256]{};
    double m2_by_E[256]{};
    const uint64_t total = 1ULL << N;
    for (uint64_t mask = 0; mask < total; ++mask) {
        int r[5]{};
        int E = 0, M = 0;
        for (int y = 0; y < L; ++y) {
            r[y] = static_cast<int>((mask >> (y * L)) & static_cast<uint64_t>(row_mask));
            E += e_row[r[y]];
            M += m_row[r[y]];
        }
        for (int y = 0; y < L; ++y) {
            E += e_vert[r[y]][r[(y + 1) % L]];
        }
        int idx = E + 2 * N;
        w_by_E[idx] += 1.0;
        m2_by_E[idx] += static_cast<double>(M * M);
    }

    double Z = 0.0, sumE = 0.0, sumM2 = 0.0;
    for (int E = -2 * N; E <= 2 * N; E += 2) {
        int idx = E + 2 * N;
        if (w_by_E[idx] == 0.0) continue;
        double boltz = std::exp(-beta * static_cast<double>(E));
        Z += w_by_E[idx] * boltz;
        sumE += static_cast<double>(E) * w_by_E[idx] * boltz;
        sumM2 += m2_by_E[idx] * boltz;
    }

    return EnumResult{beta, Z, (sumE / Z) / static_cast<double>(N), (sumM2 / Z) / (static_cast<double>(N) * static_cast<double>(N))};
}

double onsager_energy(double beta) {
    const double s2 = std::sinh(2.0 * beta);
    const double c2 = std::cosh(2.0 * beta);
    const double t2 = std::tanh(2.0 * beta);
    const double k1 = (2.0 * s2) / (c2 * c2);
    const double kp = std::sqrt(std::max(0.0, 1.0 - k1 * k1));
    const double coth2b = c2 / s2;
    if (kp < 1e-14) {
        return -coth2b;
    }
    double a = 1.0, b = kp;
    for (int it = 0; it < 30; ++it) {
        double an = 0.5 * (a + b);
        double bn = std::sqrt(a * b);
        a = an;
        b = bn;
        if (std::fabs(a - b) < 1e-16) break;
    }
    constexpr double PI = 3.14159265358979323846;
    const double K_k1 = PI / (2.0 * a);
    return -coth2b * (1.0 + (2.0 / PI) * (2.0 * t2 * t2 - 1.0) * K_k1);
}

double yang_mag(double beta) {
    const double s2 = std::sinh(2.0 * beta);
    const double val = 1.0 - std::pow(s2, -4.0);
    return (val > 0.0) ? std::pow(val, 0.125) : 0.0;
}

}  // namespace ising
