#include "ising.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void ising_init(IsingLattice *lat, int L) {
    lat->L = L;
    lat->N = L * L;
    lat->spins = (int8_t *)malloc((size_t)lat->N * sizeof(int8_t));
    lat->nbr = (int *)malloc((size_t)lat->N * 4 * sizeof(int));
    lat->uf_parent = (int *)malloc((size_t)lat->N * sizeof(int));
    lat->stack = (int *)malloc((size_t)lat->N * sizeof(int));
    lat->flip = (uint8_t *)malloc((size_t)lat->N * sizeof(uint8_t));

    for (int y = 0; y < L; ++y) {
        int yp = (y + 1 == L) ? 0 : y + 1;
        int ym = (y == 0) ? L - 1 : y - 1;
        for (int x = 0; x < L; ++x) {
            int xp = (x + 1 == L) ? 0 : x + 1;
            int xm = (x == 0) ? L - 1 : x - 1;
            int i = y * L + x;
            lat->nbr[i * 4 + 0] = y * L + xp;   /* right */
            lat->nbr[i * 4 + 1] = yp * L + x;  /* up */
            lat->nbr[i * 4 + 2] = y * L + xm;  /* left */
            lat->nbr[i * 4 + 3] = ym * L + x;  /* down */
        }
    }
    ising_set_cold(lat);
}

void ising_free(IsingLattice *lat) {
    free(lat->spins);
    free(lat->nbr);
    free(lat->uf_parent);
    free(lat->stack);
    free(lat->flip);
    memset(lat, 0, sizeof(*lat));
}

void ising_set_cold(IsingLattice *lat) {
    memset(lat->spins, 1, (size_t)lat->N);
}

void ising_set_hot(IsingLattice *lat, Xoshiro256pp *rng) {
    for (int i = 0; i < lat->N; ++i) {
        lat->spins[i] = (rng_u01(rng) < 0.5) ? 1 : -1;
    }
}

void ising_energy_mag(const IsingLattice *lat, int64_t *E_out, int64_t *M_out) {
    int64_t E = 0;
    int64_t M = 0;
    const int N = lat->N;
    const int8_t *s = lat->spins;
    const int *nbr = lat->nbr;
    for (int i = 0; i < N; ++i) {
        int8_t si = s[i];
        M += si;
        E -= (int64_t)si * (int64_t)(s[nbr[i * 4 + 0]] + s[nbr[i * 4 + 1]]);
    }
    *E_out = E;
    *M_out = M;
}

void ising_metropolis_sweep_table(IsingLattice *lat, Xoshiro256pp *rng, const double p_table[5]) {
    const int N = lat->N;
    int8_t *s = lat->spins;
    const int *nbr = lat->nbr;
    for (int i = 0; i < N; ++i) {
        const int base = i * 4;
        const int sum_nn = s[nbr[base]] + s[nbr[base + 1]] + s[nbr[base + 2]] + s[nbr[base + 3]];
        const int dE = 2 * s[i] * sum_nn;
        const double u = rng_u01(rng);
        const int idx = (dE + 8) >> 2;
        const int accept = (u < p_table[idx]);
        s[i] = accept ? (int8_t)(-s[i]) : s[i];
    }
}

void ising_metropolis_sweep_computed(IsingLattice *lat, Xoshiro256pp *rng, const double p_table[5]) {
    const int L = lat->L;
    int8_t *s = lat->spins;
    for (int y = 0; y < L; ++y) {
        const int yp = (y + 1 == L) ? 0 : y + 1;
        const int ym = (y == 0) ? L - 1 : y - 1;
        for (int x = 0; x < L; ++x) {
            const int xp = (x + 1 == L) ? 0 : x + 1;
            const int xm = (x == 0) ? L - 1 : x - 1;
            const int i = y * L + x;
            const int sum_nn = s[y * L + xp] + s[yp * L + x] + s[y * L + xm] + s[ym * L + x];
            const int dE = 2 * s[i] * sum_nn;
            const double u = rng_u01(rng);
            if (dE <= 0 || u < p_table[(dE + 8) >> 2]) {
                s[i] = (int8_t)(-s[i]);
            }
        }
    }
}

void ising_wolff_step(IsingLattice *lat, Xoshiro256pp *rng, double p_add) {
    const int N = lat->N;
    int8_t *s = lat->spins;
    const int *nbr = lat->nbr;
    int *stack = lat->stack;

    int i0 = rng_int(rng, N);
    const int8_t old_spin = s[i0];
    const int8_t new_spin = (int8_t)(-old_spin);
    s[i0] = new_spin;
    int sp = 0;
    stack[sp++] = i0;

    while (sp > 0) {
        int curr = stack[--sp];
        const int base = curr * 4;
        for (int d = 0; d < 4; ++d) {
            int nb = nbr[base + d];
            if (s[nb] == old_spin) {
                if (rng_u01(rng) < p_add) {
                    s[nb] = new_spin;
                    stack[sp++] = nb;
                }
            }
        }
    }
}

int uf_find(int *parent, int x) {
    while (parent[x] >= 0) {
        if (parent[parent[x]] >= 0) {
            parent[x] = parent[parent[x]];
        }
        x = parent[x];
    }
    return x;
}

void uf_unite(int *parent, int a, int b) {
    int ra = uf_find(parent, a);
    int rb = uf_find(parent, b);
    if (ra == rb) return;
    if (-parent[ra] < -parent[rb]) {
        int tmp = ra;
        ra = rb;
        rb = tmp;
    }
    parent[ra] += parent[rb];
    parent[rb] = ra;
}

void ising_sw_step(IsingLattice *lat, Xoshiro256pp *rng, double p_add) {
    const int N = lat->N;
    int8_t *s = lat->spins;
    const int *nbr = lat->nbr;
    int *parent = lat->uf_parent;
    uint8_t *flip = lat->flip;

    for (int i = 0; i < N; ++i) {
        parent[i] = -1;
    }

    for (int i = 0; i < N; ++i) {
        const int8_t si = s[i];
        const int r_nb = nbr[i * 4 + 0];
        if (si == s[r_nb] && rng_u01(rng) < p_add) {
            uf_unite(parent, i, r_nb);
        }
        const int u_nb = nbr[i * 4 + 1];
        if (si == s[u_nb] && rng_u01(rng) < p_add) {
            uf_unite(parent, i, u_nb);
        }
    }

    for (int i = 0; i < N; ++i) {
        if (parent[i] < 0) {
            flip[i] = (rng_u01(rng) < 0.5) ? 1 : 0;
        }
    }

    for (int i = 0; i < N; ++i) {
        int r = uf_find(parent, i);
        if (flip[r]) {
            s[i] = (int8_t)(-s[i]);
        }
    }
}

EnumResult ising_enumerate_exact(int L, double beta) {
    const int N = L * L;
    const int n_rows = 1 << L;
    const int row_mask = n_rows - 1;
    int e_row[32];
    int m_row[32];
    int e_vert[32][32];

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

    /* Bucket by (E) since E is in [-2N, 2N] and even */
    double w_by_E[256];
    double m2_by_E[256];
    memset(w_by_E, 0, sizeof(w_by_E));
    memset(m2_by_E, 0, sizeof(m2_by_E));

    const uint64_t total = 1ULL << N;
    for (uint64_t mask = 0; mask < total; ++mask) {
        int r[5];
        int E = 0, M = 0;
        for (int y = 0; y < L; ++y) {
            r[y] = (int)((mask >> (y * L)) & (uint64_t)row_mask);
            E += e_row[r[y]];
            M += m_row[r[y]];
        }
        for (int y = 0; y < L; ++y) {
            E += e_vert[r[y]][r[(y + 1) % L]];
        }
        int idx = E + 2 * N;
        w_by_E[idx] += 1.0;
        m2_by_E[idx] += (double)(M * M);
    }

    double Z = 0.0, sumE = 0.0, sumM2 = 0.0;
    for (int E = -2 * N; E <= 2 * N; E += 2) {
        int idx = E + 2 * N;
        if (w_by_E[idx] == 0.0) continue;
        double boltz = exp(-beta * (double)E);
        Z += w_by_E[idx] * boltz;
        sumE += (double)E * w_by_E[idx] * boltz;
        sumM2 += m2_by_E[idx] * boltz;
    }

    EnumResult res;
    res.beta = beta;
    res.Z = Z;
    res.e_mean = (sumE / Z) / (double)N;
    res.m2_mean = (sumM2 / Z) / ((double)N * (double)N);
    return res;
}

double ising_onsager_energy(double beta) {
    const double s2 = sinh(2.0 * beta);
    const double c2 = cosh(2.0 * beta);
    const double t2 = tanh(2.0 * beta);
    const double k1 = (2.0 * s2) / (c2 * c2);
    const double kp = sqrt(fmax(0.0, 1.0 - k1 * k1));
    const double coth2b = c2 / s2;
    if (kp < 1e-14) {
        return -coth2b;
    }
    double a = 1.0, b = kp;
    for (int it = 0; it < 30; ++it) {
        double an = 0.5 * (a + b);
        double bn = sqrt(a * b);
        a = an;
        b = bn;
        if (fabs(a - b) < 1e-16) break;
    }
    const double K_k1 = (3.14159265358979323846) / (2.0 * a);
    return -coth2b * (1.0 + (2.0 / 3.14159265358979323846) * (2.0 * t2 * t2 - 1.0) * K_k1);
}

double ising_yang_mag(double beta) {
    const double s2 = sinh(2.0 * beta);
    const double val = 1.0 - pow(s2, -4.0);
    return (val > 0.0) ? pow(val, 0.125) : 0.0;
}

void ising_multispin_metropolis(int L, double beta, int sweeps, uint64_t seed,
                                double *e_mean_out, double *m2_mean_out) {
    /* 64 independent replicas packed into 64 bits per site */
    const int N = L * L;
    uint64_t *spins = (uint64_t *)malloc((size_t)N * sizeof(uint64_t));
    int *nbr = (int *)malloc((size_t)N * 4 * sizeof(int));
    for (int y = 0; y < L; ++y) {
        int yp = (y + 1 == L) ? 0 : y + 1;
        int ym = (y == 0) ? L - 1 : y - 1;
        for (int x = 0; x < L; ++x) {
            int xp = (x + 1 == L) ? 0 : x + 1;
            int xm = (x == 0) ? L - 1 : x - 1;
            int i = y * L + x;
            nbr[i * 4 + 0] = y * L + xp;
            nbr[i * 4 + 1] = yp * L + x;
            nbr[i * 4 + 2] = y * L + xm;
            nbr[i * 4 + 3] = ym * L + x;
        }
    }
    Xoshiro256pp rng;
    xoshiro_seed(&rng, seed);
    for (int i = 0; i < N; ++i) {
        spins[i] = xoshiro_next(&rng);
    }

    const double p4 = exp(-4.0 * beta);
    const double p8 = exp(-8.0 * beta);
    const int burn = sweeps / 4;
    double sum_e = 0.0, sum_m2 = 0.0;
    int meas = 0;

    for (int sw = 0; sw < sweeps; ++sw) {
        /* Checkerboard (2-sublattice) update for parallel 64-bit spin-words */
        for (int color = 0; color < 2; ++color) {
            for (int y = 0; y < L; ++y) {
                for (int x = 0; x < L; ++x) {
                    if (((x + y) & 1) != color) continue;
                    int i = y * L + x;
                    uint64_t si = spins[i];
                    uint64_t b0 = ~(si ^ spins[nbr[i * 4 + 0]]) & 1ULL;
                    (void)b0;
                    uint64_t d0 = si ^ spins[nbr[i * 4 + 0]];
                    uint64_t d1 = si ^ spins[nbr[i * 4 + 1]];
                    uint64_t d2 = si ^ spins[nbr[i * 4 + 2]];
                    uint64_t d3 = si ^ spins[nbr[i * 4 + 3]];
                    /* same-spin = 1, diff-spin = 0 */
                    uint64_t s0 = ~d0, s1 = ~d1, s2 = ~d2, s3 = ~d3;
                    /* Count same-spin neighbors (0..4): c4=4, c3=3, c2=2 */
                    uint64_t sum01_l = s0 ^ s1;
                    uint64_t sum01_h = s0 & s1;
                    uint64_t sum23_l = s2 ^ s3;
                    uint64_t sum23_h = s2 & s3;
                    uint64_t c0 = sum01_l ^ sum23_l;
                    uint64_t carry0 = sum01_l & sum23_l;
                    uint64_t c1 = sum01_h ^ sum23_h ^ carry0;
                    uint64_t c2 = (sum01_h & sum23_h) | ((sum01_h ^ sum23_h) & carry0);

                    /* dE <= 0 when count <= 2: (c2 == 0) & !(c1 & ~c0) -> actually count<=2 iff ~c2 & ~(c1 & c0) */
                    uint64_t le2 = (~c2) & (~(c1 & c0));
                    uint64_t eq3 = (~c2) & c1 & c0;
                    uint64_t eq4 = c2;

                    uint64_t r4 = 0, r8 = 0;
                    for (int b = 0; b < 64; ++b) {
                        if ((eq3 >> b) & 1ULL) {
                            if (rng_u01(&rng) < p4) r4 |= (1ULL << b);
                        }
                        if ((eq4 >> b) & 1ULL) {
                            if (rng_u01(&rng) < p8) r8 |= (1ULL << b);
                        }
                    }
                    uint64_t flip_mask = le2 | r4 | r8;
                    spins[i] ^= flip_mask;
                }
            }
        }
        if (sw >= burn) {
            int64_t E_tot = 0;
            int64_t M_acc[64];
            memset(M_acc, 0, sizeof(M_acc));
            for (int i = 0; i < N; ++i) {
                uint64_t si = spins[i];
                uint64_t sr = spins[nbr[i * 4 + 0]];
                uint64_t su = spins[nbr[i * 4 + 1]];
                int diff_r = __builtin_popcountll(si ^ sr);
                int diff_u = __builtin_popcountll(si ^ su);
                E_tot += (2 * (diff_r + diff_u) - 2 * 64);
                for (int b = 0; b < 64; ++b) {
                    M_acc[b] += ((si >> b) & 1ULL) ? 1 : -1;
                }
            }
            sum_e += (double)E_tot / (64.0 * (double)N);
            double m2s = 0.0;
            for (int b = 0; b < 64; ++b) {
                double m = (double)M_acc[b] / (double)N;
                m2s += m * m;
            }
            sum_m2 += m2s / 64.0;
            meas++;
        }
    }
    *e_mean_out = sum_e / (double)meas;
    *m2_mean_out = sum_m2 / (double)meas;
    free(spins);
    free(nbr);
}
