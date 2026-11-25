#ifndef DAY008A_ISING_H
#define DAY008A_ISING_H

#include "rng.h"
#include <stdint.h>

typedef struct {
    int L;
    int N;
    int8_t *spins;
    int *nbr;       /* N * 4: [right, up, left, down] */
    int *uf_parent; /* N for Swendsen-Wang */
    int *stack;     /* N for Wolff LIFO */
    uint8_t *flip;  /* N for Swendsen-Wang */
} IsingLattice;

typedef struct {
    double beta;
    double Z;
    double e_mean;
    double m2_mean;
} EnumResult;

typedef struct {
    double e_mean;
    double e_err;
    double abs_m_mean;
    double abs_m_err;
    double m2_mean;
    double m2_err;
    double m4_mean;
    double m4_err;
    double binder;
    double binder_err;
    double chi;
    double cv;
    double tau_int_e;
    double tau_int_m;
    uint64_t final_hash;
} SimStats;

void ising_init(IsingLattice *lat, int L);
void ising_free(IsingLattice *lat);
void ising_set_cold(IsingLattice *lat);
void ising_set_hot(IsingLattice *lat, Xoshiro256pp *rng);

void ising_energy_mag(const IsingLattice *lat, int64_t *E_out, int64_t *M_out);

void ising_metropolis_sweep_table(IsingLattice *lat, Xoshiro256pp *rng, const double p_table[5]);
void ising_metropolis_sweep_computed(IsingLattice *lat, Xoshiro256pp *rng, const double p_table[5]);
void ising_wolff_step(IsingLattice *lat, Xoshiro256pp *rng, double p_add);
void ising_sw_step(IsingLattice *lat, Xoshiro256pp *rng, double p_add);

int uf_find(int *parent, int x);
void uf_unite(int *parent, int a, int b);

EnumResult ising_enumerate_exact(int L, double beta);
double ising_onsager_energy(double beta);
double ising_yang_mag(double beta);

void ising_multispin_metropolis(int L, double beta, int sweeps, uint64_t seed,
                                double *e_mean_out, double *m2_mean_out);

#endif
