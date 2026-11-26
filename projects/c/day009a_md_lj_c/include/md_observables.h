#ifndef MD_OBSERVABLES_H
#define MD_OBSERVABLES_H

#include "md_types.h"

/* Basic instantaneous thermodynamics */
double md_kinetic_energy(const md_system *sys);
void md_kinetic_components(const md_system *sys, double *kx, double *ky, double *kz);
int md_dof(const md_system *sys);
double md_temperature(const md_system *sys);
double md_pressure(const md_system *sys);

/* Full 3x3 symmetric pressure/stress tensor P_ab = (1/V)(sum m v_a v_b + W_ab)
 * and 5 orthogonal traceless shear stress components for Green-Kubo viscosity */
typedef struct {
    double pxx, pyy, pzz;
    double pxy, pxz, pyz;
    double shear5[5];
} md_stress_tensor;

md_stress_tensor md_compute_stress(const md_system *sys);

/* Radial distribution function g(r) accumulator */
typedef struct {
    int n_bins;
    double r_max;
    double dr;
    double *hist;
    int n_samples;
    double rho_sum;
    int dim;
    int n_particles;
} md_gr_accum;

void md_gr_init(md_gr_accum *gr, int n_bins, double r_max, int dim);
void md_gr_free(md_gr_accum *gr);
void md_gr_sample(md_gr_accum *gr, const md_system *sys);
void md_gr_finalize(const md_gr_accum *gr, double *r_centers, double *g_vals, double *coord_num);

/* Energy and pressure integrals from g(r) */
double md_gr_energy_route(const md_potential *pot, double rho, int n_bins,
                          const double *r_centers, const double *g_vals, double dr);
double md_gr_pressure_route(const md_potential *pot, double rho, double temp, int n_bins,
                            const double *r_centers, const double *g_vals, double dr);

/* Static structure factor S(k) on reciprocal lattice */
typedef struct {
    int n_bins;
    double k_max;
    double dk;
    int max_n; /* Max integer wavevector index per axis */
    double *sk_sum;
    double *k_sum;
    int *k_count;
    int n_samples;
} md_sk_accum;

void md_sk_init(md_sk_accum *sk, int n_bins, double k_max, int max_n);
void md_sk_free(md_sk_accum *sk);
void md_sk_sample(md_sk_accum *sk, const md_system *sys);
void md_sk_finalize(const md_sk_accum *sk, double *k_vals, double *sk_vals);

/* Multi-origin MSD and VACF correlator */
typedef struct {
    int n_particles;
    int max_lag;
    int max_origins;
    int n_stored;
    int sample_stride;
    double dt_sample;
    md_vec3 *pos_hist; /* [max_origins * n_particles] */
    md_vec3 *vel_hist; /* [max_origins * n_particles] */
    double *msd_sum;   /* [max_lag] */
    double *vacf_sum;  /* [max_lag] */
    int *counts;       /* [max_lag] */
} md_dyn_corr;

void md_dyn_init(md_dyn_corr *dc, int n_particles, int max_lag, int max_origins, double dt_sample);
void md_dyn_free(md_dyn_corr *dc);
void md_dyn_record(md_dyn_corr *dc, const md_system *sys);
void md_dyn_finalize(const md_dyn_corr *dc, double *t_vals, double *msd, double *vacf, double *gk_int);

/* Steinhardt global Q6 and per-particle local q6(i) bond-orientational order */
double md_compute_q6(const md_system *sys, double r_bond, double *local_q6_out);

/* Lindemann parameter relative to reference lattice positions (COM-drift corrected) */
double md_compute_lindemann(const md_system *sys, const md_vec3 *ref_unwrapped, double r_nn);

/* Numerical Mayer second virial coefficient B2(T) via composite Simpson quadrature */
double md_compute_b2_mayer(const md_potential *pot, double temp, int n_quad);

#endif /* MD_OBSERVABLES_H */
