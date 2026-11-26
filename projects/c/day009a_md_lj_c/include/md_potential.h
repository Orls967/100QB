#ifndef MD_POTENTIAL_H
#define MD_POTENTIAL_H

#include "md_types.h"

/* Initialize standard parameter sets for each potential kind */
void md_potential_init(md_potential *pot, md_pot_kind kind, double rc);
void md_potential_free(md_potential *pot);

const char *md_pot_name(md_pot_kind kind);

/* Exact analytical pair evaluation given r2 = r^2 < rc^2.
 * Returns pair energy u(r) and sets *f_div_r = -(1/r) dU/dr so F_ij = (*f_div_r) * r_ij. */
double md_pot_eval_exact(const md_potential *pot, double r2, double *f_div_r);

/* Evaluation using either lookup table (if enabled) or analytical formula */
double md_pot_eval(const md_potential *pot, double r2, double *f_div_r);

/* Build lookup table on r^2 in [rmin^2, rc^2] */
void md_pot_build_table(md_potential *pot, int n_pts, double rmin, bool cubic);

/* Numerical central finite difference -(U(r+h) - U(r-h)) / (2h) for validation */
double md_pot_force_fd(const md_potential *pot, double r, double h);

/* Analytical LJ long-range tail corrections (3D) for energy per particle and pressure */
void md_pot_tail_corrections(const md_potential *pot, int n, double vol,
                             double *u_tail_total, double *p_tail);

/* Attractive-only part of LJ (r > 2^(1/6)) for two-timescale r-RESPA splitting */
double md_pot_lj_attr_eval(const md_potential *pot, double r2, double *f_div_r);

#endif /* MD_POTENTIAL_H */
