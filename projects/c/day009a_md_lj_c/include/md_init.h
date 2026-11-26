#ifndef MD_INIT_H
#define MD_INIT_H

#include "md_types.h"

/* Build 3D cubic FCC lattice with 4 * nc^3 particles at number density rho */
void md_init_fcc(md_system *sys, int nc, double rho);

/* Build 2D triangular lattice with nx * ny particles at areal density rho_2d */
void md_init_triangular_2d(md_system *sys, int nx, int ny, double rho_2d);

/* Build non-overlapping random configuration of N particles in 3D cubic box at density rho */
int md_init_random_nonoverlap(md_system *sys, int n, double rho, double min_dist, md_rng *rng);

/* Assign Maxwell-Boltzmann velocities at temperature T, zero total COM momentum,
 * and optionally exact-rescale kinetic energy to 3(N-1)T/2 (or 2(N-1)T/2 in 2D) */
void md_init_velocities(md_system *sys, double temp, bool exact_scale, md_rng *rng);

/* Remove net center-of-mass momentum exactly */
void md_zero_momentum(md_system *sys);

/* Compute total momentum vector */
md_vec3 md_total_momentum(const md_system *sys);

/* Small analytical fixtures (dimer, equilateral triangle, regular tetrahedron) */
void md_init_fixture_dimer(md_system *sys, double r_sep, double box_l);
void md_init_fixture_triangle(md_system *sys, double side, double box_l);
void md_init_fixture_tetrahedron(md_system *sys, double side, double box_l);

#endif /* MD_INIT_H */
