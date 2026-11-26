#include "md_init.h"
#include "md_observables.h"
#include "md_rng.h"
#include "md_vec.h"
#include <math.h>

void md_init_fcc(md_system *sys, int nc, double rho) {
    int n = 4 * nc * nc * nc;
    double l = cbrt((double)n / rho);
    md_system_alloc(sys, n, md_box_cubic(l));

    double a = l / (double)nc;
    static const double basis[4][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5},
        {0.5, 0.5, 0.0}
    };
    int idx = 0;
    for (int ix = 0; ix < nc; ++ix) {
        for (int iy = 0; iy < nc; ++iy) {
            for (int iz = 0; iz < nc; ++iz) {
                for (int b = 0; b < 4; ++b) {
                    md_vec3 p;
                    p.x = ((double)ix + basis[b][0] + 0.25) * a;
                    p.y = ((double)iy + basis[b][1] + 0.25) * a;
                    p.z = ((double)iz + basis[b][2] + 0.25) * a;
                    p = md_box_wrap(&sys->box, p);
                    sys->pos[idx] = p;
                    sys->unwrapped[idx] = p;
                    idx++;
                }
            }
        }
    }
}

void md_init_triangular_2d(md_system *sys, int nx, int ny, double rho_2d) {
    /* Rectangular periodic cell hosting an exact 2D triangular lattice:
     * nx columns spaced by a, ny rows (ny even) spaced by a * sqrt(3)/2.
     * Area = nx * a * ny * a * sqrt(3)/2 = N / rho_2d. */
    if (ny % 2 != 0) ny += 1;
    int n = nx * ny;
    double s3_2 = 0.5 * sqrt(3.0);
    double a = sqrt((double)n / (rho_2d * (double)nx * (double)ny * s3_2));
    double lx = (double)nx * a;
    double ly = (double)ny * a * s3_2;
    md_system_alloc(sys, n, md_box_ortho(lx, ly, 1.0, 2));

    int idx = 0;
    for (int iy = 0; iy < ny; ++iy) {
        double shift = (iy % 2 == 1) ? (0.5 * a) : 0.0;
        for (int ix = 0; ix < nx; ++ix) {
            md_vec3 p;
            p.x = ((double)ix + 0.25) * a + shift;
            p.y = ((double)iy + 0.25) * a * s3_2;
            p.z = 0.0;
            p = md_box_wrap(&sys->box, p);
            sys->pos[idx] = p;
            sys->unwrapped[idx] = p;
            idx++;
        }
    }
}

int md_init_random_nonoverlap(md_system *sys, int n, double rho, double min_dist, md_rng *rng) {
    double l = cbrt((double)n / rho);
    md_system_alloc(sys, n, md_box_cubic(l));
    double min_d2 = min_dist * min_dist;
    int placed = 0;
    int max_attempts = 20000;

    for (int i = 0; i < n; ++i) {
        bool ok = false;
        for (int att = 0; att < max_attempts; ++att) {
            md_vec3 cand = v3_make(md_rng_uniform(rng) * l,
                                   md_rng_uniform(rng) * l,
                                   md_rng_uniform(rng) * l);
            ok = true;
            for (int j = 0; j < placed; ++j) {
                md_vec3 d = md_box_mic(&sys->box, cand, sys->pos[j]);
                if (v3_norm2(d) < min_d2) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                sys->pos[placed] = cand;
                sys->unwrapped[placed] = cand;
                placed++;
                break;
            }
        }
        if (!ok) {
            /* Fallback: slightly relaxed placement if dense packing rejected */
            md_vec3 cand = v3_make(md_rng_uniform(rng) * l,
                                   md_rng_uniform(rng) * l,
                                   md_rng_uniform(rng) * l);
            sys->pos[placed] = cand;
            sys->unwrapped[placed] = cand;
            placed++;
        }
    }
    return placed;
}

void md_zero_momentum(md_system *sys) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        sx += sys->vel[i].x;
        sy += sys->vel[i].y;
        sz += sys->vel[i].z;
    }
    double inv_n = 1.0 / (double)sys->n;
    double vx_cm = sx * inv_n;
    double vy_cm = sy * inv_n;
    double vz_cm = (sys->box.dim == 3) ? (sz * inv_n) : 0.0;
    for (int i = 0; i < sys->n; ++i) {
        sys->vel[i].x -= vx_cm;
        sys->vel[i].y -= vy_cm;
        sys->vel[i].z = (sys->box.dim == 3) ? (sys->vel[i].z - vz_cm) : 0.0;
    }
}

md_vec3 md_total_momentum(const md_system *sys) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        sx += sys->vel[i].x;
        sy += sys->vel[i].y;
        sz += sys->vel[i].z;
    }
    return v3_make(sx, sy, sz);
}

void md_init_velocities(md_system *sys, double temp, bool exact_scale, md_rng *rng) {
    double sigma = sqrt(temp);
    for (int i = 0; i < sys->n; ++i) {
        sys->vel[i].x = sigma * md_rng_normal(rng);
        sys->vel[i].y = sigma * md_rng_normal(rng);
        sys->vel[i].z = (sys->box.dim == 3) ? (sigma * md_rng_normal(rng)) : 0.0;
    }
    md_zero_momentum(sys);
    if (exact_scale && sys->n > 1) {
        double ke = md_kinetic_energy(sys);
        int dof = md_dof(sys);
        double target_ke = 0.5 * (double)dof * temp;
        if (ke > 1e-30) {
            double sc = sqrt(target_ke / ke);
            for (int i = 0; i < sys->n; ++i) {
                sys->vel[i].x *= sc;
                sys->vel[i].y *= sc;
                sys->vel[i].z *= sc;
            }
        }
    }
}

void md_init_fixture_dimer(md_system *sys, double r_sep, double box_l) {
    md_system_alloc(sys, 2, md_box_cubic(box_l));
    double cx = 0.5 * box_l;
    sys->pos[0] = v3_make(cx - 0.5 * r_sep, cx, cx);
    sys->pos[1] = v3_make(cx + 0.5 * r_sep, cx, cx);
    sys->unwrapped[0] = sys->pos[0];
    sys->unwrapped[1] = sys->pos[1];
}

void md_init_fixture_triangle(md_system *sys, double side, double box_l) {
    md_system_alloc(sys, 3, md_box_cubic(box_l));
    double cx = 0.5 * box_l;
    double h = 0.5 * sqrt(3.0) * side;
    sys->pos[0] = v3_make(cx - 0.5 * side, cx - h / 3.0, cx);
    sys->pos[1] = v3_make(cx + 0.5 * side, cx - h / 3.0, cx);
    sys->pos[2] = v3_make(cx, cx + 2.0 * h / 3.0, cx);
    for (int i = 0; i < 3; ++i) {
        sys->unwrapped[i] = sys->pos[i];
    }
}

void md_init_fixture_tetrahedron(md_system *sys, double side, double box_l) {
    md_system_alloc(sys, 4, md_box_cubic(box_l));
    /* Regular tetrahedron vertices scaled to exact edge length = side */
    double s = side / (2.0 * sqrt(2.0));
    double cx = 0.5 * box_l;
    sys->pos[0] = v3_make(cx + s, cx + s, cx + s);
    sys->pos[1] = v3_make(cx + s, cx - s, cx - s);
    sys->pos[2] = v3_make(cx - s, cx + s, cx - s);
    sys->pos[3] = v3_make(cx - s, cx - s, cx + s);
    for (int i = 0; i < 4; ++i) {
        sys->unwrapped[i] = sys->pos[i];
    }
}
