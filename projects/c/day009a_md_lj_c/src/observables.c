#include "md_observables.h"
#include "md_potential.h"
#include "md_vec.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double md_kinetic_energy(const md_system *sys) {
    double sum = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        sum += v3_norm2(sys->vel[i]);
    }
    return 0.5 * sum;
}

void md_kinetic_components(const md_system *sys, double *kx, double *ky, double *kz) {
    double sx = 0.0, sy = 0.0, sz = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        sx += sys->vel[i].x * sys->vel[i].x;
        sy += sys->vel[i].y * sys->vel[i].y;
        sz += sys->vel[i].z * sys->vel[i].z;
    }
    *kx = 0.5 * sx;
    *ky = 0.5 * sy;
    *kz = 0.5 * sz;
}

int md_dof(const md_system *sys) {
    int d = sys->box.dim;
    return d * ((sys->n > 1) ? (sys->n - 1) : 1);
}

double md_temperature(const md_system *sys) {
    return (2.0 * md_kinetic_energy(sys)) / (double)md_dof(sys);
}

double md_pressure(const md_system *sys) {
    double vol = md_box_volume(&sys->box);
    double ke = md_kinetic_energy(sys);
    double d = (double)sys->box.dim;
    /* P = (2*K + W) / (d * V) + P_tail (note 2K/d = (N-1)*T, exact virial theorem) */
    return (2.0 * ke + sys->virial) / (d * vol) + sys->p_tail;
}

md_stress_tensor md_compute_stress(const md_system *sys) {
    md_stress_tensor st;
    double kxx = 0.0, kyy = 0.0, kzz = 0.0;
    double kxy = 0.0, kxz = 0.0, kyz = 0.0;
    for (int i = 0; i < sys->n; ++i) {
        double vx = sys->vel[i].x, vy = sys->vel[i].y, vz = sys->vel[i].z;
        kxx += vx * vx;
        kyy += vy * vy;
        kzz += vz * vz;
        kxy += vx * vy;
        kxz += vx * vz;
        kyz += vy * vz;
    }
    double inv_v = 1.0 / md_box_volume(&sys->box);
    st.pxx = (kxx + sys->w_xx) * inv_v;
    st.pyy = (kyy + sys->w_yy) * inv_v;
    st.pzz = (kzz + sys->w_zz) * inv_v;
    st.pxy = (kxy + sys->w_xy) * inv_v;
    st.pxz = (kxz + sys->w_xz) * inv_v;
    st.pyz = (kyz + sys->w_yz) * inv_v;

    /* 5 orthogonal traceless shear stress components with identical Green-Kubo expectation */
    st.shear5[0] = st.pxy;
    st.shear5[1] = st.pxz;
    st.shear5[2] = st.pyz;
    st.shear5[3] = 0.5 * (st.pxx - st.pyy);
    st.shear5[4] = (st.pxx + st.pyy - 2.0 * st.pzz) / (2.0 * sqrt(3.0));
    return st;
}

void md_gr_init(md_gr_accum *gr, int n_bins, double r_max, int dim) {
    memset(gr, 0, sizeof(*gr));
    gr->n_bins = n_bins;
    gr->r_max = r_max;
    gr->dr = r_max / (double)n_bins;
    gr->dim = dim;
    gr->hist = (double *)calloc((size_t)n_bins, sizeof(double));
}

void md_gr_free(md_gr_accum *gr) {
    free(gr->hist);
    memset(gr, 0, sizeof(*gr));
}

void md_gr_sample(md_gr_accum *gr, const md_system *sys) {
    int n = sys->n;
    double r_max2 = gr->r_max * gr->r_max;
    double inv_dr = 1.0 / gr->dr;
    for (int i = 0; i < n; ++i) {
        md_vec3 pi = sys->pos[i];
        for (int j = i + 1; j < n; ++j) {
            md_vec3 d = md_box_mic(&sys->box, pi, sys->pos[j]);
            double r2 = v3_norm2(d);
            if (r2 < r_max2) {
                double r = sqrt(r2);
                int bin = (int)(r * inv_dr);
                if (bin >= 0 && bin < gr->n_bins) {
                    gr->hist[bin] += 2.0;
                }
            }
        }
    }
    gr->n_samples++;
    gr->n_particles = n;
    gr->rho_sum += (double)n / md_box_volume(&sys->box);
}

void md_gr_finalize(const md_gr_accum *gr, double *r_centers, double *g_vals, double *coord_num) {
    if (gr->n_samples == 0 || gr->n_particles == 0) return;
    double rho = gr->rho_sum / (double)gr->n_samples;
    double running = 0.0;
    for (int b = 0; b < gr->n_bins; ++b) {
        double r_lo = (double)b * gr->dr;
        double r_hi = r_lo + gr->dr;
        double r_mid = 0.5 * (r_lo + r_hi);
        double shell_vol;
        if (gr->dim == 3) {
            shell_vol = (4.0 / 3.0) * MD_PI * (r_hi * r_hi * r_hi - r_lo * r_lo * r_lo);
        } else {
            shell_vol = MD_PI * (r_hi * r_hi - r_lo * r_lo);
        }
        double norm = (double)gr->n_samples * (double)gr->n_particles * rho * shell_vol;
        double g = (norm > 0.0) ? (gr->hist[b] / norm) : 0.0;
        r_centers[b] = r_mid;
        g_vals[b] = g;
        if (coord_num) {
            running += g * rho * shell_vol;
            coord_num[b] = running;
        }
    }
}

double md_gr_energy_route(const md_potential *pot, double rho, int n_bins,
                          const double *r_centers, const double *g_vals, double dr) {
    /* U/N = 2 * pi * rho * integral_0^{rc} u(r) g(r) r^2 dr (+ U_tail/N if enabled) */
    double integral = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        double r = r_centers[b];
        if (r >= pot->rc || r <= 1e-9) continue;
        double fdr = 0.0;
        double u = md_pot_eval_exact(pot, r * r, &fdr);
        integral += u * g_vals[b] * r * r * dr;
    }
    double u_per_n = 2.0 * MD_PI * rho * integral;
    if (pot->use_tail) {
        double ut = 0.0, pt = 0.0;
        md_pot_tail_corrections(pot, 1, 1.0 / rho, &ut, &pt);
        u_per_n += ut;
    }
    return u_per_n;
}

double md_gr_pressure_route(const md_potential *pot, double rho, double temp, int n_bins,
                            const double *r_centers, const double *g_vals, double dr) {
    /* P = rho*T + (2*pi/3)*rho^2 * integral_0^{rc} (r^2 * f_div_r) * g(r) * r^2 dr (+ P_tail)
     * Note: -r * dU/dr = r^2 * f_div_r */
    double integral = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        double r = r_centers[b];
        if (r >= pot->rc || r <= 1e-9) continue;
        double fdr = 0.0;
        (void)md_pot_eval_exact(pot, r * r, &fdr);
        double w_r = r * r * fdr; /* -r dU/dr */
        integral += w_r * g_vals[b] * r * r * dr;
    }
    double p = rho * temp + (2.0 * MD_PI / 3.0) * rho * rho * integral;
    if (pot->use_tail) {
        double ut = 0.0, pt = 0.0;
        md_pot_tail_corrections(pot, 1, 1.0 / rho, &ut, &pt);
        p += pt;
    }
    return p;
}

void md_sk_init(md_sk_accum *sk, int n_bins, double k_max, int max_n) {
    memset(sk, 0, sizeof(*sk));
    sk->n_bins = n_bins;
    sk->k_max = k_max;
    sk->dk = k_max / (double)n_bins;
    sk->max_n = max_n;
    sk->sk_sum = (double *)calloc((size_t)n_bins, sizeof(double));
    sk->k_sum = (double *)calloc((size_t)n_bins, sizeof(double));
    sk->k_count = (int *)calloc((size_t)n_bins, sizeof(int));
}

void md_sk_free(md_sk_accum *sk) {
    free(sk->sk_sum);
    free(sk->k_sum);
    free(sk->k_count);
    memset(sk, 0, sizeof(*sk));
}

void md_sk_sample(md_sk_accum *sk, const md_system *sys) {
    /* S(k) = (1/N) * ( (sum cos(k.r))^2 + (sum sin(k.r))^2 ) for reciprocal k = 2*pi*(nx/Lx, ...) */
    int n = sys->n;
    int mn = sk->max_n;
    double dkx = MD_TWO_PI * sys->box.inv_lx;
    double dky = MD_TWO_PI * sys->box.inv_ly;
    double dkz = (sys->box.dim == 3) ? (MD_TWO_PI * sys->box.inv_lz) : 0.0;
    double k_max2 = sk->k_max * sk->k_max;
    double inv_dk = 1.0 / sk->dk;
    int nz_lo = (sys->box.dim == 3) ? -mn : 0;
    int nz_hi = (sys->box.dim == 3) ? mn : 0;

    for (int nx = 0; nx <= mn; ++nx) {
        double kx = (double)nx * dkx;
        int ny_lo = (nx == 0) ? 0 : -mn;
        for (int ny = ny_lo; ny <= mn; ++ny) {
            double ky = (double)ny * dky;
            int nz_start = (nx == 0 && ny == 0) ? 1 : nz_lo;
            for (int nz = nz_start; nz <= nz_hi; ++nz) {
                double kz = (double)nz * dkz;
                double k2 = kx * kx + ky * ky + kz * kz;
                if (k2 <= 1e-12 || k2 >= k_max2) continue;
                double kmag = sqrt(k2);
                int bin = (int)(kmag * inv_dk);
                if (bin < 0 || bin >= sk->n_bins) continue;

                double sc = 0.0, ss = 0.0;
                for (int i = 0; i < n; ++i) {
                    double kr = kx * sys->pos[i].x + ky * sys->pos[i].y + kz * sys->pos[i].z;
                    sc += cos(kr);
                    ss += sin(kr);
                }
                double val = (sc * sc + ss * ss) / (double)n;
                sk->sk_sum[bin] += val;
                sk->k_sum[bin] += kmag;
                sk->k_count[bin] += 1;
            }
        }
    }
    sk->n_samples++;
}

void md_sk_finalize(const md_sk_accum *sk, double *k_vals, double *sk_vals) {
    for (int b = 0; b < sk->n_bins; ++b) {
        if (sk->k_count[b] > 0) {
            k_vals[b] = sk->k_sum[b] / (double)sk->k_count[b];
            sk_vals[b] = sk->sk_sum[b] / (double)sk->k_count[b];
        } else {
            k_vals[b] = ((double)b + 0.5) * sk->dk;
            sk_vals[b] = -1.0; /* Empty reciprocal-lattice bin sentinel */
        }
    }
}

void md_dyn_init(md_dyn_corr *dc, int n_particles, int max_lag, int max_origins, double dt_sample) {
    memset(dc, 0, sizeof(*dc));
    dc->n_particles = n_particles;
    dc->max_lag = max_lag;
    dc->max_origins = max_origins;
    dc->dt_sample = dt_sample;
    dc->pos_hist = (md_vec3 *)calloc((size_t)max_origins * (size_t)n_particles, sizeof(md_vec3));
    dc->vel_hist = (md_vec3 *)calloc((size_t)max_origins * (size_t)n_particles, sizeof(md_vec3));
    dc->msd_sum = (double *)calloc((size_t)max_lag, sizeof(double));
    dc->vacf_sum = (double *)calloc((size_t)max_lag, sizeof(double));
    dc->counts = (int *)calloc((size_t)max_lag, sizeof(int));
}

void md_dyn_free(md_dyn_corr *dc) {
    free(dc->pos_hist);
    free(dc->vel_hist);
    free(dc->msd_sum);
    free(dc->vacf_sum);
    free(dc->counts);
    memset(dc, 0, sizeof(*dc));
}

void md_dyn_record(md_dyn_corr *dc, const md_system *sys) {
    int n = dc->n_particles;
    int slot = dc->n_stored % dc->max_origins;
    md_vec3 *p_dst = dc->pos_hist + (size_t)slot * (size_t)n;
    md_vec3 *v_dst = dc->vel_hist + (size_t)slot * (size_t)n;
    memcpy(p_dst, sys->unwrapped, (size_t)n * sizeof(md_vec3));
    memcpy(v_dst, sys->vel, (size_t)n * sizeof(md_vec3));

    int max_back = dc->n_stored;
    if (max_back >= dc->max_lag) max_back = dc->max_lag - 1;
    if (max_back >= dc->max_origins) max_back = dc->max_origins - 1;

    double dim_d = (double)sys->box.dim;
    for (int lag = 0; lag <= max_back; ++lag) {
        int prev_slot = (dc->n_stored - lag) % dc->max_origins;
        const md_vec3 *p0 = dc->pos_hist + (size_t)prev_slot * (size_t)n;
        const md_vec3 *v0 = dc->vel_hist + (size_t)prev_slot * (size_t)n;
        double msd = 0.0;
        double vacf = 0.0;
        for (int i = 0; i < n; ++i) {
            md_vec3 dr = v3_sub(p_dst[i], p0[i]);
            msd += v3_norm2(dr);
            vacf += v3_dot(v_dst[i], v0[i]);
        }
        dc->msd_sum[lag] += msd / (double)n;
        dc->vacf_sum[lag] += vacf / (dim_d * (double)n);
        dc->counts[lag] += 1;
    }
    dc->n_stored++;
}

void md_dyn_finalize(const md_dyn_corr *dc, double *t_vals, double *msd, double *vacf, double *gk_int) {
    double running = 0.0;
    for (int lag = 0; lag < dc->max_lag; ++lag) {
        t_vals[lag] = (double)lag * dc->dt_sample;
        if (dc->counts[lag] > 0) {
            msd[lag] = dc->msd_sum[lag] / (double)dc->counts[lag];
            vacf[lag] = dc->vacf_sum[lag] / (double)dc->counts[lag];
        } else {
            msd[lag] = 0.0;
            vacf[lag] = 0.0;
        }
        if (lag > 0) {
            running += 0.5 * (vacf[lag - 1] + vacf[lag]) * dc->dt_sample;
        }
        if (gk_int) gk_int[lag] = running;
    }
}

static void eval_y6_all(md_vec3 u, double re[7], double im[7]) {
    /* Exact spherical harmonics Y_{6,m}(theta, phi) for m = 0..6, scaled by sqrt(4*pi/13)
     * so Q6 = sqrt(|q0|^2 + 2 * sum_{m=1..6} |qm|^2) directly! */
    static const double norm_fac[7] = {
        1.0,
        0.15430334996209191,   /* 1/sqrt(42) */
        0.02439750182371333,   /* 1/sqrt(1680) */
        0.004066250303952222,  /* 1/sqrt(60480) */
        7.423923415794149e-4,  /* 1/sqrt(1814400) */
        1.5827756178436998e-4, /* 1/sqrt(39916800) */
        4.569079789715103e-5   /* 1/sqrt(479001600) */
    };
    double x = u.z;
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x2 * x2;
    double x5 = x4 * x;
    double x6 = x4 * x2;

    double rho = sqrt(u.x * u.x + u.y * u.y);
    double s = rho;
    double s2 = s * s;
    double s3 = s2 * s;
    double s4 = s2 * s2;
    double s5 = s4 * s;
    double s6 = s4 * s2;

    double p[7];
    p[0] = (231.0 * x6 - 315.0 * x4 + 105.0 * x2 - 5.0) / 16.0;
    p[1] = -(21.0 / 8.0) * s * (33.0 * x5 - 30.0 * x3 + 5.0 * x);
    p[2] = (105.0 / 8.0) * s2 * (33.0 * x4 - 18.0 * x2 + 1.0);
    p[3] = -(315.0 / 2.0) * s3 * (11.0 * x3 - 3.0 * x);
    p[4] = (945.0 / 2.0) * s4 * (11.0 * x2 - 1.0);
    p[5] = -10395.0 * s5 * x;
    p[6] = 10395.0 * s6;

    double cphi = (rho > 1e-14) ? (u.x / rho) : 1.0;
    double sphi = (rho > 1e-14) ? (u.y / rho) : 0.0;

    double cm = 1.0, sm = 0.0;
    for (int m = 0; m <= 6; ++m) {
        double amp = norm_fac[m] * p[m];
        re[m] = amp * cm;
        im[m] = amp * sm;
        double ncm = cm * cphi - sm * sphi;
        double nsm = sm * cphi + cm * sphi;
        cm = ncm;
        sm = nsm;
    }
}

double md_compute_q6(const md_system *sys, double r_bond, double *local_q6_out) {
    int n = sys->n;
    double rb2 = r_bond * r_bond;
    double g_re[7] = {0}, g_im[7] = {0};
    int total_bonds = 0;

    for (int i = 0; i < n; ++i) {
        double l_re[7] = {0}, l_im[7] = {0};
        int nb = 0;
        md_vec3 pi = sys->pos[i];
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            md_vec3 d = md_box_mic(&sys->box, sys->pos[j], pi);
            double r2 = v3_norm2(d);
            if (r2 < rb2 && r2 > 1e-16) {
                double inv_r = 1.0 / sqrt(r2);
                md_vec3 u = v3_scale(d, inv_r);
                double re[7], im[7];
                eval_y6_all(u, re, im);
                for (int m = 0; m <= 6; ++m) {
                    l_re[m] += re[m];
                    l_im[m] += im[m];
                    g_re[m] += re[m];
                    g_im[m] += im[m];
                }
                nb++;
                total_bonds++;
            }
        }
        if (local_q6_out) {
            if (nb > 0) {
                double inv_nb = 1.0 / (double)nb;
                double s2 = (l_re[0] * inv_nb) * (l_re[0] * inv_nb) +
                            (l_im[0] * inv_nb) * (l_im[0] * inv_nb);
                for (int m = 1; m <= 6; ++m) {
                    double ar = l_re[m] * inv_nb;
                    double ai = l_im[m] * inv_nb;
                    s2 += 2.0 * (ar * ar + ai * ai);
                }
                local_q6_out[i] = sqrt(s2);
            } else {
                local_q6_out[i] = 0.0;
            }
        }
    }
    if (total_bonds == 0) return 0.0;
    double inv_tb = 1.0 / (double)total_bonds;
    double sum2 = (g_re[0] * inv_tb) * (g_re[0] * inv_tb) +
                  (g_im[0] * inv_tb) * (g_im[0] * inv_tb);
    for (int m = 1; m <= 6; ++m) {
        double ar = g_re[m] * inv_tb;
        double ai = g_im[m] * inv_tb;
        sum2 += 2.0 * (ar * ar + ai * ai);
    }
    return sqrt(sum2);
}

double md_compute_lindemann(const md_system *sys, const md_vec3 *ref_unwrapped, double r_nn) {
    int n = sys->n;
    double cm_x = 0.0, cm_y = 0.0, cm_z = 0.0;
    for (int i = 0; i < n; ++i) {
        cm_x += sys->unwrapped[i].x - ref_unwrapped[i].x;
        cm_y += sys->unwrapped[i].y - ref_unwrapped[i].y;
        cm_z += sys->unwrapped[i].z - ref_unwrapped[i].z;
    }
    double inv_n = 1.0 / (double)n;
    cm_x *= inv_n;
    cm_y *= inv_n;
    cm_z *= inv_n;

    double msd = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = (sys->unwrapped[i].x - ref_unwrapped[i].x) - cm_x;
        double dy = (sys->unwrapped[i].y - ref_unwrapped[i].y) - cm_y;
        double dz = (sys->unwrapped[i].z - ref_unwrapped[i].z) - cm_z;
        msd += dx * dx + dy * dy + dz * dz;
    }
    return sqrt(msd * inv_n) / r_nn;
}

double md_compute_b2_mayer(const md_potential *pot, double temp, int n_quad) {
    /* B2(T) = -2*pi * integral_0^{rc} (exp(-u(r)/T) - 1) * r^2 dr
     * Split at r_core = 0.6 where exp(-u/T) == 0 to machine precision for steep cores:
     * integral_0^{r_core} (-1) r^2 dr = -r_core^3 / 3 -> contributes +2*pi*r_core^3/3. */
    if (n_quad % 2 != 0) n_quad++;
    double r0 = (pot->kind == POT_YUKAWA) ? 1e-4 : 0.55;
    double r1 = pot->rc;
    double core_int = -(r0 * r0 * r0) / 3.0; /* f(r) = -1 on [0, r0] */
    if (pot->kind == POT_YUKAWA) core_int = 0.0;

    double h = (r1 - r0) / (double)n_quad;
    double sum = 0.0;
    for (int i = 0; i <= n_quad; ++i) {
        double r = r0 + (double)i * h;
        double fdr = 0.0;
        double u = (i == n_quad) ? 0.0 : md_pot_eval_exact(pot, r * r, &fdr);
        double arg = -u / temp;
        double mayer = (arg < -700.0) ? -1.0 : (exp(arg) - 1.0);
        double w = (i == 0 || i == n_quad) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
        sum += w * mayer * r * r;
    }
    double tail_int = sum * h / 3.0;
    return -2.0 * MD_PI * (core_int + tail_int);
}
