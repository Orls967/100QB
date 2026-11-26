#include "md_potential.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

const char *md_pot_name(md_pot_kind kind) {
    switch (kind) {
        case POT_LJ_SHIFT: return "lj_shift";
        case POT_LJ_TAIL:  return "lj_tail";
        case POT_WCA:      return "wca";
        case POT_MORSE:    return "morse";
        case POT_YUKAWA:   return "yukawa";
        default:           return "unknown";
    }
}

static double pot_raw_eval(const md_potential *pot, double r2, double *f_div_r) {
    switch (pot->kind) {
        case POT_LJ_SHIFT:
        case POT_LJ_TAIL:
        case POT_WCA: {
            double inv_r2 = (pot->sigma * pot->sigma) / r2;
            double inv_r6 = inv_r2 * inv_r2 * inv_r2;
            double inv_r12 = inv_r6 * inv_r6;
            double u = 4.0 * pot->eps * (inv_r12 - inv_r6);
            *f_div_r = 24.0 * pot->eps * (2.0 * inv_r12 - inv_r6) / r2;
            return u;
        }
        case POT_MORSE: {
            double r = sqrt(r2);
            double e = exp(-pot->morse_alpha * (r - pot->morse_r0));
            double u = pot->morse_de * (e * e - 2.0 * e);
            /* dU/dr = 2 * D_e * alpha * (-e^2 + e) => -(1/r) dU/dr = 2*D_e*alpha*(e^2 - e)/r */
            *f_div_r = 2.0 * pot->morse_de * pot->morse_alpha * (e * e - e) / r;
            return u;
        }
        case POT_YUKAWA: {
            double r = sqrt(r2);
            double e = exp(-pot->yukawa_kappa * r);
            double u = pot->yukawa_a * e / r;
            /* -dU/dr = A * e^(-kr) * (k*r + 1) / r^2 => -(1/r) dU/dr = A * e * (k*r + 1) / (r^3) */
            *f_div_r = pot->yukawa_a * e * (pot->yukawa_kappa * r + 1.0) / (r2 * r);
            return u;
        }
        default:
            *f_div_r = 0.0;
            return 0.0;
    }
}

void md_potential_init(md_potential *pot, md_pot_kind kind, double rc) {
    memset(pot, 0, sizeof(*pot));
    pot->kind = kind;
    pot->eps = 1.0;
    pot->sigma = 1.0;
    pot->morse_de = 1.0;
    pot->morse_alpha = 6.0;
    pot->morse_r0 = MD_WCA_RC;
    pot->yukawa_a = 1.0;
    pot->yukawa_kappa = 2.0;

    if (kind == POT_WCA) {
        pot->rc = MD_WCA_RC;
    } else {
        pot->rc = (rc > 0.0) ? rc : 2.5;
    }
    pot->rc2 = pot->rc * pot->rc;
    pot->use_tail = (kind == POT_LJ_TAIL);

    double dummy_f = 0.0;
    if (kind == POT_LJ_TAIL) {
        /* Truncated LJ without energy shift; tail corrections added analytically */
        pot->u_shift = 0.0;
    } else {
        pot->u_shift = pot_raw_eval(pot, pot->rc2, &dummy_f);
    }
}

void md_potential_free(md_potential *pot) {
    free(pot->tab_u);
    free(pot->tab_du);
    free(pot->tab_fdivr);
    free(pot->tab_dfdivr);
    pot->tab_u = NULL;
    pot->tab_du = NULL;
    pot->tab_fdivr = NULL;
    pot->tab_dfdivr = NULL;
    pot->use_table = false;
    pot->tab_n = 0;
}

double md_pot_eval_exact(const md_potential *pot, double r2, double *f_div_r) {
    if (r2 >= pot->rc2) {
        *f_div_r = 0.0;
        return 0.0;
    }
    double u = pot_raw_eval(pot, r2, f_div_r) - pot->u_shift;
    return u;
}

void md_pot_build_table(md_potential *pot, int n_pts, double rmin, bool cubic) {
    md_potential_free(pot);
    pot->tab_n = n_pts;
    pot->tab_rmin2 = rmin * rmin;
    pot->tab_rmax2 = pot->rc2;
    pot->tab_dr2 = (pot->tab_rmax2 - pot->tab_rmin2) / (double)(n_pts - 1);
    pot->tab_inv_dr2 = 1.0 / pot->tab_dr2;
    pot->use_table = true;
    pot->use_cubic_table = cubic;

    pot->tab_u = (double *)malloc((size_t)n_pts * sizeof(double));
    pot->tab_du = (double *)malloc((size_t)n_pts * sizeof(double));
    pot->tab_fdivr = (double *)malloc((size_t)n_pts * sizeof(double));
    pot->tab_dfdivr = (double *)malloc((size_t)n_pts * sizeof(double));

    for (int i = 0; i < n_pts; ++i) {
        double r2 = pot->tab_rmin2 + (double)i * pot->tab_dr2;
        if (i == n_pts - 1) r2 = pot->tab_rmax2;
        double fdr = 0.0;
        double u = pot_raw_eval(pot, r2, &fdr) - pot->u_shift;
        pot->tab_u[i] = u;
        /* dU/d(r^2) = (1/(2r)) dU/dr = -0.5 * f_div_r */
        pot->tab_du[i] = -0.5 * fdr;
        pot->tab_fdivr[i] = fdr;
    }
    /* Compute d(f_div_r)/d(r^2) via central differences on table nodes for cubic Hermite */
    for (int i = 0; i < n_pts; ++i) {
        double r2 = pot->tab_rmin2 + (double)i * pot->tab_dr2;
        double h = 1e-6 * r2;
        if (h < 1e-9) h = 1e-9;
        double fp = 0.0, fm = 0.0;
        (void)pot_raw_eval(pot, r2 + h, &fp);
        (void)pot_raw_eval(pot, r2 - h, &fm);
        pot->tab_dfdivr[i] = (fp - fm) / (2.0 * h);
    }
}

double md_pot_eval(const md_potential *pot, double r2, double *f_div_r) {
    if (r2 >= pot->rc2) {
        *f_div_r = 0.0;
        return 0.0;
    }
    if (!pot->use_table || r2 < pot->tab_rmin2) {
        return md_pot_eval_exact(pot, r2, f_div_r);
    }
    double pos = (r2 - pot->tab_rmin2) * pot->tab_inv_dr2;
    int idx = (int)pos;
    if (idx < 0) idx = 0;
    if (idx >= pot->tab_n - 1) idx = pot->tab_n - 2;
    double t = pos - (double)idx;

    if (!pot->use_cubic_table) {
        double u = pot->tab_u[idx] + t * (pot->tab_u[idx + 1] - pot->tab_u[idx]);
        *f_div_r = pot->tab_fdivr[idx] + t * (pot->tab_fdivr[idx + 1] - pot->tab_fdivr[idx]);
        return u;
    }

    /* Cubic Hermite interpolation in r^2 for both U and f_div_r */
    double h = pot->tab_dr2;
    double t2 = t * t;
    double t3 = t2 * t;
    double h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    double h10 = t3 - 2.0 * t2 + t;
    double h01 = -2.0 * t3 + 3.0 * t2;
    double h11 = t3 - t2;

    double u = h00 * pot->tab_u[idx] + h10 * h * pot->tab_du[idx] +
               h01 * pot->tab_u[idx + 1] + h11 * h * pot->tab_du[idx + 1];
    *f_div_r = h00 * pot->tab_fdivr[idx] + h10 * h * pot->tab_dfdivr[idx] +
               h01 * pot->tab_fdivr[idx + 1] + h11 * h * pot->tab_dfdivr[idx + 1];
    return u;
}

double md_pot_force_fd(const md_potential *pot, double r, double h) {
    double fp = 0.0, fm = 0.0;
    double rp = r + h;
    double rm = r - h;
    double up = md_pot_eval_exact(pot, rp * rp, &fp);
    double um = md_pot_eval_exact(pot, rm * rm, &fm);
    return -(up - um) / (2.0 * h);
}

void md_pot_tail_corrections(const md_potential *pot, int n, double vol,
                             double *u_tail_total, double *p_tail) {
    if (!pot->use_tail || pot->kind != POT_LJ_TAIL) {
        *u_tail_total = 0.0;
        *p_tail = 0.0;
        return;
    }
    double rho = (double)n / vol;
    double s_rc = pot->sigma / pot->rc;
    double sr3 = s_rc * s_rc * s_rc;
    double sr9 = sr3 * sr3 * sr3;
    double sig3 = pot->sigma * pot->sigma * pot->sigma;
    /* U_tail / N = (8/3) * pi * rho * eps * sigma^3 * ((1/3)*(sigma/rc)^9 - (sigma/rc)^3) */
    double u_per_n = (8.0 / 3.0) * MD_PI * rho * pot->eps * sig3 * ((1.0 / 3.0) * sr9 - sr3);
    /* P_tail = (16/3) * pi * rho^2 * eps * sigma^3 * ((2/3)*(sigma/rc)^9 - (sigma/rc)^3) */
    double pt = (16.0 / 3.0) * MD_PI * rho * rho * pot->eps * sig3 * ((2.0 / 3.0) * sr9 - sr3);
    *u_tail_total = (double)n * u_per_n;
    *p_tail = pt;
}

double md_pot_lj_attr_eval(const md_potential *pot, double r2, double *f_div_r) {
    /* Smooth WCA decomposition:
     * U_fast(r) = U_LJ(r) + eps  for r < 2^(1/6), 0 otherwise
     * U_slow(r) = -eps - u_shift for r < 2^(1/6) (zero force!),
     *           = U_LJ(r) - u_shift for 2^(1/6) <= r < rc */
    if (r2 >= pot->rc2) {
        *f_div_r = 0.0;
        return 0.0;
    }
    double rwca2 = MD_WCA_RC * MD_WCA_RC;
    if (r2 < rwca2) {
        *f_div_r = 0.0;
        return -pot->eps - pot->u_shift;
    }
    double inv_r2 = 1.0 / r2;
    double inv_r6 = inv_r2 * inv_r2 * inv_r2;
    double inv_r12 = inv_r6 * inv_r6;
    *f_div_r = 24.0 * pot->eps * (2.0 * inv_r12 - inv_r6) / r2;
    return 4.0 * pot->eps * (inv_r12 - inv_r6) - pot->u_shift;
}
