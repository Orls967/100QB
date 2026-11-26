#include "md_forces.h"
#include "md_init.h"
#include "md_integrator.h"
#include "md_io.h"
#include "md_neighbor.h"
#include "md_observables.h"
#include "md_potential.h"
#include "md_rng.h"
#include "md_stages.h"
#include "md_stats.h"
#include "md_svg.h"
#include "md_vec.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Group F: Structure g(r), S(k), Coordination, MSD/VACF Diffusion, Viscosity & 2D Dynamics (S45..S54) */

static void stage_s45_s46_s47_structure(md_run_ctx *ctx) {
    double t45 = md_wall_time_sec();
    /* Three canonical phases: gas (rho=0.1, T=2.0), dense liquid (rho=0.8, T=1.0), solid (rho=1.0, T=0.5) */
    double ph_rho[3] = {0.1, 0.8, 1.0};
    double ph_t[3] = {2.0, 1.0, 0.5};
    const char *ph_name[3] = {"gas", "liquid", "solid"};
    int ph_nc[3] = {4, 5, 4}; /* N=500 for dense liquid so L/2 = 4.27 allows clean g(r->L/2) */
    int n_eq = ctx->small_scale ? 600 : 2200;
    int n_prod = ctx->small_scale ? 1500 : 6000;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    int n_gbins = 160;
    double r_max_ph[3] = {4.20, 4.20, 3.15}; /* Liquid L/2 = 4.275 so r->L/2 is genuinely near L/2 */
    double r_cent_ph[3][160], gr_ph[3][160], cn_ph[3][160];
    double k_cent[60], sk_ph[3][60];

    FILE *fd = fopen("out/s45_gr_diagnostics.csv", "w");
    fprintf(fd, "phase,rho,target_t,n_particles,box_l,gr_core_max_lt08,gr_first_peak,gr_peak_r,gr_tail_mean,tail_dev_from_1\n");

    for (int p = 0; p < 3; ++p) {
        md_system sys;
        md_init_fcc(&sys, ph_nc[p], ph_rho[p]);
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 45, (uint64_t)p);
        md_init_velocities(&sys, (p == 1) ? 2.0 : ph_t[p], true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, (p == 1) ? 2.0 : ph_t[p], sys.n, 3);
        if (p == 1) {
            /* Melt liquid briefly at T=2.0 before equilibrating at T=1.0 */
            for (int s = 0; s < n_eq / 2; ++s) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
            }
            th.target_t = ph_t[p];
        }
        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
        }

        md_gr_accum gr;
        md_gr_init(&gr, n_gbins, r_max_ph[p], 3);
        md_sk_accum sk;
        md_sk_init(&sk, 60, 18.0, 7);

        for (int s = 1; s <= n_prod; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
            if (s % 15 == 0) md_gr_sample(&gr, &sys);
            if (s % 60 == 0) md_sk_sample(&sk, &sys);
        }
        md_gr_finalize(&gr, r_cent_ph[p], gr_ph[p], cn_ph[p]);
        md_sk_finalize(&sk, k_cent, sk_ph[p]);

        /* Diagnostics: check tail in [r_max - 0.45, r_max] -> genuinely near L/2 */
        double core_max = 0.0, peak_v = 0.0, peak_r = 0.0;
        double tail_sum = 0.0;
        int tail_cnt = 0;
        double tail_lo = r_max_ph[p] - 0.45;
        for (int b = 0; b < n_gbins; ++b) {
            double rv = r_cent_ph[p][b];
            if (rv < 0.8 && gr_ph[p][b] > core_max) core_max = gr_ph[p][b];
            if (gr_ph[p][b] > peak_v) { peak_v = gr_ph[p][b]; peak_r = rv; }
            if (rv >= tail_lo) {
                tail_sum += gr_ph[p][b];
                tail_cnt++;
            }
        }
        double tail_m = (tail_cnt > 0) ? (tail_sum / (double)tail_cnt) : 1.0;
        double tail_dev = fabs(tail_m - 1.0);
        if (p == 1) {
            ctx->c12_core_max = core_max;
            ctx->c12_tail_dev = tail_dev;
        }
        fprintf(fd, "%s,%.2f,%.2f,%d,%.4f,%.6e,%.6f,%.4f,%.6f,%.6f\n",
                ph_name[p], ph_rho[p], ph_t[p], sys.n, sys.box.lx,
                core_max, peak_v, peak_r, tail_m, tail_dev);

        md_sk_free(&sk);
        md_gr_free(&gr);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fd);

    FILE *fg = fopen("out/s45_gr_phases.csv", "w");
    fprintf(fg, "r_liq,gr_gas,gr_liq,r_sol,gr_sol\n");
    for (int b = 0; b < n_gbins; ++b) {
        fprintf(fg, "%.4f,%.6f,%.6f,%.4f,%.6f\n",
                r_cent_ph[1][b], gr_ph[0][b], gr_ph[1][b], r_cent_ph[2][b], gr_ph[2][b]);
    }
    fclose(fg);

    md_register_csv(ctx, "S45", "out/s45_gr_phases.csv");
    md_register_csv(ctx, "S45", "out/s45_gr_diagnostics.csv");

    md_svg_series sg[3] = {
        {"Gas (rho=0.1, T=2.0)", "#16a34a", r_cent_ph[0], gr_ph[0], NULL, n_gbins, true, false, false},
        {"Liquid (rho=0.8, T=1.0)", "#2563eb", r_cent_ph[1], gr_ph[1], NULL, n_gbins, true, false, false},
        {"Solid (rho=1.0, T=0.5)", "#dc2626", r_cent_ph[2], gr_ph[2], NULL, n_gbins, true, false, false}
    };
    md_svg_opts og = {"S45: Radial Distribution Function g(r) Across Gas, Liquid & Solid",
                      "Excluded core g(r<0.8)=0; liquid relaxes to 1.0 as r->L/2; solid Bragg peaks",
                      "Distance r", "g(r)", 760, 460, false, false, true, 1.0, "Ideal gas g=1"};
    const char *svg45 = "../../../outputs/figures/day009a_36_gr_phases.svg";
    md_svg_write_xy(svg45, &og, sg, 3);
    md_register_svg(ctx, "S45", svg45);
    md_ctx_record_stage(ctx, "S45", "gr_phases", md_wall_time_sec() - t45, 3);

    /* S46: S(k) phases */
    double t46 = md_wall_time_sec();
    FILE *fsk = fopen("out/s46_sk_phases.csv", "w");
    fprintf(fsk, "k,sk_gas,sk_liq,sk_sol\n");
    double pk_k[3] = {0}, pk_s[3] = {0};
    double plot_k[60], plot_sg[60], plot_sl[60], plot_ss[60];
    int n_pk = 0;
    for (int b = 0; b < 60; ++b) {
        double g0 = (sk_ph[0][b] > 0.0) ? sk_ph[0][b] : 0.0;
        double g1 = (sk_ph[1][b] > 0.0) ? sk_ph[1][b] : 0.0;
        double g2 = (sk_ph[2][b] > 0.0) ? sk_ph[2][b] : 0.0;
        if (g0 > 0.0 || g1 > 0.0 || g2 > 0.0) {
            fprintf(fsk, "%.4f,%.6f,%.6f,%.6f\n", k_cent[b], g0, g1, g2);
            plot_k[n_pk] = k_cent[b];
            plot_sg[n_pk] = g0;
            plot_sl[n_pk] = g1;
            plot_ss[n_pk] = (g2 > 15.0) ? 15.0 : g2; /* Cap sharp Bragg spike for visual clarity */
            n_pk++;
            if (g0 > pk_s[0]) { pk_s[0] = g0; pk_k[0] = k_cent[b]; }
            if (g1 > pk_s[1]) { pk_s[1] = g1; pk_k[1] = k_cent[b]; }
            if (g2 > pk_s[2]) { pk_s[2] = g2; pk_k[2] = k_cent[b]; }
        }
    }
    fclose(fsk);

    FILE *fskp = fopen("out/s46_sk_peaks.csv", "w");
    fprintf(fskp, "phase,k_main_peak,sk_main_peak\n");
    for (int p = 0; p < 3; ++p) fprintf(fskp, "%s,%.4f,%.6f\n", ph_name[p], pk_k[p], pk_s[p]);
    fclose(fskp);

    md_register_csv(ctx, "S46", "out/s46_sk_phases.csv");
    md_register_csv(ctx, "S46", "out/s46_sk_peaks.csv");

    md_svg_series ssk[3] = {
        {"Gas", "#16a34a", plot_k, plot_sg, NULL, n_pk, true, false, false},
        {"Liquid (peak k~6.8)", "#2563eb", plot_k, plot_sl, NULL, n_pk, true, false, false},
        {"Solid (capped at 15)", "#dc2626", plot_k, plot_ss, NULL, n_pk, true, false, false}
    };
    md_svg_opts osk = {"S46: Static Structure Factor S(k) Across Phases",
                       "Reciprocal-lattice S(k) = (1/N) |sum exp(-i k.r)|^2",
                       "Wavevector magnitude k", "S(k)", 760, 460, false, false, true, 1.0, "S=1"};
    const char *svg46 = "../../../outputs/figures/day009a_37_sk_phases.svg";
    md_svg_write_xy(svg46, &osk, ssk, 3);
    md_register_svg(ctx, "S46", svg46);
    md_ctx_record_stage(ctx, "S46", "sk_phases", md_wall_time_sec() - t46, 3);

    /* S47: Running coordination number + 2D liquid g(r) */
    double t47 = md_wall_time_sec();
    FILE *fcn = fopen("out/s47_coord_number.csv", "w");
    fprintf(fcn, "r_liq,n_gas,n_liq,r_sol,n_sol\n");
    for (int b = 0; b < n_gbins; ++b) {
        fprintf(fcn, "%.4f,%.6f,%.6f,%.4f,%.6f\n",
                r_cent_ph[1][b], cn_ph[0][b], cn_ph[1][b], r_cent_ph[2][b], cn_ph[2][b]);
    }
    fclose(fcn);

    /* 2D liquid g(r) at N=1024, rho2d=0.6, T=1.5 */
    md_system s2d;
    md_init_triangular_2d(&s2d, 32, 32, 0.6);
    md_rng rng2;
    md_rng_seed_stage(&rng2, MD_MASTER_SEED, 47, 0);
    md_init_velocities(&s2d, 1.5, true, &rng2);
    md_verlet_list vl2;
    md_verlet_init(&vl2, s2d.n, pot.rc, 0.3, false);
    md_forces_verlet(&s2d, &pot, &vl2);
    md_thermostat th2;
    md_thermostat_init(&th2, THERM_LANGEVIN, 1.5, s2d.n, 2);
    int eq2 = ctx->small_scale ? 300 : 900;
    int pr2 = ctx->small_scale ? 400 : 1200;
    for (int s = 0; s < eq2; ++s) {
        md_step_vverlet(&s2d, &pot, NBR_VERLET, NULL, &vl2, NULL, 1, 0.0025, &th2, &rng2);
    }
    md_gr_accum gr2;
    md_gr_init(&gr2, 100, 4.0, 2);
    for (int s = 1; s <= pr2; ++s) {
        md_step_vverlet(&s2d, &pot, NBR_VERLET, NULL, &vl2, NULL, 1, 0.0025, &th2, &rng2);
        if (s % 25 == 0) md_gr_sample(&gr2, &s2d);
    }
    double r2[100], g2[100], cn2[100];
    md_gr_finalize(&gr2, r2, g2, cn2);
    FILE *fg2 = fopen("out/s47_gr_2d.csv", "w");
    fprintf(fg2, "r,gr_2d,coord_2d\n");
    for (int b = 0; b < 100; ++b) fprintf(fg2, "%.4f,%.6f,%.6f\n", r2[b], g2[b], cn2[b]);
    fclose(fg2);
    md_gr_free(&gr2);
    md_verlet_free(&vl2);
    md_system_free(&s2d);

    md_register_csv(ctx, "S47", "out/s47_coord_number.csv");
    md_register_csv(ctx, "S47", "out/s47_gr_2d.csv");

    md_svg_series scn[3] = {
        {"Gas", "#16a34a", r_cent_ph[0], cn_ph[0], NULL, 90, true, false, false},
        {"Liquid", "#2563eb", r_cent_ph[1], cn_ph[1], NULL, 90, true, false, false},
        {"Solid (1st shell ~12)", "#dc2626", r_cent_ph[2], cn_ph[2], NULL, 90, true, false, false}
    };
    md_svg_opts ocn = {"S47: Running Coordination Number n(r) = 4*pi*rho * int g(r) r^2 dr",
                       "First-shell coordination plateau near 12 neighbors in FCC solid & dense liquid",
                       "Distance r", "Cumulative neighbors n(r)", 760, 460, false, false, true, 12.0, "FCC 1st shell = 12"};
    const char *svg47 = "../../../outputs/figures/day009a_38_coord_number.svg";
    md_svg_write_xy(svg47, &ocn, scn, 3);
    md_register_svg(ctx, "S47", svg47);
    md_ctx_record_stage(ctx, "S47", "coord_number", md_wall_time_sec() - t47, 3);
}

typedef struct {
    double d_einstein, d_einstein_se;
    double d_gk, d_gk_se;
    double b_de[10], b_dgk[10];
} diff_result;

static diff_result measure_diffusion_nve(md_run_ctx *ctx, int nc, double rho, double temp,
                                         int n_eq, int steps_per_block, uint64_t seed_id,
                                         double *t_out, double *msd_out, double *vacf_out,
                                         double *gk_out, int max_lag) {
    diff_result res;
    memset(&res, 0, sizeof(res));
    double dt = 0.0025;
    int sample_stride = 4; /* dt_sample = 0.01 */
    double dt_s = dt * (double)sample_stride;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_system sys;
    md_init_fcc(&sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, seed_id, (uint64_t)llround(rho * 1000.0 + temp * 100.0));
    md_init_velocities(&sys, temp, true, &rng);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);
    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, temp, sys.n, 3);
    for (int s = 0; s < n_eq; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
    }

    int n_blocks = 10;
    double *msd_acc = (double *)calloc((size_t)max_lag, sizeof(double));
    double *vacf_acc = (double *)calloc((size_t)max_lag, sizeof(double));
    double *gk_acc = (double *)calloc((size_t)max_lag, sizeof(double));
    double *tv = (double *)malloc((size_t)max_lag * sizeof(double));
    double *mv = (double *)malloc((size_t)max_lag * sizeof(double));
    double *vv = (double *)malloc((size_t)max_lag * sizeof(double));
    double *gv = (double *)malloc((size_t)max_lag * sizeof(double));

    for (int b = 0; b < n_blocks; ++b) {
        /* Gentle Langevin re-thermalization between NVE blocks keeps exact target T */
        for (int s = 0; s < 150; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
        }
        md_dyn_corr dc;
        md_dyn_init(&dc, sys.n, max_lag, max_lag + 4, dt_s);
        for (int s = 0; s < steps_per_block; ++s) {
            if (s % sample_stride == 0) md_dyn_record(&dc, &sys);
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, NULL, NULL);
        }
        md_dyn_finalize(&dc, tv, mv, vv, gv);

        /* Einstein fit on linear diffusive tail [0.40 * max_lag, 0.95 * max_lag] */
        int i0 = (int)(0.40 * (double)max_lag);
        int i1 = (int)(0.95 * (double)max_lag);
        md_ols_result lr = md_ols_fit(tv + i0, mv + i0, i1 - i0);
        res.b_de[b] = lr.slope / 6.0;
        /* Green-Kubo integral evaluated at 0.70 * max_lag */
        int igk = (int)(0.70 * (double)max_lag);
        res.b_dgk[b] = gv[igk];

        for (int l = 0; l < max_lag; ++l) {
            msd_acc[l] += mv[l] / (double)n_blocks;
            vacf_acc[l] += vv[l] / (double)n_blocks;
            gk_acc[l] += gv[l] / (double)n_blocks;
        }
        md_dyn_free(&dc);
    }

    res.d_einstein = md_mean(res.b_de, n_blocks);
    res.d_einstein_se = sqrt(md_variance(res.b_de, n_blocks) / (double)n_blocks);
    res.d_gk = md_mean(res.b_dgk, n_blocks);
    res.d_gk_se = sqrt(md_variance(res.b_dgk, n_blocks) / (double)n_blocks);

    if (t_out) memcpy(t_out, tv, (size_t)max_lag * sizeof(double));
    if (msd_out) memcpy(msd_out, msd_acc, (size_t)max_lag * sizeof(double));
    if (vacf_out) memcpy(vacf_out, vacf_acc, (size_t)max_lag * sizeof(double));
    if (gk_out) memcpy(gk_out, gk_acc, (size_t)max_lag * sizeof(double));

    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;

    free(msd_acc); free(vacf_acc); free(gk_acc);
    free(tv); free(mv); free(vv); free(gv);
    md_verlet_free(&vl);
    md_system_free(&sys);
    return res;
}

static void stage_s48_s49_s50_diffusion(md_run_ctx *ctx) {
    double t48 = md_wall_time_sec();
    int max_lag = 300; /* t_max = 3.0 */
    int n_eq = ctx->small_scale ? 600 : 2000;
    int spb = ctx->small_scale ? 1400 : 2800; /* 28,000 total NVE steps across 10 blocks */

    double t_v[300], msd_v[300], vacf_v[300], gk_v[300];
    diff_result dr = measure_diffusion_nve(ctx, 4, 0.6, 1.5, n_eq, spb, 48,
                                           t_v, msd_v, vacf_v, gk_v, max_lag);

    /* S48: MSD & Einstein */
    FILE *fm = fopen("out/s48_msd_curve.csv", "w");
    fprintf(fm, "t,msd,einstein_line_6dt\n");
    double fit_line[300];
    for (int l = 0; l < max_lag; ++l) {
        fit_line[l] = 6.0 * dr.d_einstein * t_v[l];
        fprintf(fm, "%.4f,%.6f,%.6f\n", t_v[l], msd_v[l], fit_line[l]);
    }
    fclose(fm);

    FILE *fmb = fopen("out/s48_msd_blocks.csv", "w");
    fprintf(fmb, "block,d_einstein\n");
    for (int b = 0; b < 10; ++b) fprintf(fmb, "%d,%.6f\n", b, dr.b_de[b]);
    fclose(fmb);
    md_register_csv(ctx, "S48", "out/s48_msd_curve.csv");
    md_register_csv(ctx, "S48", "out/s48_msd_blocks.csv");

    char sub48[128];
    snprintf(sub48, sizeof(sub48), "rho=0.6, T=1.5, N=256: D_Einstein = %.4f +- %.4f",
             dr.d_einstein, dr.d_einstein_se);
    md_svg_series sm[2] = {
        {"MD MSD(t)", "#2563eb", t_v, msd_v, NULL, max_lag, true, false, false},
        {"Linear 6*D_E*t", "#dc2626", t_v, fit_line, NULL, max_lag, true, false, true}
    };
    md_svg_opts om = {"S48: Mean Squared Displacement MSD(t) & Einstein Diffusion", sub48,
                      "Time lag t", "MSD(t)", 760, 460, false, false, false, 0.0, NULL};
    const char *svg48 = "../../../outputs/figures/day009a_39_msd_einstein.svg";
    md_svg_write_xy(svg48, &om, sm, 2);
    md_register_svg(ctx, "S48", svg48);
    md_ctx_record_stage(ctx, "S48", "msd_einstein", md_wall_time_sec() - t48, 3);

    /* S49: VACF & Green-Kubo + dense liquid backscattering cage */
    double t49 = md_wall_time_sec();
    double t_d[300], msd_d[300], vacf_d[300], gk_d[300];
    diff_result dr_dense = measure_diffusion_nve(ctx, 4, 0.85, 0.85, n_eq, spb / 2, 49,
                                                 t_d, msd_d, vacf_d, gk_d, max_lag);

    FILE *fv = fopen("out/s49_vacf_curve.csv", "w");
    fprintf(fv, "t,vacf_rho06_t15,gk_d_rho06_t15,vacf_dense_cage\n");
    for (int l = 0; l < max_lag; ++l) {
        fprintf(fv, "%.4f,%.6f,%.6f,%.6f\n", t_v[l], vacf_v[l], gk_v[l], vacf_d[l]);
    }
    fclose(fv);
    FILE *fvb = fopen("out/s49_vacf_blocks.csv", "w");
    fprintf(fvb, "block,d_gk\n");
    for (int b = 0; b < 10; ++b) fprintf(fvb, "%d,%.6f\n", b, dr.b_dgk[b]);
    fclose(fvb);

    md_register_csv(ctx, "S49", "out/s49_vacf_curve.csv");
    md_register_csv(ctx, "S49", "out/s49_vacf_blocks.csv");

    md_svg_series sv[3] = {
        {"VACF (rho=0.6, T=1.5)", "#2563eb", t_v, vacf_v, NULL, 200, true, false, false},
        {"VACF dense cage (rho=0.85)", "#dc2626", t_d, vacf_d, NULL, 200, true, false, false},
        {"Running int D_GK(t)", "#16a34a", t_v, gk_v, NULL, 200, true, false, true}
    };
    md_svg_opts ov = {"S49: Velocity Autocorrelation VACF(t) & Green-Kubo Integral",
                      "Dense liquid shows negative cage-backscattering lobe; D_GK plateaus cleanly",
                      "Time lag t", "VACF(t) / D_GK(t)", 760, 460, false, false, true, 0.0, "Zero"};
    const char *svg49 = "../../../outputs/figures/day009a_40_vacf_gk.svg";
    md_svg_write_xy(svg49, &ov, sv, 3);
    md_register_svg(ctx, "S49", svg49);
    md_ctx_record_stage(ctx, "S49", "vacf_gk", md_wall_time_sec() - t49, 3);

    /* S50: Direct consistency check Einstein vs Green-Kubo (c11 part 1) */
    double t50 = md_wall_time_sec();
    double rel_primary = fabs(dr.d_einstein - dr.d_gk) / (fabs(dr.d_einstein) + 1e-12);
    double rel_dense = fabs(dr_dense.d_einstein - dr_dense.d_gk) / (fabs(dr_dense.d_einstein) + 1e-12);
    ctx->c11_diff_rel_err = rel_primary;

    FILE *fc = fopen("out/s50_diff_consistency.csv", "w");
    fprintf(fc, "state,rho,target_t,d_einstein,d_einstein_se,d_gk,d_gk_se,rel_diff\n");
    fprintf(fc, "primary_fluid,0.60,1.50,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            dr.d_einstein, dr.d_einstein_se, dr.d_gk, dr.d_gk_se, rel_primary);
    fprintf(fc, "dense_liquid,0.85,0.85,%.6f,%.6f,%.6f,%.6f,%.6f\n",
            dr_dense.d_einstein, dr_dense.d_einstein_se, dr_dense.d_gk, dr_dense.d_gk_se, rel_dense);
    fclose(fc);

    FILE *fcb = fopen("out/s50_diff_blocks.csv", "w");
    fprintf(fcb, "block,d_einstein,d_gk,diff_block\n");
    for (int b = 0; b < 10; ++b) {
        fprintf(fcb, "%d,%.6f,%.6f,%.6f\n", b, dr.b_de[b], dr.b_dgk[b], dr.b_de[b] - dr.b_dgk[b]);
    }
    fclose(fcb);

    md_register_csv(ctx, "S50", "out/s50_diff_consistency.csv");
    md_register_csv(ctx, "S50", "out/s50_diff_blocks.csv");
    md_ctx_record_stage(ctx, "S50", "diff_consistency", md_wall_time_sec() - t50, 2);
}

static void stage_s51_diff_scan(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rhos[4] = {0.2, 0.4, 0.6, 0.8};
    double ts[4] = {1.0, 1.5, 2.0, 2.5};
    double d_rho[4], d_rho_se[4], d_t[4], d_t_se[4];
    int n_eq = ctx->small_scale ? 500 : 1400;
    int spb = ctx->small_scale ? 1000 : 1800;

    FILE *fr = fopen("out/s51_diff_rho.csv", "w");
    fprintf(fr, "rho,target_t,d_einstein,d_einstein_se,d_gk,d_gk_se\n");
    for (int i = 0; i < 4; ++i) {
        diff_result r = measure_diffusion_nve(ctx, 4, rhos[i], 1.5, n_eq, spb, (uint64_t)(510 + i),
                                              NULL, NULL, NULL, NULL, 200);
        d_rho[i] = r.d_einstein;
        d_rho_se[i] = r.d_einstein_se;
        fprintf(fr, "%.2f,1.50,%.6f,%.6f,%.6f,%.6f\n",
                rhos[i], r.d_einstein, r.d_einstein_se, r.d_gk, r.d_gk_se);
    }
    fclose(fr);

    bool mono = true;
    double min_z = 1e9;
    for (int i = 0; i < 3; ++i) {
        double gap = d_rho[i] - d_rho[i + 1];
        double se_g = sqrt(d_rho_se[i] * d_rho_se[i] + d_rho_se[i + 1] * d_rho_se[i + 1]);
        double z = (se_g > 1e-12) ? (gap / se_g) : 0.0;
        if (z < min_z) min_z = z;
        if (z < 3.0) mono = false;
    }
    ctx->c11_mono_3se = mono;
    ctx->c11_min_step_z = min_z;

    FILE *ft = fopen("out/s51_diff_t.csv", "w");
    fprintf(ft, "rho,target_t,d_einstein,d_einstein_se,d_gk,d_gk_se\n");
    for (int i = 0; i < 4; ++i) {
        diff_result r = measure_diffusion_nve(ctx, 4, 0.6, ts[i], n_eq, spb, (uint64_t)(515 + i),
                                              NULL, NULL, NULL, NULL, 200);
        d_t[i] = r.d_einstein;
        d_t_se[i] = r.d_einstein_se;
        fprintf(ft, "0.60,%.2f,%.6f,%.6f,%.6f,%.6f\n",
                ts[i], r.d_einstein, r.d_einstein_se, r.d_gk, r.d_gk_se);
    }
    fclose(ft);

    md_register_csv(ctx, "S51", "out/s51_diff_rho.csv");
    md_register_csv(ctx, "S51", "out/s51_diff_t.csv");

    md_svg_series ser[1] = {
        {"D(rho) at T=1.5", "#2563eb", rhos, d_rho, d_rho_se, 4, true, true, false}
    };
    md_svg_opts o = {"S51: Self-Diffusion Coefficient D vs Density rho (T=1.5)",
                     "Monotonic decrease with density (>= 3 SE between consecutive points)",
                     "Density rho", "Diffusion coefficient D", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_41_diff_rho_t.svg";
    md_svg_write_xy(svgp, &o, ser, 1);
    md_register_svg(ctx, "S51", svgp);

    md_svg_series sert[1] = {
        {"D(T) at rho=0.6", "#dc2626", ts, d_t, d_t_se, 4, true, true, false}
    };
    md_svg_opts ot = {"S51: Self-Diffusion Coefficient D vs Temperature T (rho=0.6)",
                      "Thermal enhancement of fluid self-diffusion with 10-block SE bars",
                      "Temperature T", "Diffusion coefficient D", 760, 460, false, false, false, 0.0, NULL};
    const char *svgt = "../../../outputs/figures/day009a_41b_diff_temp.svg";
    md_svg_write_xy(svgt, &ot, sert, 1);
    md_register_svg(ctx, "S51", svgt);
    md_ctx_record_stage(ctx, "S51", "diff_rho_t_scan", md_wall_time_sec() - t0, 4);
}

static void stage_s52_viscosity_gk(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    /* Triple-point LJ liquid: rho = 0.8442, T = 0.722; literature memory ~ 3.2, a priori [2.0, 4.5] */
    double rho = 0.8442, T = 0.722;
    double dt = 0.0025;
    int sample_stride = 2; /* dt_s = 0.005 */
    double dt_s = dt * (double)sample_stride;
    int max_lag = 440;     /* t_max = 2.2 */
    int tc_lag = 400;      /* A priori cutoff t_c = 2.0 */
    int n_blocks = 10;
    int steps_per_block = ctx->small_scale ? 2400 : 8000; /* 80,000 NVE steps in full scale */

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_system sys;
    md_init_fcc(&sys, 4, rho); /* N = 256 */
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 52, 0);
    md_init_velocities(&sys, 2.5, true, &rng);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);

    /* Melt thoroughly at T = 2.5 then cool to triple point T = 0.722 so state is genuine liquid */
    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, 2.5, sys.n, 3);
    int n_melt = ctx->small_scale ? 600 : 1800;
    int n_cool = ctx->small_scale ? 800 : 2500;
    for (int s = 0; s < n_melt; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
    }
    th.target_t = T;
    for (int s = 0; s < n_cool; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
    }

    double vol = md_box_volume(&sys.box);
    double *acf_avg = (double *)calloc((size_t)max_lag, sizeof(double));
    double *eta_run = (double *)calloc((size_t)max_lag, sizeof(double));
    double t_lags[440];
    double block_eta[10], block_t[10], block_q6[10];

    int n_samp = steps_per_block / sample_stride;
    double *sh_buf = (double *)malloc((size_t)(5 * n_samp) * sizeof(double));

    for (int b = 0; b < n_blocks; ++b) {
        /* Brief NVT re-equilibration between blocks to anchor T = 0.722 */
        for (int s = 0; s < 300; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
        }
        double t_sum = 0.0;
        int si = 0;
        for (int s = 0; s < steps_per_block; ++s) {
            if (s % sample_stride == 0 && si < n_samp) {
                md_stress_tensor st = md_compute_stress(&sys);
                for (int c = 0; c < 5; ++c) sh_buf[c * n_samp + si] = st.shear5[c];
                t_sum += md_temperature(&sys);
                si++;
            }
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, NULL, NULL);
        }
        double tb = t_sum / (double)si;
        block_t[b] = tb;
        block_q6[b] = md_compute_q6(&sys, 1.38, NULL);

        double *bacf = (double *)calloc((size_t)max_lag, sizeof(double));
        for (int l = 0; l < max_lag; ++l) {
            int cnt = si - l;
            double sum = 0.0;
            for (int c = 0; c < 5; ++c) {
                const double *arr = sh_buf + c * n_samp;
                for (int i = 0; i < cnt; ++i) {
                    sum += arr[i] * arr[i + l];
                }
            }
            bacf[l] = sum / (5.0 * (double)cnt);
            acf_avg[l] += bacf[l] / (double)n_blocks;
        }
        double integ = 0.0;
        for (int l = 1; l <= tc_lag; ++l) {
            integ += 0.5 * (bacf[l - 1] + bacf[l]) * dt_s;
        }
        block_eta[b] = (vol / tb) * integ;
        free(bacf);
    }

    double eta_m = md_mean(block_eta, n_blocks);
    double eta_se = sqrt(md_variance(block_eta, n_blocks) / (double)n_blocks);
    double t_avg = md_mean(block_t, n_blocks);
    ctx->c13_eta = eta_m;
    ctx->c13_eta_se = eta_se;

    double run_i = 0.0;
    for (int l = 0; l < max_lag; ++l) {
        t_lags[l] = (double)l * dt_s;
        if (l > 0) run_i += 0.5 * (acf_avg[l - 1] + acf_avg[l]) * dt_s;
        eta_run[l] = (vol / t_avg) * run_i;
    }

    FILE *fa = fopen("out/s52_stress_acf.csv", "w");
    fprintf(fa, "t,stress_acf,eta_running\n");
    for (int l = 0; l < max_lag; l += 2) {
        fprintf(fa, "%.4f,%.8f,%.6f\n", t_lags[l], acf_avg[l], eta_run[l]);
    }
    fclose(fa);

    FILE *fb = fopen("out/s52_viscosity_blocks.csv", "w");
    fprintf(fb, "block,t_block,q6_liquid_check,eta_block\n");
    for (int b = 0; b < n_blocks; ++b) {
        fprintf(fb, "%d,%.6f,%.6f,%.6f\n", b, block_t[b], block_q6[b], block_eta[b]);
    }
    fclose(fb);

    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;

    md_register_csv(ctx, "S52", "out/s52_stress_acf.csv");
    md_register_csv(ctx, "S52", "out/s52_viscosity_blocks.csv");

    char sub[128];
    snprintf(sub, sizeof(sub), "rho=0.8442, T=0.722: eta(t_c=2.0) = %.3f +- %.3f (lit memory ~3.2)",
             eta_m, eta_se);
    md_svg_series ser[1] = {
        {"Running Green-Kubo eta(t)", "#2563eb", t_lags, eta_run, NULL, max_lag, true, false, false}
    };
    md_svg_opts o = {"S52: Green-Kubo Shear Viscosity Running Integral (Triple Point)", sub,
                     "Upper integration limit t", "Viscosity eta(t)", 760, 460, false, false, true, 3.2, "Lit memory ~3.2"};
    const char *svgp = "../../../outputs/figures/day009a_42_viscosity_gk.svg";
    md_svg_write_xy(svgp, &o, ser, 1);
    md_register_svg(ctx, "S52", svgp);

    free(sh_buf); free(acf_avg); free(eta_run);
    md_verlet_free(&vl);
    md_system_free(&sys);
    md_ctx_record_stage(ctx, "S52", "viscosity_gk", md_wall_time_sec() - t0, 3);
}

static void stage_s53_diff_finitesize(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int ncs[4] = {3, 4, 5, 6}; /* N = 108, 256, 500, 864 */
    double inv_l[4], d_vals[4], d_ses[4];
    int n_eq = ctx->small_scale ? 400 : 1200;
    int spb = ctx->small_scale ? 800 : 1600;

    FILE *fp = fopen("out/s53_diff_finitesize.csv", "w");
    fprintf(fp, "n_particles,box_l,inv_l,d_einstein,d_se,d_gk\n");
    FILE *fb = fopen("out/s53_diff_fs_blocks.csv", "w");
    fprintf(fb, "n_particles,block,d_einstein\n");

    for (int i = 0; i < 4; ++i) {
        int n = 4 * ncs[i] * ncs[i] * ncs[i];
        double l = cbrt((double)n / 0.6);
        inv_l[i] = 1.0 / l;
        diff_result r = measure_diffusion_nve(ctx, ncs[i], 0.6, 1.5, n_eq, spb, (uint64_t)(530 + i),
                                              NULL, NULL, NULL, NULL, 180);
        d_vals[i] = r.d_einstein;
        d_ses[i] = r.d_einstein_se;
        fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,%.6f\n", n, l, inv_l[i], r.d_einstein, r.d_einstein_se, r.d_gk);
        for (int b = 0; b < 10; ++b) fprintf(fb, "%d,%d,%.6f\n", n, b, r.b_de[b]);
    }
    fclose(fp);
    fclose(fb);

    md_ols_result lr = md_ols_fit(inv_l, d_vals, 4);
    double fit_y[4];
    for (int i = 0; i < 4; ++i) fit_y[i] = lr.intercept + lr.slope * inv_l[i];

    md_register_csv(ctx, "S53", "out/s53_diff_finitesize.csv");
    md_register_csv(ctx, "S53", "out/s53_diff_fs_blocks.csv");

    char sub[128];
    snprintf(sub, sizeof(sub), "OLS vs 1/L: D_inf=%.4f, slope=%.4f +- %.4f (R2=%.3f)",
             lr.intercept, lr.slope, lr.slope_se, lr.r2);
    md_svg_series ser[2] = {
        {"MD D(L) (N=108..864)", "#2563eb", inv_l, d_vals, d_ses, 4, false, true, false},
        {"Linear fit in 1/L", "#dc2626", inv_l, fit_y, NULL, 4, true, false, true}
    };
    md_svg_opts o = {"S53: Finite-Size Scaling of Diffusion D(L) vs 1/L (rho=0.6, T=1.5)", sub,
                     "Inverse box length 1/L", "Diffusion coefficient D(L)", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_43_diff_finitesize.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S53", svgp);
    md_ctx_record_stage(ctx, "S53", "diff_finitesize", md_wall_time_sec() - t0, 3);
}

static void stage_s54_dynamics_2d(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double dt = 0.0025;
    int max_lag = 150;
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_system sys;
    md_init_triangular_2d(&sys, 32, 32, 0.6); /* N = 1024 */
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 54, 0);
    md_init_velocities(&sys, 1.5, true, &rng);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);
    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 2);

    int n_eq = ctx->small_scale ? 300 : 1000;
    int n_prod = ctx->small_scale ? 600 : 1800;
    for (int s = 0; s < n_eq; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
    }
    md_dyn_corr dc;
    md_dyn_init(&dc, sys.n, max_lag, max_lag + 4, dt * 4.0);
    double u_sum = 0.0, p_sum = 0.0;
    int cnt = 0;
    for (int s = 0; s < n_prod; ++s) {
        if (s % 4 == 0) md_dyn_record(&dc, &sys);
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, NULL, NULL);
        if (s % 20 == 0) {
            u_sum += sys.pe / (double)sys.n;
            p_sum += md_pressure(&sys);
            cnt++;
        }
    }
    double tv[150], mv[150], vv[150], gv[150];
    md_dyn_finalize(&dc, tv, mv, vv, gv);
    md_ols_result lr = md_ols_fit(tv + 60, mv + 60, 80);

    FILE *fs = fopen("out/s54_dynamics_2d.csv", "w");
    fprintf(fs, "dim,n_particles,rho_2d,target_t,u_per_n,p_2d,d_einstein_2d,d_gk_2d\n");
    fprintf(fs, "2,%d,0.60,1.50,%.6f,%.6f,%.6f,%.6f\n",
            sys.n, u_sum / (double)cnt, p_sum / (double)cnt, lr.slope / 4.0, gv[120]);
    fclose(fs);

    FILE *fm = fopen("out/s54_msd_2d.csv", "w");
    fprintf(fm, "t,msd_2d,vacf_2d,gk_2d\n");
    for (int l = 0; l < max_lag; ++l) {
        fprintf(fm, "%.4f,%.6f,%.6f,%.6f\n", tv[l], mv[l], vv[l], gv[l]);
    }
    fclose(fm);

    md_dyn_free(&dc);
    md_verlet_free(&vl);
    md_system_free(&sys);
    md_register_csv(ctx, "S54", "out/s54_dynamics_2d.csv");
    md_register_csv(ctx, "S54", "out/s54_msd_2d.csv");
    md_ctx_record_stage(ctx, "S54", "dynamics_2d", md_wall_time_sec() - t0, 2);
}

void md_run_stages_f(md_run_ctx *ctx) {
    stage_s45_s46_s47_structure(ctx);
    stage_s48_s49_s50_diffusion(ctx);
    stage_s51_diff_scan(ctx);
    stage_s52_viscosity_gk(ctx);
    stage_s53_diff_finitesize(ctx);
    stage_s54_dynamics_2d(ctx);
}
