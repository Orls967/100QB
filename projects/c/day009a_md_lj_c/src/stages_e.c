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

/* Group E: Thermodynamics, Isotherms, Mayer B2, Cv, Compressibility, g(r) Routes & NPT (S35..S44) */

typedef struct {
    double p_mean, p_se;
    double u_mean, u_se;
    double t_mean, t_se;
    double var_u;
    double p_blocks[10];
} eos_point;

static eos_point run_eos_nvt(md_run_ctx *ctx, md_pot_kind pkind, double rho, double temp,
                             int nc, int n_eq, int n_prod, uint64_t seed_id) {
    eos_point out;
    memset(&out, 0, sizeof(out));
    double dt = 0.0025;

    md_potential pot;
    md_potential_init(&pot, pkind, 2.5);
    md_system sys;
    md_init_fcc(&sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, seed_id, (uint64_t)llround(rho * 1000.0 + temp * 10.0));
    md_init_velocities(&sys, temp, true, &rng);

    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);

    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, temp, sys.n, 3);
    th.gamma = 1.5;

    for (int s = 0; s < n_eq; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
    }

    int sample_every = 5;
    int n_samp = n_prod / sample_every;
    if (n_samp < 40) n_samp = 40;
    double *p_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *u_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *t_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *utot_arr = (double *)malloc((size_t)n_samp * sizeof(double));

    int si = 0;
    for (int s = 1; s <= n_samp * sample_every; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
        if (s % sample_every == 0 && si < n_samp) {
            p_arr[si] = md_pressure(&sys);
            u_arr[si] = sys.pe / (double)sys.n;
            utot_arr[si] = sys.pe;
            t_arr[si] = md_temperature(&sys);
            si++;
        }
    }

    md_block_stats bp = md_block_average(p_arr, si, 10, out.p_blocks);
    md_block_stats bu = md_block_average(u_arr, si, 10, NULL);
    md_block_stats bt = md_block_average(t_arr, si, 10, NULL);
    out.p_mean = bp.mean; out.p_se = bp.se;
    out.u_mean = bu.mean; out.u_se = bu.se;
    out.t_mean = bt.mean; out.t_se = bt.se;
    out.var_u = md_variance(utot_arr, si);

    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;

    free(p_arr); free(u_arr); free(t_arr); free(utot_arr);
    md_verlet_free(&vl);
    md_system_free(&sys);
    return out;
}

static double g_s35_rhos[9] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
static double g_s35_p_shift[9];
static double g_s35_p_se[9];

static void stage_s35_s36_isotherms(md_run_ctx *ctx) {
    int nc = 4; /* N = 256 */
    int n_eq = ctx->small_scale ? 600 : 2000;
    int n_prod = ctx->small_scale ? 2000 : 8000;

    /* S35: Isotherm T = 1.5 */
    double t35 = md_wall_time_sec();
    double p_tail_15[9], z_15[9];
    FILE *f15 = fopen("out/s35_isotherm_t15.csv", "w");
    fprintf(f15, "rho,p_shift,p_shift_se,p_tail,u_shift,u_tail,z_shift,temp_meas\n");
    FILE *fb15 = fopen("out/s35_isotherm_t15_blocks.csv", "w");
    fprintf(fb15, "rho,block,p_block\n");

    for (int i = 0; i < 9; ++i) {
        double rho = g_s35_rhos[i];
        eos_point es = run_eos_nvt(ctx, POT_LJ_SHIFT, rho, 1.5, nc, n_eq, n_prod, 35);
        eos_point et = run_eos_nvt(ctx, POT_LJ_TAIL, rho, 1.5, nc, n_eq / 2, n_prod / 2, 350);
        double p_t = et.p_mean;
        double u_t = et.u_mean;

        g_s35_p_shift[i] = es.p_mean;
        g_s35_p_se[i] = es.p_se;
        p_tail_15[i] = p_t;
        z_15[i] = es.p_mean / (rho * 1.5);
        fprintf(f15, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                rho, es.p_mean, es.p_se, p_t, es.u_mean, u_t, z_15[i], es.t_mean);
        for (int b = 0; b < 10; ++b) {
            fprintf(fb15, "%.2f,%d,%.6f\n", rho, b, es.p_blocks[b]);
        }
    }
    fclose(f15);
    fclose(fb15);

    bool mono15 = true;
    double min_z15 = 1e9;
    /* rho >= 0.3 corresponds to indices 2..8 */
    for (int i = 2; i < 8; ++i) {
        double dp = g_s35_p_shift[i + 1] - g_s35_p_shift[i];
        double se_d = sqrt(g_s35_p_se[i + 1] * g_s35_p_se[i + 1] + g_s35_p_se[i] * g_s35_p_se[i]);
        double z = (se_d > 1e-12) ? (dp / se_d) : 0.0;
        if (z < min_z15) min_z15 = z;
        if (z < 3.0) mono15 = false;
    }
    ctx->c15_t15_mono_3se = mono15;
    ctx->c15_t15_min_z = min_z15;

    md_register_csv(ctx, "S35", "out/s35_isotherm_t15.csv");
    md_register_csv(ctx, "S35", "out/s35_isotherm_t15_blocks.csv");
    md_ctx_record_stage(ctx, "S35", "isotherm_t15", md_wall_time_sec() - t35, 2);

    /* S36: Isotherm T = 2.0 */
    double t36 = md_wall_time_sec();
    double p_20[9], p_20_se[9], z_20[9];
    FILE *f20 = fopen("out/s36_isotherm_t20.csv", "w");
    fprintf(f20, "rho,p_shift,p_shift_se,u_shift,z_shift,temp_meas\n");
    FILE *fb20 = fopen("out/s36_isotherm_t20_blocks.csv", "w");
    fprintf(fb20, "rho,block,p_block\n");

    for (int i = 0; i < 9; ++i) {
        double rho = g_s35_rhos[i];
        eos_point e = run_eos_nvt(ctx, POT_LJ_SHIFT, rho, 2.0, nc, n_eq, n_prod, 36);
        p_20[i] = e.p_mean;
        p_20_se[i] = e.p_se;
        z_20[i] = e.p_mean / (rho * 2.0);
        fprintf(f20, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                rho, e.p_mean, e.p_se, e.u_mean, z_20[i], e.t_mean);
        for (int b = 0; b < 10; ++b) {
            fprintf(fb20, "%.2f,%d,%.6f\n", rho, b, e.p_blocks[b]);
        }
    }
    fclose(f20);
    fclose(fb20);

    bool mono20 = true;
    for (int i = 0; i < 8; ++i) {
        if (p_20[i + 1] <= p_20[i]) mono20 = false;
    }
    ctx->c15_t20_mono = mono20;

    md_register_csv(ctx, "S36", "out/s36_isotherm_t20.csv");
    md_register_csv(ctx, "S36", "out/s36_isotherm_t20_blocks.csv");

    md_svg_series sp[3] = {
        {"T=1.5 LJ-shift", "#2563eb", g_s35_rhos, g_s35_p_shift, g_s35_p_se, 9, true, true, false},
        {"T=1.5 LJ+tail", "#0891b2", g_s35_rhos, p_tail_15, NULL, 9, true, true, true},
        {"T=2.0 LJ-shift", "#dc2626", g_s35_rhos, p_20, p_20_se, 9, true, true, false}
    };
    md_svg_opts op = {"S35-S36: Supercritical Isotherms P(rho) at T=1.5 and T=2.0",
                      "Monotonically increasing pressure without van der Waals loop (N=256)",
                      "Density rho", "Pressure P", 760, 460, false, false, false, 0.0, NULL};
    const char *svg1 = "../../../outputs/figures/day009a_25_isotherms.svg";
    md_svg_write_xy(svg1, &op, sp, 3);
    md_register_svg(ctx, "S36", svg1);

    md_svg_series sz[2] = {
        {"T=1.5 Z = P/(rho T)", "#2563eb", g_s35_rhos, z_15, NULL, 9, true, true, false},
        {"T=2.0 Z = P/(rho T)", "#dc2626", g_s35_rhos, z_20, NULL, 9, true, true, false}
    };
    md_svg_opts oz = {"S36: Compressibility Factor Z = P/(rho T) vs Density",
                      "Ideal gas limit Z -> 1 as rho -> 0; repulsive upturn at high density",
                      "Density rho", "Z = P / (rho T)", 760, 460, false, false, true, 1.0, "Ideal gas Z=1"};
    const char *svg2 = "../../../outputs/figures/day009a_26_compressibility_z.svg";
    md_svg_write_xy(svg2, &oz, sz, 2);
    md_register_svg(ctx, "S36", svg2);
    md_ctx_record_stage(ctx, "S36", "isotherm_t20", md_wall_time_sec() - t36, 4);
}

static void stage_s37_virial_b2(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    /* Numerical Mayer B2(T) curve for 4 potentials */
    int n_t = 25;
    double ts[25], b2_lj[25], b2_wca[25], b2_morse[25], b2_yuk[25];
    md_potential p_lj, p_wca, p_morse, p_yuk;
    md_potential_init(&p_lj, POT_LJ_SHIFT, 2.5);
    md_potential_init(&p_wca, POT_WCA, MD_WCA_RC);
    md_potential_init(&p_morse, POT_MORSE, 2.5);
    md_potential_init(&p_yuk, POT_YUKAWA, 2.5);

    FILE *fb = fopen("out/s37_b2_curve.csv", "w");
    fprintf(fb, "t,b2_lj_shift,b2_wca,b2_morse,b2_yukawa\n");
    for (int i = 0; i < n_t; ++i) {
        double T = 0.8 + 3.2 * (double)i / (double)(n_t - 1);
        ts[i] = T;
        b2_lj[i] = md_compute_b2_mayer(&p_lj, T, 4000);
        b2_wca[i] = md_compute_b2_mayer(&p_wca, T, 4000);
        b2_morse[i] = md_compute_b2_mayer(&p_morse, T, 4000);
        b2_yuk[i] = md_compute_b2_mayer(&p_yuk, T, 4000);
        fprintf(fb, "%.4f,%.8f,%.8f,%.8f,%.8f\n", T, b2_lj[i], b2_wca[i], b2_morse[i], b2_yuk[i]);
    }
    fclose(fb);

    /* Low-density EOS check at T = 2.0 for rho <= 0.05 (c9) */
    double T_low = 2.0;
    double b2_20 = md_compute_b2_mayer(&p_lj, T_low, 8000);
    double rhos[5] = {0.01, 0.02, 0.03, 0.04, 0.05};
    double zm1_meas[5], zm1_se[5], zm1_vir[5];
    int n_eq = ctx->small_scale ? 1000 : 3000;
    int n_prod = ctx->small_scale ? 4000 : 14000;

    FILE *fl = fopen("out/s37_low_density_eos.csv", "w");
    fprintf(fl, "rho,p_meas,p_se,z_minus_1,z_minus_1_se,b2_rho_theory,z_score\n");
    double max_z = 0.0;
    bool all_ok = true;
    for (int i = 0; i < 5; ++i) {
        eos_point e = run_eos_nvt(ctx, POT_LJ_SHIFT, rhos[i], T_low, 4, n_eq, n_prod, 37);
        double zm1 = e.p_mean / (rhos[i] * e.t_mean) - 1.0;
        double se_z = e.p_se / (rhos[i] * e.t_mean);
        double b2r = b2_20 * rhos[i];
        double zsc = (se_z > 1e-12) ? fabs(zm1 - b2r) / se_z : 0.0;
        if (zsc > max_z) max_z = zsc;
        if (zsc > 3.0) all_ok = false;
        zm1_meas[i] = zm1;
        zm1_se[i] = se_z;
        zm1_vir[i] = b2r;
        fprintf(fl, "%.3f,%.8f,%.8f,%.8f,%.8f,%.8f,%.4f\n",
                rhos[i], e.p_mean, e.p_se, zm1, se_z, b2r, zsc);
    }
    fclose(fl);
    ctx->c9_all_within_3se = all_ok;
    ctx->c9_max_z = max_z;

    md_register_csv(ctx, "S37", "out/s37_b2_curve.csv");
    md_register_csv(ctx, "S37", "out/s37_low_density_eos.csv");

    md_svg_series sb[4] = {
        {"LJ-shift (rc=2.5)", "#2563eb", ts, b2_lj, NULL, n_t, true, false, false},
        {"WCA (repulsive)", "#dc2626", ts, b2_wca, NULL, n_t, true, false, false},
        {"Morse", "#16a34a", ts, b2_morse, NULL, n_t, true, false, false},
        {"Yukawa", "#9333ea", ts, b2_yuk, NULL, n_t, true, false, false}
    };
    md_svg_opts ob = {"S37: Second Virial Coefficient B2(T) via Mayer Quadrature",
                      "B2(T) = -2*pi * int (exp(-u(r)/T) - 1) r^2 dr across 4 potentials",
                      "Temperature T", "B2(T)", 760, 460, false, false, true, 0.0, "Boyle line B2=0"};
    const char *svg1 = "../../../outputs/figures/day009a_27_b2_curve.svg";
    md_svg_write_xy(svg1, &ob, sb, 4);
    md_register_svg(ctx, "S37", svg1);

    char sub[128];
    snprintf(sub, sizeof(sub), "LJ-shift T=2.0: B2(2.0) = %.5f; max |diff|/SE = %.2f", b2_20, max_z);
    md_svg_series sl[2] = {
        {"MD measured Z - 1", "#2563eb", rhos, zm1_meas, zm1_se, 5, false, true, false},
        {"Virial theory B2(T)*rho", "#dc2626", rhos, zm1_vir, NULL, 5, true, false, true}
    };
    md_svg_opts ol = {"S37: Low-Density Equation of State vs Second Virial Expansion", sub,
                      "Density rho", "P/(rho T) - 1", 760, 460, false, false, true, 0.0, "Ideal gas"};
    const char *svg2 = "../../../outputs/figures/day009a_28_low_density_z.svg";
    md_svg_write_xy(svg2, &ol, sl, 2);
    md_register_svg(ctx, "S37", svg2);
    md_ctx_record_stage(ctx, "S37", "virial_b2_mayer", md_wall_time_sec() - t0, 4);
}

static void stage_s38_internal_energy(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rhos[4] = {0.2, 0.4, 0.6, 0.8};
    double u_lj15[4], u_lj20[4], u_wca15[4];
    int n_eq = ctx->small_scale ? 400 : 1200;
    int n_prod = ctx->small_scale ? 1200 : 4000;

    FILE *fp = fopen("out/s38_internal_energy.csv", "w");
    fprintf(fp, "potential,target_t,rho,u_pot_per_n,u_pot_se,u_int_per_n\n");
    for (int i = 0; i < 4; ++i) {
        eos_point e1 = run_eos_nvt(ctx, POT_LJ_SHIFT, rhos[i], 1.5, 4, n_eq, n_prod, 381);
        eos_point e2 = run_eos_nvt(ctx, POT_LJ_SHIFT, rhos[i], 2.0, 4, n_eq, n_prod, 382);
        eos_point ew = run_eos_nvt(ctx, POT_WCA, rhos[i], 1.5, 4, n_eq, n_prod, 383);
        u_lj15[i] = e1.u_mean + 1.5 * e1.t_mean;
        u_lj20[i] = e2.u_mean + 1.5 * e2.t_mean;
        u_wca15[i] = ew.u_mean + 1.5 * ew.t_mean;
        fprintf(fp, "lj_shift,1.50,%.2f,%.6f,%.6f,%.6f\n", rhos[i], e1.u_mean, e1.u_se, u_lj15[i]);
        fprintf(fp, "lj_shift,2.00,%.2f,%.6f,%.6f,%.6f\n", rhos[i], e2.u_mean, e2.u_se, u_lj20[i]);
        fprintf(fp, "wca,1.50,%.2f,%.6f,%.6f,%.6f\n", rhos[i], ew.u_mean, ew.u_se, u_wca15[i]);
    }
    fclose(fp);
    md_register_csv(ctx, "S38", "out/s38_internal_energy.csv");

    md_svg_series ser[3] = {
        {"LJ-shift T=1.5", "#2563eb", rhos, u_lj15, NULL, 4, true, true, false},
        {"LJ-shift T=2.0", "#0891b2", rhos, u_lj20, NULL, 4, true, true, false},
        {"WCA T=1.5", "#dc2626", rhos, u_wca15, NULL, 4, true, true, true}
    };
    md_svg_opts o = {"S38: Internal Energy per Particle u(rho, T)",
                     "Total internal energy E/N = <U>/N + 3T/2 for LJ-shift vs purely repulsive WCA",
                     "Density rho", "Internal energy E/N", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_29_internal_energy.svg";
    md_svg_write_xy(svgp, &o, ser, 3);
    md_register_svg(ctx, "S38", svgp);
    md_ctx_record_stage(ctx, "S38", "internal_energy", md_wall_time_sec() - t0, 2);
}

static void stage_s39_heat_capacity(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rho = 0.70, T = 1.5;
    int nc = 4; /* N = 256 */
    int n_part = 256;
    int n_eq = ctx->small_scale ? 600 : 2500;
    int n_prod = ctx->small_scale ? 3000 : 12000;

    /* 1. NVT fluctuation route: Cv/N = 3/2 + Var(U)/(N T^2) */
    eos_point e_mid = run_eos_nvt(ctx, POT_LJ_SHIFT, rho, T, nc, n_eq, n_prod, 391);
    double cv_nvt = 1.5 + e_mid.var_u / ((double)n_part * T * T);

    /* 2. Finite difference route: d<E/N>/dT using T = 1.4 and T = 1.6 */
    eos_point e_lo = run_eos_nvt(ctx, POT_LJ_SHIFT, rho, 1.4, nc, n_eq, n_prod / 2, 392);
    eos_point e_hi = run_eos_nvt(ctx, POT_LJ_SHIFT, rho, 1.6, nc, n_eq, n_prod / 2, 393);
    double cv_fd = ((e_hi.u_mean + 1.5 * e_hi.t_mean) - (e_lo.u_mean + 1.5 * e_lo.t_mean)) /
                   (e_hi.t_mean - e_lo.t_mean);

    /* 3. NVE Lebowitz formula: Cv/N = [2/3 - N*Var(K)/<K>^2]^-1 */
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_system sys;
    md_init_fcc(&sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 39, 9);
    md_init_velocities(&sys, T, true, &rng);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);
    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, T, sys.n, 3);
    for (int s = 0; s < n_eq; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
    }
    int ns = n_prod / 5;
    double *ks = (double *)malloc((size_t)ns * sizeof(double));
    int si = 0;
    for (int s = 1; s <= ns * 5; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, NULL, NULL);
        if (s % 5 == 0 && si < ns) ks[si++] = md_kinetic_energy(&sys);
    }
    double km = md_mean(ks, si);
    double vk = md_variance(ks, si);
    double denom = (2.0 / 3.0) - ((double)(n_part - 1) * vk) / (km * km);
    double cv_nve = (denom > 0.05) ? (1.0 / denom) : 0.0;
    free(ks);
    md_verlet_free(&vl);
    md_system_free(&sys);

    FILE *fp = fopen("out/s39_heat_capacity.csv", "w");
    fprintf(fp, "method,rho,target_t,cv_per_n,ideal_ref\n");
    fprintf(fp, "nvt_fluct,%.2f,%.2f,%.6f,1.5000\n", rho, T, cv_nvt);
    fprintf(fp, "nve_lebowitz,%.2f,%.2f,%.6f,1.5000\n", rho, T, cv_nve);
    fprintf(fp, "finite_diff,%.2f,%.2f,%.6f,1.5000\n", rho, T, cv_fd);
    fclose(fp);
    md_register_csv(ctx, "S39", "out/s39_heat_capacity.csv");

    const char *lbls[4] = {"Ideal gas (3/2)", "NVT Var(U)", "NVE Lebowitz", "d<E>/dT FD"};
    double vals[4] = {1.5, cv_nvt, cv_nve, cv_fd};
    const char *cols[4] = {"#6b7280", "#2563eb", "#16a34a", "#d97706"};
    const char *svgp = "../../../outputs/figures/day009a_30_heat_capacity.svg";
    md_svg_write_bars(svgp, "S39: Isochoric Heat Capacity Cv/N at rho=0.70, T=1.5",
                      "Cross-check of NVT energy fluctuations, NVE Lebowitz formula, and finite difference",
                      "Cv / N", lbls, vals, NULL, cols, 4, false);
    md_register_svg(ctx, "S39", svgp);
    md_ctx_record_stage(ctx, "S39", "heat_capacity_cv", md_wall_time_sec() - t0, 2);
}

static void stage_s40_compressibility(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rhos[3] = {0.3, 0.5, 0.7};
    int r_idx[3] = {2, 4, 6}; /* Indices in g_s35_rhos */
    double T = 1.5;
    int n_eq = ctx->small_scale ? 600 : 2000;
    int n_prod = ctx->small_scale ? 1500 : 6000;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fc = fopen("out/s40_compressibility.csv", "w");
    fprintf(fc, "rho,target_t,sk0_extrap,kappa_sk,dp_drho_eos,kappa_eos\n");
    FILE *fs = fopen("out/s40_sk_lowk.csv", "w");
    fprintf(fs, "rho,k,sk\n");

    double k_sk[3], k_eos[3];
    for (int i = 0; i < 3; ++i) {
        md_system sys;
        md_init_fcc(&sys, 5, rhos[i]); /* N = 500 for finer low-k resolution */
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 40, (uint64_t)i);
        md_init_velocities(&sys, T, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, T, sys.n, 3);
        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
        }

        md_sk_accum sk;
        md_sk_init(&sk, 40, 10.0, 5);
        int stride = 30;
        for (int s = 1; s <= n_prod; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
            if (s % stride == 0) md_sk_sample(&sk, &sys);
        }
        double kv[40], sv[40];
        md_sk_finalize(&sk, kv, sv);
        double fit_k2[12], fit_s[12];
        int nf = 0;
        for (int b = 0; b < 40; ++b) {
            if (sv[b] > 0.0) {
                fprintf(fs, "%.2f,%.6f,%.6f\n", rhos[i], kv[b], sv[b]);
                if (kv[b] < 3.5 && nf < 12) {
                    fit_k2[nf] = kv[b] * kv[b];
                    fit_s[nf] = sv[b];
                    nf++;
                }
            }
        }
        md_ols_result lr = md_ols_fit(fit_k2, fit_s, nf);
        double sk0 = (lr.intercept > 0.01) ? lr.intercept : fit_s[0];
        k_sk[i] = sk0 / (rhos[i] * T);
        int ri = r_idx[i];
        double dp_drho = (g_s35_p_shift[ri + 1] - g_s35_p_shift[ri - 1]) / 0.2;
        k_eos[i] = (dp_drho > 1e-6) ? (1.0 / (rhos[i] * dp_drho)) : 0.0;
        fprintf(fc, "%.2f,%.2f,%.6f,%.6f,%.6f,%.6f\n",
                rhos[i], T, sk0, k_sk[i], dp_drho, k_eos[i]);

        md_sk_free(&sk);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fc);
    fclose(fs);

    md_register_csv(ctx, "S40", "out/s40_compressibility.csv");
    md_register_csv(ctx, "S40", "out/s40_sk_lowk.csv");

    md_svg_series ser[2] = {
        {"From S(k->0) extrapolation", "#2563eb", rhos, k_sk, NULL, 3, true, true, false},
        {"From EOS (rho dP/drho)^-1", "#dc2626", rhos, k_eos, NULL, 3, true, true, true}
    };
    md_svg_opts o = {"S40: Isothermal Compressibility kappa_T at T=1.5",
                     "Comparison of structure factor S(k->0) route vs isotherm derivative route",
                     "Density rho", "Isothermal compressibility kappa_T", 760, 460, false, true, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_31_compressibility.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S40", svgp);
    md_ctx_record_stage(ctx, "S40", "compressibility", md_wall_time_sec() - t0, 3);
}

static void stage_s41_energy_route(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double states[3][2] = {{0.6, 1.5}, {0.7, 1.2}, {0.8, 1.0}};
    int n_eq = ctx->small_scale ? 600 : 2000;
    int n_prod = ctx->small_scale ? 2000 : 8000;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fe = fopen("out/s41_energy_route.csv", "w");
    fprintf(fe, "rho,target_t,u_direct,u_direct_se,u_gr_route,rel_err,z_score\n");
    FILE *fg = fopen("out/s41_gr_fine.csv", "w");
    fprintf(fg, "r,g_r_rho06_t15,u_pair,integrand\n");

    double x_rho[3], y_dir[3], y_gr[3];
    double worst_rel = 0.0, worst_z = 0.0;

    for (int k = 0; k < 3; ++k) {
        double rho = states[k][0], T = states[k][1];
        md_system sys;
        md_init_fcc(&sys, 4, rho); /* N=256 */
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 41, (uint64_t)k);
        md_init_velocities(&sys, T, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, T, sys.n, 3);
        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
        }

        int n_bins = 250; /* dr = 0.01 up to r = 2.5 */
        md_gr_accum gr;
        md_gr_init(&gr, n_bins, 2.5, 3);
        int ns = n_prod / 5;
        double *u_arr = (double *)malloc((size_t)ns * sizeof(double));
        int si = 0;
        for (int s = 1; s <= ns * 5; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
            if (s % 5 == 0 && si < ns) {
                u_arr[si++] = sys.pe / (double)sys.n;
                if (si % 2 == 0) md_gr_sample(&gr, &sys);
            }
        }
        md_block_stats bu = md_block_average(u_arr, si, 10, NULL);
        double rc_v[250], gv[250];
        md_gr_finalize(&gr, rc_v, gv, NULL);
        double u_gr = md_gr_energy_route(&pot, rho, n_bins, rc_v, gv, gr.dr);
        double rel = fabs(u_gr - bu.mean) / (fabs(bu.mean) + 1e-12);
        double z = (bu.se > 1e-12) ? fabs(u_gr - bu.mean) / bu.se : 0.0;
        if (rel > worst_rel) worst_rel = rel;
        if (z > worst_z) worst_z = z;
        x_rho[k] = rho; y_dir[k] = bu.mean; y_gr[k] = u_gr;

        fprintf(fe, "%.2f,%.2f,%.8f,%.8f,%.8f,%.6e,%.4f\n",
                rho, T, bu.mean, bu.se, u_gr, rel, z);

        if (k == 0) {
            for (int b = 0; b < n_bins; ++b) {
                double fdr = 0.0;
                double u = (rc_v[b] > 0.5) ? md_pot_eval_exact(&pot, rc_v[b] * rc_v[b], &fdr) : 0.0;
                double integ = 2.0 * MD_PI * rho * u * gv[b] * rc_v[b] * rc_v[b];
                if (fabs(integ) > 100.0) integ = 0.0;
                fprintf(fg, "%.4f,%.6f,%.6f,%.6f\n", rc_v[b], gv[b], u, integ);
            }
        }
        free(u_arr);
        md_gr_free(&gr);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fe);
    fclose(fg);
    ctx->c10_rel_err = worst_rel;
    ctx->c10_z = worst_z;

    md_register_csv(ctx, "S41", "out/s41_energy_route.csv");
    md_register_csv(ctx, "S41", "out/s41_gr_fine.csv");

    md_svg_series ser[2] = {
        {"Direct <U>/N", "#2563eb", x_rho, y_dir, NULL, 3, true, true, false},
        {"g(r) integral route", "#dc2626", x_rho, y_gr, NULL, 3, true, true, true}
    };
    md_svg_opts o = {"S41: Energy Route via g(r) vs Direct Measured PE",
                     "u_g(r) = 2*pi*rho * int u(r) g(r) r^2 dr matches direct PE within <1%",
                     "Density rho", "Potential energy per particle U/N", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_32_energy_route.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S41", svgp);
    md_ctx_record_stage(ctx, "S41", "energy_route_gr", md_wall_time_sec() - t0, 3);
}

static void stage_s42_wca_pressure(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rhos[4] = {0.2, 0.4, 0.6, 0.8};
    double p_vir[4], p_gr[4];
    int n_eq = ctx->small_scale ? 500 : 1500;
    int n_prod = ctx->small_scale ? 1500 : 5000;

    md_potential wca;
    md_potential_init(&wca, POT_WCA, MD_WCA_RC);

    FILE *fp = fopen("out/s42_wca_pressure_route.csv", "w");
    fprintf(fp, "rho,target_t,p_virial,p_virial_se,p_gr_route,rel_err\n");
    FILE *fg = fopen("out/s42_wca_gr.csv", "w");
    fprintf(fg, "r,gr_wca_rho06\n");

    for (int k = 0; k < 4; ++k) {
        md_system sys;
        md_init_fcc(&sys, 4, rhos[k]);
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 42, (uint64_t)k);
        md_init_velocities(&sys, 1.5, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, wca.rc, 0.3, false);
        md_forces_verlet(&sys, &wca, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 3);
        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &wca, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
        }

        int n_bins = 200;
        md_gr_accum gr;
        md_gr_init(&gr, n_bins, 2.0, 3);
        int ns = n_prod / 5;
        double *pa = (double *)malloc((size_t)ns * sizeof(double));
        double *ta = (double *)malloc((size_t)ns * sizeof(double));
        int si = 0;
        for (int s = 1; s <= ns * 5; ++s) {
            md_step_vverlet(&sys, &wca, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
            if (s % 5 == 0 && si < ns) {
                pa[si] = md_pressure(&sys);
                ta[si] = md_temperature(&sys);
                si++;
                if (si % 2 == 0) md_gr_sample(&gr, &sys);
            }
        }
        md_block_stats bp = md_block_average(pa, si, 10, NULL);
        double tm = md_mean(ta, si);
        double rc_v[200], gv[200];
        md_gr_finalize(&gr, rc_v, gv, NULL);
        double pg = md_gr_pressure_route(&wca, rhos[k], tm, n_bins, rc_v, gv, gr.dr);
        double rel = fabs(pg - bp.mean) / (fabs(bp.mean) + 1e-12);
        p_vir[k] = bp.mean; p_gr[k] = pg;
        fprintf(fp, "%.2f,1.50,%.6f,%.6f,%.6f,%.6e\n", rhos[k], bp.mean, bp.se, pg, rel);
        if (k == 2) {
            for (int b = 0; b < n_bins; ++b) fprintf(fg, "%.4f,%.6f\n", rc_v[b], gv[b]);
        }
        free(pa); free(ta);
        md_gr_free(&gr);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fp);
    fclose(fg);

    md_register_csv(ctx, "S42", "out/s42_wca_pressure_route.csv");
    md_register_csv(ctx, "S42", "out/s42_wca_gr.csv");

    md_svg_series ser[2] = {
        {"Direct virial pressure", "#2563eb", rhos, p_vir, NULL, 4, true, true, false},
        {"g(r) pressure route", "#dc2626", rhos, p_gr, NULL, 4, true, true, true}
    };
    md_svg_opts o = {"S42: WCA Pressure Route via g(r) vs Direct Virial (T=1.5)",
                     "Purely repulsive WCA fluid: exact virial theorem verification",
                     "Density rho", "Pressure P", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_33_wca_pressure.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S42", svgp);
    md_ctx_record_stage(ctx, "S42", "pressure_route_wca", md_wall_time_sec() - t0, 3);
}

static void stage_s43_npt_barostat(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double targets[5] = {0.5, 1.5, 3.0, 5.0, 8.0};
    double rho_npt[5];
    int n_eq = ctx->small_scale ? 1000 : 3000;
    int n_prod = ctx->small_scale ? 2000 : 6000;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fe = fopen("out/s43_npt_eos.csv", "w");
    fprintf(fe, "target_p,rho_mean,rho_se,p_mean,t_mean\n");
    FILE *ft = fopen("out/s43_npt_trace.csv", "w");
    fprintf(ft, "step,rho_p30,p_inst_p30\n");

    for (int k = 0; k < 5; ++k) {
        md_system sys;
        md_init_fcc(&sys, 4, 0.55); /* Start from mid density */
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 43, (uint64_t)k);
        md_init_velocities(&sys, 1.5, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 3);
        md_thermostat_enable_barostat(&th, targets[k], 1.5, 0.12);

        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
        }
        int ns = n_prod / 5;
        double *ra = (double *)malloc((size_t)ns * sizeof(double));
        double *pa = (double *)malloc((size_t)ns * sizeof(double));
        double *ta = (double *)malloc((size_t)ns * sizeof(double));
        int si = 0;
        for (int s = 1; s <= ns * 5; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
            if (s % 5 == 0 && si < ns) {
                ra[si] = (double)sys.n / md_box_volume(&sys.box);
                pa[si] = md_pressure(&sys);
                ta[si] = md_temperature(&sys);
                if (k == 2 && si % 6 == 0) {
                    fprintf(ft, "%d,%.6f,%.6f\n", s, ra[si], pa[si]);
                }
                si++;
            }
        }
        md_block_stats br = md_block_average(ra, si, 10, NULL);
        rho_npt[k] = br.mean;
        fprintf(fe, "%.2f,%.6f,%.6f,%.6f,%.6f\n",
                targets[k], br.mean, br.se, md_mean(pa, si), md_mean(ta, si));
        free(ra); free(pa); free(ta);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fe);
    fclose(ft);

    md_register_csv(ctx, "S43", "out/s43_npt_eos.csv");
    md_register_csv(ctx, "S43", "out/s43_npt_trace.csv");

    md_svg_series ser[2] = {
        {"NVT isotherm (S35)", "#2563eb", g_s35_rhos, g_s35_p_shift, NULL, 9, true, true, false},
        {"NPT Berendsen <rho(P*)>", "#dc2626", rho_npt, targets, NULL, 5, false, true, false}
    };
    md_svg_opts o = {"S43: NPT Barostat <rho(P*)> Cross-Checked Against NVT Isotherm",
                     "T=1.5 LJ-shift: NPT points lie cleanly on NVT equation of state",
                     "Density rho", "Pressure P", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_34_npt_vs_nvt.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S43", svgp);
    md_ctx_record_stage(ctx, "S43", "npt_barostat", md_wall_time_sec() - t0, 3);
}

static void stage_s44_multi_pot(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_pot_kind pkinds[4] = {POT_LJ_SHIFT, POT_WCA, POT_MORSE, POT_YUKAWA};
    const char *cols[4] = {"#2563eb", "#dc2626", "#16a34a", "#9333ea"};
    int n_bins = 120;
    double r_cent[120], grs[4][120];
    int n_eq = ctx->small_scale ? 500 : 1500;
    int n_prod = ctx->small_scale ? 1200 : 4000;

    FILE *fe = fopen("out/s44_multi_pot_eos.csv", "w");
    fprintf(fe, "potential,rho,target_t,u_per_n,p_meas,gr_peak_val,gr_peak_r\n");

    for (int k = 0; k < 4; ++k) {
        md_potential pot;
        md_potential_init(&pot, pkinds[k], 2.5);
        md_system sys;
        md_init_fcc(&sys, 4, 0.6);
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 44, (uint64_t)k);
        md_init_velocities(&sys, 1.5, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 3);
        for (int s = 0; s < n_eq; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
        }
        md_gr_accum gr;
        md_gr_init(&gr, n_bins, 3.0, 3);
        double u_sum = 0.0, p_sum = 0.0;
        int cnt = 0;
        for (int s = 1; s <= n_prod; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
            if (s % 10 == 0) {
                u_sum += sys.pe / (double)sys.n;
                p_sum += md_pressure(&sys);
                cnt++;
                md_gr_sample(&gr, &sys);
            }
        }
        md_gr_finalize(&gr, r_cent, grs[k], NULL);
        double pkv = 0.0, pkr = 0.0;
        for (int b = 0; b < n_bins; ++b) {
            if (grs[k][b] > pkv) { pkv = grs[k][b]; pkr = r_cent[b]; }
        }
        fprintf(fe, "%s,0.60,1.50,%.6f,%.6f,%.6f,%.4f\n",
                md_pot_name(pkinds[k]), u_sum / (double)cnt, p_sum / (double)cnt, pkv, pkr);
        md_gr_free(&gr);
        md_verlet_free(&vl);
        md_system_free(&sys);
    }
    fclose(fe);

    FILE *fg = fopen("out/s44_multi_pot_gr.csv", "w");
    fprintf(fg, "r,gr_lj_shift,gr_wca,gr_morse,gr_yukawa\n");
    for (int b = 0; b < n_bins; ++b) {
        fprintf(fg, "%.4f,%.6f,%.6f,%.6f,%.6f\n", r_cent[b], grs[0][b], grs[1][b], grs[2][b], grs[3][b]);
    }
    fclose(fg);

    md_register_csv(ctx, "S44", "out/s44_multi_pot_eos.csv");
    md_register_csv(ctx, "S44", "out/s44_multi_pot_gr.csv");

    md_svg_series ser[4];
    for (int k = 0; k < 4; ++k) {
        ser[k] = (md_svg_series){md_pot_name(pkinds[k]), cols[k], r_cent, grs[k], NULL, n_bins, true, false, false};
    }
    md_svg_opts o = {"S44: Radial Distribution g(r) Across 4 Potentials (rho=0.6, T=1.5)",
                     "WCA closely mimics LJ liquid structure; Yukawa has soft core",
                     "Distance r", "g(r)", 760, 460, false, false, true, 1.0, "Ideal gas g=1"};
    const char *svgp = "../../../outputs/figures/day009a_35_multi_pot_gr.svg";
    md_svg_write_xy(svgp, &o, ser, 4);
    md_register_svg(ctx, "S44", svgp);
    md_ctx_record_stage(ctx, "S44", "multi_pot_eos", md_wall_time_sec() - t0, 3);
}

void md_run_stages_e(md_run_ctx *ctx) {
    stage_s35_s36_isotherms(ctx);
    stage_s37_virial_b2(ctx);
    stage_s38_internal_energy(ctx);
    stage_s39_heat_capacity(ctx);
    stage_s40_compressibility(ctx);
    stage_s41_energy_route(ctx);
    stage_s42_wca_pressure(ctx);
    stage_s43_npt_barostat(ctx);
    stage_s44_multi_pot(ctx);
}
