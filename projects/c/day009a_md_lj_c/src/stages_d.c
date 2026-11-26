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

/* Group D: Thermostats, Equipartition, Maxwell KS & Kinetic Energy Variance (S27..S34) */

typedef struct {
    double t_mean, t_se, k_mean, k_se, var_k_mean, var_k_se;
    double ks_d, ks_p, kx_mean, ky_mean, kz_mean, aux_val;
} therm_eval;

static therm_eval run_ensemble_point(md_run_ctx *ctx, md_thermostat_kind kind, double target_t,
                                     double tau_override, int nc, double rho, int n_eq, int n_prod,
                                     uint64_t seed_id, double *speed_hist_out, int n_vbins,
                                     double vmax, double *block_k_out) {
    therm_eval res;
    memset(&res, 0, sizeof(res));
    double dt = 0.0025;
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_system sys;
    md_init_fcc(&sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, seed_id, (uint64_t)llround(target_t * 100.0));
    md_init_velocities(&sys, target_t, true, &rng);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);

    md_thermostat th_eq;
    md_thermostat_init(&th_eq, THERM_LANGEVIN, target_t, sys.n, 3);
    for (int s = 0; s < n_eq; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th_eq, &rng);
    }

    md_thermostat th;
    md_thermostat_init(&th, kind, target_t, sys.n, 3);
    if (tau_override > 0.0 && kind == THERM_BERENDSEN) th.tau_t = tau_override;
    else if (kind == THERM_BERENDSEN) th.tau_t = 0.05;
    if (kind == THERM_LANGEVIN) th.gamma = 1.5;
    if (kind == THERM_ANDERSEN) th.nu_coll = 4.0;
    if (kind == THERM_NONE) {
        double cur_t = md_temperature(&sys);
        if (cur_t > 1e-6) {
            double sc = sqrt(target_t / cur_t);
            for (int i = 0; i < sys.n; ++i) sys.vel[i] = v3_scale(sys.vel[i], sc);
        }
    }
    for (int s = 0; s < n_eq / 2; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt,
                        (kind == THERM_NONE) ? NULL : &th, &rng);
    }

    int n_blocks = 10, sample_every = 5;
    int b_len = n_prod / n_blocks;
    if (b_len < 50) b_len = 50;
    int total_steps = n_blocks * b_len;
    int n_samp = total_steps / sample_every;

    double *k_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *kx_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *ky_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *kz_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *t_arr = (double *)malloc((size_t)n_samp * sizeof(double));
    double *h_arr = (double *)malloc((size_t)n_samp * sizeof(double));

    int n_ks_snaps = 8;
    int ks_stride = total_steps / n_ks_snaps;
    if (ks_stride < 1) ks_stride = 1;
    double *speeds = (double *)malloc((size_t)(n_ks_snaps * sys.n) * sizeof(double));
    int n_sp = 0, si = 0;

    for (int s = 1; s <= total_steps; ++s) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt,
                        (kind == THERM_NONE) ? NULL : &th, &rng);
        if (s % sample_every == 0 && si < n_samp) {
            double ke = md_kinetic_energy(&sys);
            md_kinetic_components(&sys, &kx_arr[si], &ky_arr[si], &kz_arr[si]);
            k_arr[si] = ke;
            t_arr[si] = md_temperature(&sys);
            h_arr[si] = (kind == THERM_NOSEHOOVER) ? md_nosehoover_conserved(&sys, &th) : (ke + sys.pe);
            si++;
        }
        if (s % ks_stride == 0 && n_sp + sys.n <= n_ks_snaps * sys.n) {
            for (int i = 0; i < sys.n; ++i) speeds[n_sp++] = v3_norm(sys.vel[i]);
        }
    }

    double b_k[10], b_var_jk[10];
    md_block_stats bs_k = md_block_average(k_arr, si, n_blocks, b_k);
    md_block_stats bs_t = md_block_average(t_arr, si, n_blocks, NULL);
    md_block_stats bs_kx = md_block_average(kx_arr, si, n_blocks, NULL);
    md_block_stats bs_ky = md_block_average(ky_arr, si, n_blocks, NULL);
    md_block_stats bs_kz = md_block_average(kz_arr, si, n_blocks, NULL);

    /* Unbiased full-run Var(K) + 10-block delete-1 jackknife SE */
    int sb_len = si / n_blocks;
    int n_use = sb_len * n_blocks;
    double full_var_k = md_variance(k_arr, n_use);
    double *jk_tmp = (double *)malloc((size_t)(n_use - sb_len) * sizeof(double));
    for (int b = 0; b < n_blocks; ++b) {
        int dst = 0;
        for (int j = 0; j < n_blocks; ++j) {
            if (j == b) continue;
            memcpy(jk_tmp + dst, k_arr + j * sb_len, (size_t)sb_len * sizeof(double));
            dst += sb_len;
        }
        b_var_jk[b] = md_variance(jk_tmp, dst);
        if (block_k_out) block_k_out[b] = b_k[b];
    }
    free(jk_tmp);
    double jk_m = md_mean(b_var_jk, n_blocks), jk_ss = 0.0;
    for (int b = 0; b < n_blocks; ++b) jk_ss += (b_var_jk[b] - jk_m) * (b_var_jk[b] - jk_m);

    res.k_mean = bs_k.mean; res.k_se = bs_k.se;
    res.t_mean = bs_t.mean; res.t_se = bs_t.se;
    res.var_k_mean = full_var_k;
    res.var_k_se = sqrt(((double)(n_blocks - 1) / (double)n_blocks) * jk_ss);
    res.kx_mean = bs_kx.mean; res.ky_mean = bs_ky.mean; res.kz_mean = bs_kz.mean;
    if (kind == THERM_NONE) {
        double k3 = bs_k.mean / 3.0;
        double zx = (bs_kx.se > 1e-12) ? fabs(bs_kx.mean - k3) / bs_kx.se : 0.0;
        double zy = (bs_ky.se > 1e-12) ? fabs(bs_ky.mean - k3) / bs_ky.se : 0.0;
        double zz = (bs_kz.se > 1e-12) ? fabs(bs_kz.mean - k3) / bs_kz.se : 0.0;
        res.aux_val = (zx > zy ? (zx > zz ? zx : zz) : (zy > zz ? zy : zz));
    } else {
        res.aux_val = sqrt(md_variance(h_arr, si));
    }
    md_ks_maxwell(speeds, n_sp, target_t, &res.ks_d, &res.ks_p);

    if (speed_hist_out && n_vbins > 0) {
        double dv = vmax / (double)n_vbins;
        memset(speed_hist_out, 0, (size_t)n_vbins * sizeof(double));
        for (int i = 0; i < n_sp; ++i) {
            int b = (int)(speeds[i] / dv);
            if (b >= 0 && b < n_vbins) speed_hist_out[b] += 1.0 / ((double)n_sp * dv);
        }
    }
    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;
    free(speeds); free(k_arr); free(kx_arr); free(ky_arr); free(kz_arr); free(t_arr); free(h_arr);
    md_verlet_free(&vl);
    md_system_free(&sys);
    return res;
}

static void run_canon_stage(md_run_ctx *ctx, const char *sid, const char *sname,
                            md_thermostat_kind kind, uint64_t seed_id,
                            const char *csv_main, const char *csv_sp, const char *sp_col,
                            const double temps[3], int nc, double rho, int n_eq, int n_prod,
                            const double v_cent[24], double vmax,
                            therm_eval out_arr[3], double hist_out[24]) {
    double t0 = md_wall_time_sec();
    int n_part = 4 * nc * nc * nc;
    FILE *fp = fopen(csv_main, "w");
    fprintf(fp, "target_t,t_mean,t_se,k_mean,k_se,k_eq_ref,var_k,var_k_se,var_k_canon,ks_d,ks_p,h_sd\n");
    for (int ti = 0; ti < 3; ++ti) {
        therm_eval e = run_ensemble_point(ctx, kind, temps[ti], 0.0, nc, rho, n_eq, n_prod,
                                          seed_id, (ti == 1) ? hist_out : NULL, 24, vmax, NULL);
        out_arr[ti] = e;
        double keq = 1.5 * (double)(n_part - 1) * temps[ti];
        double vkc = 1.5 * (double)n_part * temps[ti] * temps[ti];
        fprintf(fp, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6e\n",
                temps[ti], e.t_mean, e.t_se, e.k_mean, e.k_se, keq,
                e.var_k_mean, e.var_k_se, vkc, e.ks_d, e.ks_p, e.aux_val);
    }
    fclose(fp);
    FILE *fs = fopen(csv_sp, "w");
    fprintf(fs, "v,%s,pdf_maxwell_theory\n", sp_col);
    for (int b = 0; b < 24; ++b) {
        double v = v_cent[b];
        fprintf(fs, "%.4f,%.6f,%.6f\n", v, hist_out[b], sqrt(2.0 / MD_PI) * v * v * exp(-0.5 * v * v));
    }
    fclose(fs);
    md_register_csv(ctx, sid, csv_main);
    md_register_csv(ctx, sid, csv_sp);
    md_ctx_record_stage(ctx, sid, sname, md_wall_time_sec() - t0, 2);
}

void md_run_stages_d(md_run_ctx *ctx) {
    double temps[3] = {0.7, 1.0, 1.5};
    int nc = 4, n_part = 256;
    double rho = 0.80, vmax = 4.8, dv = 4.8 / 24.0;
    int n_eq = ctx->small_scale ? 1000 : 2500;
    int n_prod = ctx->small_scale ? 4000 : 18000;

    therm_eval res_all[5][3];
    double hist_lang[24] = {0}, hist_nh[24] = {0}, hist_and[24] = {0}, v_cent[24];
    for (int b = 0; b < 24; ++b) v_cent[b] = ((double)b + 0.5) * dv;

    /* S27: NVE */
    {
        double t0 = md_wall_time_sec();
        FILE *fp = fopen("out/s27_nve_stats.csv", "w");
        fprintf(fp, "target_t,t_mean,t_se,k_mean,k_se,kx_mean,ky_mean,kz_mean,dir_eq_z,var_k,var_k_se,var_k_canon\n");
        FILE *fb = fopen("out/s27_nve_blocks.csv", "w");
        fprintf(fb, "target_t,block,k_block_mean\n");
        for (int ti = 0; ti < 3; ++ti) {
            double bk[10];
            therm_eval e = run_ensemble_point(ctx, THERM_NONE, temps[ti], 0.0, nc, rho,
                                              n_eq, n_prod, 27, NULL, 0, vmax, bk);
            res_all[0][ti] = e;
            double vkc = 1.5 * (double)n_part * temps[ti] * temps[ti];
            fprintf(fp, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.4f,%.6f,%.6f,%.6f\n",
                    temps[ti], e.t_mean, e.t_se, e.k_mean, e.k_se,
                    e.kx_mean, e.ky_mean, e.kz_mean, e.aux_val, e.var_k_mean, e.var_k_se, vkc);
            for (int b = 0; b < 10; ++b) fprintf(fb, "%.2f,%d,%.6f\n", temps[ti], b, bk[b]);
        }
        fclose(fp); fclose(fb);
        md_register_csv(ctx, "S27", "out/s27_nve_stats.csv");
        md_register_csv(ctx, "S27", "out/s27_nve_blocks.csv");
        md_ctx_record_stage(ctx, "S27", "nve_equipart", md_wall_time_sec() - t0, 2);
    }

    /* S28, S29, S31: Canonical thermostats */
    run_canon_stage(ctx, "S28", "therm_langevin", THERM_LANGEVIN, 28,
                    "out/s28_langevin.csv", "out/s28_langevin_speeds.csv", "pdf_langevin",
                    temps, nc, rho, n_eq, n_prod, v_cent, vmax, res_all[1], hist_lang);
    run_canon_stage(ctx, "S29", "therm_nosehoover", THERM_NOSEHOOVER, 29,
                    "out/s29_nosehoover.csv", "out/s29_nh_speeds.csv", "pdf_nosehoover",
                    temps, nc, rho, n_eq, n_prod, v_cent, vmax, res_all[2], hist_nh);

    /* S30: Berendsen */
    {
        double t0 = md_wall_time_sec();
        FILE *fp = fopen("out/s30_berendsen.csv", "w");
        fprintf(fp, "target_t,t_mean,t_se,k_mean,k_se,k_eq_ref,var_k,var_k_se,var_k_canon,var_ratio\n");
        for (int ti = 0; ti < 3; ++ti) {
            therm_eval e = run_ensemble_point(ctx, THERM_BERENDSEN, temps[ti], 0.05, nc, rho,
                                              n_eq, n_prod, 30, NULL, 0, vmax, NULL);
            res_all[3][ti] = e;
            double keq = 1.5 * (double)(n_part - 1) * temps[ti];
            double vkc = 1.5 * (double)n_part * temps[ti] * temps[ti];
            fprintf(fp, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    temps[ti], e.t_mean, e.t_se, e.k_mean, e.k_se, keq,
                    e.var_k_mean, e.var_k_se, vkc, e.var_k_mean / vkc);
        }
        fclose(fp);
        double taus[3] = {0.02, 0.1, 0.5};
        FILE *ft = fopen("out/s30_berendsen_tau.csv", "w");
        fprintf(ft, "tau_t,t_mean,var_k,var_k_se,var_ratio_canon\n");
        for (int k = 0; k < 3; ++k) {
            therm_eval e = run_ensemble_point(ctx, THERM_BERENDSEN, 1.0, taus[k], nc, rho,
                                              n_eq / 2, n_prod / 2, (uint64_t)(300 + k), NULL, 0, vmax, NULL);
            fprintf(ft, "%.3f,%.6f,%.6f,%.6f,%.6f\n",
                    taus[k], e.t_mean, e.var_k_mean, e.var_k_se, e.var_k_mean / (1.5 * (double)n_part));
        }
        fclose(ft);
        md_register_csv(ctx, "S30", "out/s30_berendsen.csv");
        md_register_csv(ctx, "S30", "out/s30_berendsen_tau.csv");
        md_ctx_record_stage(ctx, "S30", "therm_berendsen", md_wall_time_sec() - t0, 2);
    }

    run_canon_stage(ctx, "S31", "therm_andersen", THERM_ANDERSEN, 31,
                    "out/s31_andersen.csv", "out/s31_andersen_speeds.csv", "pdf_andersen",
                    temps, nc, rho, n_eq, n_prod, v_cent, vmax, res_all[4], hist_and);

    /* S32: Var(K) & Equipartition Summary */
    {
        double t0 = md_wall_time_sec();
        const char *enames[5] = {"NVE", "Langevin", "NoseHoover", "Berendsen", "Andersen"};
        const char *ecols[5] = {"#6b7280", "#2563eb", "#16a34a", "#dc2626", "#9333ea"};
        FILE *fp = fopen("out/s32_var_k_summary.csv", "w");
        fprintf(fp, "ensemble,target_t,k_z_score,var_k_ratio,var_k_ratio_se,z_below_canon\n");

        double max_z_c6 = 0.0, max_canon_rel = 0.0, min_noncanon_z = 1e9;
        bool all_c6 = true, canon_ok = true;
        double bar_vals[5], bar_errs[5];

        for (int e = 0; e < 5; ++e) {
            for (int ti = 0; ti < 3; ++ti) {
                double T = temps[ti];
                double keq = 1.5 * (double)(n_part - 1) * T;
                double vkc = 1.5 * (double)n_part * T * T;
                therm_eval r = res_all[e][ti];
                double zk = (e == 0) ? r.aux_val : ((r.k_se > 1e-12) ? fabs(r.k_mean - keq) / r.k_se : 0.0);
                double ratio = r.var_k_mean / vkc;
                double ratio_se = r.var_k_se / vkc;
                double z_below = (ratio_se > 1e-12) ? (1.0 - ratio) / ratio_se : 0.0;
                fprintf(fp, "%s,%.2f,%.4f,%.6f,%.6f,%.4f\n", enames[e], T, zk, ratio, ratio_se, z_below);

                if (ti == 1) {
                    bar_vals[e] = ratio; bar_errs[e] = ratio_se;
                    if (zk > max_z_c6) max_z_c6 = zk;
                    if (zk > 3.0) all_c6 = false;
                    if (e == 1 || e == 2 || e == 4) {
                        double rel = fabs(ratio - 1.0);
                        if (rel > max_canon_rel) max_canon_rel = rel;
                        if (rel > 0.10) canon_ok = false;
                    } else {
                        if (e == 0) ctx->c8_nve_z = z_below;
                        if (e == 3) ctx->c8_berendsen_z = z_below;
                        if (z_below < min_noncanon_z) min_noncanon_z = z_below;
                    }
                }
            }
        }
        fclose(fp);
        ctx->c6_all_within_3se = all_c6; ctx->c6_max_z = max_z_c6;
        ctx->c8_canon_within_10pct = canon_ok; ctx->c8_canon_max_rel = max_canon_rel;
        ctx->c8_noncanon_smaller_3se = (min_noncanon_z >= 3.0);

        md_register_csv(ctx, "S32", "out/s32_var_k_summary.csv");
        const char *svgp = "../../../outputs/figures/day009a_21_var_k_ensembles.svg";
        md_svg_write_bars(svgp, "S32: Normalized Kinetic Energy Variance Var(K) / (3N T^2 / 2)",
                          "T=1.0, rho=0.80, N=256: Canonical thermostats ~1.0; NVE & Berendsen suppressed",
                          "Var(K) / Var_canon(K)", enames, bar_vals, bar_errs, ecols, 5, false);
        md_register_svg(ctx, "S32", svgp);
        md_ctx_record_stage(ctx, "S32", "therm_var_k_compare", md_wall_time_sec() - t0, 2);
    }

    /* S33: KS Maxwell summary */
    {
        double t0 = md_wall_time_sec();
        const char *cnames[3] = {"Langevin", "NoseHoover", "Andersen"};
        int eidx[3] = {1, 2, 4};
        double min_p = 1.0;
        bool all_p = true;
        FILE *fp = fopen("out/s33_ks_summary.csv", "w");
        fprintf(fp, "ensemble,target_t,ks_d,ks_p,pass_p001\n");
        for (int c = 0; c < 3; ++c) {
            for (int ti = 0; ti < 3; ++ti) {
                therm_eval r = res_all[eidx[c]][ti];
                if (ti == 1) {
                    if (r.ks_p < min_p) min_p = r.ks_p;
                    if (r.ks_p <= 0.01) all_p = false;
                }
                fprintf(fp, "%s,%.2f,%.6f,%.6f,%d\n", cnames[c], temps[ti], r.ks_d, r.ks_p, r.ks_p > 0.01 ? 1 : 0);
            }
        }
        fclose(fp);
        ctx->c7_min_ks_p = min_p; ctx->c7_all_pass = all_p;
        md_register_csv(ctx, "S33", "out/s33_ks_summary.csv");

        double cdf_l[24], cdf_n[24], cdf_a[24], cdf_th[24], pdf_th[24];
        double cl = 0.0, cn = 0.0, ca = 0.0;
        for (int b = 0; b < 24; ++b) {
            cl += hist_lang[b] * dv; cn += hist_nh[b] * dv; ca += hist_and[b] * dv;
            cdf_l[b] = cl; cdf_n[b] = cn; cdf_a[b] = ca;
            cdf_th[b] = md_maxwell_speed_cdf(v_cent[b], 1.0);
            pdf_th[b] = sqrt(2.0 / MD_PI) * v_cent[b] * v_cent[b] * exp(-0.5 * v_cent[b] * v_cent[b]);
        }
        md_svg_series scdf[4] = {
            {"Maxwell exact CDF", "#111827", v_cent, cdf_th, NULL, 24, true, false, true},
            {"Langevin BAOAB", "#2563eb", v_cent, cdf_l, NULL, 24, true, false, false},
            {"Nose-Hoover", "#16a34a", v_cent, cdf_n, NULL, 24, true, false, false},
            {"Andersen pair", "#9333ea", v_cent, cdf_a, NULL, 24, true, false, false}
        };
        md_svg_opts oc = {"S33: Speed CDF vs Exact Maxwell-Boltzmann (T=1.0)",
                          "Canonical thermostats match exact 3D Maxwell speed distribution",
                          "Speed v", "CDF F(v)", 760, 460, false, false, false, 0.0, NULL};
        const char *svg1 = "../../../outputs/figures/day009a_22_maxwell_cdf.svg";
        md_svg_write_xy(svg1, &oc, scdf, 4);
        md_register_svg(ctx, "S33", svg1);
        const char *svg2 = "../../../outputs/figures/day009a_23_speed_hist.svg";
        md_svg_write_hist(svg2, "S33: Langevin BAOAB Speed Distribution (T=1.0)",
                          "Empirical speed PDF vs exact Maxwell-Boltzmann curve",
                          "Speed v", "f(v)", v_cent, hist_lang, 24, pdf_th, "Maxwell theory T=1.0");
        md_register_svg(ctx, "S33", svg2);
        md_ctx_record_stage(ctx, "S33", "therm_ks_maxwell", md_wall_time_sec() - t0, 3);
    }

    /* S34: Thermal relaxation T0=2.5 -> 1.0 */
    {
        double t0 = md_wall_time_sec();
        md_thermostat_kind kinds[4] = {THERM_LANGEVIN, THERM_NOSEHOOVER, THERM_BERENDSEN, THERM_ANDERSEN};
        const char *names[4] = {"Langevin", "NoseHoover", "Berendsen", "Andersen"};
        const char *cols[4] = {"#2563eb", "#16a34a", "#dc2626", "#9333ea"};
        double tx[120], tr[4][120];
        md_potential pot;
        md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
        FILE *ff = fopen("out/s34_relaxation_fit.csv", "w");
        fprintf(ff, "ensemble,t_init,t_final,tau_half_est\n");
        for (int k = 0; k < 4; ++k) {
            md_system sys;
            md_init_fcc(&sys, 4, 0.8);
            md_rng rng;
            md_rng_seed_stage(&rng, MD_MASTER_SEED, 34, (uint64_t)k);
            md_init_velocities(&sys, 2.5, true, &rng);
            md_verlet_list vl;
            md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
            md_forces_verlet(&sys, &pot, &vl);
            md_thermostat th;
            md_thermostat_init(&th, kinds[k], 1.0, sys.n, 3);
            th.tau_t = 0.2; th.gamma = 2.0; th.nu_coll = 4.0;
            double t_half = -1.0;
            int ti = 0;
            for (int st = 1; st <= 1200; ++st) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.0025, &th, &rng);
                if (st % 10 == 0 && ti < 120) {
                    tx[ti] = (double)st * 0.0025;
                    tr[k][ti] = md_temperature(&sys);
                    if (t_half < 0.0 && tr[k][ti] <= 1.75) t_half = tx[ti];
                    ti++;
                }
            }
            fprintf(ff, "%s,2.5000,%.6f,%.4f\n", names[k], tr[k][119], t_half < 0.0 ? tx[119] : t_half);
            md_verlet_free(&vl);
            md_system_free(&sys);
        }
        fclose(ff);
        FILE *ft = fopen("out/s34_relaxation_trace.csv", "w");
        fprintf(ft, "t,temp_langevin,temp_nosehoover,temp_berendsen,temp_andersen\n");
        for (int i = 0; i < 120; ++i) fprintf(ft, "%.4f,%.6f,%.6f,%.6f,%.6f\n", tx[i], tr[0][i], tr[1][i], tr[2][i], tr[3][i]);
        fclose(ft);
        md_register_csv(ctx, "S34", "out/s34_relaxation_trace.csv");
        md_register_csv(ctx, "S34", "out/s34_relaxation_fit.csv");
        md_svg_series ser[4];
        for (int k = 0; k < 4; ++k) ser[k] = (md_svg_series){names[k], cols[k], tx, tr[k], NULL, 120, true, false, false};
        md_svg_opts o = {"S34: Thermal Relaxation T0=2.5 -> T_target=1.0",
                         "Berendsen decays monotonically; Nose-Hoover exhibits inertial overshoot",
                         "Time t", "Instantaneous temperature T(t)", 760, 460, false, false, true, 1.0, "T_target = 1.0"};
        const char *svgp = "../../../outputs/figures/day009a_24_therm_relax.svg";
        md_svg_write_xy(svgp, &o, ser, 4);
        md_register_svg(ctx, "S34", svgp);
        md_ctx_record_stage(ctx, "S34", "therm_relaxation", md_wall_time_sec() - t0, 3);
    }
}
