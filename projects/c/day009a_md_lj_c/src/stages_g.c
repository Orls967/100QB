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

/* Group G: Phase Transition, Melting (Q6, Lindemann), Hysteresis & Cold-Solid Negative Control (S55..S58) */

void md_run_stages_g(md_run_ctx *ctx) {
    double temps[9] = {0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2};
    double rho = 1.0;
    int nc = 4; /* N = 256 */
    double a_fcc = cbrt(4.0 / rho);
    double r_nn = a_fcc * sqrt(0.5); /* 1.1225 */
    double r_bond = 1.38;            /* Between 1st (1.12) and 2nd (1.59) FCC shells */
    double dt = 0.002;
    int n_eq = ctx->small_scale ? 600 : 2000;
    int n_prod = ctx->small_scale ? 1200 : 4000;

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    double heat_q6[9], heat_q6_se[9], heat_lind[9], heat_u[9];
    double cool_q6[9], cool_u[9];
    double local_q6_solid[256], local_q6_liquid[256];

    /* S55: Heating FCC solid at rho = 1.0 */
    {
        double t0 = md_wall_time_sec();
        md_system sys;
        md_init_fcc(&sys, nc, rho);
        md_vec3 *ref_pos = (md_vec3 *)malloc((size_t)sys.n * sizeof(md_vec3));
        memcpy(ref_pos, sys.unwrapped, (size_t)sys.n * sizeof(md_vec3));

        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 55, 0);
        md_init_velocities(&sys, temps[0], true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, temps[0], sys.n, 3);

        FILE *fh = fopen("out/s55_heating.csv", "w");
        fprintf(fh, "target_t,q6_mean,q6_se,lindemann_mean,u_per_n,u_se,p_mean\n");
        FILE *fhb = fopen("out/s55_heating_blocks.csv", "w");
        fprintf(fhb, "target_t,block,q6_block\n");

        for (int ti = 0; ti < 9; ++ti) {
            th.target_t = temps[ti];
            for (int s = 0; s < n_eq; ++s) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
            }
            int stride = 20;
            int ns = n_prod / stride;
            double *q6_a = (double *)malloc((size_t)ns * sizeof(double));
            double *li_a = (double *)malloc((size_t)ns * sizeof(double));
            double *u_a = (double *)malloc((size_t)ns * sizeof(double));
            double *p_a = (double *)malloc((size_t)ns * sizeof(double));
            int si = 0;
            for (int s = 1; s <= ns * stride; ++s) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
                if (s % stride == 0 && si < ns) {
                    q6_a[si] = md_compute_q6(&sys, r_bond, NULL);
                    li_a[si] = md_compute_lindemann(&sys, ref_pos, r_nn);
                    u_a[si] = sys.pe / (double)sys.n;
                    p_a[si] = md_pressure(&sys);
                    si++;
                }
            }
            if (ti == 0) (void)md_compute_q6(&sys, r_bond, local_q6_solid);
            if (ti == 8) (void)md_compute_q6(&sys, r_bond, local_q6_liquid);

            double bq[10];
            md_block_stats bq6 = md_block_average(q6_a, si, 10, bq);
            md_block_stats bu = md_block_average(u_a, si, 10, NULL);
            heat_q6[ti] = bq6.mean;
            heat_q6_se[ti] = bq6.se;
            heat_lind[ti] = md_mean(li_a, si);
            heat_u[ti] = bu.mean;

            fprintf(fh, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    temps[ti], bq6.mean, bq6.se, heat_lind[ti], bu.mean, bu.se, md_mean(p_a, si));
            for (int b = 0; b < 10; ++b) fprintf(fhb, "%.2f,%d,%.6f\n", temps[ti], b, bq[b]);
            free(q6_a); free(li_a); free(u_a); free(p_a);
        }
        fclose(fh);
        fclose(fhb);

        /* Evaluate c14 Q6 drop between solid T=0.6 (idx 0) and melted liquid T=2.2 (idx 8, >= 1.4) */
        ctx->c14_q6_solid = heat_q6[0];
        ctx->c14_q6_liquid = heat_q6[8];
        double se_d = sqrt(heat_q6_se[0] * heat_q6_se[0] + heat_q6_se[8] * heat_q6_se[8]);
        ctx->c14_q6_drop_z = (se_d > 1e-12) ? ((heat_q6[0] - heat_q6[8]) / se_d) : 0.0;

        ctx->total_nbr_overflows += vl.overflow_count;
        ctx->total_late_rebuilds += vl.late_rebuild_count;

        md_register_csv(ctx, "S55", "out/s55_heating.csv");
        md_register_csv(ctx, "S55", "out/s55_heating_blocks.csv");

        /* Cap Lindemann at 0.8 for clean co-plot with Q6 */
        double lind_plot[9];
        for (int i = 0; i < 9; ++i) lind_plot[i] = (heat_lind[i] > 0.8) ? 0.8 : heat_lind[i];
        md_svg_series ser[2] = {
            {"Global Q6 order", "#2563eb", temps, heat_q6, heat_q6_se, 9, true, true, false},
            {"Lindemann delta_L (cap 0.8)", "#dc2626", temps, lind_plot, NULL, 9, true, true, true}
        };
        md_svg_opts o = {"S55: Heating FCC Solid at rho=1.0 — Melting Detection",
                         "Sharp drop in Q6 and surge in Lindemann parameter signal crystal melting",
                         "Temperature T", "Order parameter Q6 / Lindemann", 760, 460, false, false, true, 0.15, "Lindemann ~0.15"};
        const char *svgp = "../../../outputs/figures/day009a_44_melting_q6_lind.svg";
        md_svg_write_xy(svgp, &o, ser, 2);
        md_register_svg(ctx, "S55", svgp);

        /* Continue directly from melted T=2.2 state into S56 cooling! */
        md_ctx_record_stage(ctx, "S55", "melt_heating", md_wall_time_sec() - t0, 3);

        /* S56: Cooling liquid T=2.2 -> 0.6 hysteresis */
        double t56 = md_wall_time_sec();
        FILE *fc = fopen("out/s56_cooling.csv", "w");
        fprintf(fc, "target_t,q6_cool,q6_cool_se,u_cool_per_n,p_cool\n");
        for (int ti = 8; ti >= 0; --ti) {
            th.target_t = temps[ti];
            for (int s = 0; s < n_eq; ++s) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
            }
            int stride = 20;
            int ns = n_prod / stride;
            double *q6_a = (double *)malloc((size_t)ns * sizeof(double));
            double *u_a = (double *)malloc((size_t)ns * sizeof(double));
            double *p_a = (double *)malloc((size_t)ns * sizeof(double));
            int si = 0;
            for (int s = 1; s <= ns * stride; ++s) {
                md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
                if (s % stride == 0 && si < ns) {
                    q6_a[si] = md_compute_q6(&sys, r_bond, NULL);
                    u_a[si] = sys.pe / (double)sys.n;
                    p_a[si] = md_pressure(&sys);
                    si++;
                }
            }
            md_block_stats bq6 = md_block_average(q6_a, si, 10, NULL);
            cool_q6[ti] = bq6.mean;
            cool_u[ti] = md_mean(u_a, si);
            fprintf(fc, "%.2f,%.6f,%.6f,%.6f,%.6f\n",
                    temps[ti], bq6.mean, bq6.se, cool_u[ti], md_mean(p_a, si));
            free(q6_a); free(u_a); free(p_a);
        }
        fclose(fc);

        FILE *fhy = fopen("out/s56_hysteresis.csv", "w");
        fprintf(fhy, "target_t,u_heat,u_cool,u_hyst_diff,q6_heat,q6_cool\n");
        for (int ti = 0; ti < 9; ++ti) {
            fprintf(fhy, "%.2f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                    temps[ti], heat_u[ti], cool_u[ti], cool_u[ti] - heat_u[ti], heat_q6[ti], cool_q6[ti]);
        }
        fclose(fhy);

        md_register_csv(ctx, "S56", "out/s56_cooling.csv");
        md_register_csv(ctx, "S56", "out/s56_hysteresis.csv");

        md_svg_series sh[2] = {
            {"Heating branch (solid->liq)", "#dc2626", temps, heat_u, NULL, 9, true, true, false},
            {"Cooling branch (liq->supercool)", "#2563eb", temps, cool_u, NULL, 9, true, true, true}
        };
        md_svg_opts oh = {"S56: Heating-Cooling Hysteresis Loop at rho=1.0 (N=256)",
                          "Superheating and supercooling create clear first-order hysteresis loop",
                          "Temperature T", "Potential energy U/N", 760, 460, false, false, false, 0.0, NULL};
        const char *svg56 = "../../../outputs/figures/day009a_45_hysteresis.svg";
        md_svg_write_xy(svg56, &oh, sh, 2);
        md_register_svg(ctx, "S56", svg56);

        free(ref_pos);
        md_verlet_free(&vl);
        md_system_free(&sys);
        md_ctx_record_stage(ctx, "S56", "melt_cooling_hyst", md_wall_time_sec() - t56, 3);
    }

    /* S57: Negative control — cold solid at rho = 1.0, T = 0.4 does NOT melt */
    {
        double t0 = md_wall_time_sec();
        md_system sys;
        md_init_fcc(&sys, nc, rho);
        md_vec3 *ref_pos = (md_vec3 *)malloc((size_t)sys.n * sizeof(md_vec3));
        memcpy(ref_pos, sys.unwrapped, (size_t)sys.n * sizeof(md_vec3));
        md_rng rng;
        md_rng_seed_stage(&rng, MD_MASTER_SEED, 57, 0);
        md_init_velocities(&sys, 0.4, true, &rng);
        md_verlet_list vl;
        md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th;
        md_thermostat_init(&th, THERM_LANGEVIN, 0.4, sys.n, 3);

        int n_steps = ctx->small_scale ? 2000 : 8000;
        int stride = n_steps / 80;
        double *q6s = (double *)malloc(80 * sizeof(double));
        double *lds = (double *)malloc(80 * sizeof(double));
        FILE *ft = fopen("out/s57_cold_solid_trace.csv", "w");
        fprintf(ft, "step,t,q6,lindemann,u_per_n,temp\n");
        int si = 0;
        for (int s = 1; s <= n_steps; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
            if (s % stride == 0 && si < 80) {
                q6s[si] = md_compute_q6(&sys, r_bond, NULL);
                lds[si] = md_compute_lindemann(&sys, ref_pos, r_nn);
                fprintf(ft, "%d,%.4f,%.6f,%.6f,%.6f,%.6f\n",
                        s, (double)s * dt, q6s[si], lds[si], sys.pe / (double)sys.n, md_temperature(&sys));
                si++;
            }
        }
        fclose(ft);
        double q6_m = md_mean(q6s, si);
        double ld_m = md_mean(lds, si);
        ctx->c14_cold_q6 = q6_m;

        FILE *fs = fopen("out/s57_cold_solid_summary.csv", "w");
        fprintf(fs, "target_t,rho,n_steps,q6_mean,q6_min,lindemann_mean,stayed_crystalline\n");
        double qmin = q6s[0];
        for (int i = 1; i < si; ++i) if (q6s[i] < qmin) qmin = q6s[i];
        fprintf(fs, "0.40,1.00,%d,%.6f,%.6f,%.6f,%d\n", n_steps, q6_m, qmin, ld_m, (q6_m > 0.4) ? 1 : 0);
        fclose(fs);

        free(q6s); free(lds); free(ref_pos);
        md_verlet_free(&vl);
        md_system_free(&sys);
        md_register_csv(ctx, "S57", "out/s57_cold_solid_trace.csv");
        md_register_csv(ctx, "S57", "out/s57_cold_solid_summary.csv");
        md_ctx_record_stage(ctx, "S57", "melt_negctrl_cold", md_wall_time_sec() - t0, 2);
    }

    /* S58: Local per-particle q6(i) distribution solid vs liquid */
    {
        double t0 = md_wall_time_sec();
        int nb = 25;
        double qmin = 0.0, qmax = 0.75;
        double dq = (qmax - qmin) / (double)nb;
        double qc[25], hs[25] = {0}, hl[25] = {0};
        for (int b = 0; b < nb; ++b) qc[b] = qmin + ((double)b + 0.5) * dq;
        for (int i = 0; i < 256; ++i) {
            int bs = (int)(local_q6_solid[i] / dq);
            if (bs >= 0 && bs < nb) hs[bs] += 1.0 / (256.0 * dq);
            int bl = (int)(local_q6_liquid[i] / dq);
            if (bl >= 0 && bl < nb) hl[bl] += 1.0 / (256.0 * dq);
        }
        FILE *fp = fopen("out/s58_q6_hist.csv", "w");
        fprintf(fp, "q6_local,pdf_solid_t06,pdf_liquid_t22\n");
        for (int b = 0; b < nb; ++b) {
            fprintf(fp, "%.4f,%.6f,%.6f\n", qc[b], hs[b], hl[b]);
        }
        fclose(fp);
        md_register_csv(ctx, "S58", "out/s58_q6_hist.csv");

        md_svg_series ser[2] = {
            {"FCC solid (T=0.6)", "#2563eb", qc, hs, NULL, nb, true, true, false},
            {"Melted liquid (T=2.2)", "#dc2626", qc, hl, NULL, nb, true, true, false}
        };
        md_svg_opts o = {"S58: Local Bond-Order Distribution q6(i) — Solid vs Liquid",
                         "Clean separation between crystalline FCC peak (~0.52) and disordered liquid (~0.34)",
                         "Local Steinhardt q6(i)", "Probability density", 760, 460, false, false, false, 0.0, NULL};
        const char *svgp = "../../../outputs/figures/day009a_46_q6_hist.svg";
        md_svg_write_xy(svgp, &o, ser, 2);
        md_register_svg(ctx, "S58", svgp);
        md_ctx_record_stage(ctx, "S58", "melt_q6_dist", md_wall_time_sec() - t0, 2);
    }
}
