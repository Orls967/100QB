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

/* Group C: Symplectic Integration, Reversibility, Order, Momentum, Chaos & r-RESPA (S19..S26) */

static void prep_liquid_state(md_system *sys, int nc, double rho, double temp, uint64_t stage_id) {
    md_init_fcc(sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, stage_id, 101);
    md_init_velocities(sys, temp, true, &rng);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_verlet_list vl;
    md_verlet_init(&vl, sys->n, pot.rc, 0.3, false);
    md_forces_verlet(sys, &pot, &vl);
    md_thermostat th;
    md_thermostat_init(&th, THERM_LANGEVIN, temp, sys->n, 3);
    for (int s = 0; s < 400; ++s) {
        md_step_vverlet(sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
    }
    md_zero_momentum(sys);
    md_forces_brute(sys, &pot);
    md_verlet_free(&vl);
}

static void run_reversibility_leg(const md_system *base, const md_potential *pot, bool use_euler,
                                  double dt, int n_steps, double *step_x, double *trace_out,
                                  double *max_p_out, double *max_v_out, double *de_out) {
    md_system s;
    memset(&s, 0, sizeof(s));
    md_system_copy(&s, base);
    double e0 = md_kinetic_energy(&s) + s.pe;
    md_vec3 *fwd = trace_out ? (md_vec3 *)malloc((size_t)(100 * s.n) * sizeof(md_vec3)) : NULL;

    for (int st = 0; st < n_steps; ++st) {
        if (fwd && st % 5 == 0) {
            memcpy(fwd + (size_t)(st / 5) * (size_t)s.n, s.unwrapped, (size_t)s.n * sizeof(md_vec3));
        }
        if (use_euler) md_step_euler(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt);
        else md_step_vverlet(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt, NULL, NULL);
    }
    for (int i = 0; i < s.n; ++i) s.vel[i] = v3_scale(s.vel[i], -1.0);
    for (int st = 0; st < n_steps; ++st) {
        if (use_euler) md_step_euler(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt);
        else md_step_vverlet(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt, NULL, NULL);
        int rem = n_steps - 1 - st;
        if (fwd && rem % 5 == 0) {
            int idx = rem / 5;
            double mpe = 0.0;
            const md_vec3 *ref = fwd + (size_t)idx * (size_t)s.n;
            for (int i = 0; i < s.n; ++i) {
                double d = v3_norm(v3_sub(s.unwrapped[i], ref[i]));
                if (d > mpe) mpe = d;
            }
            trace_out[99 - idx] = (mpe > 1e-16) ? mpe : 1e-16;
            if (step_x) step_x[99 - idx] = (double)(st + 1);
        }
    }
    for (int i = 0; i < s.n; ++i) s.vel[i] = v3_scale(s.vel[i], -1.0);

    double max_p = 0.0, max_v = 0.0;
    for (int i = 0; i < s.n; ++i) {
        double dp = v3_norm(v3_sub(s.unwrapped[i], base->unwrapped[i]));
        double dv = v3_norm(v3_sub(s.vel[i], base->vel[i]));
        if (!isfinite(dp) || dp > 1e8) dp = 1e8;
        if (!isfinite(dv) || dv > 1e8) dv = 1e8;
        if (dp > max_p) max_p = dp;
        if (dv > max_v) max_v = dv;
    }
    *max_p_out = max_p;
    *max_v_out = max_v;
    double de = fabs((md_kinetic_energy(&s) + s.pe) - e0);
    *de_out = (!isfinite(de) || de > 1e8) ? 1e8 : de;
    free(fwd);
    md_system_free(&s);
}

static void stage_s19_reversibility(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_system base;
    prep_liquid_state(&base, 3, 0.7, 1.2, 19);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    double dts[3] = {0.0005, 0.001, 0.002};
    double step_x[100], vv_tr[100], eu_tr[100];
    FILE *fp = fopen("out/s19_reversibility.csv", "w");
    fprintf(fp, "integrator,dt,steps_fwd,max_pos_err,max_vel_err,energy_err\n");
    FILE *ft = fopen("out/s19_rev_trace.csv", "w");
    fprintf(ft, "return_step,vv_pos_err_dt002,euler_pos_err_dt002\n");

    for (int k = 0; k < 3; ++k) {
        double mp, mv, de;
        run_reversibility_leg(&base, &pot, false, dts[k], 500, (k == 2) ? step_x : NULL,
                              (k == 2) ? vv_tr : NULL, &mp, &mv, &de);
        fprintf(fp, "vverlet,%.4f,500,%.6e,%.6e,%.6e\n", dts[k], mp, mv, de);
        if (k == 2) ctx->c4_rev_max_pos_err = mp;

        run_reversibility_leg(&base, &pot, true, dts[k], 500, NULL,
                              (k == 2) ? eu_tr : NULL, &mp, &mv, &de);
        fprintf(fp, "euler,%.4f,500,%.6e,%.6e,%.6e\n", dts[k], mp, mv, de);
    }
    for (int i = 0; i < 100; ++i) fprintf(ft, "%.0f,%.6e,%.6e\n", step_x[i], vv_tr[i], eu_tr[i]);
    fclose(fp); fclose(ft);
    md_system_free(&base);

    md_register_csv(ctx, "S19", "out/s19_reversibility.csv");
    md_register_csv(ctx, "S19", "out/s19_rev_trace.csv");
    md_svg_series ser[2] = {
        {"Velocity Verlet (reversible)", "#2563eb", step_x, vv_tr, NULL, 100, true, false, false},
        {"Explicit Euler (neg control)", "#dc2626", step_x, eu_tr, NULL, 100, true, false, false}
    };
    md_svg_opts o = {"S19: Time Reversibility (N=108, 500 Fwd + 500 Bwd, dt=0.002)",
                     "Max position return error along reverse leg (c4 threshold 1e-9)",
                     "Backward step", "Max position error", 760, 460, false, true, true, 1e-9, "c4 limit (1e-9)"};
    const char *svgp = "../../../outputs/figures/day009a_14_reversibility.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S19", svgp);
    md_ctx_record_stage(ctx, "S19", "time_reversibility", md_wall_time_sec() - t0, 3);
}

static md_ols_result run_dt_sweep(const md_system *base, const md_potential *pot, md_integrator_kind kind,
                                  const double dts[5], double t_total, const char *csv_dt, const char *csv_tr,
                                  const char *col_tr, const double *ref_sig, double sig_out[5],
                                  double *max_de_002_out, double trace_t[160], double trace_de[160]) {
    FILE *fd = fopen(csv_dt, "w");
    if (ref_sig) fprintf(fd, "dt,steps,e0,sigma_e,max_abs_de,end_drift,ratio_vs_vv\n");
    else fprintf(fd, "dt,steps,e0,sigma_e,max_abs_de,end_drift\n");
    FILE *ft = fopen(csv_tr, "w");
    fprintf(ft, "t,%s\n", col_tr);

    double log_dt[5], log_s[5];
    for (int k = 0; k < 5; ++k) {
        double dt = dts[k];
        int steps = (int)llround(t_total / dt);
        md_system s;
        memset(&s, 0, sizeof(s));
        md_system_copy(&s, base);
        double e0 = md_kinetic_energy(&s) + s.pe;
        double *es = (double *)malloc((size_t)steps * sizeof(double));
        double max_de = 0.0;
        int stride = steps / 160;
        if (stride < 1) stride = 1;
        int ti = 0;
        for (int st = 0; st < steps; ++st) {
            if (kind == INT_VVERLET) md_step_vverlet(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt, NULL, NULL);
            else if (kind == INT_LEAPFROG) md_step_leapfrog(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt);
            else md_step_euler(&s, pot, NBR_BRUTE, NULL, NULL, NULL, 1, dt);
            double e = md_kinetic_energy(&s) + s.pe;
            if (!isfinite(e) || fabs(e) > 1e8) e = 1e8;
            es[st] = e;
            double de = fabs(e - e0);
            if (de > max_de) max_de = de;
            if (k == 2 && (st % stride == 0) && ti < 160) {
                double tv = (double)(st + 1) * dt;
                if (trace_t) trace_t[ti] = tv;
                trace_de[ti] = e - e0;
                fprintf(ft, "%.4f,%.8e\n", tv, trace_de[ti]);
                ti++;
            }
        }
        double sig = sqrt(md_variance(es, steps));
        double drift = fabs(es[steps - 1] - e0);
        sig_out[k] = sig;
        log_dt[k] = log(dt);
        log_s[k] = log(sig);
        if (k == 2 && max_de_002_out) *max_de_002_out = max_de;
        if (ref_sig) {
            fprintf(fd, "%.4f,%d,%.8f,%.6e,%.6e,%.6e,%.2f\n",
                    dt, steps, e0, sig, max_de, drift, sig / (ref_sig[k] + 1e-30));
        } else {
            fprintf(fd, "%.4f,%d,%.8f,%.6e,%.6e,%.6e\n", dt, steps, e0, sig, max_de, drift);
        }
        free(es);
        md_system_free(&s);
    }
    fclose(fd);
    fclose(ft);
    return md_ols_fit(log_dt, log_s, 5);
}

static void stage_s20_s21_s22_integrator_order(md_run_ctx *ctx) {
    md_system base;
    prep_liquid_state(&base, 3, 0.7, 1.2, 20);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    double dts[5] = {0.0005, 0.001, 0.002, 0.004, 0.008};
    double sig_vv[5], sig_lf[5], sig_eu[5];
    double drift_vv = 0.0, drift_eu = 0.0;
    double trace_t[160], tr_vv[160], tr_lf[160], tr_eu[160];

    double t20 = md_wall_time_sec();
    md_ols_result ovv = run_dt_sweep(&base, &pot, INT_VVERLET, dts, 1.6,
                                     "out/s20_vv_dt.csv", "out/s20_vv_trace.csv", "de_vv_dt002",
                                     NULL, sig_vv, &drift_vv, trace_t, tr_vv);
    ctx->c5_vv_slope = ovv.slope;
    ctx->c5_vv_slope_se = ovv.slope_se;
    md_register_csv(ctx, "S20", "out/s20_vv_dt.csv");
    md_register_csv(ctx, "S20", "out/s20_vv_trace.csv");
    md_ctx_record_stage(ctx, "S20", "energy_dt_vv", md_wall_time_sec() - t20, 2);

    double t21 = md_wall_time_sec();
    md_ols_result olf = run_dt_sweep(&base, &pot, INT_LEAPFROG, dts, 1.6,
                                     "out/s21_leapfrog_dt.csv", "out/s21_leapfrog_trace.csv", "de_leapfrog_dt002",
                                     NULL, sig_lf, NULL, NULL, tr_lf);
    md_register_csv(ctx, "S21", "out/s21_leapfrog_dt.csv");
    md_register_csv(ctx, "S21", "out/s21_leapfrog_trace.csv");
    md_ctx_record_stage(ctx, "S21", "energy_dt_leapfrog", md_wall_time_sec() - t21, 2);

    double t22 = md_wall_time_sec();
    (void)run_dt_sweep(&base, &pot, INT_EULER, dts, 1.6,
                       "out/s22_euler_negctrl.csv", "out/s22_euler_trace.csv", "de_euler_dt002",
                       sig_vv, sig_eu, &drift_eu, NULL, tr_eu);
    md_system_free(&base);
    ctx->c5_euler_ratio = drift_eu / (drift_vv + 1e-30);
    md_register_csv(ctx, "S22", "out/s22_euler_negctrl.csv");
    md_register_csv(ctx, "S22", "out/s22_euler_trace.csv");

    char sub[128];
    snprintf(sub, sizeof(sub), "Log-log slope: VV=%.3f+-%.3f, Leapfrog=%.3f+-%.3f; Euler/VV=%.0fx",
             ovv.slope, ovv.slope_se, olf.slope, olf.slope_se, ctx->c5_euler_ratio);
    md_svg_series s_ord[3] = {
        {"Velocity Verlet (O(dt^2))", "#2563eb", dts, sig_vv, NULL, 5, true, true, false},
        {"Leapfrog DKD (O(dt^2))", "#16a34a", dts, sig_lf, NULL, 5, true, true, true},
        {"Explicit Euler (neg ctrl)", "#dc2626", dts, sig_eu, NULL, 5, true, true, false}
    };
    md_svg_opts o1 = {"S20-S22: Energy Fluctuation sigma_E vs Timestep dt (log-log)", sub,
                      "Timestep dt", "Energy std dev sigma_E", 760, 460, true, true, false, 0.0, NULL};
    const char *svg1 = "../../../outputs/figures/day009a_15_integrator_order.svg";
    md_svg_write_xy(svg1, &o1, s_ord, 3);
    md_register_svg(ctx, "S22", svg1);

    md_svg_series s_tr[3] = {
        {"Velocity Verlet", "#2563eb", trace_t, tr_vv, NULL, 160, true, false, false},
        {"Leapfrog DKD", "#16a34a", trace_t, tr_lf, NULL, 160, true, false, true},
        {"Explicit Euler", "#dc2626", trace_t, tr_eu, NULL, 160, true, false, false}
    };
    md_svg_opts o2 = {"S22: Energy Drift E(t) - E(0) at dt = 0.002",
                      "Symplectic integrators remain bounded while Euler drifts secularly",
                      "Time t", "Delta E(t)", 760, 460, false, false, true, 0.0, "Exact conservation"};
    const char *svg2 = "../../../outputs/figures/day009a_16_euler_vs_vv.svg";
    md_svg_write_xy(svg2, &o2, s_tr, 3);
    md_register_svg(ctx, "S22", svg2);
    md_ctx_record_stage(ctx, "S22", "energy_dt_euler", md_wall_time_sec() - t22, 4);
}

static void stage_s23_momentum_nve(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int n_steps = ctx->small_scale ? 10000 : 100000;
    double dt = 0.002;
    md_system sys;
    prep_liquid_state(&sys, 3, 0.7, 1.2, 23);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);

    double e0 = md_kinetic_energy(&sys) + sys.pe;
    int stride = n_steps / 200;
    double tx[200], py[200];
    int ti = 0;
    FILE *ft = fopen("out/s23_nve_long_trace.csv", "w");
    fprintf(ft, "step,t,p_norm,e_tot,de_rel,temp\n");
    double max_p = 0.0;
    for (int st = 1; st <= n_steps; ++st) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, NULL, NULL);
        if (st % stride == 0 && ti < 200) {
            double pn = v3_norm(md_total_momentum(&sys));
            if (pn > max_p) max_p = pn;
            double e = md_kinetic_energy(&sys) + sys.pe;
            tx[ti] = (double)st * dt;
            py[ti] = (pn > 1e-17) ? pn : 1e-17;
            fprintf(ft, "%d,%.4f,%.6e,%.8f,%.6e,%.6f\n",
                    st, tx[ti], pn, e, fabs(e - e0) / fabs(e0), md_temperature(&sys));
            ti++;
        }
    }
    fclose(ft);
    double final_p = v3_norm(md_total_momentum(&sys));
    if (final_p > max_p) max_p = final_p;
    ctx->c3_p_tot_1e5 = max_p;
    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;

    FILE *fs = fopen("out/s23_momentum_nve.csv", "w");
    fprintf(fs, "n_particles,n_steps,dt,final_p_norm,max_p_norm,e0,e_final,rebuilds\n");
    fprintf(fs, "%d,%d,%.4f,%.6e,%.6e,%.8f,%.8f,%llu\n",
            sys.n, n_steps, dt, final_p, max_p, e0, md_kinetic_energy(&sys) + sys.pe,
            (unsigned long long)vl.rebuild_count);
    fclose(fs);

    md_register_csv(ctx, "S23", "out/s23_momentum_nve.csv");
    md_register_csv(ctx, "S23", "out/s23_nve_long_trace.csv");
    md_svg_series ser[1] = {{"||P_tot(t)||", "#2563eb", tx, py, NULL, ti, true, false, false}};
    md_svg_opts o = {"S23: Total Momentum Conservation Over 100,000 NVE Steps",
                     "N=108, dt=0.002, Velocity Verlet + Newton's 3rd law (c3 limit 1e-10)",
                     "Time t", "||P_tot(t)||", 760, 460, false, true, true, 1e-10, "c3 limit (1e-10)"};
    const char *svgp = "../../../outputs/figures/day009a_17_momentum_drift.svg";
    md_svg_write_xy(svgp, &o, ser, 1);
    md_register_svg(ctx, "S23", svgp);
    md_verlet_free(&vl);
    md_system_free(&sys);
    md_ctx_record_stage(ctx, "S23", "momentum_nve", md_wall_time_sec() - t0, 3);
}

static void stage_s24_long_drift(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int n_steps = ctx->small_scale ? 4000 : 20000;
    double dts[2] = {0.002, 0.005};
    double t_pts[200], de1[200], de2[200];
    md_system base;
    prep_liquid_state(&base, 4, 0.75, 1.1, 24);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fs = fopen("out/s24_long_drift_summary.csv", "w");
    fprintf(fs, "dt,n_steps,e0,mean_de,sigma_de,slope_de_per_t\n");
    for (int k = 0; k < 2; ++k) {
        md_system s;
        memset(&s, 0, sizeof(s));
        md_system_copy(&s, &base);
        md_verlet_list vl;
        md_verlet_init(&vl, s.n, pot.rc, 0.3, false);
        md_forces_verlet(&s, &pot, &vl);
        double e0 = md_kinetic_energy(&s) + s.pe;
        int stride = n_steps / 200, ti = 0;
        double arr_t[200], arr_de[200];
        for (int st = 1; st <= n_steps; ++st) {
            md_step_vverlet(&s, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dts[k], NULL, NULL);
            if (st % stride == 0 && ti < 200) {
                arr_t[ti] = (double)st * dts[k];
                arr_de[ti] = (md_kinetic_energy(&s) + s.pe) - e0;
                if (k == 0) { t_pts[ti] = arr_t[ti]; de1[ti] = arr_de[ti]; }
                else de2[ti] = arr_de[ti];
                ti++;
            }
        }
        md_ols_result lr = md_ols_fit(arr_t, arr_de, ti);
        fprintf(fs, "%.4f,%d,%.8f,%.6e,%.6e,%.6e\n",
                dts[k], n_steps, e0, md_mean(arr_de, ti), sqrt(md_variance(arr_de, ti)), lr.slope);
        ctx->total_nbr_overflows += vl.overflow_count;
        ctx->total_late_rebuilds += vl.late_rebuild_count;
        md_verlet_free(&vl);
        md_system_free(&s);
    }
    fclose(fs);
    md_system_free(&base);

    FILE *ft = fopen("out/s24_long_drift.csv", "w");
    fprintf(ft, "idx,t_dt002,de_dt002,de_dt005\n");
    for (int i = 0; i < 200; ++i) fprintf(ft, "%d,%.4f,%.6e,%.6e\n", i, t_pts[i], de1[i], de2[i]);
    fclose(ft);
    md_register_csv(ctx, "S24", "out/s24_long_drift.csv");
    md_register_csv(ctx, "S24", "out/s24_long_drift_summary.csv");

    md_svg_series ser[2] = {
        {"dt = 0.002", "#2563eb", t_pts, de1, NULL, 200, true, false, false},
        {"dt = 0.005", "#dc2626", t_pts, de2, NULL, 200, true, false, false}
    };
    md_svg_opts o = {"S24: Long-Time NVE Energy Stability (N=256, 20,000 Steps)",
                     "Bounded symplectic shadow-Hamiltonian oscillations without blow-up",
                     "Time t (dt=0.002 scale)", "E(t) - E(0)", 760, 460, false, false, true, 0.0, "Delta E = 0"};
    const char *svgp = "../../../outputs/figures/day009a_18_long_energy.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S24", svgp);
    md_ctx_record_stage(ctx, "S24", "long_energy_drift", md_wall_time_sec() - t0, 3);
}

static void stage_s25_lyapunov(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double states[2][2] = {{0.6, 1.5}, {0.8, 1.0}};
    double tx[125], ly1[125], ly2[125];
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *ff = fopen("out/s25_lyapunov_fit.csv", "w");
    fprintf(ff, "rho,temp,lambda_rough,lambda_se,r2_linear_window\n");
    FILE *ft = fopen("out/s25_lyapunov_trace.csv", "w");
    fprintf(ft, "t,dist_rho06_t15,dist_rho08_t10\n");

    for (int si = 0; si < 2; ++si) {
        md_system s1, s2;
        prep_liquid_state(&s1, 3, states[si][0], states[si][1], (uint64_t)(25 + si));
        memset(&s2, 0, sizeof(s2));
        md_system_copy(&s2, &s1);
        s2.pos[0].x += 1e-9; s2.unwrapped[0].x += 1e-9;
        md_forces_brute(&s1, &pot); md_forces_brute(&s2, &pot);

        double fit_t[60], fit_y[60];
        int nf = 0, ti = 0;
        for (int st = 1; st <= 2500; ++st) {
            md_step_vverlet(&s1, &pot, NBR_BRUTE, NULL, NULL, NULL, 1, 0.002, NULL, NULL);
            md_step_vverlet(&s2, &pot, NBR_BRUTE, NULL, NULL, NULL, 1, 0.002, NULL, NULL);
            if (st % 20 == 0 && ti < 125) {
                double sum2 = 0.0;
                for (int i = 0; i < s1.n; ++i) sum2 += v3_norm2(v3_sub(s2.unwrapped[i], s1.unwrapped[i]));
                double dist = sqrt(sum2);
                tx[ti] = (double)st * 0.002;
                if (si == 0) ly1[ti] = dist; else ly2[ti] = dist;
                if (dist > 1e-8 && dist < 1e-2 && nf < 60) {
                    fit_t[nf] = tx[ti]; fit_y[nf] = log(dist); nf++;
                }
                ti++;
            }
        }
        md_ols_result lr = md_ols_fit(fit_t, fit_y, nf);
        fprintf(ff, "%.2f,%.2f,%.6f,%.6f,%.6f\n", states[si][0], states[si][1], lr.slope, lr.slope_se, lr.r2);
        md_system_free(&s1); md_system_free(&s2);
    }
    for (int i = 0; i < 125; ++i) fprintf(ft, "%.4f,%.6e,%.6e\n", tx[i], ly1[i], ly2[i]);
    fclose(ff); fclose(ft);

    md_register_csv(ctx, "S25", "out/s25_lyapunov_trace.csv");
    md_register_csv(ctx, "S25", "out/s25_lyapunov_fit.csv");
    md_svg_series ser[2] = {
        {"rho=0.6, T=1.5", "#2563eb", tx, ly1, NULL, 125, true, false, false},
        {"rho=0.8, T=1.0", "#dc2626", tx, ly2, NULL, 125, true, false, false}
    };
    md_svg_opts o = {"S25: Twin-Trajectory Divergence (Rough Lyapunov Report)",
                     "Initial perturbation delta_r0 = 1e-9 (N=108 NVE; diagnostic report only)",
                     "Time t", "Euclidean separation ||Delta R(t)||", 760, 460, false, true, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_19_lyapunov.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S25", svgp);
    md_ctx_record_stage(ctx, "S25", "lyapunov_rough", md_wall_time_sec() - t0, 3);
}

static void stage_s26_respa(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_system base;
    prep_liquid_state(&base, 3, 0.7, 1.2, 26);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    int ms[3] = {1, 2, 4};
    double tx[100], de_m1[100], de_m4[100];

    FILE *fp = fopen("out/s26_respa.csv", "w");
    fprintf(fp, "m_inner,dt_outer,dt_inner,sigma_e,max_abs_de\n");
    FILE *ft = fopen("out/s26_respa_trace.csv", "w");
    fprintf(ft, "t,de_respa_m1,de_respa_m4\n");

    for (int k = 0; k < 3; ++k) {
        md_system s;
        memset(&s, 0, sizeof(s));
        md_system_copy(&s, &base);
        md_verlet_list vl;
        md_verlet_init(&vl, s.n, pot.rc, 0.3, false);
        md_forces_verlet(&s, &pot, &vl);
        double e0 = md_kinetic_energy(&s) + s.pe;
        double es[400], max_de = 0.0;
        for (int st = 0; st < 400; ++st) {
            md_step_respa(&s, &pot, &vl, 0.006, ms[k]);
            double e = md_kinetic_energy(&s) + s.pe;
            es[st] = e;
            double de = fabs(e - e0);
            if (de > max_de) max_de = de;
            if (st % 4 == 0) {
                int idx = st / 4;
                tx[idx] = (double)(st + 1) * 0.006;
                if (k == 0) de_m1[idx] = e - e0;
                if (k == 2) de_m4[idx] = e - e0;
            }
        }
        fprintf(fp, "%d,0.0060,%.5f,%.6e,%.6e\n", ms[k], 0.006 / (double)ms[k], sqrt(md_variance(es, 400)), max_de);
        md_verlet_free(&vl);
        md_system_free(&s);
    }
    for (int i = 0; i < 100; ++i) fprintf(ft, "%.4f,%.6e,%.6e\n", tx[i], de_m1[i], de_m4[i]);
    fclose(fp); fclose(ft);
    md_system_free(&base);

    md_register_csv(ctx, "S26", "out/s26_respa.csv");
    md_register_csv(ctx, "S26", "out/s26_respa_trace.csv");
    md_svg_series ser[2] = {
        {"r-RESPA m=1 (dt_in=0.006)", "#dc2626", tx, de_m1, NULL, 100, true, false, false},
        {"r-RESPA m=4 (dt_in=0.0015)", "#2563eb", tx, de_m4, NULL, 100, true, false, false}
    };
    md_svg_opts o = {"S26: Two-Timescale r-RESPA Energy Conservation (dt_outer=0.006)",
                     "Sub-stepping stiff WCA core (m=4) sharply reduces energy fluctuations",
                     "Time t", "Delta E(t)", 760, 460, false, false, true, 0.0, "Delta E = 0"};
    const char *svgp = "../../../outputs/figures/day009a_20_respa.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S26", svgp);
    md_ctx_record_stage(ctx, "S26", "respa_twoscale", md_wall_time_sec() - t0, 3);
}

void md_run_stages_c(md_run_ctx *ctx) {
    stage_s19_reversibility(ctx);
    stage_s20_s21_s22_integrator_order(ctx);
    stage_s23_momentum_nve(ctx);
    stage_s24_long_drift(ctx);
    stage_s25_lyapunov(ctx);
    stage_s26_respa(ctx);
}
