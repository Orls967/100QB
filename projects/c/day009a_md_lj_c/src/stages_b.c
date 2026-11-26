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

/* Group B: Neighbor Structures, Scaling & Optimization Ladder (S11..S18) */

static double max_rel_force_diff(const md_vec3 *a, const md_vec3 *b, int n) {
    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        double rel = v3_norm(v3_sub(a[i], b[i])) / (v3_norm(a[i]) + 1.0);
        if (rel > max_err) max_err = rel;
    }
    return max_err;
}

static void init_jittered_fcc(md_system *sys, int nc, double rho, double amp, uint64_t sid, uint64_t sub) {
    md_init_fcc(sys, nc, rho);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, sid, sub);
    for (int i = 0; i < sys->n; ++i) {
        sys->pos[i].x += md_rng_range(&rng, -amp, amp);
        sys->pos[i].y += md_rng_range(&rng, -amp, amp);
        sys->pos[i].z += md_rng_range(&rng, -amp, amp);
        sys->pos[i] = md_box_wrap(&sys->box, sys->pos[i]);
        sys->unwrapped[i] = sys->pos[i];
    }
}

static void stage_s11_nbr_equiv(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double rhos[3] = {0.2, 0.5, 0.85};
    int ncs[2] = {5, 6}; /* N=500, N=864 */
    double worst_err = 0.0;
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fp = fopen("out/s11_nbr_equiv.csv", "w");
    fprintf(fp, "n_particles,rho,err_cell,err_verlet,err_soa,err_threaded4,err_cellsorted,pe_brute,pe_verlet_diff\n");

    for (int ni = 0; ni < 2; ++ni) {
        for (int ri = 0; ri < 3; ++ri) {
            md_system sys;
            init_jittered_fcc(&sys, ncs[ni], rhos[ri], 0.10, 11, (uint64_t)(ni * 3 + ri));
            md_vec3 *f_ref = (md_vec3 *)malloc((size_t)sys.n * sizeof(md_vec3));
            md_forces_brute(&sys, &pot);
            double pe_ref = sys.pe;
            memcpy(f_ref, sys.force, (size_t)sys.n * sizeof(md_vec3));

            md_cell_list cl; md_cell_init(&cl, sys.n);
            md_forces_cell(&sys, &pot, &cl);
            double ec = max_rel_force_diff(f_ref, sys.force, sys.n);

            md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, 0.3, true);
            md_forces_verlet(&sys, &pot, &vl);
            double ev = max_rel_force_diff(f_ref, sys.force, sys.n);
            double pe_vd = fabs(sys.pe - pe_ref);

            md_soa soa; md_soa_alloc(&soa, sys.n);
            md_forces_soa_verlet(&sys, &soa, &pot, &vl);
            double es = max_rel_force_diff(f_ref, sys.force, sys.n);

            md_forces_threaded(&sys, &pot, &vl, 4);
            double et = max_rel_force_diff(f_ref, sys.force, sys.n);

            md_forces_cell_sorted_soa(&sys, &soa, &pot, &vl, 4);
            double ecs = max_rel_force_diff(f_ref, sys.force, sys.n);

            double errs[5] = {ec, ev, es, et, ecs};
            for (int k = 0; k < 5; ++k) if (errs[k] > worst_err) worst_err = errs[k];
            ctx->total_nbr_overflows += vl.overflow_count;
            ctx->total_late_rebuilds += vl.late_rebuild_count;

            fprintf(fp, "%d,%.2f,%.6e,%.6e,%.6e,%.6e,%.6e,%.8f,%.6e\n",
                    sys.n, rhos[ri], ec, ev, es, et, ecs, pe_ref, pe_vd);
            free(f_ref); md_soa_free(&soa); md_verlet_free(&vl); md_cell_free(&cl); md_system_free(&sys);
        }
    }
    fclose(fp);
    ctx->c1_max_rel_err = worst_err;

    FILE *f2 = fopen("out/s11_nbr_per_pot.csv", "w");
    fprintf(f2, "potential,n_particles,err_cell,err_verlet,pe_brute\n");
    md_system sp;
    init_jittered_fcc(&sp, 5, 0.7, 0.09, 11, 99);
    md_vec3 *fr = (md_vec3 *)malloc((size_t)sp.n * sizeof(md_vec3));
    for (int k = 0; k < POT_COUNT; ++k) {
        md_potential pk; md_potential_init(&pk, (md_pot_kind)k, 2.5);
        md_forces_brute(&sp, &pk);
        double peb = sp.pe;
        memcpy(fr, sp.force, (size_t)sp.n * sizeof(md_vec3));
        md_cell_list cl; md_cell_init(&cl, sp.n);
        md_forces_cell(&sp, &pk, &cl);
        double ec = max_rel_force_diff(fr, sp.force, sp.n);
        md_verlet_list vl; md_verlet_init(&vl, sp.n, pk.rc, 0.3, false);
        md_forces_verlet(&sp, &pk, &vl);
        double ev = max_rel_force_diff(fr, sp.force, sp.n);
        fprintf(f2, "%s,%d,%.6e,%.6e,%.8f\n", md_pot_name((md_pot_kind)k), sp.n, ec, ev, peb);
        md_verlet_free(&vl); md_cell_free(&cl); md_potential_free(&pk);
    }
    free(fr); md_system_free(&sp); fclose(f2);

    md_register_csv(ctx, "S11", "out/s11_nbr_equiv.csv");
    md_register_csv(ctx, "S11", "out/s11_nbr_per_pot.csv");
    md_ctx_record_stage(ctx, "S11", "nbr_equivalence", md_wall_time_sec() - t0, 2);
}

static void stage_s12_verlet_trigger(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int n_steps = ctx->small_scale ? 600 : 2000;
    double dt = 0.002, skin = 0.3;
    md_system sys; md_init_fcc(&sys, 4, 0.7);
    md_rng rng; md_rng_seed_stage(&rng, MD_MASTER_SEED, 12, 0);
    md_init_velocities(&sys, 1.5, true, &rng);
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, skin, false);
    md_forces_verlet(&sys, &pot, &vl);
    md_thermostat th; md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 3);

    int stride = n_steps / 200;
    if (stride < 1) stride = 1;
    double tx[202], td[202];
    int nt = 0;
    double max_f_audit = 0.0;
    FILE *ft = fopen("out/s12_verlet_disp_trace.csv", "w");
    fprintf(ft, "step,t,max_disp,r_skin_half,rebuilds_cum\n");

    for (int step = 1; step <= n_steps; ++step) {
        md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, dt, &th, &rng);
        if (step % stride == 0 && nt < 200) {
            tx[nt] = (double)step * dt;
            td[nt] = vl.last_max_disp;
            fprintf(ft, "%d,%.4f,%.8f,%.4f,%llu\n",
                    step, tx[nt], td[nt], 0.5 * skin, (unsigned long long)vl.rebuild_count);
            nt++;
        }
        if (step % 250 == 0) {
            md_system copy; memset(&copy, 0, sizeof(copy));
            md_system_copy(&copy, &sys);
            md_forces_brute(&copy, &pot);
            double e = max_rel_force_diff(copy.force, sys.force, sys.n);
            if (e > max_f_audit) max_f_audit = e;
            md_system_free(&copy);
        }
    }
    fclose(ft);

    FILE *fs = fopen("out/s12_verlet_trigger.csv", "w");
    fprintf(fs, "n_steps,r_skin,r_skin_half,rebuild_count,avg_steps_per_rebuild,max_disp_before_rebuild,overflow_count,late_rebuild_count,audit_max_ferr\n");
    fprintf(fs, "%d,%.3f,%.3f,%llu,%.4f,%.8f,%llu,%llu,%.6e\n",
            n_steps, skin, 0.5 * skin, (unsigned long long)vl.rebuild_count,
            (double)n_steps / (double)(vl.rebuild_count > 0 ? vl.rebuild_count : 1),
            vl.max_disp_before_rebuild, (unsigned long long)vl.overflow_count,
            (unsigned long long)vl.late_rebuild_count, max_f_audit);
    fclose(fs);
    ctx->total_nbr_overflows += vl.overflow_count;
    ctx->total_late_rebuilds += vl.late_rebuild_count;

    md_register_csv(ctx, "S12", "out/s12_verlet_trigger.csv");
    md_register_csv(ctx, "S12", "out/s12_verlet_disp_trace.csv");
    md_svg_series ser[1] = {{"max_i ||Delta r_i||", "#2563eb", tx, td, NULL, nt, true, false, false}};
    md_svg_opts o = {"S12: Verlet List Displacement Trigger Audit",
                     "Sawtooth max displacement resets cleanly below r_skin/2 = 0.15 (0 late rebuilds)",
                     "Time t", "Max displacement since rebuild", 760, 460, false, false, true, 0.5 * skin, "Safety limit r_skin/2"};
    const char *svgp = "../../../outputs/figures/day009a_09_verlet_disp.svg";
    md_svg_write_xy(svgp, &o, ser, 1);
    md_register_svg(ctx, "S12", svgp);
    md_verlet_free(&vl); md_system_free(&sys);
    md_ctx_record_stage(ctx, "S12", "verlet_trigger", md_wall_time_sec() - t0, 3);
}

static void stage_s13_skin_sweep(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double skins[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
    double reb_rate[5], us_step[5];
    int n_steps = ctx->small_scale ? 250 : 800;
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fp = fopen("out/s13_skin_sweep.csv", "w");
    fprintf(fp, "r_skin,rebuilds,avg_nbrs_per_atom,us_per_step,overflow_count,late_rebuild_count\n");
    FILE *fr = fopen("out/s13_skin_rebuilds.csv", "w");
    fprintf(fr, "r_skin,steps_per_rebuild,rebuild_freq_per_100\n");

    for (int k = 0; k < 5; ++k) {
        md_system sys; md_init_fcc(&sys, 5, 0.75);
        md_rng rng; md_rng_seed_stage(&rng, MD_MASTER_SEED, 13, 0);
        md_init_velocities(&sys, 1.5, true, &rng);
        md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, skins[k], false);
        md_forces_verlet(&sys, &pot, &vl);
        md_thermostat th; md_thermostat_init(&th, THERM_LANGEVIN, 1.5, sys.n, 3);

        double ts = md_wall_time_sec();
        for (int s = 0; s < n_steps; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, &th, &rng);
        }
        us_step[k] = (md_wall_time_sec() - ts) * 1e6 / (double)n_steps;
        reb_rate[k] = 100.0 * (double)vl.rebuild_count / (double)n_steps;
        ctx->total_nbr_overflows += vl.overflow_count;
        ctx->total_late_rebuilds += vl.late_rebuild_count;

        fprintf(fp, "%.2f,%llu,%.4f,%.4f,%llu,%llu\n", skins[k], (unsigned long long)vl.rebuild_count,
                2.0 * (double)vl.half_count / (double)sys.n, us_step[k],
                (unsigned long long)vl.overflow_count, (unsigned long long)vl.late_rebuild_count);
        fprintf(fr, "%.2f,%.4f,%.4f\n", skins[k],
                (double)n_steps / (double)(vl.rebuild_count > 0 ? vl.rebuild_count : 1), reb_rate[k]);
        md_verlet_free(&vl); md_system_free(&sys);
    }
    fclose(fp); fclose(fr);
    md_register_csv(ctx, "S13", "out/s13_skin_sweep.csv");
    md_register_csv(ctx, "S13", "out/s13_skin_rebuilds.csv");

    md_svg_series ser[2] = {
        {"Time per step (us)", "#2563eb", skins, us_step, NULL, 5, true, true, false},
        {"Rebuilds per 100 steps", "#dc2626", skins, reb_rate, NULL, 5, true, true, true}
    };
    md_svg_opts o = {"S13: Verlet Skin Trade-off (N=500, rho=0.75, T=1.5)",
                     "Larger skin cuts rebuild frequency but enlarges neighbor list",
                     "Verlet skin r_skin", "Time (us/step) / Rebuilds per 100 steps",
                     760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_10_skin_sweep.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S13", svgp);
    md_ctx_record_stage(ctx, "S13", "verlet_skin_sweep", md_wall_time_sec() - t0, 3);
}

static void stage_s14_s15_scaling(md_run_ctx *ctx) {
    double t14 = md_wall_time_sec();
    int ncs_b[5] = {5, 6, 7, 8, 9}; /* N = 500..2916 */
    int n_rep = ctx->small_scale ? 3 : 5;
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    double nb_x[5], nb_y[5], log_nb[5], log_tb[5];

    FILE *fb = fopen("out/s14_scaling_brute.csv", "w");
    fprintf(fb, "nc,n_particles,median_us,mad_us\n");
    FILE *fbr = fopen("out/s14_scaling_brute_repeats.csv", "w");
    fprintf(fbr, "n_particles,repeat,time_us\n");

    for (int k = 0; k < 5; ++k) {
        md_system sys; init_jittered_fcc(&sys, ncs_b[k], 0.8, 0.05, 14, (uint64_t)k);
        int calls = (sys.n <= 864) ? 6 : 2;
        double vals[8], med = 0.0, mad = 0.0;
        for (int r = 0; r < n_rep; ++r) {
            double ts = md_wall_time_sec();
            for (int c = 0; c < calls; ++c) md_forces_brute(&sys, &pot);
            vals[r] = (md_wall_time_sec() - ts) * 1e6 / (double)calls;
            fprintf(fbr, "%d,%d,%.4f\n", sys.n, r, vals[r]);
        }
        md_median_mad(vals, n_rep, &med, &mad);
        fprintf(fb, "%d,%d,%.4f,%.4f\n", ncs_b[k], sys.n, med, mad);
        nb_x[k] = (double)sys.n; nb_y[k] = med;
        log_nb[k] = log((double)sys.n); log_tb[k] = log(med);
        md_system_free(&sys);
    }
    fclose(fb); fclose(fbr);
    md_ols_result ob = md_ols_fit(log_nb, log_tb, 5);
    ctx->c17_brute_slope = ob.slope; ctx->c17_brute_se = ob.slope_se;
    md_register_csv(ctx, "S14", "out/s14_scaling_brute.csv");
    md_register_csv(ctx, "S14", "out/s14_scaling_brute_repeats.csv");
    md_ctx_record_stage(ctx, "S14", "scaling_brute", md_wall_time_sec() - t14, 2);

    double t15 = md_wall_time_sec();
    int ncs_f[6] = {5, 6, 7, 8, 9, 10}; /* N = 500..4000 */
    double nf_x[6], cell_y[6], verl_y[6], log_nf[6], log_tc[6], log_tv[6];
    FILE *ff = fopen("out/s15_scaling_fast.csv", "w");
    fprintf(ff, "nc,n_particles,cell_median_us,cell_mad_us,verlet_median_us,verlet_mad_us\n");
    FILE *ffr = fopen("out/s15_scaling_fast_repeats.csv", "w");
    fprintf(ffr, "n_particles,repeat,cell_us,verlet_us\n");

    for (int k = 0; k < 6; ++k) {
        md_system sys; init_jittered_fcc(&sys, ncs_f[k], 0.8, 0.05, 15, (uint64_t)k);
        md_cell_list cl; md_cell_init(&cl, sys.n);
        md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
        md_verlet_build(&vl, &sys);
        int calls = (sys.n <= 1372) ? 12 : 6;
        double vc[8], vv[8], mc = 0.0, madc = 0.0, mv = 0.0, madv = 0.0;
        for (int r = 0; r < n_rep; ++r) {
            double tc = md_wall_time_sec();
            for (int c = 0; c < calls; ++c) md_forces_cell(&sys, &pot, &cl);
            vc[r] = (md_wall_time_sec() - tc) * 1e6 / (double)calls;
            double tv = md_wall_time_sec();
            for (int c = 0; c < calls; ++c) md_forces_verlet(&sys, &pot, &vl);
            vv[r] = (md_wall_time_sec() - tv) * 1e6 / (double)calls;
            fprintf(ffr, "%d,%d,%.4f,%.4f\n", sys.n, r, vc[r], vv[r]);
        }
        md_median_mad(vc, n_rep, &mc, &madc);
        md_median_mad(vv, n_rep, &mv, &madv);
        fprintf(ff, "%d,%d,%.4f,%.4f,%.4f,%.4f\n", ncs_f[k], sys.n, mc, madc, mv, madv);
        nf_x[k] = (double)sys.n; cell_y[k] = mc; verl_y[k] = mv;
        log_nf[k] = log((double)sys.n); log_tc[k] = log(mc); log_tv[k] = log(mv);
        md_verlet_free(&vl); md_cell_free(&cl); md_system_free(&sys);
    }
    fclose(ff); fclose(ffr);
    md_ols_result oc = md_ols_fit(log_nf, log_tc, 6);
    md_ols_result ov = md_ols_fit(log_nf, log_tv, 6);
    ctx->c17_cell_slope = oc.slope; ctx->c17_cell_se = oc.slope_se;
    ctx->c17_verlet_slope = ov.slope; ctx->c17_verlet_se = ov.slope_se;

    md_register_csv(ctx, "S15", "out/s15_scaling_fast.csv");
    md_register_csv(ctx, "S15", "out/s15_scaling_fast_repeats.csv");
    char sub[128];
    snprintf(sub, sizeof(sub), "Slopes (N>=500): Brute=%.2f+-%.2f, Cell=%.2f+-%.2f, Verlet=%.2f+-%.2f",
             ob.slope, ob.slope_se, oc.slope, oc.slope_se, ov.slope, ov.slope_se);
    md_svg_series ser[3] = {
        {"Brute force O(N^2)", "#dc2626", nb_x, nb_y, NULL, 5, true, true, false},
        {"Cell list O(N)", "#d97706", nf_x, cell_y, NULL, 6, true, true, false},
        {"Verlet list O(N)", "#2563eb", nf_x, verl_y, NULL, 6, true, true, false}
    };
    md_svg_opts o = {"S14-S15: Force Evaluation Time Scaling vs N (log-log)", sub,
                     "Number of particles N", "Time per force eval (us)",
                     760, 460, true, true, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_11_scaling_loglog.svg";
    md_svg_write_xy(svgp, &o, ser, 3);
    md_register_svg(ctx, "S15", svgp);
    md_ctx_record_stage(ctx, "S15", "scaling_cell_verlet", md_wall_time_sec() - t15, 3);
}

static void stage_s16_opt_ladder(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int nc = ctx->small_scale ? 6 : 7;
    int n_rep = ctx->small_scale ? 5 : 7;
    int steps = ctx->small_scale ? 8 : 16;
    md_system base; init_jittered_fcc(&base, nc, 0.8, 0.06, 16, 0);
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    const char *vnames[6] = {"V0_brute", "V1_cell", "V2_verlet", "V3_soa", "V4_thr8", "V5_sort_soa8"};
    const char *vcols[6] = {"#dc2626", "#ea580c", "#d97706", "#0891b2", "#2563eb", "#16a34a"};
    md_nbr_kind kinds[6] = {NBR_BRUTE, NBR_CELL, NBR_VERLET, NBR_SOA_VERLET, NBR_THREADED, NBR_CELL_SORTED_SOA};
    int threads_v[6] = {1, 1, 1, 1, 8, 8};

    FILE *fr = fopen("out/s16_ladder_repeats.csv", "w");
    fprintf(fr, "version,repeat,n_particles,us_per_eval\n");
    FILE *fs = fopen("out/s16_ladder_summary.csv", "w");
    fprintf(fs, "version,n_particles,median_us,mad_us,speedup_vs_v0\n");

    double meds[6], mads[6];
    for (int v = 0; v < 6; ++v) {
        md_system sys; memset(&sys, 0, sizeof(sys)); md_system_copy(&sys, &base);
        md_cell_list cl; md_cell_init(&cl, sys.n);
        md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, 0.3, v >= 4);
        md_verlet_build(&vl, &sys);
        md_soa soa; md_soa_alloc(&soa, sys.n);
        if (v == 5) md_soa_sort_by_cell(&soa, &sys, vl.r_list);

        int iters = (v == 0) ? 3 : steps;
        double times[10];
        md_forces_compute(&sys, &pot, kinds[v], &cl, &vl, &soa, threads_v[v]);
        for (int r = 0; r < n_rep; ++r) {
            double ts = md_wall_time_sec();
            for (int it = 0; it < iters; ++it) md_forces_compute(&sys, &pot, kinds[v], &cl, &vl, &soa, threads_v[v]);
            times[r] = (md_wall_time_sec() - ts) * 1e6 / (double)iters;
            fprintf(fr, "%s,%d,%d,%.4f\n", vnames[v], r, sys.n, times[r]);
        }
        md_median_mad(times, n_rep, &meds[v], &mads[v]);
        ctx->c18_medians_us[v] = meds[v]; ctx->c18_mads_us[v] = mads[v];
        md_soa_free(&soa); md_verlet_free(&vl); md_cell_free(&cl); md_system_free(&sys);
    }

    bool mono = true;
    for (int v = 0; v < 6; ++v) {
        if (v > 0 && meds[v] >= meds[v - 1]) mono = false;
        fprintf(fs, "%s,%d,%.4f,%.4f,%.4f\n", vnames[v], base.n, meds[v], mads[v], meds[0] / meds[v]);
    }
    ctx->c18_monotone_faster = mono;
    fclose(fr); fclose(fs); md_system_free(&base);
    md_register_csv(ctx, "S16", "out/s16_ladder_repeats.csv");
    md_register_csv(ctx, "S16", "out/s16_ladder_summary.csv");

    const char *svgp = "../../../outputs/figures/day009a_12_opt_ladder.svg";
    md_svg_write_bars(svgp, "S16: Optimization Ladder V0..V5 (Median Time per Force Eval)",
                      "N=1372, rho=0.8, 7 repeats (log y-axis; lower is faster)",
                      "Median time per force eval (us)", vnames, meds, NULL, vcols, 6, true);
    md_register_svg(ctx, "S16", svgp);
    md_ctx_record_stage(ctx, "S16", "opt_ladder_raw", md_wall_time_sec() - t0, 3);
}

static void stage_s17_thread_scaling(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int thr_list[4] = {1, 2, 4, 8}, ncs[2] = {6, 8};
    double tx[4] = {1.0, 2.0, 4.0, 8.0}, sp_864[4], sp_2048[4];
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    int n_rep = ctx->small_scale ? 3 : 5;

    FILE *fs = fopen("out/s17_thread_scaling.csv", "w");
    fprintf(fs, "n_particles,threads,median_us,mad_us,speedup,efficiency,force_hash\n");
    FILE *fr = fopen("out/s17_thread_repeats.csv", "w");
    fprintf(fr, "n_particles,threads,repeat,time_us\n");

    for (int ni = 0; ni < 2; ++ni) {
        md_system sys; init_jittered_fcc(&sys, ncs[ni], 0.8, 0.05, 17, (uint64_t)ni);
        md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, 0.3, true);
        md_verlet_build(&vl, &sys);
        double base_t = 1.0;
        for (int ti = 0; ti < 4; ++ti) {
            int nt = thr_list[ti], iters = (sys.n <= 864) ? 20 : 12;
            double times[8], med = 0.0, mad = 0.0;
            md_forces_threaded(&sys, &pot, &vl, nt);
            uint64_t h = md_fnv1a_system_state(&sys);
            for (int r = 0; r < n_rep; ++r) {
                double ts = md_wall_time_sec();
                for (int it = 0; it < iters; ++it) md_forces_threaded(&sys, &pot, &vl, nt);
                times[r] = (md_wall_time_sec() - ts) * 1e6 / (double)iters;
                fprintf(fr, "%d,%d,%d,%.4f\n", sys.n, nt, r, times[r]);
            }
            md_median_mad(times, n_rep, &med, &mad);
            if (ti == 0) base_t = med;
            double sp = base_t / med;
            if (ni == 0) sp_864[ti] = sp; else sp_2048[ti] = sp;
            fprintf(fs, "%d,%d,%.4f,%.4f,%.4f,%.4f,%016llx\n",
                    sys.n, nt, med, mad, sp, sp / (double)nt, (unsigned long long)h);
        }
        md_verlet_free(&vl); md_system_free(&sys);
    }
    fclose(fs); fclose(fr);
    md_register_csv(ctx, "S17", "out/s17_thread_scaling.csv");
    md_register_csv(ctx, "S17", "out/s17_thread_repeats.csv");

    md_svg_series ser[3] = {
        {"N=864", "#2563eb", tx, sp_864, NULL, 4, true, true, false},
        {"N=2048", "#16a34a", tx, sp_2048, NULL, 4, true, true, false},
        {"Ideal linear", "#6b7280", tx, tx, NULL, 4, true, false, true}
    };
    md_svg_opts o = {"S17: Deterministic pthread Strong Scaling (1/2/4/8 Threads)",
                     "Bitwise-identical force hashes across thread counts",
                     "Threads", "Speedup vs 1 thread", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_13_thread_scaling.svg";
    md_svg_write_xy(svgp, &o, ser, 3);
    md_register_svg(ctx, "S17", svgp);
    md_ctx_record_stage(ctx, "S17", "thread_scaling", md_wall_time_sec() - t0, 3);
}

static void stage_s18_neighbor_2d(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_system sys; md_init_triangular_2d(&sys, 32, 32, 0.7);
    md_rng rng; md_rng_seed_stage(&rng, MD_MASTER_SEED, 18, 0);
    for (int i = 0; i < sys.n; ++i) {
        sys.pos[i].x += md_rng_range(&rng, -0.07, 0.07);
        sys.pos[i].y += md_rng_range(&rng, -0.07, 0.07);
        sys.pos[i] = md_box_wrap(&sys.box, sys.pos[i]);
        sys.unwrapped[i] = sys.pos[i];
    }
    md_potential pot; md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_vec3 *f_ref = (md_vec3 *)malloc((size_t)sys.n * sizeof(md_vec3));
    double tb = md_wall_time_sec();
    md_forces_brute(&sys, &pot);
    double dt_b = (md_wall_time_sec() - tb) * 1e6, pe_b = sys.pe;
    memcpy(f_ref, sys.force, (size_t)sys.n * sizeof(md_vec3));

    md_cell_list cl; md_cell_init(&cl, sys.n);
    double tc = md_wall_time_sec();
    for (int r = 0; r < 10; ++r) md_forces_cell(&sys, &pot, &cl);
    double dt_c = (md_wall_time_sec() - tc) * 1e6 / 10.0;
    double ec = max_rel_force_diff(f_ref, sys.force, sys.n);

    md_verlet_list vl; md_verlet_init(&vl, sys.n, pot.rc, 0.3, true);
    md_verlet_build(&vl, &sys);
    double tv = md_wall_time_sec();
    for (int r = 0; r < 10; ++r) md_forces_verlet(&sys, &pot, &vl);
    double dt_v = (md_wall_time_sec() - tv) * 1e6 / 10.0;
    double ev = max_rel_force_diff(f_ref, sys.force, sys.n);
    md_forces_threaded(&sys, &pot, &vl, 4);
    double et = max_rel_force_diff(f_ref, sys.force, sys.n);

    FILE *fp = fopen("out/s18_neighbor_2d.csv", "w");
    fprintf(fp, "method,n_particles,max_rel_err,pe_meas\n");
    fprintf(fp, "brute,%d,0.0,%.8f\ncell,%d,%.6e,%.8f\nverlet,%d,%.6e,%.8f\nthreaded4,%d,%.6e,%.8f\n",
            sys.n, pe_b, sys.n, ec, sys.pe, sys.n, ev, sys.pe, sys.n, et, sys.pe);
    fclose(fp);
    FILE *ft = fopen("out/s18_neighbor_2d_timing.csv", "w");
    fprintf(ft, "method,n_particles,time_us,speedup_vs_brute\n");
    fprintf(ft, "brute,%d,%.4f,1.0000\ncell,%d,%.4f,%.4f\nverlet,%d,%.4f,%.4f\n",
            sys.n, dt_b, sys.n, dt_c, dt_b / dt_c, sys.n, dt_v, dt_b / dt_v);
    fclose(ft);
    free(f_ref); md_verlet_free(&vl); md_cell_free(&cl); md_system_free(&sys);
    md_register_csv(ctx, "S18", "out/s18_neighbor_2d.csv");
    md_register_csv(ctx, "S18", "out/s18_neighbor_2d_timing.csv");
    md_ctx_record_stage(ctx, "S18", "neighbor_2d", md_wall_time_sec() - t0, 2);
}

void md_run_stages_b(md_run_ctx *ctx) {
    stage_s11_nbr_equiv(ctx);
    stage_s12_verlet_trigger(ctx);
    stage_s13_skin_sweep(ctx);
    stage_s14_s15_scaling(ctx);
    stage_s16_opt_ladder(ctx);
    stage_s17_thread_scaling(ctx);
    stage_s18_neighbor_2d(ctx);
}
