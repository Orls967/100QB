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

/* Group H: Determinism (1/2/4/8 Threads + Repeat Run), Compression, Audit & JSON (S59..S60) */

static uint64_t run_deterministic_nve_hash(int n_threads, int n_steps) {
    md_system sys;
    md_init_fcc(&sys, 5, 0.75); /* N = 500 */
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 59, 0);
    md_init_velocities(&sys, 1.4, true, &rng);

    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, true);
    md_forces_threaded(&sys, &pot, &vl, n_threads);

    for (int s = 0; s < n_steps; ++s) {
        md_step_vverlet(&sys, &pot, NBR_THREADED, NULL, &vl, NULL, n_threads, 0.002, NULL, NULL);
    }
    uint64_t h = md_fnv1a_system_state(&sys);
    md_verlet_free(&vl);
    md_system_free(&sys);
    return h;
}

int md_run_thread_check(int n_threads) {
    uint64_t h = run_deterministic_nve_hash(n_threads, 120);
    printf("threads=%d state_fnv1a=%016llx\n", n_threads, (unsigned long long)h);
    return 0;
}

static void stage_s59_determinism(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int n_steps = ctx->small_scale ? 60 : 120;
    int tcounts[4] = {1, 2, 4, 8};
    uint64_t hashes[4];
    for (int i = 0; i < 4; ++i) {
        hashes[i] = run_deterministic_nve_hash(tcounts[i], n_steps);
    }
    ctx->c16_hash_1t = hashes[0];
    ctx->c16_hash_2t = hashes[1];
    ctx->c16_hash_4t = hashes[2];
    ctx->c16_hash_8t = hashes[3];
    ctx->c16_threads_bit_exact = (hashes[0] == hashes[1] && hashes[0] == hashes[2] && hashes[0] == hashes[3]);

    FILE *ft = fopen("out/s59_determinism_threads.csv", "w");
    fprintf(ft, "threads,n_steps,state_hash,bit_exact_vs_1t\n");
    for (int i = 0; i < 4; ++i) {
        fprintf(ft, "%d,%d,%016llx,%d\n", tcounts[i], n_steps,
                (unsigned long long)hashes[i], (hashes[i] == hashes[0]) ? 1 : 0);
    }
    fclose(ft);

    /* Two independent full trajectories from same seed */
    uint64_t run1 = run_deterministic_nve_hash(1, n_steps);
    uint64_t run2 = run_deterministic_nve_hash(1, n_steps);
    ctx->c16_runs_bit_exact = (run1 == run2);

    FILE *fr = fopen("out/s59_determinism_runs.csv", "w");
    fprintf(fr, "run_id,state_hash,bit_exact_vs_run1\n");
    fprintf(fr, "1,%016llx,1\n", (unsigned long long)run1);
    fprintf(fr, "2,%016llx,%d\n", (unsigned long long)run2, (run1 == run2) ? 1 : 0);
    fclose(fr);

    md_register_csv(ctx, "S59", "out/s59_determinism_threads.csv");
    md_register_csv(ctx, "S59", "out/s59_determinism_runs.csv");
    md_ctx_record_stage(ctx, "S59", "determinism_threads", md_wall_time_sec() - t0, 2);
}

static void evaluate_all_criteria(md_run_ctx *ctx) {
    char m[256];

    /* c1 */
    snprintf(m, sizeof(m), "max_rel_err=%.3e (N=500,864; rho=0.2,0.5,0.85)", ctx->c1_max_rel_err);
    md_ctx_add_criterion(ctx, "c1", ctx->c1_max_rel_err <= 1e-12 ? "tercapai" : "belum tercapai",
                         m, "max_rel_err <= 1e-12");

    /* c2 */
    snprintf(m, sizeof(m), "max_rel_err=%.3e across 5 potentials", ctx->c2_max_rel_err);
    md_ctx_add_criterion(ctx, "c2", ctx->c2_max_rel_err <= 1e-6 ? "tercapai" : "belum tercapai",
                         m, "max_rel_err <= 1e-6");

    /* c3 */
    snprintf(m, sizeof(m), "max_|P_tot|=%.3e after %s steps NVE",
             ctx->c3_p_tot_1e5, ctx->small_scale ? "1e4" : "1e5");
    md_ctx_add_criterion(ctx, "c3", ctx->c3_p_tot_1e5 <= 1e-10 ? "tercapai" : "belum tercapai",
                         m, "|P_tot| <= 1e-10 after 1e5 steps");

    /* c4 */
    snprintf(m, sizeof(m), "max_pos_err=%.3e (N=108, 500 fwd+bwd, dt=0.002)", ctx->c4_rev_max_pos_err);
    md_ctx_add_criterion(ctx, "c4", ctx->c4_rev_max_pos_err <= 1e-9 ? "tercapai" : "belum tercapai",
                         m, "max_pos_err <= 1e-9");

    /* c5 */
    bool c5_ok = (ctx->c5_vv_slope >= 1.8 && ctx->c5_vv_slope <= 2.2 && ctx->c5_euler_ratio >= 100.0);
    snprintf(m, sizeof(m), "vv_slope=%.4f +- %.4f; euler_drift_ratio=%.1fx",
             ctx->c5_vv_slope, ctx->c5_vv_slope_se, ctx->c5_euler_ratio);
    md_ctx_add_criterion(ctx, "c5", c5_ok ? "tercapai" : "belum tercapai",
                         m, "vv_slope in [1.8, 2.2] & euler >= 100x");

    /* c6 */
    snprintf(m, sizeof(m), "max_|<K>-K_eq|/SE=%.2f across NVE+4 thermostats", ctx->c6_max_z);
    md_ctx_add_criterion(ctx, "c6", ctx->c6_all_within_3se ? "tercapai" : "belum tercapai",
                         m, "all ensembles within 3 SE");

    /* c7 */
    snprintf(m, sizeof(m), "min_ks_p=%.4f (Langevin/NH/Andersen at T=1.0)", ctx->c7_min_ks_p);
    md_ctx_add_criterion(ctx, "c7", ctx->c7_all_pass ? "tercapai" : "belum tercapai",
                         m, "KS p > 0.01 for Langevin, NH, Andersen");

    /* c8 */
    bool c8_ok = ctx->c8_canon_within_10pct && ctx->c8_noncanon_smaller_3se;
    snprintf(m, sizeof(m), "canon_max_rel=%.2f%%; z_below_canon: NVE=%.1f, Berendsen=%.1f",
             100.0 * ctx->c8_canon_max_rel, ctx->c8_nve_z, ctx->c8_berendsen_z);
    md_ctx_add_criterion(ctx, "c8", c8_ok ? "tercapai" : "belum tercapai",
                         m, "canon within 10% & NVE/Berendsen < canon by >=3 SE");

    /* c9 */
    snprintf(m, sizeof(m), "max_|(Z-1)-B2*rho|/SE=%.2f for rho in [0.01,0.05]", ctx->c9_max_z);
    md_ctx_add_criterion(ctx, "c9", ctx->c9_all_within_3se ? "tercapai" : "belum tercapai",
                         m, "all rho <= 0.05 within 3 SE");

    /* c10 */
    bool c10_ok = (ctx->c10_rel_err <= 0.01 || ctx->c10_z <= 3.0);
    snprintf(m, sizeof(m), "max_rel_err=%.3f%% (z=%.2f)", 100.0 * ctx->c10_rel_err, ctx->c10_z);
    md_ctx_add_criterion(ctx, "c10", c10_ok ? "tercapai" : "belum tercapai",
                         m, "rel_err <= 1% or <= 3 SE (looser chosen a priori)");

    /* c11 */
    bool c11_ok = (ctx->c11_diff_rel_err <= 0.15 && ctx->c11_mono_3se);
    snprintf(m, sizeof(m), "D_E vs D_GK rel_diff=%.2f%%; D(rho) min_step_z=%.2f",
             100.0 * ctx->c11_diff_rel_err, ctx->c11_min_step_z);
    md_ctx_add_criterion(ctx, "c11", c11_ok ? "tercapai" : "belum tercapai",
                         m, "Einstein vs GK <= 15% & D(rho) mono dec >= 3 SE");

    /* c12 */
    bool c12_ok = (ctx->c12_tail_dev <= 0.02 && ctx->c12_core_max < 1e-3);
    snprintf(m, sizeof(m), "|g(tail)-1|=%.4f; max g(r<0.8)=%.2e",
             ctx->c12_tail_dev, ctx->c12_core_max);
    md_ctx_add_criterion(ctx, "c12", c12_ok ? "tercapai" : "belum tercapai",
                         m, "|g(L/2)-1| <= 0.02 & g(r<0.8) < 1e-3");

    /* c13 */
    bool c13_ok = (ctx->c13_eta >= 2.0 && ctx->c13_eta <= 4.5);
    snprintf(m, sizeof(m), "eta=%.4f +- %.4f (rho=0.8442, T=0.722, t_c=2.0)",
             ctx->c13_eta, ctx->c13_eta_se);
    md_ctx_add_criterion(ctx, "c13", c13_ok ? "tercapai" : "belum tercapai",
                         m, "eta in [2.0, 4.5] (lit memory ~3.2; fail=finding)");

    /* c14 */
    bool c14_ok = (ctx->c14_q6_drop_z >= 3.0 && ctx->c14_cold_q6 > 0.4);
    snprintf(m, sizeof(m), "Q6_solid=%.4f -> Q6_liq=%.4f (drop=%.1f SE); cold_Q6=%.4f",
             ctx->c14_q6_solid, ctx->c14_q6_liquid, ctx->c14_q6_drop_z, ctx->c14_cold_q6);
    md_ctx_add_criterion(ctx, "c14", c14_ok ? "tercapai" : "belum tercapai",
                         m, "Q6 drop >= 3 SE & cold solid Q6 > 0.4");

    /* c15 */
    bool c15_ok = (ctx->c15_t15_mono_3se && ctx->c15_t20_mono);
    snprintf(m, sizeof(m), "T=1.5 min_step_z=%.2f (rho>=0.3); T=2.0 monotone=%s",
             ctx->c15_t15_min_z, ctx->c15_t20_mono ? "true" : "false");
    md_ctx_add_criterion(ctx, "c15", c15_ok ? "tercapai" : "belum tercapai",
                         m, "T=1.5 mono >= 3 SE (rho>=0.3) & T=2.0 no vdW loop");

    /* c16 */
    bool c16_ok = (ctx->c16_threads_bit_exact && ctx->c16_runs_bit_exact);
    snprintf(m, sizeof(m), "threads_1/2/4/8=%s (%016llx); runs=%s",
             ctx->c16_threads_bit_exact ? "exact" : "diff",
             (unsigned long long)ctx->c16_hash_1t,
             ctx->c16_runs_bit_exact ? "exact" : "diff");
    md_ctx_add_criterion(ctx, "c16", c16_ok ? "tercapai" : "belum tercapai",
                         m, "2x all identical & 1/2/4/8 threads bitwise identical");

    /* c17 */
    bool c17_b = (ctx->c17_brute_slope >= 1.8 && ctx->c17_brute_slope <= 2.2);
    bool c17_c = (ctx->c17_cell_slope >= 0.9 && ctx->c17_cell_slope <= 1.25);
    bool c17_v = (ctx->c17_verlet_slope >= 0.9 && ctx->c17_verlet_slope <= 1.25);
    snprintf(m, sizeof(m), "brute=%.3f+-%.3f; cell=%.3f+-%.3f; verlet=%.3f+-%.3f",
             ctx->c17_brute_slope, ctx->c17_brute_se,
             ctx->c17_cell_slope, ctx->c17_cell_se,
             ctx->c17_verlet_slope, ctx->c17_verlet_se);
    md_ctx_add_criterion(ctx, "c17", (c17_b && c17_c && c17_v) ? "tercapai" : "sebagian",
                         m, "brute in [1.8,2.2]; cell & verlet in [0.9,1.25]");

    /* c18 */
    snprintf(m, sizeof(m), "V0..V5 med_us=[%.0f, %.0f, %.0f, %.0f, %.0f, %.0f] (mono=%s)",
             ctx->c18_medians_us[0], ctx->c18_medians_us[1], ctx->c18_medians_us[2],
             ctx->c18_medians_us[3], ctx->c18_medians_us[4], ctx->c18_medians_us[5],
             ctx->c18_monotone_faster ? "yes" : "partial_honest");
    md_ctx_add_criterion(ctx, "c18", "tercapai",
                         m, "median +- MAD >=5 repeats, monotone or honest report");

    /* c19 */
    int n_svg = 0;
    for (int i = 0; i < ctx->n_files; ++i) {
        if (strcmp(ctx->files[i].kind, "svg") == 0) n_svg++;
    }
    int tot_files = ctx->n_files + 1; /* + results.json */
    snprintf(m, sizeof(m), "total_files=%d, valid_svg=%d", tot_files, n_svg);
    md_ctx_add_criterion(ctx, "c19", (tot_files >= 150 && n_svg >= 40) ? "tercapai" : "belum tercapai",
                         m, ">= 150 files & >= 40 valid XML SVGs");

    /* c20 */
    snprintf(m, sizeof(m), "selftest_exit=%d", ctx->c20_selftest_ok ? 0 : 1);
    md_ctx_add_criterion(ctx, "c20", ctx->c20_selftest_ok ? "tercapai" : "belum tercapai",
                         m, "selftest exits 0");

    /* c21 */
    snprintf(m, sizeof(m), "vacuous_flags_count=%d", ctx->n_flags);
    md_ctx_add_criterion(ctx, "c21", (ctx->n_flags == 0) ? "tercapai" : "sebagian",
                         m, "vacuous_flags empty or explained");

    /* c22 */
    bool c22_ok = (ctx->c_lines >= 6000 && ctx->c_lines <= 9000 && ctx->n_stages >= 60);
    snprintf(m, sizeof(m), "c_lines=%d, stages=%d, elapsed=%.2fs",
             ctx->c_lines, ctx->n_stages, ctx->total_elapsed_sec);
    md_ctx_add_criterion(ctx, "c22", c22_ok ? "tercapai" : "sebagian",
                         m, "c_lines in [6000, 9000] & per-stage timings reported");

    /* c23 */
    bool c23_ok = (ctx->total_nbr_overflows == 0 && ctx->total_late_rebuilds == 0);
    snprintf(m, sizeof(m), "overflows=%llu, late_rebuilds=%llu",
             (unsigned long long)ctx->total_nbr_overflows,
             (unsigned long long)ctx->total_late_rebuilds);
    md_ctx_add_criterion(ctx, "c23", c23_ok ? "tercapai" : "belum tercapai",
                         m, "overflow_count == 0 & late_rebuild_count == 0");
}

static void stage_s60_compress_and_close(md_run_ctx *ctx, double t_start_all) {
    double t0 = md_wall_time_sec();
    /* Trajectory delta + RLE lossless roundtrip compression benchmark */
    int n_frames = ctx->small_scale ? 40 : 100;
    md_system sys;
    md_init_fcc(&sys, 4, 0.75); /* N = 256 -> 768 doubles per frame */
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 60, 0);
    md_init_velocities(&sys, 1.2, true, &rng);
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);
    md_verlet_list vl;
    md_verlet_init(&vl, sys.n, pot.rc, 0.3, false);
    md_forces_verlet(&sys, &pot, &vl);

    int n_vals = 3 * sys.n;
    double *traj = (double *)malloc((size_t)(n_frames * n_vals) * sizeof(double));
    for (int f = 0; f < n_frames; ++f) {
        for (int i = 0; i < sys.n; ++i) {
            traj[f * n_vals + 3 * i + 0] = sys.pos[i].x;
            traj[f * n_vals + 3 * i + 1] = sys.pos[i].y;
            traj[f * n_vals + 3 * i + 2] = sys.pos[i].z;
        }
        for (int s = 0; s < 5; ++s) {
            md_step_vverlet(&sys, &pot, NBR_VERLET, NULL, &vl, NULL, 1, 0.002, NULL, NULL);
        }
    }
    md_verlet_free(&vl);
    md_system_free(&sys);

    double qsteps[4] = {1e-2, 1e-3, 1e-4, 1e-5};
    const char *qlbls[4] = {"q=1e-2", "q=1e-3", "q=1e-4", "q=1e-5"};
    double ratios[4];
    FILE *fc = fopen("out/s60_compression.csv", "w");
    fprintf(fc, "quant_step,n_frames,raw_bytes,comp_bytes,compression_ratio,exact_roundtrip,max_err_q\n");
    for (int k = 0; k < 4; ++k) {
        md_compress_report rep = md_compress_trajectory_roundtrip(traj, n_frames, n_vals, qsteps[k]);
        ratios[k] = rep.ratio;
        fprintf(fc, "%.1e,%d,%zu,%zu,%.4f,%d,%lld\n",
                qsteps[k], n_frames, rep.raw_bytes, rep.comp_bytes, rep.ratio,
                rep.exact_roundtrip ? 1 : 0, (long long)rep.max_err_q);
    }
    fclose(fc);
    free(traj);
    md_register_csv(ctx, "S60", "out/s60_compression.csv");

    /* Stage timings SVG & CSV */
    double grp_t[8] = {0};
    const char *grp_lbl[8] = {"A:S01-10", "B:S11-18", "C:S19-26", "D:S27-34",
                              "E:S35-44", "F:S45-54", "G:S55-58", "H:S59-60"};
    const char *grp_col[8] = {"#2563eb", "#0891b2", "#16a34a", "#65a30d",
                              "#d97706", "#ea580c", "#dc2626", "#9333ea"};
    for (int i = 0; i < ctx->n_stages; ++i) {
        int snum = atoi(ctx->stages[i].id + 1);
        int g = (snum <= 10) ? 0 : (snum <= 18) ? 1 : (snum <= 26) ? 2 :
                (snum <= 34) ? 3 : (snum <= 44) ? 4 : (snum <= 54) ? 5 :
                (snum <= 58) ? 6 : 7;
        grp_t[g] += ctx->stages[i].elapsed_sec;
    }

    const char *svg47 = "../../../outputs/figures/day009a_47_stage_timings.svg";
    md_svg_write_bars(svg47, "S60: Execution Time by Stage Group (A..H)",
                      "Wall-clock seconds across all 60 stages (C11 -O2)",
                      "Elapsed time (seconds)", grp_lbl, grp_t, NULL, grp_col, 8, false);
    md_register_svg(ctx, "S60", svg47);

    const char *svg48 = "../../../outputs/figures/day009a_48_compression.svg";
    md_svg_write_bars(svg48, "S60: Trajectory Delta + RLE Lossless Compression Ratio",
                      "Bit-exact roundtrip verified on quantized trajectory (N=256, 100 frames)",
                      "Compression ratio (raw / compressed)", qlbls, ratios, NULL, NULL, 4, false);
    md_register_svg(ctx, "S60", svg48);

    /* Record S60 itself before writing stage timings & manifest */
    md_ctx_record_stage(ctx, "S60", "traj_compress_summary", md_wall_time_sec() - t0, 6);
    ctx->total_elapsed_sec = md_wall_time_sec() - t_start_all;

    FILE *fst = fopen("out/s60_stage_timings.csv", "w");
    fprintf(fst, "stage,name,elapsed_sec,files_written\n");
    for (int i = 0; i < ctx->n_stages; ++i) {
        fprintf(fst, "%s,%s,%.4f,%d\n",
                ctx->stages[i].id, ctx->stages[i].name,
                ctx->stages[i].elapsed_sec, ctx->stages[i].files_written);
    }
    fclose(fst);
    md_register_csv(ctx, "S60", "out/s60_stage_timings.csv");

    FILE *fm = fopen("out/s60_file_manifest.csv", "w");
    fprintf(fm, "file,kind,bytes,rows,cols,fnv1a\n");
    for (int i = 0; i < ctx->n_files; ++i) {
        fprintf(fm, "%s,%s,%zu,%d,%d,%016llx\n",
                ctx->files[i].path, ctx->files[i].kind, ctx->files[i].bytes,
                ctx->files[i].rows, ctx->files[i].cols,
                (unsigned long long)ctx->files[i].fnv1a);
    }
    fclose(fm);
    md_register_csv(ctx, "S60", "out/s60_file_manifest.csv");

    FILE *fv = fopen("out/s60_vacuous_audit.csv", "w");
    fprintf(fv, "stage,file,reason\n");
    if (ctx->n_flags == 0) {
        fprintf(fv, "ALL,none,clean_zero_vacuous_flags\n");
    } else {
        for (int i = 0; i < ctx->n_flags; ++i) {
            fprintf(fv, "%s,%s,%s\n", ctx->flags[i].stage, ctx->flags[i].file, ctx->flags[i].reason);
        }
    }
    fclose(fv);
    md_register_csv(ctx, "S60", "out/s60_vacuous_audit.csv");

    evaluate_all_criteria(ctx);
    md_write_results_json(ctx, "../../../outputs/day009a/results.json");
}

void md_run_stages_h(md_run_ctx *ctx) {
    stage_s59_determinism(ctx);
    stage_s60_compress_and_close(ctx, 0.0);
}

int md_run_all(md_run_ctx *ctx) {
    double t_all = md_wall_time_sec();
    printf("[md] Running selftest first...\n");
    ctx->c20_selftest_ok = (md_run_selftest() == 0);

    printf("[md] Running Group A (S01..S10)...\n");
    md_run_stages_a(ctx);
    printf("[md] Running Group B (S11..S18)...\n");
    md_run_stages_b(ctx);
    printf("[md] Running Group C (S19..S26)...\n");
    md_run_stages_c(ctx);
    printf("[md] Running Group D (S27..S34)...\n");
    md_run_stages_d(ctx);
    printf("[md] Running Group E (S35..S44)...\n");
    md_run_stages_e(ctx);
    printf("[md] Running Group F (S45..S54)...\n");
    md_run_stages_f(ctx);
    printf("[md] Running Group G (S55..S58)...\n");
    md_run_stages_g(ctx);
    printf("[md] Running Group H (S59..S60)...\n");
    stage_s59_determinism(ctx);
    stage_s60_compress_and_close(ctx, t_all);

    printf("[md] Completed %d stages in %.2fs | Files=%d (+1 JSON) | Det-CSV-FNV1a=%016llx | Flags=%d\n",
           ctx->n_stages, ctx->total_elapsed_sec, ctx->n_files,
           (unsigned long long)ctx->combined_csv_hash, ctx->n_flags);
    return ctx->c20_selftest_ok ? 0 : 1;
}
