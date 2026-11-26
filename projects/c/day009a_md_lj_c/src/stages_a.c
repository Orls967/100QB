#include "md_forces.h"
#include "md_init.h"
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

/* Group A: Foundation stages S01..S10 */

static void stage_s01_rng(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int n_u = ctx->small_scale ? 20000 : 100000;
    double *u = (double *)malloc((size_t)n_u * sizeof(double));
    double *g = (double *)malloc((size_t)n_u * sizeof(double));

    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 1, 0);
    int n_bins = 20;
    double bins[20] = {0};
    double bin_c[20], theory_u[20];
    for (int i = 0; i < n_u; ++i) {
        u[i] = md_rng_uniform(&rng);
        int b = (int)(u[i] * (double)n_bins);
        if (b >= n_bins) b = n_bins - 1;
        bins[b] += 1.0;
        g[i] = md_rng_normal(&rng);
    }

    double u_mean = md_mean(u, n_u);
    double u_var = md_variance(u, n_u);
    double g_mean = md_mean(g, n_u);
    double g_var = md_variance(g, n_u);
    double acf_u[6] = {0};
    (void)md_autocorr_time(u, n_u, 6, acf_u);

    double expected = (double)n_u / (double)n_bins;
    double chi2 = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        double d = bins[b] - expected;
        chi2 += (d * d) / expected;
        bin_c[b] = ((double)b + 0.5) / (double)n_bins;
        theory_u[b] = 1.0;
        bins[b] /= ((double)n_u / (double)n_bins);
    }

    /* Check cross-stream independence across 16 stages */
    FILE *fp = fopen("out/s01_rng_moments.csv", "w");
    fprintf(fp, "stream_id,u_mean,u_var,g_mean,g_var,lag1_acf,chi2_20\n");
    fprintf(fp, "0,%.8f,%.8f,%.8f,%.8f,%.8f,%.6f\n",
            u_mean, u_var, g_mean, g_var, acf_u[1], chi2);
    int n_sub = ctx->small_scale ? 4000 : 12000;
    double *su = (double *)malloc((size_t)n_sub * sizeof(double));
    double *sg = (double *)malloc((size_t)n_sub * sizeof(double));
    for (int st = 1; st <= 15; ++st) {
        md_rng sr;
        md_rng_seed_stage(&sr, MD_MASTER_SEED, (uint64_t)(st + 1), 0);
        double sb[20] = {0};
        for (int i = 0; i < n_sub; ++i) {
            su[i] = md_rng_uniform(&sr);
            int b = (int)(su[i] * 20.0);
            if (b >= 20) b = 19;
            sb[b] += 1.0;
            sg[i] = md_rng_normal(&sr);
        }
        double sa[3] = {0};
        (void)md_autocorr_time(su, n_sub, 3, sa);
        double c2 = 0.0, ex = (double)n_sub / 20.0;
        for (int b = 0; b < 20; ++b) c2 += (sb[b] - ex) * (sb[b] - ex) / ex;
        fprintf(fp, "%d,%.8f,%.8f,%.8f,%.8f,%.8f,%.6f\n",
                st, md_mean(su, n_sub), md_variance(su, n_sub),
                md_mean(sg, n_sub), md_variance(sg, n_sub), sa[1], c2);
    }
    fclose(fp);
    free(su);
    free(sg);

    FILE *fb = fopen("out/s01_rng_bins.csv", "w");
    fprintf(fb, "bin,center,density,theory_density\n");
    for (int b = 0; b < n_bins; ++b) {
        fprintf(fb, "%d,%.4f,%.6f,%.6f\n", b, bin_c[b], bins[b], theory_u[b]);
    }
    fclose(fb);

    md_register_csv(ctx, "S01", "out/s01_rng_moments.csv");
    md_register_csv(ctx, "S01", "out/s01_rng_bins.csv");

    const char *svgp = "../../../outputs/figures/day009a_01_rng_hist.svg";
    md_svg_write_hist(svgp, "S01: xoshiro256++ Uniformity Audit",
                      "SplitMix64 -> xoshiro256++ (20 bins, normalized PDF)",
                      "u in [0, 1)", "Probability density", bin_c, bins, n_bins,
                      theory_u, "Uniform theory = 1.0");
    md_register_svg(ctx, "S01", svgp);

    free(u);
    free(g);
    md_ctx_record_stage(ctx, "S01", "rng_selftest", md_wall_time_sec() - t0, 3);
}

static void stage_s02_lattice(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    FILE *fp = fopen("out/s02_lattice_geom.csv", "w");
    fprintf(fp, "label,dim,n_particles,rho_target,rho_actual,box_l,min_pair_dist\n");

    int ncs[5] = {2, 3, 4, 5, 6};
    for (int k = 0; k < 5; ++k) {
        md_system sys;
        md_init_fcc(&sys, ncs[k], 0.8);
        double dmin = 1e30;
        for (int i = 0; i < sys.n; ++i) {
            for (int j = i + 1; j < sys.n; ++j) {
                double r = v3_norm(md_box_mic(&sys.box, sys.pos[i], sys.pos[j]));
                if (r < dmin) dmin = r;
            }
        }
        double rho_act = (double)sys.n / md_box_volume(&sys.box);
        fprintf(fp, "fcc_nc%d,3,%d,0.8000,%.8f,%.6f,%.8f\n",
                ncs[k], sys.n, rho_act, sys.box.lx, dmin);
        md_system_free(&sys);
    }

    /* Non-overlapping random placement for N = 108*k */
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 2, 0);
    int rns[3] = {108, 216, 324};
    for (int k = 0; k < 3; ++k) {
        md_system sys;
        (void)md_init_random_nonoverlap(&sys, rns[k], 0.3, 0.95, &rng);
        double dmin = 1e30;
        for (int i = 0; i < sys.n; ++i) {
            for (int j = i + 1; j < sys.n; ++j) {
                double r = v3_norm(md_box_mic(&sys.box, sys.pos[i], sys.pos[j]));
                if (r < dmin) dmin = r;
            }
        }
        fprintf(fp, "rand_%d,3,%d,0.3000,%.8f,%.6f,%.8f\n",
                rns[k], sys.n, (double)sys.n / md_box_volume(&sys.box), sys.box.lx, dmin);
        md_system_free(&sys);
    }

    /* 2D triangular lattice N=1024 */
    md_system sys2d;
    md_init_triangular_2d(&sys2d, 32, 32, 0.8);
    double dmin2 = 1e30;
    for (int i = 0; i < sys2d.n; ++i) {
        for (int j = i + 1; j < sys2d.n; ++j) {
            double r = v3_norm(md_box_mic(&sys2d.box, sys2d.pos[i], sys2d.pos[j]));
            if (r < dmin2) dmin2 = r;
        }
    }
    fprintf(fp, "tri2d_32x32,2,%d,0.8000,%.8f,%.6f,%.8f\n",
            sys2d.n, (double)sys2d.n / md_box_volume(&sys2d.box), sys2d.box.lx, dmin2);
    fclose(fp);

    FILE *f2 = fopen("out/s02_lattice_2d_coords.csv", "w");
    fprintf(f2, "idx,x,y\n");
    for (int i = 0; i < 128; ++i) {
        fprintf(f2, "%d,%.6f,%.6f\n", i, sys2d.pos[i].x, sys2d.pos[i].y);
    }
    fclose(f2);

    const char *svgp = "../../../outputs/figures/day009a_02_lattice_2d.svg";
    md_svg_write_particles_2d(svgp, "S02: 2D Triangular Lattice (N=1024, rho=0.8)",
                              "Periodic rectangular cell matching exact 6-fold triangular geometry",
                              &sys2d, NULL);
    md_system_free(&sys2d);

    md_register_csv(ctx, "S02", "out/s02_lattice_geom.csv");
    md_register_csv(ctx, "S02", "out/s02_lattice_2d_coords.csv");
    md_register_svg(ctx, "S02", svgp);
    md_ctx_record_stage(ctx, "S02", "lattice_gen", md_wall_time_sec() - t0, 3);
}

static void stage_s03_zero_momentum(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double temps[4] = {0.5, 1.0, 1.5, 2.0};
    FILE *fp = fopen("out/s03_init_velocities.csv", "w");
    fprintf(fp, "target_t,n_particles,p_norm,t_meas,ke_meas,ks_d,ks_p\n");

    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 3, 0);

    int n_bins = 24;
    double v_max = 5.0;
    double dv = v_max / (double)n_bins;
    double b_cent[24], b_pdf[24], b_theo[24];
    for (int b = 0; b < n_bins; ++b) {
        b_cent[b] = ((double)b + 0.5) * dv;
        b_pdf[b] = 0.0;
        double v = b_cent[b];
        b_theo[b] = sqrt(2.0 / MD_PI) * v * v * exp(-0.5 * v * v);
    }

    for (int k = 0; k < 4; ++k) {
        md_system sys;
        md_init_fcc(&sys, 6, 0.8); /* N=864 */
        md_init_velocities(&sys, temps[k], false, &rng);
        double p_norm = v3_norm(md_total_momentum(&sys));
        double t_m = md_temperature(&sys);
        double ke = md_kinetic_energy(&sys);
        double *sp = (double *)malloc((size_t)sys.n * sizeof(double));
        for (int i = 0; i < sys.n; ++i) {
            sp[i] = v3_norm(sys.vel[i]);
            if (k == 1) {
                int b = (int)(sp[i] / dv);
                if (b >= 0 && b < n_bins) b_pdf[b] += 1.0 / ((double)sys.n * dv);
            }
        }
        double kd = 0.0, kp = 0.0;
        md_ks_maxwell(sp, sys.n, temps[k], &kd, &kp);
        fprintf(fp, "%.2f,%d,%.6e,%.8f,%.6f,%.6f,%.6f\n",
                temps[k], sys.n, p_norm, t_m, ke, kd, kp);
        free(sp);
        md_system_free(&sys);
    }
    fclose(fp);

    FILE *fh = fopen("out/s03_vel_hist.csv", "w");
    fprintf(fh, "v,empirical_pdf,maxwell_theory_pdf\n");
    for (int b = 0; b < n_bins; ++b) {
        fprintf(fh, "%.4f,%.6f,%.6f\n", b_cent[b], b_pdf[b], b_theo[b]);
    }
    fclose(fh);

    md_register_csv(ctx, "S03", "out/s03_init_velocities.csv");
    md_register_csv(ctx, "S03", "out/s03_vel_hist.csv");

    const char *svgp = "../../../outputs/figures/day009a_03_init_maxwell.svg";
    md_svg_write_hist(svgp, "S03: Initial Velocity Distribution (T=1.0, N=864)",
                      "Zero-COM Maxwell-Boltzmann initialization vs exact 3D Maxwell PDF",
                      "Speed v", "f(v)", b_cent, b_pdf, n_bins, b_theo, "Maxwell-Boltzmann T=1.0");
    md_register_svg(ctx, "S03", svgp);
    md_ctx_record_stage(ctx, "S03", "zero_momentum", md_wall_time_sec() - t0, 3);
}

static void stage_s04_pbc_mic(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 4, 0);
    md_box b3 = md_box_cubic(8.0);
    md_box b2 = md_box_ortho(10.0, 9.0, 1.0, 2);

    FILE *f3 = fopen("out/s04_pbc_mic.csv", "w");
    fprintf(f3, "idx,shift_x,shift_y,shift_z,mic_dist,antisym_err,wrap_err\n");
    for (int i = 0; i < 32; ++i) {
        md_vec3 p1 = v3_make(md_rng_range(&rng, -25.0, 25.0),
                             md_rng_range(&rng, -25.0, 25.0),
                             md_rng_range(&rng, -25.0, 25.0));
        md_vec3 p2 = v3_make(md_rng_range(&rng, -25.0, 25.0),
                             md_rng_range(&rng, -25.0, 25.0),
                             md_rng_range(&rng, -25.0, 25.0));
        md_vec3 w1 = md_box_wrap(&b3, p1);
        md_vec3 w2 = md_box_wrap(&b3, p2);
        md_vec3 d12 = md_box_mic(&b3, w1, w2);
        md_vec3 d21 = md_box_mic(&b3, w2, w1);
        md_vec3 d_raw = md_box_mic(&b3, p1, p2);
        double anti = v3_norm(v3_add(d12, d21));
        double werr = v3_norm(v3_sub(d12, d_raw));
        fprintf(f3, "%d,%.4f,%.4f,%.4f,%.8f,%.6e,%.6e\n",
                i, p1.x - w1.x, p1.y - w1.y, p1.z - w1.z, v3_norm(d12), anti, werr);
    }
    fclose(f3);

    FILE *f2 = fopen("out/s04_pbc_2d.csv", "w");
    fprintf(f2, "idx,mic_dist_2d,antisym_err,wrap_err\n");
    for (int i = 0; i < 24; ++i) {
        md_vec3 p1 = v3_make(md_rng_range(&rng, -20.0, 20.0), md_rng_range(&rng, -20.0, 20.0), 0.0);
        md_vec3 p2 = v3_make(md_rng_range(&rng, -20.0, 20.0), md_rng_range(&rng, -20.0, 20.0), 0.0);
        md_vec3 w1 = md_box_wrap(&b2, p1);
        md_vec3 w2 = md_box_wrap(&b2, p2);
        md_vec3 d12 = md_box_mic(&b2, w1, w2);
        md_vec3 d21 = md_box_mic(&b2, w2, w1);
        md_vec3 d_raw = md_box_mic(&b2, p1, p2);
        fprintf(f2, "%d,%.8f,%.6e,%.6e\n",
                i, v3_norm(d12), v3_norm(v3_add(d12, d21)), v3_norm(v3_sub(d12, d_raw)));
    }
    fclose(f2);

    md_register_csv(ctx, "S04", "out/s04_pbc_mic.csv");
    md_register_csv(ctx, "S04", "out/s04_pbc_2d.csv");
    md_ctx_record_stage(ctx, "S04", "pbc_mic", md_wall_time_sec() - t0, 2);
}

static void stage_s05_force_fd(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    double h = 1e-6;
    int n_r = 60;
    double r_arr[POT_COUNT][60], u_arr[POT_COUNT][60], err_arr[POT_COUNT][60];
    double max_rel_all = 0.0;

    FILE *fs = fopen("out/s05_force_fd_scan.csv", "w");
    fprintf(fs, "potential,r,u_exact,f_analytic,f_fd,rel_err\n");
    FILE *fsum = fopen("out/s05_force_fd_summary.csv", "w");
    fprintf(fsum, "potential,rc,max_rel_err,mean_rel_err\n");

    for (int k = 0; k < POT_COUNT; ++k) {
        md_potential pot;
        md_potential_init(&pot, (md_pot_kind)k, 2.5);
        double r_lo = 0.90;
        double r_hi = pot.rc - 0.02;
        double max_e = 0.0, sum_e = 0.0;
        for (int i = 0; i < n_r; ++i) {
            double r = r_lo + (r_hi - r_lo) * (double)i / (double)(n_r - 1);
            double fdr = 0.0;
            double u = md_pot_eval_exact(&pot, r * r, &fdr);
            double f_ana = fdr * r;
            double f_fd = md_pot_force_fd(&pot, r, h);
            double rel = fabs(f_ana - f_fd) / (fabs(f_ana) + 1e-6);
            if (rel > max_e) max_e = rel;
            sum_e += rel;
            r_arr[k][i] = r;
            u_arr[k][i] = (u > 5.0) ? 5.0 : u;
            err_arr[k][i] = (rel > 1e-16) ? rel : 1e-16;
            fprintf(fs, "%s,%.6f,%.8f,%.8f,%.8f,%.6e\n",
                    md_pot_name((md_pot_kind)k), r, u, f_ana, f_fd, rel);
        }
        if (max_e > max_rel_all) max_rel_all = max_e;
        fprintf(fsum, "%s,%.6f,%.6e,%.6e\n",
                md_pot_name((md_pot_kind)k), pot.rc, max_e, sum_e / (double)n_r);
        md_potential_free(&pot);
    }
    fclose(fs);
    fclose(fsum);
    ctx->c2_max_rel_err = max_rel_all;

    md_register_csv(ctx, "S05", "out/s05_force_fd_scan.csv");
    md_register_csv(ctx, "S05", "out/s05_force_fd_summary.csv");

    const char *cols[POT_COUNT] = {"#2563eb", "#0891b2", "#dc2626", "#16a34a", "#9333ea"};
    md_svg_series ser_u[POT_COUNT], ser_e[POT_COUNT];
    for (int k = 0; k < POT_COUNT; ++k) {
        ser_u[k] = (md_svg_series){md_pot_name((md_pot_kind)k), cols[k], r_arr[k], u_arr[k], NULL, n_r, true, false, k == 1};
        ser_e[k] = (md_svg_series){md_pot_name((md_pot_kind)k), cols[k], r_arr[k], err_arr[k], NULL, n_r, true, false, false};
    }

    md_svg_opts o1 = {"S05: Pair Potentials u(r)", "All 5 implemented potentials (capped at u=5 for display)",
                      "Distance r", "u(r)", 760, 460, false, false, true, 0.0, "u = 0"};
    const char *svg1 = "../../../outputs/figures/day009a_04_potentials.svg";
    md_svg_write_xy(svg1, &o1, ser_u, POT_COUNT);
    md_register_svg(ctx, "S05", svg1);

    md_svg_opts o2 = {"S05: Analytic Force vs Central Finite Difference -dU/dr",
                      "Relative error across all 5 potentials (h = 1e-6; threshold 1e-6)",
                      "Distance r", "Relative error |F_ana - F_fd| / (|F_ana| + 1e-6)",
                      760, 460, false, true, true, 1e-6, "c2 threshold (1e-6)"};
    const char *svg2 = "../../../outputs/figures/day009a_05_force_fd_err.svg";
    md_svg_write_xy(svg2, &o2, ser_e, POT_COUNT);
    md_register_svg(ctx, "S05", svg2);

    md_ctx_record_stage(ctx, "S05", "force_fd", md_wall_time_sec() - t0, 4);
}

static void stage_s06_newton3(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 6, 0);

    md_system s3;
    md_init_fcc(&s3, 4, 0.8);
    for (int i = 0; i < s3.n; ++i) {
        s3.pos[i].x += md_rng_range(&rng, -0.08, 0.08);
        s3.pos[i].y += md_rng_range(&rng, -0.08, 0.08);
        s3.pos[i].z += md_rng_range(&rng, -0.08, 0.08);
        s3.pos[i] = md_box_wrap(&s3.box, s3.pos[i]);
    }

    FILE *f3 = fopen("out/s06_newton3.csv", "w");
    fprintf(f3, "potential,dim,n_particles,net_force_norm,max_pair_asym,pe_per_n\n");
    for (int k = 0; k < POT_COUNT; ++k) {
        md_potential pot;
        md_potential_init(&pot, (md_pot_kind)k, 2.5);
        md_forces_brute(&s3, &pot);
        md_vec3 sumf = v3_zero();
        for (int i = 0; i < s3.n; ++i) sumf = v3_add(sumf, s3.force[i]);
        /* Direct pair antisymmetry check on first 30 pairs */
        double max_asym = 0.0;
        for (int i = 0; i < 15; ++i) {
            for (int j = i + 1; j < 20; ++j) {
                md_vec3 dij = md_box_mic(&s3.box, s3.pos[i], s3.pos[j]);
                md_vec3 dji = md_box_mic(&s3.box, s3.pos[j], s3.pos[i]);
                double fij = 0.0, fji = 0.0;
                (void)md_pot_eval_exact(&pot, v3_norm2(dij), &fij);
                (void)md_pot_eval_exact(&pot, v3_norm2(dji), &fji);
                md_vec3 vsum = v3_add(v3_scale(dij, fij), v3_scale(dji, fji));
                double e = v3_norm(vsum);
                if (e > max_asym) max_asym = e;
            }
        }
        fprintf(f3, "%s,3,%d,%.6e,%.6e,%.8f\n",
                md_pot_name((md_pot_kind)k), s3.n, v3_norm(sumf), max_asym, s3.pe / (double)s3.n);
        md_potential_free(&pot);
    }
    fclose(f3);
    md_system_free(&s3);

    md_system s2;
    md_init_triangular_2d(&s2, 16, 16, 0.7);
    for (int i = 0; i < s2.n; ++i) {
        s2.pos[i].x += md_rng_range(&rng, -0.07, 0.07);
        s2.pos[i].y += md_rng_range(&rng, -0.07, 0.07);
        s2.pos[i] = md_box_wrap(&s2.box, s2.pos[i]);
    }
    FILE *f2 = fopen("out/s06_newton3_2d.csv", "w");
    fprintf(f2, "potential,dim,n_particles,net_force_norm,pe_per_n\n");
    for (int k = 0; k < POT_COUNT; ++k) {
        md_potential pot;
        md_potential_init(&pot, (md_pot_kind)k, 2.5);
        md_forces_brute(&s2, &pot);
        md_vec3 sumf = v3_zero();
        for (int i = 0; i < s2.n; ++i) sumf = v3_add(sumf, s2.force[i]);
        fprintf(f2, "%s,2,%d,%.6e,%.8f\n",
                md_pot_name((md_pot_kind)k), s2.n, v3_norm(sumf), s2.pe / (double)s2.n);
        md_potential_free(&pot);
    }
    fclose(f2);
    md_system_free(&s2);

    md_register_csv(ctx, "S06", "out/s06_newton3.csv");
    md_register_csv(ctx, "S06", "out/s06_newton3_2d.csv");
    md_ctx_record_stage(ctx, "S06", "newton3", md_wall_time_sec() - t0, 2);
}

static void stage_s07_pot_tables(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    int ns[4] = {256, 1024, 4096, 16384};
    double nx[4], lin_err[4], cub_err[4];

    FILE *fa = fopen("out/s07_table_accuracy.csv", "w");
    fprintf(fa, "potential,n_pts,mode,max_u_err,max_f_err\n");

    for (int k = 0; k < POT_COUNT; ++k) {
        for (int m = 0; m < 4; ++m) {
            for (int cub = 0; cub <= 1; ++cub) {
                md_potential pot;
                md_potential_init(&pot, (md_pot_kind)k, 2.5);
                md_pot_build_table(&pot, ns[m], 0.70, cub == 1);
                double max_u = 0.0, max_f = 0.0;
                for (int i = 0; i < 200; ++i) {
                    double r = 0.88 + (pot.rc - 0.90) * ((double)i + 0.41) / 200.0;
                    double fe = 0.0, ft = 0.0;
                    double ue = md_pot_eval_exact(&pot, r * r, &fe);
                    double ut = md_pot_eval(&pot, r * r, &ft);
                    double du = fabs(ut - ue);
                    double df = fabs(ft - fe);
                    if (du > max_u) max_u = du;
                    if (df > max_f) max_f = df;
                }
                fprintf(fa, "%s,%d,%s,%.6e,%.6e\n",
                        md_pot_name((md_pot_kind)k), ns[m], cub ? "cubic" : "linear", max_u, max_f);
                if (k == POT_LJ_SHIFT) {
                    nx[m] = (double)ns[m];
                    if (cub == 0) lin_err[m] = max_f;
                    else cub_err[m] = (max_f > 1e-16) ? max_f : 1e-16;
                }
                md_potential_free(&pot);
            }
        }
    }
    fclose(fa);

    FILE *fp = fopen("out/s07_table_profile.csv", "w");
    fprintf(fp, "r,u_exact,u_cubic_4096,f_exact,f_cubic_4096,abs_f_err\n");
    md_potential lj;
    md_potential_init(&lj, POT_LJ_SHIFT, 2.5);
    md_pot_build_table(&lj, 4096, 0.70, true);
    for (int i = 0; i < 50; ++i) {
        double r = 0.90 + 1.55 * (double)i / 49.0;
        double fe = 0.0, ft = 0.0;
        double ue = md_pot_eval_exact(&lj, r * r, &fe);
        double ut = md_pot_eval(&lj, r * r, &ft);
        fprintf(fp, "%.6f,%.8f,%.8f,%.8f,%.8f,%.6e\n", r, ue, ut, fe * r, ft * r, fabs(ft - fe) * r);
    }
    md_potential_free(&lj);
    fclose(fp);

    md_register_csv(ctx, "S07", "out/s07_table_accuracy.csv");
    md_register_csv(ctx, "S07", "out/s07_table_profile.csv");

    md_svg_series ser[2] = {
        {"Linear table (O(h^2))", "#dc2626", nx, lin_err, NULL, 4, true, true, false},
        {"Cubic Hermite (O(h^4))", "#2563eb", nx, cub_err, NULL, 4, true, true, false}
    };
    md_svg_opts o = {"S07: Potential Lookup Table Convergence (LJ-shift)",
                     "Max force error |f/r_tab - f/r_exact| vs table size M",
                     "Table grid points M", "Max force error", 760, 460, true, true, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_06_table_convergence.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S07", svgp);

    md_ctx_record_stage(ctx, "S07", "pot_tables", md_wall_time_sec() - t0, 3);
}

static void stage_s08_small_fixtures(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 2.5);

    FILE *fp = fopen("out/s08_small_fixtures.csv", "w");
    fprintf(fp, "label,n_particles,r_sep,u_analytic,u_meas,u_err,vir_analytic,vir_meas,vir_err\n");
    FILE *ff = fopen("out/s08_fixture_forces.csv", "w");
    fprintf(ff, "label,particle,fx,fy,fz,f_norm\n");

    double r_test = 1.25;
    double fdr = 0.0;
    double u_pair = md_pot_eval_exact(&pot, r_test * r_test, &fdr);
    double w_pair = r_test * r_test * fdr;

    /* 2-particle dimer (1 pair) */
    md_system s2;
    md_init_fixture_dimer(&s2, r_test, 10.0);
    md_forces_brute(&s2, &pot);
    fprintf(fp, "dimer,2,%.4f,%.10f,%.10f,%.4e,%.10f,%.10f,%.4e\n",
            r_test, u_pair, s2.pe, fabs(s2.pe - u_pair), w_pair, s2.virial, fabs(s2.virial - w_pair));
    for (int i = 0; i < 2; ++i)
        fprintf(ff, "dimer,%d,%.8f,%.8f,%.8f,%.8f\n", i, s2.force[i].x, s2.force[i].y, s2.force[i].z, v3_norm(s2.force[i]));
    md_system_free(&s2);

    /* 3-particle equilateral triangle (3 pairs) */
    md_system s3;
    md_init_fixture_triangle(&s3, r_test, 10.0);
    md_forces_brute(&s3, &pot);
    fprintf(fp, "triangle,3,%.4f,%.10f,%.10f,%.4e,%.10f,%.10f,%.4e\n",
            r_test, 3.0 * u_pair, s3.pe, fabs(s3.pe - 3.0 * u_pair),
            3.0 * w_pair, s3.virial, fabs(s3.virial - 3.0 * w_pair));
    for (int i = 0; i < 3; ++i)
        fprintf(ff, "triangle,%d,%.8f,%.8f,%.8f,%.8f\n", i, s3.force[i].x, s3.force[i].y, s3.force[i].z, v3_norm(s3.force[i]));
    md_system_free(&s3);

    /* 4-particle regular tetrahedron (6 pairs) */
    md_system s4;
    md_init_fixture_tetrahedron(&s4, r_test, 10.0);
    md_forces_brute(&s4, &pot);
    fprintf(fp, "tetrahedron,4,%.4f,%.10f,%.10f,%.4e,%.10f,%.10f,%.4e\n",
            r_test, 6.0 * u_pair, s4.pe, fabs(s4.pe - 6.0 * u_pair),
            6.0 * w_pair, s4.virial, fabs(s4.virial - 6.0 * w_pair));
    for (int i = 0; i < 4; ++i)
        fprintf(ff, "tetrahedron,%d,%.8f,%.8f,%.8f,%.8f\n", i, s4.force[i].x, s4.force[i].y, s4.force[i].z, v3_norm(s4.force[i]));
    md_system_free(&s4);

    fclose(fp);
    fclose(ff);
    md_register_csv(ctx, "S08", "out/s08_small_fixtures.csv");
    md_register_csv(ctx, "S08", "out/s08_fixture_forces.csv");
    md_ctx_record_stage(ctx, "S08", "small_fixtures", md_wall_time_sec() - t0, 2);
}

static void stage_s09_fcc_shells(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    /* Exact FCC coordination shells: r/a = sqrt(m/2) for m = 1,2,3,4,5,6 */
    static const int expected_cnt[6] = {12, 6, 24, 12, 24, 8};
    double a = cbrt(4.0 / 1.0); /* rho = 1.0 */

    md_system sys;
    md_init_fcc(&sys, 5, 1.0); /* N=500, L=7.937 > 2*r_shell6 */

    FILE *fs = fopen("out/s09_fcc_shells.csv", "w");
    fprintf(fs, "shell,r_expected,expected_count,measured_count,u_shell_contrib\n");
    md_potential pot;
    md_potential_init(&pot, POT_LJ_SHIFT, 3.0);

    double sh_r[6], sh_cnt[6];
    const char *sh_lbl[6] = {"S1 (12)", "S2 (6)", "S3 (24)", "S4 (12)", "S5 (24)", "S6 (8)"};
    for (int s = 0; s < 6; ++s) {
        double r_exp = a * sqrt(0.5 * (double)(s + 1));
        int cnt = 0;
        for (int j = 1; j < sys.n; ++j) {
            double r = v3_norm(md_box_mic(&sys.box, sys.pos[0], sys.pos[j]));
            if (fabs(r - r_exp) < 1e-4) cnt++;
        }
        double fdr = 0.0;
        double up = md_pot_eval_exact(&pot, r_exp * r_exp, &fdr);
        fprintf(fs, "%d,%.6f,%d,%d,%.8f\n", s + 1, r_exp, expected_cnt[s], cnt, 0.5 * (double)cnt * up);
        sh_r[s] = r_exp;
        sh_cnt[s] = (double)cnt;
    }
    (void)sh_r;
    fclose(fs);
    md_system_free(&sys);

    /* FCC lattice energy vs N at rho = 1.0 */
    FILE *fe = fopen("out/s09_fcc_energy.csv", "w");
    fprintf(fe, "nc,n_particles,u_fcc_formula,u_meas_per_n,abs_err\n");
    md_potential lj;
    md_potential_init(&lj, POT_LJ_SHIFT, 2.5);
    double u_form = 0.0;
    for (int s = 0; s < 6; ++s) {
        double r_exp = a * sqrt(0.5 * (double)(s + 1));
        if (r_exp < lj.rc) {
            double fdr = 0.0;
            u_form += 0.5 * (double)expected_cnt[s] * md_pot_eval_exact(&lj, r_exp * r_exp, &fdr);
        }
    }
    int ncs[4] = {3, 4, 5, 6};
    for (int k = 0; k < 4; ++k) {
        md_system fsys;
        md_init_fcc(&fsys, ncs[k], 1.0);
        md_forces_brute(&fsys, &lj);
        double um = fsys.pe / (double)fsys.n;
        fprintf(fe, "%d,%d,%.10f,%.10f,%.4e\n", ncs[k], fsys.n, u_form, um, fabs(um - u_form));
        md_system_free(&fsys);
    }
    fclose(fe);

    md_register_csv(ctx, "S09", "out/s09_fcc_shells.csv");
    md_register_csv(ctx, "S09", "out/s09_fcc_energy.csv");

    const char *svgp = "../../../outputs/figures/day009a_07_fcc_shells.svg";
    md_svg_write_bars(svgp, "S09: FCC Neighbor Shell Multiplicities (N=500, rho=1.0)",
                      "Exact coordination numbers (12, 6, 24, 12, 24, 8) verified in periodic box",
                      "Neighbors in shell", sh_lbl, sh_cnt, NULL, NULL, 6, false);
    md_register_svg(ctx, "S09", svgp);
    md_ctx_record_stage(ctx, "S09", "fcc_shells", md_wall_time_sec() - t0, 3);
}

static void stage_s10_tail_corr(md_run_ctx *ctx) {
    double t0 = md_wall_time_sec();
    /* Compare truncated vs tail-corrected energy on an FCC crystal (N=864, L=10.26) as rc increases */
    md_system sys;
    md_init_fcc(&sys, 6, 0.8);
    md_rng rng;
    md_rng_seed_stage(&rng, MD_MASTER_SEED, 10, 0);
    for (int i = 0; i < sys.n; ++i) {
        sys.pos[i].x += md_rng_range(&rng, -0.05, 0.05);
        sys.pos[i].y += md_rng_range(&rng, -0.05, 0.05);
        sys.pos[i].z += md_rng_range(&rng, -0.05, 0.05);
        sys.pos[i] = md_box_wrap(&sys.box, sys.pos[i]);
    }

    double rcs[6] = {2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
    double u_shift[6], u_tail[6], p_tail_arr[6];
    FILE *fc = fopen("out/s10_tail_convergence.csv", "w");
    fprintf(fc, "rc,u_shift_per_n,u_trunc_per_n,u_tail_corr_per_n,u_total_tail_per_n,p_tail_corr\n");
    for (int k = 0; k < 6; ++k) {
        md_potential ps, pt;
        md_potential_init(&ps, POT_LJ_SHIFT, rcs[k]);
        md_potential_init(&pt, POT_LJ_TAIL, rcs[k]);
        md_forces_brute(&sys, &ps);
        double us = sys.pe / (double)sys.n;
        md_forces_brute(&sys, &pt);
        double ut_tot = sys.pe / (double)sys.n;
        double ut_corr = sys.pe_tail / (double)sys.n;
        double ut_raw = ut_tot - ut_corr;
        u_shift[k] = us;
        u_tail[k] = ut_tot;
        p_tail_arr[k] = sys.p_tail;
        fprintf(fc, "%.2f,%.8f,%.8f,%.8f,%.8f,%.8f\n",
                rcs[k], us, ut_raw, ut_corr, ut_tot, p_tail_arr[k]);
    }
    fclose(fc);
    md_system_free(&sys);

    /* Analytical tail corrections vs density at rc = 2.5 */
    FILE *fd = fopen("out/s10_tail_density.csv", "w");
    fprintf(fd, "rho,u_tail_per_n,p_tail\n");
    md_potential p25;
    md_potential_init(&p25, POT_LJ_TAIL, 2.5);
    for (int i = 1; i <= 10; ++i) {
        double rho = 0.1 * (double)i;
        double ut = 0.0, pt = 0.0;
        md_pot_tail_corrections(&p25, 1, 1.0 / rho, &ut, &pt);
        fprintf(fd, "%.2f,%.8f,%.8f\n", rho, ut, pt);
    }
    fclose(fd);

    md_register_csv(ctx, "S10", "out/s10_tail_convergence.csv");
    md_register_csv(ctx, "S10", "out/s10_tail_density.csv");

    md_svg_series ser[2] = {
        {"LJ truncated & shifted", "#dc2626", rcs, u_shift, NULL, 6, true, true, false},
        {"LJ truncated + analytical tail", "#2563eb", rcs, u_tail, NULL, 6, true, true, false}
    };
    md_svg_opts o = {"S10: Cutoff Convergence — Shifted vs Tail-Corrected LJ",
                     "Perturbed FCC (N=864, rho=0.8): tail correction accelerates convergence with r_c",
                     "Cutoff radius r_c", "Potential energy U/N", 760, 460, false, false, false, 0.0, NULL};
    const char *svgp = "../../../outputs/figures/day009a_08_tail_corr.svg";
    md_svg_write_xy(svgp, &o, ser, 2);
    md_register_svg(ctx, "S10", svgp);
    md_ctx_record_stage(ctx, "S10", "pot_tail_corr", md_wall_time_sec() - t0, 3);
}

void md_run_stages_a(md_run_ctx *ctx) {
    stage_s01_rng(ctx);
    stage_s02_lattice(ctx);
    stage_s03_zero_momentum(ctx);
    stage_s04_pbc_mic(ctx);
    stage_s05_force_fd(ctx);
    stage_s06_newton3(ctx);
    stage_s07_pot_tables(ctx);
    stage_s08_small_fixtures(ctx);
    stage_s09_fcc_shells(ctx);
    stage_s10_tail_corr(ctx);
}
