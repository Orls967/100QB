#define _POSIX_C_SOURCE 200809L
#include "sandpile.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void ols_fit(const double *x, const double *y, int n, double *slope, double *se) {
    double sx = 0, sy = 0;
    for (int i = 0; i < n; i++) { sx += x[i]; sy += y[i]; }
    double mx = sx / n, my = sy / n, sxx = 0, sxy = 0, ssr = 0;
    for (int i = 0; i < n; i++) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); }
    *slope = (sxx > 0) ? (sxy / sxx) : 0.0;
    double a = my - (*slope) * mx;
    for (int i = 0; i < n; i++) { double e = y[i] - (a + (*slope) * x[i]); ssr += e * e; }
    *se = (n > 2 && sxx > 0) ? sqrt((ssr / (n - 2)) / sxx) : 0.0;
}

static void sort5(double a[5]) {
    for (int i = 0; i < 5; i++)
        for (int j = i + 1; j < 5; j++)
            if (a[j] < a[i]) { double t = a[i]; a[i] = a[j]; a[j] = t; }
}

int cmd_radius(Metrics *m, const char *out_dir) {
    const int64_t Ns[4] = {1000, 10000, 100000, 1000000};
    double lx[4], ly[4];
    char path[128]; snprintf(path, sizeof(path), "%s/radius.csv", out_dir);
    FILE *fp = fopen(path, "w");
    if (fp) fprintf(fp, "N,grid_n,R,area,sum_z,sink,mass_diff\n");
    for (int k = 0; k < 4; k++) {
        int64_t N = Ns[k], area = 0, sum_z = 0;
        int n = 2 * ((int)(sqrt((double)N) / 1.4) + 8) + 1, c0 = n / 2;
        Grid g; grid_init(&g, n);
        g.z[c0 * n + c0] = (int32_t)N;
        stabilize_det(&g);
        double max_r2 = 0;
        for (int i = 0; i < n * n; i++)
            if (g.z[i] > 0) {
                area++; sum_z += g.z[i];
                double dr = (double)(i / n - c0), dc = (double)(i % n - c0), r2 = dr * dr + dc * dc;
                if (r2 > max_r2) max_r2 = r2;
            }
        double R = sqrt(max_r2);
        int64_t md = llabs((sum_z + g.sink) - N);
        if (md > m->mass_max_diff) m->mass_max_diff = md;
        lx[k] = log((double)N); ly[k] = log(R);
        if (fp) fprintf(fp, "%lld,%d,%.6f,%lld,%lld,%lld,%lld\n",
                        (long long)N, n, R, (long long)area, (long long)sum_z, (long long)g.sink, (long long)md);
        grid_free(&g);
    }
    if (fp) fclose(fp);
    ols_fit(lx, ly, 4, &m->rad_slope, &m->rad_se);
    return 0;
}

static void fit_hist_window(const int64_t *cnt, const double *sum_s, int64_t tot,
                            int b_lo, int b_hi, double *tau, double *se) {
    double xs[24], ys[24]; int k = 0;
    for (int b = b_lo; b <= b_hi; b++) {
        if (cnt[b] < 20) continue;
        double ms = sum_s[b] / (double)cnt[b], ps = (double)cnt[b] / ((double)tot * (double)(1ULL << b));
        xs[k] = log(ms); ys[k] = -log(ps); k++;
    }
    ols_fit(xs, ys, k, tau, se);
}

static void run_soc_one(int n, uint64_t seed_tag, const char *out_dir, Metrics *m, int32_t *save_128) {
    int sz = n * n, stable_streak = 0, blk_snaps = 0;
    Grid g; grid_init(&g, n);
    Rng rng; rng_seed(&rng, MASTER_SEED ^ seed_tag);
    int64_t added = 0, w_drops = 0, hist_cnt[24] = {0}, target_aval = 2000000, aval_done = 0;
    double prev_mean = -1.0, hist_s[24] = {0}, hist_a[24] = {0}, hist_d[24] = {0};
    double blk_p[20][4] = {{0}}, blk_mean[20] = {0}, full_sum = 0;
    while (stable_streak < 3 || w_drops < 4LL * sz) {
        for (int step = 0; step < sz; step++) {
            int p = rng_int(&rng, 0, sz - 1);
            added++;
            if (++g.z[p] >= 4) stabilize_det(&g);
        }
        w_drops += sz;
        int64_t s = 0;
        for (int i = 0; i < sz; i++) s += g.z[i];
        double cur = (double)s / (double)sz;
        stable_streak = (prev_mean > 0.5 && fabs(cur - prev_mean) / prev_mean < 0.01) ? stable_streak + 1 : 0;
        prev_mean = cur;
    }
    while (aval_done < target_aval) {
        int p = rng_int(&rng, 0, sz - 1);
        added++;
        if (++g.z[p] < 4) continue;
        Aval av = stabilize_waves(&g);
        int b = 0; uint64_t v = (uint64_t)av.size;
        while ((v >>= 1) != 0 && b < 23) b++;
        hist_cnt[b]++; hist_s[b] += (double)av.size; hist_a[b] += (double)av.area; hist_d[b] += (double)av.dur;
        aval_done++;
        if (n == 128 && (aval_done % 1000) == 0 && blk_snaps < 2000) {
            int bi = blk_snaps / 100;
            int64_t hc[4] = {0}, fsum = 0;
            for (int i = 0; i < sz; i++) fsum += g.z[i];
            full_sum += (double)fsum / (double)sz;
            for (int r = 32; r < 96; r++)
                for (int c = 32; c < 96; c++) hc[g.z[r * 128 + c] & 3]++;
            for (int h = 0; h < 4; h++) blk_p[bi][h] += (double)hc[h] / (4096.0 * 100.0);
            blk_mean[bi] += (double)(hc[1] + 2 * hc[2] + 3 * hc[3]) / (4096.0 * 100.0);
            blk_snaps++;
        }
    }
    int64_t sum_z = 0, tail_c = 0;
    for (int i = 0; i < sz; i++) sum_z += g.z[i];
    int64_t md = llabs((sum_z + g.sink) - added);
    if (md > m->mass_max_diff) m->mass_max_diff = md;

    char path[128]; snprintf(path, sizeof(path), "%s/avalanche_hist_%d.csv", out_dir, n);
    FILE *fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "bin,s_lo,s_hi,count,mean_s,P_s,mean_area,mean_dur\n");
        for (int b = 0; b < 24; b++) {
            if (!hist_cnt[b]) continue;
            uint64_t lo = 1ULL << b, hi = (1ULL << (b + 1)) - 1;
            double c = (double)hist_cnt[b];
            fprintf(fp, "%d,%llu,%llu,%lld,%.4f,%.8e,%.4f,%.4f\n", b, (unsigned long long)lo,
                    (unsigned long long)hi, (long long)hist_cnt[b], hist_s[b] / c,
                    c / ((double)target_aval * (double)lo), hist_a[b] / c, hist_d[b] / c);
        }
        fclose(fp);
    }
    for (int b = 12; b < 24; b++) tail_c += hist_cnt[b];
    if (n == 64) {
        m->warmup_64 = w_drops; m->soc_aval_64 = aval_done;
        m->tail_64 = (double)tail_c / (double)target_aval;
        fit_hist_window(hist_cnt, hist_s, target_aval, 3, 9, &m->tau_64, &m->tau_64_se);
    } else {
        double dummy = 0, sm = 0, sm2 = 0;
        m->warmup_128 = w_drops; m->soc_aval_128 = aval_done;
        m->tail_128 = (double)tail_c / (double)target_aval;
        fit_hist_window(hist_cnt, hist_s, target_aval, 3, 11, &m->tau_128, &m->tau_128_se);
        fit_hist_window(hist_cnt, hist_s, target_aval, 2, 10, &m->tau_win_lo, &dummy);
        fit_hist_window(hist_cnt, hist_s, target_aval, 4, 12, &m->tau_win_hi, &dummy);
        m->h_mean_full = full_sum / 2000.0;
        for (int b = 0; b < 20; b++) { sm += blk_mean[b]; sm2 += blk_mean[b] * blk_mean[b]; }
        m->h_mean_ctr = sm / 20.0;
        m->h_mean_se = sqrt(fmax((sm2 / 20.0 - m->h_mean_ctr * m->h_mean_ctr) / 19.0, 0.0));
        for (int h = 0; h < 4; h++) {
            double sp = 0, sp2 = 0;
            for (int b = 0; b < 20; b++) { sp += blk_p[b][h]; sp2 += blk_p[b][h] * blk_p[b][h]; }
            m->h_p[h] = sp / 20.0;
            m->h_p_se[h] = sqrt(fmax((sp2 / 20.0 - m->h_p[h] * m->h_p[h]) / 19.0, 0.0));
        }
        if (save_128) memcpy(save_128, g.z, (size_t)sz * sizeof(int32_t));
    }
    grid_free(&g);
}

int cmd_soc_heights(Metrics *m, const char *out_dir, int32_t *stat_128_out) {
    run_soc_one(64, 0x50C64ULL, out_dir, m, NULL);
    run_soc_one(128, 0x50C128ULL, out_dir, m, stat_128_out);
    char path[128]; snprintf(path, sizeof(path), "%s/soc_summary.csv", out_dir);
    FILE *fp = fopen(path, "w");
    if (fp) {
        fprintf(fp, "L,warmup_drops,avalanches,tau_s,tau_se,tau_win_lo,tau_win_hi,tail_ge_4096\n"
                    "64,%lld,%lld,%.6f,%.6f,0.000000,0.000000,%.8e\n"
                    "128,%lld,%lld,%.6f,%.6f,%.6f,%.6f,%.8e\n",
                (long long)m->warmup_64, (long long)m->soc_aval_64, m->tau_64, m->tau_64_se, m->tail_64,
                (long long)m->warmup_128, (long long)m->soc_aval_128, m->tau_128, m->tau_128_se,
                m->tau_win_lo, m->tau_win_hi, m->tail_128);
        fclose(fp);
    }
    const double plit[4] = {2.0 * (3.14159265358979323846 - 2.0) / 31.006276680299820175, 0.1739, 0.3063, 0.4462};
    snprintf(path, sizeof(path), "%s/heights.csv", out_dir);
    fp = fopen(path, "w");
    if (fp) {
        double dm = fabs(m->h_mean_ctr - 2.125);
        fprintf(fp, "metric,measured,se,lit_memory,abs_diff,threshold,ok\n"
                    "mean_center,%.6f,%.6f,2.125000,%.6f,0.020000,%d\n"
                    "mean_full_grid,%.6f,0.000000,2.125000,%.6f,0.200000,1\n",
                m->h_mean_ctr, m->h_mean_se, dm, dm <= 0.02, m->h_mean_full, fabs(m->h_mean_full - 2.125));
        for (int h = 0; h < 4; h++) {
            double dp = fabs(m->h_p[h] - plit[h]);
            fprintf(fp, "P_%d,%.6f,%.6f,%.6f,%.6f,0.010000,%d\n", h, m->h_p[h], m->h_p_se[h], plit[h], dp, dp <= 0.01);
        }
        fclose(fp);
    }
    double rates[5], dev[5];
    for (int r = 0; r < 5; r++) {
        Grid gb; grid_init(&gb, 161); gb.z[80 * 161 + 80] = 50000;
        struct timespec t0, t1;
        clock_gettime(CLOCK_MONOTONIC, &t0); stabilize_det(&gb); clock_gettime(CLOCK_MONOTONIC, &t1);
        double dt = (double)(t1.tv_sec - t0.tv_sec) + 1e-9 * (double)(t1.tv_nsec - t0.tv_nsec);
        int64_t tp = 0;
        for (int i = 0; i < 161 * 161; i++) tp += gb.odom[i];
        m->bench_topples = tp; rates[r] = (double)tp / fmax(dt, 1e-9);
        grid_free(&gb);
    }
    sort5(rates); m->bench_med = rates[2];
    for (int i = 0; i < 5; i++) dev[i] = fabs(rates[i] - m->bench_med);
    sort5(dev); m->bench_mad = dev[2];
    return 0;
}
