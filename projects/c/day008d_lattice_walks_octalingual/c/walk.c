#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define MASTER_SEED 20261126ULL
#define W_WALKS 200000
#define W_RUIN 20000
#define BENCH_STEPS 10000000
#define T_EXACT_MAX 5000

typedef struct { uint64_t s[4]; } Rng;

static inline uint64_t sm64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rng_seed(Rng *r, uint64_t seed) {
    uint64_t st = seed;
    for (int i = 0; i < 4; ++i) { st += 0x9E3779B97F4A7C15ULL; r->s[i] = sm64(st); }
}
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rng_u64(Rng *r) {
    uint64_t res = rotl(r->s[0] + r->s[3], 23) + r->s[0];
    uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1]; r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t; r->s[3] = rotl(r->s[3], 45);
    return res;
}
static inline double rng_u01(Rng *r) { return (double)(rng_u64(r) >> 11) * 0x1.0p-53; }
static inline int rng_int(Rng *r, int n) {
    int k = (int)floor(rng_u01(r) * (double)n);
    return (k >= n) ? (n - 1) : k;
}
static uint64_t task_seed(int d, int st, int idx) {
    uint64_t h = 0x517CC1B727220A95ULL;
    h = sm64(h ^ (uint64_t)d); h = sm64(h ^ (uint64_t)st); h = sm64(h ^ (uint64_t)idx);
    return sm64(MASTER_SEED ^ h);
}
static uint64_t fnv1a_i32(const int32_t *v, int n) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (int i = 0; i < n; ++i) {
        uint32_t u = (uint32_t)v[i];
        for (int b = 0; b < 4; ++b) { h ^= (u >> (8 * b)) & 0xFFU; h *= 0x100000001B3ULL; }
    }
    return h;
}
static double now_s(void) {
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}
static void ensure_dirs(void) {
    mkdir("out", 0755); mkdir("out/c", 0755); mkdir("out/c/figures", 0755);
}

static double g_msd[3][13], g_msd_se[3][13], g_ret_sim[3][4], g_ret_se[3][4], g_ret_ex[3][4];
static double g_ruin_p[6], g_ruin_pe[6], g_ruin_t[6], g_ruin_te[6];
static double g_s03_1024[3][W_WALKS];
static uint8_t g_s04_10k[3][W_WALKS];
static double g_stage_sec[12];
static const int TS_MSD[13] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096};
static const int TS_RET[4] = {10,100,1000,10000};

static int cmp_dbl(const void *a, const void *b) {
    double da = *(const double *)a, db = *(const double *)b;
    return (da > db) - (da < db);
}

static void run_s01(void) {
    Rng r; rng_seed(&r, MASTER_SEED);
    uint64_t x0 = rng_u64(&r), x1 = rng_u64(&r), x2 = rng_u64(&r);
    int ok = (x0 == 0xA9920E0E7DAB3EF0ULL && x1 == 0x033B9CEB4B3BA7A0ULL && x2 == 0xE7473A6A39AFCEA6ULL);
    double w[5] = {1.0, 1.0 / 3.0, 0, 0, 0};
    for (int n = 1; n < 4; ++n) w[n+1] = ((10.0*n*n + 10.0*n + 3.0)*w[n] - (double)(n*n)*w[n-1]) / (9.0*(n+1)*(n+1));
    int enum_ok = (fabs(w[2] - 15.0/81.0) < 1e-14 && fabs(w[3] - 93.0/729.0) < 1e-14 && fabs(w[4] - 639.0/6561.0) < 1e-14);
    FILE *fp = fopen("out/c/s01_selftest.txt", "w");
    fprintf(fp, "ok,x0,x1,x2,enum_ok\n%d,%016llx,%016llx,%016llx,%d\n",
            ok && enum_ok, (unsigned long long)x0, (unsigned long long)x1, (unsigned long long)x2, enum_ok);
    fclose(fp);
    if (!(ok && enum_ok)) { fprintf(stderr, "C selftest FAILED\n"); exit(1); }
}

static void run_s02(void) {
    FILE *fp = fopen("out/c/s02_parity.csv", "w");
    fprintf(fp, "d,walk,x0,x1,x2,fnv1a_hex\n");
    for (int d = 1; d <= 3; ++d) {
        int32_t coords[24], pos_all[8][3];
        int ci = 0;
        for (int w = 0; w < 8; ++w) {
            Rng r; rng_seed(&r, task_seed(d, 2, w));
            int32_t p[3] = {0, 0, 0};
            for (int s = 0; s < 200; ++s) {
                int k = rng_int(&r, 2 * d);
                p[k >> 1] += (k & 1) ? -1 : 1;
            }
            for (int a = 0; a < 3; ++a) pos_all[w][a] = p[a];
            for (int a = 0; a < d; ++a) coords[ci++] = p[a];
        }
        uint64_t h = fnv1a_i32(coords, 8 * d);
        for (int w = 0; w < 8; ++w)
            fprintf(fp, "%d,%d,%d,%d,%d,%016llx\n", d, w, pos_all[w][0], pos_all[w][1], pos_all[w][2], (unsigned long long)h);
    }
    fclose(fp);
}

static void run_s03(void) {
    FILE *fp = fopen("out/c/s03_msd.csv", "w");
    fprintf(fp, "d,t,msd,se,exact\n");
    for (int d = 1; d <= 3; ++d) {
        double sum1[13] = {0}, sum2[13] = {0};
        for (int w = 0; w < W_WALKS; ++w) {
            Rng r; rng_seed(&r, task_seed(d, 3, w));
            int p[3] = {0, 0, 0}, ti = 0;
            for (int t = 1; t <= 4096; ++t) {
                int k = rng_int(&r, 2 * d);
                p[k >> 1] += (k & 1) ? -1 : 1;
                if (t == TS_MSD[ti]) {
                    double r2 = (double)(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                    sum1[ti] += r2; sum2[ti] += r2 * r2;
                    if (ti == 10) g_s03_1024[d - 1][w] = r2;
                    ++ti;
                }
            }
        }
        for (int ti = 0; ti < 13; ++ti) {
            double m = sum1[ti] / W_WALKS;
            double v = fmax(0.0, sum2[ti] / W_WALKS - m * m);
            double se = sqrt(v / W_WALKS);
            g_msd[d - 1][ti] = m; g_msd_se[d - 1][ti] = se;
            fprintf(fp, "%d,%d,%.8f,%.8f,%d\n", d, TS_MSD[ti], m, se, TS_MSD[ti]);
        }
    }
    fclose(fp);
}

static void run_s04(void) {
    FILE *fp = fopen("out/c/s04_return_sim.csv", "w");
    fprintf(fp, "d,T,p_hat,se\n");
    for (int d = 1; d <= 3; ++d) {
        int hits[4] = {0, 0, 0, 0};
        for (int w = 0; w < W_WALKS; ++w) {
            Rng r; rng_seed(&r, task_seed(d, 4, w));
            int p[3] = {0, 0, 0}, ret_t = 0;
            for (int t = 1; t <= 10000; ++t) {
                int k = rng_int(&r, 2 * d);
                p[k >> 1] += (k & 1) ? -1 : 1;
                if ((p[0] | p[1] | p[2]) == 0) { ret_t = t; break; }
            }
            g_s04_10k[d - 1][w] = (ret_t > 0) ? 1 : 0;
            if (ret_t > 0) {
                for (int ti = 0; ti < 4; ++ti) if (ret_t <= TS_RET[ti]) ++hits[ti];
            }
        }
        for (int ti = 0; ti < 4; ++ti) {
            double ph = (double)hits[ti] / W_WALKS;
            double se = sqrt(fmax(0.0, ph * (1.0 - ph)) / W_WALKS);
            g_ret_sim[d - 1][ti] = ph; g_ret_se[d - 1][ti] = se;
            fprintf(fp, "%d,%d,%.8f,%.8f\n", d, TS_RET[ti], ph, se);
        }
    }
    fclose(fp);
}

static void run_s05(void) {
    static double u1[T_EXACT_MAX + 1], u2[T_EXACT_MAX + 1], w3[T_EXACT_MAX + 1], u3[T_EXACT_MAX + 1], f[T_EXACT_MAX + 1];
    u1[0] = u2[0] = w3[0] = u3[0] = 1.0; w3[1] = 1.0 / 3.0;
    for (int n = 1; n < T_EXACT_MAX; ++n)
        w3[n + 1] = ((10.0*n*n + 10.0*n + 3.0)*w3[n] - (double)n*n*w3[n - 1]) / (9.0*(n + 1.0)*(n + 1.0));
    for (int n = 1; n <= T_EXACT_MAX; ++n) {
        u1[n] = u1[n - 1] * (2.0*n - 1.0) / (2.0*n);
        u2[n] = u1[n] * u1[n];
        u3[n] = u1[n] * w3[n];
    }
    double *us[3] = {u1, u2, u3};
    FILE *fp = fopen("out/c/s05_return_exact.csv", "w");
    fprintf(fp, "d,T,F_exact\n");
    for (int d = 1; d <= 3; ++d) {
        double *u = us[d - 1], acc = 0.0;
        f[0] = 0.0;
        for (int n = 1; n <= T_EXACT_MAX; ++n) {
            double s = u[n];
            for (int k = 1; k < n; ++k) s -= f[k] * u[n - k];
            f[n] = s; acc += s;
            for (int ti = 0; ti < 4; ++ti)
                if (2 * n == TS_RET[ti]) { g_ret_ex[d - 1][ti] = acc; fprintf(fp, "%d,%d,%.10f\n", d, TS_RET[ti], acc); }
        }
    }
    fclose(fp);
}

static void run_s06(void) {
    FILE *fp = fopen("out/c/s06_return_compare.csv", "w");
    fprintf(fp, "d,T,p_hat,se,F_exact,z_score,p_inf_est\n");
    double p3_inf = g_ret_sim[2][3] + (g_ret_sim[2][3] - g_ret_sim[2][2]) / (sqrt(10.0) - 1.0);
    for (int d = 1; d <= 3; ++d)
        for (int ti = 0; ti < 4; ++ti) {
            double z = (g_ret_se[d - 1][ti] > 1e-12) ? (g_ret_sim[d - 1][ti] - g_ret_ex[d - 1][ti]) / g_ret_se[d - 1][ti] : 0.0;
            fprintf(fp, "%d,%d,%.8f,%.8f,%.10f,%.4f,%.8f\n", d, TS_RET[ti], g_ret_sim[d - 1][ti], g_ret_se[d - 1][ti],
                    g_ret_ex[d - 1][ti], z, (d == 3 && ti == 3) ? p3_inf : 0.0);
        }
    fclose(fp);
}

static void run_s07(void) {
    const int Ns[6] = {10, 10, 10, 50, 50, 50}, xs[6] = {2, 5, 8, 10, 25, 40};
    FILE *fp = fopen("out/c/s07_ruin.csv", "w");
    fprintf(fp, "N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z\n");
    for (int ci = 0; ci < 6; ++ci) {
        int N = Ns[ci], x0 = xs[ci], wins = 0;
        double st1 = 0.0, st2 = 0.0;
        for (int w = 0; w < W_RUIN; ++w) {
            Rng r; rng_seed(&r, task_seed(1, 7, ci * 100000 + w));
            int x = x0, steps = 0;
            while (x > 0 && x < N) { x += (rng_int(&r, 2) == 0) ? 1 : -1; ++steps; }
            if (x == N) ++wins;
            st1 += steps; st2 += (double)steps * steps;
        }
        double ph = (double)wins / W_RUIN, pse = sqrt(fmax(0.0, ph * (1.0 - ph)) / W_RUIN), pe = (double)x0 / N;
        double th = st1 / W_RUIN, tse = sqrt(fmax(0.0, st2 / W_RUIN - th * th) / W_RUIN), te = (double)x0 * (N - x0);
        g_ruin_p[ci] = ph; g_ruin_pe[ci] = pe; g_ruin_t[ci] = th; g_ruin_te[ci] = te;
        fprintf(fp, "%d,%d,%.8f,%.8f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                N, x0, ph, pse, pe, (ph - pe) / pse, th, tse, te, (th - te) / tse);
    }
    fclose(fp);
}

static void run_s08(void) {
    FILE *fp = fopen("out/c/s08_endpoint.csv", "w");
    fprintf(fp, "d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count\n");
    int hist[40] = {0};
    for (int d = 1; d <= 3; ++d) {
        double m1 = 0, m2 = 0, m4 = 0;
        int n_ep = 20000;
        for (int w = 0; w < n_ep; ++w) {
            Rng r; rng_seed(&r, task_seed(d, 8, w));
            int p[3] = {0, 0, 0};
            for (int s = 0; s < 1024; ++s) { int k = rng_int(&r, 2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
            double x = (double)p[0];
            m1 += x; m2 += x * x; m4 += x * x * x * x;
            if (d == 1) {
                int b = (int)floor((x + 128.0) / 6.4);
                if (b < 0) b = 0;
                if (b > 39) b = 39;
                ++hist[b];
            }
        }
        m1 /= n_ep; double v = m2 / n_ep - m1 * m1;
        double kurt = (m4 / n_ep) / (v * v);
        fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,-1,0,0,%d\n", d, v, 1024.0 / d, kurt, 3.0 + (d - 3.0) / 1024.0, n_ep);
    }
    for (int b = 0; b < 40; ++b)
        fprintf(fp, "1,0,0,0,0,%d,%.2f,%.2f,%d\n", b, -128.0 + b * 6.4, -128.0 + (b + 1) * 6.4, hist[b]);
    fclose(fp);
}

static void run_s09(void) {
    FILE *fp = fopen("out/c/s09_ci.csv", "w");
    fprintf(fp, "metric,d,estimate,boot_se,ci_lo,ci_hi\n");
    double boot[200];
    for (int d = 1; d <= 3; ++d) {
        Rng r; rng_seed(&r, task_seed(d, 9, 0));
        for (int b = 0; b < 200; ++b) {
            double s = 0;
            for (int i = 0; i < W_WALKS; ++i) s += g_s03_1024[d - 1][rng_int(&r, W_WALKS)];
            boot[b] = s / W_WALKS;
        }
        qsort(boot, 200, sizeof(double), cmp_dbl);
        double bm = 0, bv = 0;
        for (int b = 0; b < 200; ++b) bm += boot[b];
        bm /= 200.0;
        for (int b = 0; b < 200; ++b) bv += (boot[b] - bm) * (boot[b] - bm);
        fprintf(fp, "msd_1024,%d,%.6f,%.6f,%.6f,%.6f\n", d, g_msd[d - 1][10], sqrt(bv / 199.0), boot[5], boot[194]);

        rng_seed(&r, task_seed(d, 9, 1));
        for (int b = 0; b < 200; ++b) {
            int s = 0;
            for (int i = 0; i < W_WALKS; ++i) s += g_s04_10k[d - 1][rng_int(&r, W_WALKS)];
            boot[b] = (double)s / W_WALKS;
        }
        qsort(boot, 200, sizeof(double), cmp_dbl);
        bm = bv = 0;
        for (int b = 0; b < 200; ++b) bm += boot[b];
        bm /= 200.0;
        for (int b = 0; b < 200; ++b) bv += (boot[b] - bm) * (boot[b] - bm);
        fprintf(fp, "return_10000,%d,%.6f,%.6f,%.6f,%.6f\n", d, g_ret_sim[d - 1][3], sqrt(bv / 199.0), boot[5], boot[194]);
    }
    fclose(fp);
}

static void run_s10(void) {
    FILE *fp = fopen("out/c/s10_fit.csv", "w");
    fprintf(fp, "d,slope,slope_se,intercept\n");
    for (int d = 1; d <= 3; ++d) {
        double mx = 0, my = 0;
        for (int i = 0; i < 13; ++i) { mx += log((double)TS_MSD[i]); my += log(g_msd[d - 1][i]); }
        mx /= 13.0; my /= 13.0;
        double num = 0, den = 0;
        for (int i = 0; i < 13; ++i) {
            double dx = log((double)TS_MSD[i]) - mx;
            num += dx * (log(g_msd[d - 1][i]) - my); den += dx * dx;
        }
        double sl = num / den, ic = my - sl * mx, rss = 0;
        for (int i = 0; i < 13; ++i) { double e = log(g_msd[d - 1][i]) - (sl * log((double)TS_MSD[i]) + ic); rss += e * e; }
        fprintf(fp, "%d,%.6f,%.6f,%.6f\n", d, sl, sqrt((rss / 11.0) / den), ic);
    }
    fclose(fp);
}

static void run_s11(void) {
    FILE *fp = fopen("out/c/s11_bench.csv", "w");
    fprintf(fp, "d,steps_per_rep,reps,ns_med,ns_mad\n");
    volatile int64_t sink = 0;
    for (int d = 1; d <= 3; ++d) {
        double vals[7], dev[7];
        for (int rep = 0; rep < 7; ++rep) {
            Rng r; rng_seed(&r, task_seed(d, 11, rep));
            int p[3] = {0, 0, 0};
            double t0 = now_s();
            for (int s = 0; s < BENCH_STEPS; ++s) { int k = rng_int(&r, 2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
            vals[rep] = (now_s() - t0) * 1e9 / (double)BENCH_STEPS;
            sink += p[0] + p[1] + p[2];
        }
        qsort(vals, 7, sizeof(double), cmp_dbl);
        double med = vals[3];
        for (int rep = 0; rep < 7; ++rep) dev[rep] = fabs(vals[rep] - med);
        qsort(dev, 7, sizeof(double), cmp_dbl);
        fprintf(fp, "%d,%d,7,%.4f,%.4f\n", d, BENCH_STEPS, med, dev[3]);
    }
    (void)sink;
    fclose(fp);
}

static void run_s12(void) {
    const char *cols[3] = {"#1f77b4", "#ff7f0e", "#2ca02c"};
    FILE *f1 = fopen("out/c/figures/msd.svg", "w");
    fprintf(f1, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for (int d = 0; d < 3; ++d) {
        fprintf(f1, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.6\" points=\"", cols[d]);
        for (int i = 0; i < 13; ++i)
            fprintf(f1, "%.1f,%.1f ", 60.0 + 500.0 * log10((double)TS_MSD[i]) / log10(4096.0), 310.0 - 280.0 * log10(fmax(1.0, g_msd[d][i])) / log10(4096.0));
        fprintf(f1, "\"/>");
    }
    fprintf(f1, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 MSD Log-Log (y=t)</text></svg>\n");
    fclose(f1);

    FILE *f2 = fopen("out/c/figures/return_prob.svg", "w");
    fprintf(f2, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int d = 0; d < 3; ++d) {
        fprintf(f2, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.6\" points=\"", cols[d]);
        for (int i = 0; i < 4; ++i) fprintf(f2, "%.1f,%.1f ", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * g_ret_ex[d][i]);
        fprintf(f2, "\"/>");
        for (int i = 0; i < 4; ++i) fprintf(f2, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4\" fill=\"%s\"/>", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * g_ret_sim[d][i], cols[d]);
    }
    fprintf(f2, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 Return Prob F_d(T) Exact vs Sim</text></svg>\n");
    fclose(f2);

    FILE *f3 = fopen("out/c/figures/ruin.svg", "w");
    fprintf(f3, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for (int i = 0; i < 6; ++i)
        fprintf(f3, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"5\" fill=\"#1f77b4\"/>", 60.0 + 500.0 * g_ruin_pe[i], 310.0 - 280.0 * g_ruin_p[i]);
    fprintf(f3, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 Gambler's Ruin Exact vs Sim</text></svg>\n");
    fclose(f3);

    double tot = 0; for (int i = 0; i < 12; ++i) tot += g_stage_sec[i];
    FILE *fj = fopen("out/c/s12_summary.json", "w");
    fprintf(fj, "{\"lang\":\"c\",\"W\":%d,\"ruin_W\":%d,\"bench_steps\":%d,\"toolchain\":\"Apple clang 16.0.0 (-std=c11)\",\"total_sec\":%.3f}\n",
            W_WALKS, W_RUIN, BENCH_STEPS, tot);
    fclose(fj);
}

int main(int argc, char **argv) {
    ensure_dirs();
    const char *cmd = (argc > 1) ? argv[1] : "all";
    typedef void (*Fn)(void);
    Fn fns[12] = {run_s01,run_s02,run_s03,run_s04,run_s05,run_s06,run_s07,run_s08,run_s09,run_s10,run_s11,run_s12};
    const char *names[12] = {"s01","s02","s03","s04","s05","s06","s07","s08","s09","s10","s11","s12"};
    if (strcmp(cmd, "all") == 0) {
        for (int i = 0; i < 12; ++i) {
            printf("[c] tahap %d/12 (%s)\n", i + 1, names[i]);
            double t0 = now_s(); fns[i](); g_stage_sec[i] = now_s() - t0;
        }
        return 0;
    }
    for (int i = 0; i < 12; ++i) if (strcmp(cmd, names[i]) == 0) { fns[i](); return 0; }
    fprintf(stderr, "Unknown cmd: %s\n", cmd);
    return 1;
}
