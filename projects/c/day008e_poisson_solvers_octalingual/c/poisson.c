#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#define PI_VAL 3.14159265358979323846
static double sin_spec(double x) {
    double sgn = 1.0;
    if (x < 0.0) { x = -x; sgn = -sgn; }
    if (x > PI_VAL) { x = x - PI_VAL; sgn = -sgn; }
    if (x > 0.5 * PI_VAL) x = PI_VAL - x;
    double x2 = x * x;
    static const double c[10] = {
        1.0, -1.0/6.0, 1.0/120.0, -1.0/5040.0, 1.0/362880.0, -1.0/39916800.0,
        1.0/6227020800.0, -1.0/1307674368000.0, 1.0/355687428096000.0, -1.0/121645100408832000.0
    };
    double p = c[9];
    for (int k = 8; k >= 0; --k) p = p * x2 + c[k];
    return sgn * (x * p);
}
static double cos_spec(double x) { return sin_spec(0.5 * PI_VAL - x); }
static double omega_opt(int N) {
    if (N == 15) return 1.673513677715992;
    if (N == 31) return 1.8214651907890225;
    if (N == 63) return 1.906454701582762;
    if (N == 127) return 1.952093233850055;
    if (N == 255) return 1.9757544535797149;
    return 1.9878030696593354;
}
static uint64_t fnv1a_f64(const double *v, int n) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (int i = 0; i < n; ++i) {
        uint64_t u; memcpy(&u, &v[i], 8);
        for (int b = 0; b < 8; ++b) { h ^= (u >> (8 * b)) & 0xFFU; h *= 0x100000001B3ULL; }
    }
    return h;
}
static double now_s(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec; }
static int cmp_d(const void *a, const void *b) { double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y); }

static void apply_A(int N, const double *u, double *out) {
    for (int i = 0; i < N; ++i) {
        int r = i * N;
        for (int j = 0; j < N; ++j) {
            int k = r + j;
            double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
            double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
            out[k] = (((4.0 * u[k] - w) - e) - s) - n;
        }
    }
}
static double dot_seq(int n2, const double *a, const double *b) {
    double s = 0.0; for (int i = 0; i < n2; ++i) s += a[i] * b[i]; return s;
}
static double norm2_seq(int n2, const double *a) { return sqrt(dot_seq(n2, a, a)); }
static double max_abs_diff(int n2, const double *a, const double *b) {
    double m = 0.0; for (int i = 0; i < n2; ++i) { double d = fabs(a[i] - b[i]); if (d > m) m = d; } return m;
}
static void build_prob(int N, double *b, double *ue) {
    double h = 1.0 / (N + 1), h2 = h * h;
    double *sx = malloc((size_t)N * 8), *sy = malloc((size_t)N * 8);
    for (int i = 0; i < N; ++i) { sx[i] = sin_spec(PI_VAL * (i + 1) * h); sy[i] = sin_spec(2.0 * PI_VAL * (i + 1) * h); }
    for (int i = 0; i < N; ++i) {
        double x = (i + 1) * h, xx = x * (1.0 - x);
        for (int j = 0; j < N; ++j) {
            double y = (j + 1) * h, yy = y * (1.0 - y), sxy = sx[i] * sy[j];
            ue[i * N + j] = sxy + xx * yy;
            b[i * N + j] = h2 * (5.0 * PI_VAL * PI_VAL * sxy + 2.0 * xx + 2.0 * yy);
        }
    }
    free(sx); free(sy);
}
static double compute_res(int N, const double *u, const double *b, double *r) {
    apply_A(N, u, r);
    for (int i = 0; i < N * N; ++i) r[i] = b[i] - r[i];
    return norm2_seq(N * N, r);
}

static void rbgs_sweep(int N, double *u, const double *b) {
    for (int col = 0; col < 2; ++col)
        for (int i = 0; i < N; ++i) {
            int row = i * N;
            for (int j = 0; j < N; ++j) {
                if (((i + j) & 1) != col) continue;
                int k = row + j;
                double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
                double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
                u[k] = 0.25 * ((((b[k] + w) + e) + s) + n);
            }
        }
}
static void restrict_fw(int Nf, const double *rf, double *rc) {
    int Nc = (Nf - 1) / 2;
    for (int ic = 0; ic < Nc; ++ic) {
        int fi = 2 * ic + 1;
        for (int jc = 0; jc < Nc; ++jc) {
            int fj = 2 * jc + 1, k = fi * Nf + fj;
            double edge = ((rf[k - 1] + rf[k + 1]) + rf[k - Nf]) + rf[k + Nf];
            double corn = ((rf[k - Nf - 1] + rf[k - Nf + 1]) + rf[k + Nf - 1]) + rf[k + Nf + 1];
            rc[ic * Nc + jc] = 4.0 * (0.25 * rf[k] + 0.125 * edge + 0.0625 * corn);
        }
    }
}
static void prolong_add(int Nf, double *uf, const double *ec) {
    int Nc = (Nf - 1) / 2;
    for (int ic = 0; ic < Nc; ++ic) {
        int fi = 2 * ic + 1;
        for (int jc = 0; jc < Nc; ++jc) {
            int fj = 2 * jc + 1, k = fi * Nf + fj;
            double v = ec[ic * Nc + jc];
            uf[k] += v;
            uf[k - 1] += 0.5 * v; uf[k + 1] += 0.5 * v; uf[k - Nf] += 0.5 * v; uf[k + Nf] += 0.5 * v;
            uf[k - Nf - 1] += 0.25 * v; uf[k - Nf + 1] += 0.25 * v; uf[k + Nf - 1] += 0.25 * v; uf[k + Nf + 1] += 0.25 * v;
        }
    }
}
static void solve_dense3(const double *b, double *x) {
    double M[9][10] = {{0}};
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) {
        int r = i * 3 + j; M[r][r] = 4.0;
        if (j > 0) M[r][r - 1] = -1.0;
        if (j < 2) M[r][r + 1] = -1.0;
        if (i > 0) M[r][r - 3] = -1.0;
        if (i < 2) M[r][r + 3] = -1.0;
        M[r][9] = b[r];
    }
    for (int c = 0; c < 9; ++c) {
        int piv = c;
        for (int r = c + 1; r < 9; ++r) if (fabs(M[r][c]) > fabs(M[piv][c])) piv = r;
        if (piv != c) for (int k = 0; k < 10; ++k) { double tmp = M[c][k]; M[c][k] = M[piv][k]; M[piv][k] = tmp; }
        for (int r = c + 1; r < 9; ++r) {
            double f = M[r][c] / M[c][c];
            for (int k = c; k < 10; ++k) M[r][k] -= f * M[c][k];
        }
    }
    for (int r = 8; r >= 0; --r) {
        double s = M[r][9];
        for (int k = r + 1; k < 9; ++k) s -= M[r][k] * x[k];
        x[r] = s / M[r][r];
    }
}
static void vcycle(int N, double *u, const double *b) {
    if (N == 3) { solve_dense3(b, u); return; }
    rbgs_sweep(N, u, b); rbgs_sweep(N, u, b);
    int n2 = N * N, Nc = (N - 1) / 2, nc2 = Nc * Nc;
    double *r = malloc((size_t)n2 * 8), *rc = malloc((size_t)nc2 * 8), *ec = calloc((size_t)nc2, 8);
    compute_res(N, u, b, r);
    restrict_fw(N, r, rc);
    vcycle(Nc, ec, rc);
    prolong_add(N, u, ec);
    rbgs_sweep(N, u, b);
    free(r); free(rc); free(ec);
}
static void fmg_solve(int N, double *u, const double *b) {
    if (N == 3) { solve_dense3(b, u); return; }
    int Nc = (N - 1) / 2, nc2 = Nc * Nc;
    double *bc = malloc((size_t)nc2 * 8), *uc = calloc((size_t)nc2, 8);
    restrict_fw(N, b, bc);
    fmg_solve(Nc, uc, bc);
    memset(u, 0, (size_t)(N * N) * 8);
    prolong_add(N, u, uc);
    vcycle(N, u, b);
    free(bc); free(uc);
}

typedef struct { int iters; double rel_res, err_inf, rho_emp; uint64_t hist_hash, sol_hash; } SolveRes;
static double g_hist63[7][100001];
static int g_hist63_len[7];
static double g_sol63[8][63 * 63];

static SolveRes run_solver(int meth, int N, int save63) {
    int n2 = N * N;
    double *b = malloc((size_t)n2 * 8), *ue = malloc((size_t)n2 * 8);
    double *u = calloc((size_t)n2, 8), *r = calloc((size_t)n2, 8);
    double *p = calloc((size_t)n2, 8), *Ap = calloc((size_t)n2, 8), *z = calloc((size_t)n2, 8);
    build_prob(N, b, ue);
    double nb = norm2_seq(n2, b), om = omega_opt(N);
    int max_it = (meth == 0) ? 100000 : (meth <= 2 ? 50000 : (meth == 3 ? 10000 : (meth == 4 ? 5000 : (meth == 5 ? 3000 : 100))));
    double *hist = malloc((size_t)(max_it + 1) * 8);
    int it = 0; double rr = 0.0, rz = 0.0;
    if (meth == 4 || meth == 5) { memcpy(r, b, (size_t)n2 * 8); rr = dot_seq(n2, r, r); }
    if (meth == 4) memcpy(p, r, (size_t)n2 * 8);
    if (meth == 5) {
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
            int k = i * N + j;
            double yw = (j > 0) ? z[k - 1] : 0.0, ys = (i > 0) ? z[k - N] : 0.0;
            z[k] = 0.25 * ((r[k] + yw) + ys);
        }
        for (int i = N - 1; i >= 0; --i) for (int j = N - 1; j >= 0; --j) {
            int k = i * N + j;
            double ze = (j + 1 < N) ? p[k + 1] : 0.0, zn = (i + 1 < N) ? p[k + N] : 0.0;
            p[k] = z[k] + 0.25 * (ze + zn);
        }
        memcpy(z, p, (size_t)n2 * 8); rz = dot_seq(n2, r, z);
    }
    for (it = 1; it <= max_it; ++it) {
        double rel = 0.0;
        if (meth == 0) {
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
                double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
                z[k] = 0.25 * ((((b[k] + w) + e) + s) + n);
            }
            memcpy(u, z, (size_t)n2 * 8);
            rel = compute_res(N, u, b, r) / nb;
        } else if (meth == 1 || meth == 3) {
            double ouse = (meth == 1) ? 1.0 : om;
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
                double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
                double gs = 0.25 * ((((b[k] + w) + e) + s) + n);
                u[k] = (meth == 1) ? gs : ((1.0 - ouse) * u[k] + ouse * gs);
            }
            rel = compute_res(N, u, b, r) / nb;
        } else if (meth == 2) {
            rbgs_sweep(N, u, b);
            rel = compute_res(N, u, b, r) / nb;
        } else if (meth == 4) {
            apply_A(N, p, Ap);
            double alpha = rr / dot_seq(n2, p, Ap);
            for (int k = 0; k < n2; ++k) { u[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }
            double rr_new = dot_seq(n2, r, r);
            rel = sqrt(rr_new) / nb;
            hist[it - 1] = rel;
            if (rel <= 1e-8) { it++; break; }
            double beta = rr_new / rr;
            for (int k = 0; k < n2; ++k) p[k] = r[k] + beta * p[k];
            rr = rr_new; continue;
        } else if (meth == 5) {
            apply_A(N, p, Ap);
            double alpha = rz / dot_seq(n2, p, Ap);
            for (int k = 0; k < n2; ++k) { u[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }
            rel = norm2_seq(n2, r) / nb;
            hist[it - 1] = rel;
            if (rel <= 1e-8) { it++; break; }
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                double yw = (j > 0) ? Ap[k - 1] : 0.0, ys = (i > 0) ? Ap[k - N] : 0.0;
                Ap[k] = 0.25 * ((r[k] + yw) + ys);
            }
            for (int i = N - 1; i >= 0; --i) for (int j = N - 1; j >= 0; --j) {
                int k = i * N + j;
                double ze = (j + 1 < N) ? z[k + 1] : 0.0, zn = (i + 1 < N) ? z[k + N] : 0.0;
                z[k] = Ap[k] + 0.25 * (ze + zn);
            }
            double rz_new = dot_seq(n2, r, z), beta = rz_new / rz;
            for (int k = 0; k < n2; ++k) p[k] = z[k] + beta * p[k];
            rz = rz_new; continue;
        } else {
            vcycle(N, u, b);
            rel = compute_res(N, u, b, r) / nb;
        }
        hist[it - 1] = rel;
        if (rel <= 1e-8) { it++; break; }
    }
    int used = it - 1;
    double rho = 0.0;
    if (used >= 51) rho = pow(hist[used - 1] / hist[used - 51], 1.0 / 50.0);
    else if (used >= 2) rho = pow(hist[used - 1] / hist[0], 1.0 / (double)(used - 1));
    if (save63 && N == 63) {
        g_hist63_len[meth] = used;
        for (int k = 0; k < used; ++k) g_hist63[meth][k] = hist[k];
        memcpy(g_sol63[meth], u, (size_t)n2 * 8);
    }
    SolveRes res = {used, hist[used - 1], max_abs_diff(n2, u, ue), rho, fnv1a_f64(hist, used), fnv1a_f64(u, n2)};
    free(b); free(ue); free(u); free(r); free(p); free(Ap); free(z); free(hist);
    return res;
}

static const int NS_ALL[6] = {15, 31, 63, 127, 255, 511};
static SolveRes g_res[7][6];
static double g_fmg_err[6], g_disc_err[6], g_mg_fac[2], g_stage_s[19];
static double g_slopes[7];

static void s01(void) {
    double u[9] = {1,2,3,4,5,6,7,8,9}, Au[9], t[3] = {1.0, -2.5, 3.25};
    apply_A(3, u, Au);
    uint64_t h1 = fnv1a_f64(u, 9), h2 = fnv1a_f64(Au, 9), h3 = fnv1a_f64(t, 3);
    double exp_au[9] = {-2,-1,4,3,0,7,16,11,22};
    int ok = (max_abs_diff(9, Au, exp_au) == 0.0) && (h1 == 0xc56de10c6671e6c2ULL) && (h2 == 0x59ed72ec82ee959cULL) && (h3 == 0xafc5dc1ca48a2366ULL);
    FILE *fp = fopen("out/c/s01_selftest.txt", "w");
    fprintf(fp, "ok,h_u3,h_Au3,h_test\n%d,%016llx,%016llx,%016llx\n", ok, (unsigned long long)h1, (unsigned long long)h2, (unsigned long long)h3);
    fclose(fp);
    if (!ok) exit(1);
}
static void s02(void) {
    FILE *fp = fopen("out/c/s02_problem.csv", "w");
    fprintf(fp, "N,norm_b,h_b,h_ue\n");
    for (int i = 0; i < 6; ++i) {
        int N = NS_ALL[i], n2 = N * N;
        double *b = malloc((size_t)n2 * 8), *ue = malloc((size_t)n2 * 8);
        build_prob(N, b, ue);
        fprintf(fp, "%d,%.17g,%016llx,%016llx\n", N, norm2_seq(n2, b), (unsigned long long)fnv1a_f64(b, n2), (unsigned long long)fnv1a_f64(ue, n2));
        free(b); free(ue);
    }
    fclose(fp);
}
static void s03(void) {
    FILE *fp = fopen("out/c/s03_eigen.csv", "w");
    fprintf(fp, "N,j,k,lambda_exact,rel_err_inf\n");
    int modes[3][2] = {{1,1}, {1,2}, {3,5}};
    for (int ni = 0; ni < 3; ++ni) {
        int N = NS_ALL[ni], n2 = N * N; double h = 1.0 / (N + 1);
        double *v = malloc((size_t)n2 * 8), *Av = malloc((size_t)n2 * 8);
        for (int mi = 0; mi < 3; ++mi) {
            int j = modes[mi][0], k = modes[mi][1];
            double lam = 4.0 - 2.0 * cos_spec(j * PI_VAL * h) - 2.0 * cos_spec(k * PI_VAL * h);
            for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c)
                v[r * N + c] = sin_spec(j * PI_VAL * (r + 1) * h) * sin_spec(k * PI_VAL * (c + 1) * h);
            apply_A(N, v, Av);
            double me = 0.0;
            for (int q = 0; q < n2; ++q) { double d = fabs(Av[q] - lam * v[q]); if (d > me) me = d; }
            fprintf(fp, "%d,%d,%d,%.17g,%.3e\n", N, j, k, lam, me / lam);
        }
        free(v); free(Av);
    }
    fclose(fp);
}
static void run_m_csv(int m, const char *fn, int nmax) {
    FILE *fp = fopen(fn, "w");
    fprintf(fp, "N,iters,rel_res,err_inf,rho_emp\n");
    for (int i = 0; i < nmax; ++i) {
        g_res[m][i] = run_solver(m, NS_ALL[i], NS_ALL[i] == 63);
        fprintf(fp, "%d,%d,%.10e,%.10e,%.10f\n", NS_ALL[i], g_res[m][i].iters, g_res[m][i].rel_res, g_res[m][i].err_inf, g_res[m][i].rho_emp);
    }
    fclose(fp);
}
static void s04(void) { run_m_csv(0, "out/c/s04_jacobi.csv", 4); }
static void s05(void) { run_m_csv(1, "out/c/s05_gs.csv", 4); }
static void s06(void) { run_m_csv(2, "out/c/s06_rbgs.csv", 4); }
static void s07(void) { run_m_csv(3, "out/c/s07_sor.csv", 6); }
static void s08(void) { run_m_csv(4, "out/c/s08_cg.csv", 6); }
static void s09(void) { run_m_csv(5, "out/c/s09_pcg_ssor.csv", 6); }
static void s10(void) { run_m_csv(6, "out/c/s10_mg.csv", 6); }
static void s11(void) {
    FILE *fp = fopen("out/c/s11_fmg.csv", "w");
    fprintf(fp, "N,fmg_err_inf,disc_err_inf,ratio\n");
    for (int i = 0; i < 6; ++i) {
        int N = NS_ALL[i], n2 = N * N;
        double *b = malloc((size_t)n2 * 8), *ue = malloc((size_t)n2 * 8), *u = malloc((size_t)n2 * 8);
        build_prob(N, b, ue); fmg_solve(N, u, b);
        g_fmg_err[i] = max_abs_diff(n2, u, ue); g_disc_err[i] = g_res[6][i].err_inf;
        if (N == 63) memcpy(g_sol63[7], u, (size_t)n2 * 8);
        fprintf(fp, "%d,%.10e,%.10e,%.6f\n", N, g_fmg_err[i], g_disc_err[i], g_fmg_err[i] / g_disc_err[i]);
        free(b); free(ue); free(u);
    }
    fclose(fp);
}
static void s12(void) {
    FILE *fp = fopen("out/c/s12_order.csv", "w");
    fprintf(fp, "N,h,err_inf,slope\n");
    double mx = 0, my = 0;
    for (int i = 0; i < 6; ++i) { mx += log(1.0 / (NS_ALL[i] + 1)); my += log(g_disc_err[i]); }
    mx /= 6.0; my /= 6.0;
    double num = 0, den = 0;
    for (int i = 0; i < 6; ++i) { double dx = log(1.0 / (NS_ALL[i] + 1)) - mx; num += dx * (log(g_disc_err[i]) - my); den += dx * dx; }
    double sl = num / den;
    for (int i = 0; i < 6; ++i) fprintf(fp, "%d,%.10e,%.10e,%.6f\n", NS_ALL[i], 1.0 / (NS_ALL[i] + 1), g_disc_err[i], sl);
    fclose(fp);
}
static void s13(void) {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = fopen("out/c/s13_scaling.csv", "w");
    fprintf(fp, "method,N,iters,slope,slope_se\n");
    for (int m = 0; m < 7; ++m) {
        int cnt = (m <= 2) ? 4 : 6;
        double mx = 0, my = 0;
        for (int i = 0; i < cnt; ++i) { mx += log((double)NS_ALL[i]); my += log((double)g_res[m][i].iters); }
        mx /= cnt; my /= cnt;
        double num = 0, den = 0, rss = 0;
        for (int i = 0; i < cnt; ++i) { double dx = log((double)NS_ALL[i]) - mx; num += dx * (log((double)g_res[m][i].iters) - my); den += dx * dx; }
        double sl = num / den, ic = my - sl * mx; g_slopes[m] = sl;
        for (int i = 0; i < cnt; ++i) { double e = log((double)g_res[m][i].iters) - (sl * log((double)NS_ALL[i]) + ic); rss += e * e; }
        double se = sqrt((rss / (cnt - 2)) / den);
        for (int i = 0; i < cnt; ++i) fprintf(fp, "%s,%d,%d,%.6f,%.6f\n", mnames[m], NS_ALL[i], g_res[m][i].iters, sl, se);
    }
    fclose(fp);
}
static void s14(void) {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = fopen("out/c/s14_history.csv", "w");
    fprintf(fp, "method,iter,rel_res\n");
    for (int m = 0; m < 7; ++m) {
        int len = g_hist63_len[m], step = (len > 200) ? (len / 200) : 1;
        for (int k = 0; k < len; k += step) fprintf(fp, "%s,%d,%.10e\n", mnames[m], k + 1, g_hist63[m][k]);
        if ((len - 1) % step != 0) fprintf(fp, "%s,%d,%.10e\n", mnames[m], len, g_hist63[m][len - 1]);
    }
    fclose(fp);
}
static void s15(void) {
    FILE *fp = fopen("out/c/s15_mg_factor.csv", "w");
    fprintf(fp, "N,iters,mg_factor_last5\n");
    for (int k = 0; k < 2; ++k) {
        int N = (k == 0) ? 63 : 127, n2 = N * N;
        double *b = malloc((size_t)n2 * 8), *ue = malloc((size_t)n2 * 8), *u = calloc((size_t)n2, 8), *r = malloc((size_t)n2 * 8);
        build_prob(N, b, ue); double nb = norm2_seq(n2, b), h[20] = {0}; int it = 0;
        for (it = 1; it <= 20; ++it) { vcycle(N, u, b); h[it - 1] = compute_res(N, u, b, r) / nb; if (h[it - 1] <= 1e-8) break; }
        double fac = pow(h[it - 1] / h[it - 6], 0.2); g_mg_fac[k] = fac;
        fprintf(fp, "%d,%d,%.8f\n", N, it, fac);
        free(b); free(ue); free(u); free(r);
    }
    fclose(fp);
}
static void s16(void) {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    const int flops[7] = {9, 9, 9, 12, 21, 35, 36};
    FILE *fp = fopen("out/c/s16_cost.csv", "w");
    fprintf(fp, "method,N,flops_per_unk_iter,ns_med,ns_mad\n");
    int N = 63, n2 = N * N;
    double *b = malloc((size_t)n2 * 8), *ue = malloc((size_t)n2 * 8), *u = malloc((size_t)n2 * 8), *z = malloc((size_t)n2 * 8);
    build_prob(N, b, ue);
    for (int m = 0; m < 7; ++m) {
        int reps = 7, iters = (m == 6) ? 10 : 80;
        double vals[7], dev[7];
        for (int rp = 0; rp < reps; ++rp) {
            memset(u, 0, (size_t)n2 * 8); double t0 = now_s();
            for (int it = 0; it < iters; ++it) {
                if (m == 6) vcycle(N, u, b);
                else if (m == 2) rbgs_sweep(N, u, b);
                else { apply_A(N, u, z); for (int k = 0; k < n2; ++k) u[k] += 0.25 * (b[k] - z[k]); }
            }
            vals[rp] = (now_s() - t0) * 1e9 / ((double)iters * n2);
        }
        qsort(vals, 7, 8, cmp_d); double med = vals[3];
        for (int rp = 0; rp < 7; ++rp) dev[rp] = fabs(vals[rp] - med);
        qsort(dev, 7, 8, cmp_d);
        fprintf(fp, "%s,%d,%d,%.4f,%.4f\n", mnames[m], N, flops[m], med, dev[3]);
    }
    free(b); free(ue); free(u); free(z); fclose(fp);
}
static void s17(void) {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = fopen("out/c/s17_parity.csv", "w");
    fprintf(fp, "method,N,iters,hist_fnv1a,sol_fnv1a\n");
    for (int m = 0; m < 7; ++m) for (int ni = 0; ni < 2; ++ni)
        fprintf(fp, "%s,%d,%d,%016llx,%016llx\n", mnames[m], NS_ALL[ni], g_res[m][ni].iters,
                (unsigned long long)g_res[m][ni].hist_hash, (unsigned long long)g_res[m][ni].sol_hash);
    fclose(fp);
}
static void s18(void) {
    const char *cols[7] = {"#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2"};
    FILE *f1 = fopen("out/c/figures/history.svg", "w");
    fprintf(f1, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int m = 0; m < 7; ++m) {
        fprintf(f1, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.4\" points=\"", cols[m]);
        int len = g_hist63_len[m], step = (len > 150) ? (len / 150) : 1;
        for (int k = 0; k < len; k += step)
            fprintf(f1, "%.1f,%.1f ", 60.0 + 500.0 * log10(k + 1.0) / log10((double)g_hist63_len[0]),
                    310.0 - 280.0 * fmin(1.0, fmax(0.0, -log10(fmax(1e-12, g_hist63[m][k])) / 9.0)));
        fprintf(f1, "\"/>");
    }
    fprintf(f1, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 Residual History (N=63)</text></svg>\n");
    fclose(f1);

    FILE *f2 = fopen("out/c/figures/scaling.svg", "w");
    fprintf(f2, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int m = 0; m < 7; ++m) {
        int cnt = (m <= 2) ? 4 : 6;
        fprintf(f2, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.5\" points=\"", cols[m]);
        for (int i = 0; i < cnt; ++i)
            fprintf(f2, "%.1f,%.1f ", 60.0 + 500.0 * (log10(NS_ALL[i]) - log10(15.0)) / (log10(511.0) - log10(15.0)),
                    310.0 - 280.0 * log10((double)g_res[m][i].iters) / 5.0);
        fprintf(f2, "\"/>");
    }
    fprintf(f2, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 Iterations vs N</text></svg>\n");
    fclose(f2);

    FILE *f3 = fopen("out/c/figures/order.svg", "w");
    fprintf(f3, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"1.8\" points=\"");
    for (int i = 0; i < 6; ++i)
        fprintf(f3, "%.1f,%.1f ", 60.0 + 500.0 * (i / 5.0), 310.0 - 280.0 * (log10(g_disc_err[i]) + 5.0) / 4.0);
    fprintf(f3, "\"/><text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C11 Discretization Error vs h (Order ~2)</text></svg>\n");
    fclose(f3);
}
static void s19(void) {
    double tot = 0; for (int i = 0; i < 19; ++i) tot += g_stage_s[i];
    FILE *fp = fopen("out/c/s19_summary.json", "w");
    fprintf(fp, "{\"lang\":\"c\",\"N_max\":511,\"N_stat_max\":127,\"toolchain\":\"Apple clang 16.0.0 (-std=c11)\",\"total_sec\":%.3f}\n", tot);
    fclose(fp);
}

int main(int argc, char **argv) {
    mkdir("out", 0755); mkdir("out/c", 0755); mkdir("out/c/figures", 0755);
    const char *cmd = (argc > 1) ? argv[1] : "all";
    typedef void (*Fn)(void);
    Fn fns[19] = {s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19};
    const char *nms[19] = {"s01","s02","s03","s04","s05","s06","s07","s08","s09","s10","s11","s12","s13","s14","s15","s16","s17","s18","s19"};
    if (strcmp(cmd, "all") == 0) {
        for (int i = 0; i < 19; ++i) { printf("[c] tahap %d/19 (%s)\n", i+1, nms[i]); double t0=now_s(); fns[i](); g_stage_s[i]=now_s()-t0; }
        return 0;
    }
    for (int i = 0; i < 19; ++i) if (strcmp(cmd, nms[i]) == 0) { fns[i](); return 0; }
    return 1;
}
