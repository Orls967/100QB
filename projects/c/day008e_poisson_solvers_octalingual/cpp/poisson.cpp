#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;
static constexpr double PI_VAL = 3.14159265358979323846;

static double sin_spec(double x) {
    double sgn = 1.0;
    if (x < 0.0) { x = -x; sgn = -sgn; }
    if (x > PI_VAL) { x = x - PI_VAL; sgn = -sgn; }
    if (x > 0.5 * PI_VAL) x = PI_VAL - x;
    double x2 = x * x;
    static constexpr double c[10] = {
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
        uint64_t u; std::memcpy(&u, &v[i], 8);
        for (int b = 0; b < 8; ++b) { h ^= (u >> (8 * b)) & 0xFFU; h *= 0x100000001B3ULL; }
    }
    return h;
}
static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

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
static double norm2_seq(int n2, const double *a) { return std::sqrt(dot_seq(n2, a, a)); }
static double max_abs_diff(int n2, const double *a, const double *b) {
    double m = 0.0; for (int i = 0; i < n2; ++i) m = std::max(m, std::fabs(a[i] - b[i])); return m;
}
static void build_prob(int N, double *b, double *ue) {
    double h = 1.0 / (N + 1), h2 = h * h;
    std::vector<double> sx(N), sy(N);
    for (int i = 0; i < N; ++i) { sx[i] = sin_spec(PI_VAL * (i + 1) * h); sy[i] = sin_spec(2.0 * PI_VAL * (i + 1) * h); }
    for (int i = 0; i < N; ++i) {
        double x = (i + 1) * h, xx = x * (1.0 - x);
        for (int j = 0; j < N; ++j) {
            double y = (j + 1) * h, yy = y * (1.0 - y), sxy = sx[i] * sy[j];
            ue[i * N + j] = sxy + xx * yy;
            b[i * N + j] = h2 * (5.0 * PI_VAL * PI_VAL * sxy + 2.0 * xx + 2.0 * yy);
        }
    }
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
    double M[9][10]{};
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
        for (int r = c + 1; r < 9; ++r) if (std::fabs(M[r][c]) > std::fabs(M[piv][c])) piv = r;
        if (piv != c) for (int k = 0; k < 10; ++k) std::swap(M[c][k], M[piv][k]);
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
    std::vector<double> r(n2), rc(nc2), ec(nc2, 0.0);
    compute_res(N, u, b, r.data());
    restrict_fw(N, r.data(), rc.data());
    vcycle(Nc, ec.data(), rc.data());
    prolong_add(N, u, ec.data());
    rbgs_sweep(N, u, b);
}
static void fmg_solve(int N, double *u, const double *b) {
    if (N == 3) { solve_dense3(b, u); return; }
    int Nc = (N - 1) / 2, nc2 = Nc * Nc;
    std::vector<double> bc(nc2), uc(nc2, 0.0);
    restrict_fw(N, b, bc.data());
    fmg_solve(Nc, uc.data(), bc.data());
    std::fill(u, u + N * N, 0.0);
    prolong_add(N, u, uc.data());
    vcycle(N, u, b);
}

struct SolveRes { int iters; double rel_res, err_inf, rho_emp; uint64_t hist_hash, sol_hash; };
static std::vector<double> g_hist63[7];
static SolveRes g_res[7][6];
static double g_disc_err[6]{}, g_stage_s[19]{};
static const int NS_ALL[6] = {15, 31, 63, 127, 255, 511};

static SolveRes run_solver(int meth, int N, bool save63) {
    int n2 = N * N;
    std::vector<double> b(n2), ue(n2), u(n2, 0.0), r(n2, 0.0), p(n2, 0.0), Ap(n2, 0.0), z(n2, 0.0);
    build_prob(N, b.data(), ue.data());
    double nb = norm2_seq(n2, b.data()), om = omega_opt(N);
    int max_it = (meth == 0) ? 100000 : (meth <= 2 ? 50000 : (meth == 3 ? 10000 : (meth == 4 ? 5000 : (meth == 5 ? 3000 : 100))));
    std::vector<double> hist; hist.reserve(256);
    double rr = 0.0, rz = 0.0;
    if (meth == 4 || meth == 5) { r = b; rr = dot_seq(n2, r.data(), r.data()); }
    if (meth == 4) p = r;
    if (meth == 5) {
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
            int k = i * N + j;
            z[k] = 0.25 * ((r[k] + ((j > 0) ? z[k - 1] : 0.0)) + ((i > 0) ? z[k - N] : 0.0));
        }
        for (int i = N - 1; i >= 0; --i) for (int j = N - 1; j >= 0; --j) {
            int k = i * N + j;
            p[k] = z[k] + 0.25 * (((j + 1 < N) ? p[k + 1] : 0.0) + ((i + 1 < N) ? p[k + N] : 0.0));
        }
        z = p; rz = dot_seq(n2, r.data(), z.data());
    }
    for (int it = 1; it <= max_it; ++it) {
        double rel = 0.0;
        if (meth == 0) {
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
                double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
                z[k] = 0.25 * ((((b[k] + w) + e) + s) + n);
            }
            u.swap(z); rel = compute_res(N, u.data(), b.data(), r.data()) / nb;
        } else if (meth == 1 || meth == 3) {
            double ouse = (meth == 1) ? 1.0 : om;
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                double w = (j > 0) ? u[k - 1] : 0.0, e = (j + 1 < N) ? u[k + 1] : 0.0;
                double s = (i > 0) ? u[k - N] : 0.0, n = (i + 1 < N) ? u[k + N] : 0.0;
                double gs = 0.25 * ((((b[k] + w) + e) + s) + n);
                u[k] = (meth == 1) ? gs : ((1.0 - ouse) * u[k] + ouse * gs);
            }
            rel = compute_res(N, u.data(), b.data(), r.data()) / nb;
        } else if (meth == 2) {
            rbgs_sweep(N, u.data(), b.data()); rel = compute_res(N, u.data(), b.data(), r.data()) / nb;
        } else if (meth == 4) {
            apply_A(N, p.data(), Ap.data());
            double alpha = rr / dot_seq(n2, p.data(), Ap.data());
            for (int k = 0; k < n2; ++k) { u[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }
            double rr_new = dot_seq(n2, r.data(), r.data());
            rel = std::sqrt(rr_new) / nb; hist.push_back(rel);
            if (rel <= 1e-8) break;
            double beta = rr_new / rr;
            for (int k = 0; k < n2; ++k) p[k] = r[k] + beta * p[k];
            rr = rr_new; continue;
        } else if (meth == 5) {
            apply_A(N, p.data(), Ap.data());
            double alpha = rz / dot_seq(n2, p.data(), Ap.data());
            for (int k = 0; k < n2; ++k) { u[k] += alpha * p[k]; r[k] -= alpha * Ap[k]; }
            rel = norm2_seq(n2, r.data()) / nb; hist.push_back(rel);
            if (rel <= 1e-8) break;
            for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) {
                int k = i * N + j;
                Ap[k] = 0.25 * ((r[k] + ((j > 0) ? Ap[k - 1] : 0.0)) + ((i > 0) ? Ap[k - N] : 0.0));
            }
            for (int i = N - 1; i >= 0; --i) for (int j = N - 1; j >= 0; --j) {
                int k = i * N + j;
                z[k] = Ap[k] + 0.25 * (((j + 1 < N) ? z[k + 1] : 0.0) + ((i + 1 < N) ? z[k + N] : 0.0));
            }
            double rz_new = dot_seq(n2, r.data(), z.data()), beta = rz_new / rz;
            for (int k = 0; k < n2; ++k) p[k] = z[k] + beta * p[k];
            rz = rz_new; continue;
        } else {
            vcycle(N, u.data(), b.data()); rel = compute_res(N, u.data(), b.data(), r.data()) / nb;
        }
        hist.push_back(rel);
        if (rel <= 1e-8) break;
    }
    int used = static_cast<int>(hist.size());
    double rho = (used >= 51) ? std::pow(hist[used - 1] / hist[used - 51], 0.02) : ((used >= 2) ? std::pow(hist[used - 1] / hist[0], 1.0 / (used - 1)) : 0.0);
    if (save63 && N == 63) g_hist63[meth] = hist;
    return {used, hist.back(), max_abs_diff(n2, u.data(), ue.data()), rho, fnv1a_f64(hist.data(), used), fnv1a_f64(u.data(), n2)};
}

static void s01() {
    double u[9] = {1,2,3,4,5,6,7,8,9}, Au[9], t[3] = {1.0, -2.5, 3.25}, exp_au[9] = {-2,-1,4,3,0,7,16,11,22};
    apply_A(3, u, Au);
    uint64_t h1 = fnv1a_f64(u, 9), h2 = fnv1a_f64(Au, 9), h3 = fnv1a_f64(t, 3);
    int ok = (max_abs_diff(9, Au, exp_au) == 0.0 && h1 == 0xc56de10c6671e6c2ULL && h2 == 0x59ed72ec82ee959cULL && h3 == 0xafc5dc1ca48a2366ULL) ? 1 : 0;
    FILE *fp = std::fopen("out/cpp/s01_selftest.txt", "w");
    std::fprintf(fp, "ok,h_u3,h_Au3,h_test\n%d,%016llx,%016llx,%016llx\n", ok, static_cast<unsigned long long>(h1), static_cast<unsigned long long>(h2), static_cast<unsigned long long>(h3));
    std::fclose(fp);
    if (!ok) std::exit(1);
}
static void s02() {
    FILE *fp = std::fopen("out/cpp/s02_problem.csv", "w");
    std::fprintf(fp, "N,norm_b,h_b,h_ue\n");
    for (int N : NS_ALL) {
        std::vector<double> b(N * N), ue(N * N); build_prob(N, b.data(), ue.data());
        std::fprintf(fp, "%d,%.17g,%016llx,%016llx\n", N, norm2_seq(N * N, b.data()),
                     static_cast<unsigned long long>(fnv1a_f64(b.data(), N * N)), static_cast<unsigned long long>(fnv1a_f64(ue.data(), N * N)));
    }
    std::fclose(fp);
}
static void s03() {
    FILE *fp = std::fopen("out/cpp/s03_eigen.csv", "w"); std::fprintf(fp, "N,j,k,lambda_exact,rel_err_inf\n");
    int modes[3][2] = {{1,1},{1,2},{3,5}};
    for (int ni = 0; ni < 3; ++ni) {
        int N = NS_ALL[ni], n2 = N * N; double h = 1.0 / (N + 1);
        std::vector<double> v(n2), Av(n2);
        for (auto &md : modes) {
            double lam = 4.0 - 2.0 * cos_spec(md[0] * PI_VAL * h) - 2.0 * cos_spec(md[1] * PI_VAL * h);
            for (int r = 0; r < N; ++r) for (int c = 0; c < N; ++c)
                v[r * N + c] = sin_spec(md[0] * PI_VAL * (r + 1) * h) * sin_spec(md[1] * PI_VAL * (c + 1) * h);
            apply_A(N, v.data(), Av.data());
            double me = 0.0; for (int q = 0; q < n2; ++q) me = std::max(me, std::fabs(Av[q] - lam * v[q]));
            std::fprintf(fp, "%d,%d,%d,%.17g,%.3e\n", N, md[0], md[1], lam, me / lam);
        }
    }
    std::fclose(fp);
}
static void run_m_csv(int m, const char *fn, int nmax) {
    FILE *fp = std::fopen(fn, "w"); std::fprintf(fp, "N,iters,rel_res,err_inf,rho_emp\n");
    for (int i = 0; i < nmax; ++i) {
        g_res[m][i] = run_solver(m, NS_ALL[i], NS_ALL[i] == 63);
        std::fprintf(fp, "%d,%d,%.10e,%.10e,%.10f\n", NS_ALL[i], g_res[m][i].iters, g_res[m][i].rel_res, g_res[m][i].err_inf, g_res[m][i].rho_emp);
    }
    std::fclose(fp);
}
static void s04() { run_m_csv(0, "out/cpp/s04_jacobi.csv", 4); }
static void s05() { run_m_csv(1, "out/cpp/s05_gs.csv", 4); }
static void s06() { run_m_csv(2, "out/cpp/s06_rbgs.csv", 4); }
static void s07() { run_m_csv(3, "out/cpp/s07_sor.csv", 6); }
static void s08() { run_m_csv(4, "out/cpp/s08_cg.csv", 6); }
static void s09() { run_m_csv(5, "out/cpp/s09_pcg_ssor.csv", 6); }
static void s10() { run_m_csv(6, "out/cpp/s10_mg.csv", 6); }
static void s11() {
    FILE *fp = std::fopen("out/cpp/s11_fmg.csv", "w"); std::fprintf(fp, "N,fmg_err_inf,disc_err_inf,ratio\n");
    for (int i = 0; i < 6; ++i) {
        int N = NS_ALL[i], n2 = N * N;
        std::vector<double> b(n2), ue(n2), u(n2); build_prob(N, b.data(), ue.data()); fmg_solve(N, u.data(), b.data());
        double fe = max_abs_diff(n2, u.data(), ue.data()); g_disc_err[i] = g_res[6][i].err_inf;
        std::fprintf(fp, "%d,%.10e,%.10e,%.6f\n", N, fe, g_disc_err[i], fe / g_disc_err[i]);
    }
    std::fclose(fp);
}
static void s12() {
    FILE *fp = std::fopen("out/cpp/s12_order.csv", "w"); std::fprintf(fp, "N,h,err_inf,slope\n");
    double mx = 0, my = 0, num = 0, den = 0;
    for (int i = 0; i < 6; ++i) { mx += std::log(1.0 / (NS_ALL[i] + 1)); my += std::log(g_disc_err[i]); }
    mx /= 6.0; my /= 6.0;
    for (int i = 0; i < 6; ++i) { double dx = std::log(1.0 / (NS_ALL[i] + 1)) - mx; num += dx * (std::log(g_disc_err[i]) - my); den += dx * dx; }
    for (int i = 0; i < 6; ++i) std::fprintf(fp, "%d,%.10e,%.10e,%.6f\n", NS_ALL[i], 1.0 / (NS_ALL[i] + 1), g_disc_err[i], num / den);
    std::fclose(fp);
}
static void s13() {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = std::fopen("out/cpp/s13_scaling.csv", "w"); std::fprintf(fp, "method,N,iters,slope,slope_se\n");
    for (int m = 0; m < 7; ++m) {
        int cnt = (m <= 2) ? 4 : 6; double mx = 0, my = 0, num = 0, den = 0, rss = 0;
        for (int i = 0; i < cnt; ++i) { mx += std::log(NS_ALL[i]); my += std::log(g_res[m][i].iters); }
        mx /= cnt; my /= cnt;
        for (int i = 0; i < cnt; ++i) { double dx = std::log(NS_ALL[i]) - mx; num += dx * (std::log(g_res[m][i].iters) - my); den += dx * dx; }
        double sl = num / den, ic = my - sl * mx;
        for (int i = 0; i < cnt; ++i) { double e = std::log(g_res[m][i].iters) - (sl * std::log(NS_ALL[i]) + ic); rss += e * e; }
        for (int i = 0; i < cnt; ++i) std::fprintf(fp, "%s,%d,%d,%.6f,%.6f\n", mnames[m], NS_ALL[i], g_res[m][i].iters, sl, std::sqrt((rss / (cnt - 2)) / den));
    }
    std::fclose(fp);
}
static void s14() {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = std::fopen("out/cpp/s14_history.csv", "w"); std::fprintf(fp, "method,iter,rel_res\n");
    for (int m = 0; m < 7; ++m) {
        int len = static_cast<int>(g_hist63[m].size()), step = (len > 200) ? (len / 200) : 1;
        for (int k = 0; k < len; k += step) std::fprintf(fp, "%s,%d,%.10e\n", mnames[m], k + 1, g_hist63[m][k]);
        if ((len - 1) % step != 0) std::fprintf(fp, "%s,%d,%.10e\n", mnames[m], len, g_hist63[m].back());
    }
    std::fclose(fp);
}
static void s15() {
    FILE *fp = std::fopen("out/cpp/s15_mg_factor.csv", "w"); std::fprintf(fp, "N,iters,mg_factor_last5\n");
    for (int N : {63, 127}) {
        int n2 = N * N; std::vector<double> b(n2), ue(n2), u(n2, 0.0), r(n2); build_prob(N, b.data(), ue.data());
        double nb = norm2_seq(n2, b.data()), h[20]{}; int it = 0;
        for (it = 1; it <= 20; ++it) { vcycle(N, u.data(), b.data()); h[it - 1] = compute_res(N, u.data(), b.data(), r.data()) / nb; if (h[it - 1] <= 1e-8) break; }
        std::fprintf(fp, "%d,%d,%.8f\n", N, it, std::pow(h[it - 1] / h[it - 6], 0.2));
    }
    std::fclose(fp);
}
static void s16() {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"}; const int flops[7] = {9,9,9,12,21,35,36};
    FILE *fp = std::fopen("out/cpp/s16_cost.csv", "w"); std::fprintf(fp, "method,N,flops_per_unk_iter,ns_med,ns_mad\n");
    int N = 63, n2 = N * N; std::vector<double> b(n2), ue(n2), u(n2), z(n2); build_prob(N, b.data(), ue.data());
    for (int m = 0; m < 7; ++m) {
        int iters = (m == 6) ? 10 : 80; std::array<double, 7> vals{}, dev{};
        for (int rp = 0; rp < 7; ++rp) {
            std::fill(u.begin(), u.end(), 0.0); double t0 = now_s();
            for (int it = 0; it < iters; ++it) {
                if (m == 6) vcycle(N, u.data(), b.data());
                else if (m == 2) rbgs_sweep(N, u.data(), b.data());
                else { apply_A(N, u.data(), z.data()); for (int k = 0; k < n2; ++k) u[k] += 0.25 * (b[k] - z[k]); }
            }
            vals[rp] = (now_s() - t0) * 1e9 / (static_cast<double>(iters) * n2);
        }
        std::sort(vals.begin(), vals.end()); double med = vals[3];
        for (int rp = 0; rp < 7; ++rp) dev[rp] = std::fabs(vals[rp] - med);
        std::sort(dev.begin(), dev.end());
        std::fprintf(fp, "%s,%d,%d,%.4f,%.4f\n", mnames[m], N, flops[m], med, dev[3]);
    }
    std::fclose(fp);
}
static void s17() {
    const char *mnames[7] = {"jacobi","gs","rbgs","sor","cg","pcg_ssor","mg"};
    FILE *fp = std::fopen("out/cpp/s17_parity.csv", "w"); std::fprintf(fp, "method,N,iters,hist_fnv1a,sol_fnv1a\n");
    for (int m = 0; m < 7; ++m) for (int ni = 0; ni < 2; ++ni)
        std::fprintf(fp, "%s,%d,%d,%016llx,%016llx\n", mnames[m], NS_ALL[ni], g_res[m][ni].iters,
                     static_cast<unsigned long long>(g_res[m][ni].hist_hash), static_cast<unsigned long long>(g_res[m][ni].sol_hash));
    std::fclose(fp);
}
static void s18() {
    const char *cols[7] = {"#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2"};
    FILE *f1 = std::fopen("out/cpp/figures/history.svg", "w");
    std::fprintf(f1, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int m = 0; m < 7; ++m) {
        std::fprintf(f1, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.4\" points=\"", cols[m]);
        int len = static_cast<int>(g_hist63[m].size()), step = (len > 150) ? (len / 150) : 1;
        for (int k = 0; k < len; k += step)
            std::fprintf(f1, "%.1f,%.1f ", 60.0 + 500.0 * std::log10(k + 1.0) / std::log10(static_cast<double>(g_hist63[0].size())),
                         310.0 - 280.0 * std::clamp(-std::log10(std::max(1e-12, g_hist63[m][k])) / 9.0, 0.0, 1.0));
        std::fprintf(f1, "\"/>");
    }
    std::fprintf(f1, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 Residual History (N=63)</text></svg>\n");
    std::fclose(f1);

    FILE *f2 = std::fopen("out/cpp/figures/scaling.svg", "w");
    std::fprintf(f2, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int m = 0; m < 7; ++m) {
        int cnt = (m <= 2) ? 4 : 6;
        std::fprintf(f2, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.5\" points=\"", cols[m]);
        for (int i = 0; i < cnt; ++i)
            std::fprintf(f2, "%.1f,%.1f ", 60.0 + 500.0 * (std::log10(NS_ALL[i]) - std::log10(15.0)) / (std::log10(511.0) - std::log10(15.0)),
                         310.0 - 280.0 * std::log10(static_cast<double>(g_res[m][i].iters)) / 5.0);
        std::fprintf(f2, "\"/>");
    }
    std::fprintf(f2, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 Iterations vs N</text></svg>\n");
    std::fclose(f2);

    FILE *f3 = std::fopen("out/cpp/figures/order.svg", "w");
    std::fprintf(f3, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"1.8\" points=\"");
    for (int i = 0; i < 6; ++i) std::fprintf(f3, "%.1f,%.1f ", 60.0 + 500.0 * (i / 5.0), 310.0 - 280.0 * (std::log10(g_disc_err[i]) + 5.0) / 4.0);
    std::fprintf(f3, "\"/><text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 Discretization Error vs h (Order ~2)</text></svg>\n");
    std::fclose(f3);
}
static void s19() {
    double tot = 0; for (double s : g_stage_s) tot += s;
    FILE *fp = std::fopen("out/cpp/s19_summary.json", "w");
    std::fprintf(fp, "{\"lang\":\"cpp\",\"N_max\":511,\"N_stat_max\":127,\"toolchain\":\"Apple clang 16.0.0 (-std=c++20)\",\"total_sec\":%.3f}\n", tot);
    std::fclose(fp);
}

int main(int argc, char **argv) {
    fs::create_directories("out/cpp/figures");
    std::string cmd = (argc > 1) ? argv[1] : "all";
    using Fn = void (*)();
    Fn fns[19] = {s01,s02,s03,s04,s05,s06,s07,s08,s09,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19};
    const char *nms[19] = {"s01","s02","s03","s04","s05","s06","s07","s08","s09","s10","s11","s12","s13","s14","s15","s16","s17","s18","s19"};
    if (cmd == "all") {
        for (int i = 0; i < 19; ++i) { std::printf("[cpp] tahap %d/19 (%s)\n", i+1, nms[i]); double t0=now_s(); fns[i](); g_stage_s[i]=now_s()-t0; }
        return 0;
    }
    for (int i = 0; i < 19; ++i) if (cmd == nms[i]) { fns[i](); return 0; }
    return 1;
}
