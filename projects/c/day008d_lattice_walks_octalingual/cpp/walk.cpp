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

static constexpr uint64_t MASTER_SEED = 20261126ULL;
static constexpr int W_WALKS = 200000;
static constexpr int W_RUIN = 20000;
static constexpr int BENCH_STEPS = 10000000;
static constexpr int T_EXACT_MAX = 5000;
static constexpr int TS_MSD[13] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096};
static constexpr int TS_RET[4] = {10,100,1000,10000};

static inline uint64_t sm64(uint64_t z) {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
struct Rng {
    uint64_t s[4]{};
    explicit Rng(uint64_t seed) {
        uint64_t st = seed;
        for (int i = 0; i < 4; ++i) { st += 0x9E3779B97F4A7C15ULL; s[i] = sm64(st); }
    }
    inline uint64_t u64() {
        auto rotl = [](uint64_t x, int k) { return (x << k) | (x >> (64 - k)); };
        uint64_t res = rotl(s[0] + s[3], 23) + s[0];
        uint64_t t = s[1] << 17;
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
        s[2] ^= t; s[3] = rotl(s[3], 45);
        return res;
    }
    inline double u01() { return static_cast<double>(u64() >> 11) * 0x1.0p-53; }
    inline int rint(int n) {
        int k = static_cast<int>(std::floor(u01() * static_cast<double>(n)));
        return (k >= n) ? (n - 1) : k;
    }
};
static uint64_t task_seed(int d, int st, int idx) {
    uint64_t h = 0x517CC1B727220A95ULL;
    h = sm64(h ^ static_cast<uint64_t>(d));
    h = sm64(h ^ static_cast<uint64_t>(st));
    h = sm64(h ^ static_cast<uint64_t>(idx));
    return sm64(MASTER_SEED ^ h);
}
static uint64_t fnv1a_i32(const int32_t *v, int n) {
    uint64_t h = 0xCBF29CE484222325ULL;
    for (int i = 0; i < n; ++i) {
        uint32_t u = static_cast<uint32_t>(v[i]);
        for (int b = 0; b < 4; ++b) { h ^= (u >> (8 * b)) & 0xFFU; h *= 0x100000001B3ULL; }
    }
    return h;
}
static double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

static double g_msd[3][13]{}, g_msd_se[3][13]{}, g_ret_sim[3][4]{}, g_ret_se[3][4]{}, g_ret_ex[3][4]{};
static double g_ruin_p[6]{}, g_ruin_pe[6]{};
static std::vector<double> g_s03_1024[3];
static std::vector<uint8_t> g_s04_10k[3];
static double g_stage_sec[12]{};

static void run_s01() {
    Rng r(MASTER_SEED);
    uint64_t x0 = r.u64(), x1 = r.u64(), x2 = r.u64();
    bool ok = (x0 == 0xA9920E0E7DAB3EF0ULL && x1 == 0x033B9CEB4B3BA7A0ULL && x2 == 0xE7473A6A39AFCEA6ULL);
    double w[5] = {1.0, 1.0 / 3.0, 0, 0, 0};
    for (int n = 1; n < 4; ++n) w[n+1] = ((10.0*n*n + 10.0*n + 3.0)*w[n] - static_cast<double>(n*n)*w[n-1]) / (9.0*(n+1)*(n+1));
    bool enum_ok = (std::fabs(w[2] - 15.0/81.0) < 1e-14 && std::fabs(w[3] - 93.0/729.0) < 1e-14 && std::fabs(w[4] - 639.0/6561.0) < 1e-14);
    FILE *fp = std::fopen("out/cpp/s01_selftest.txt", "w");
    std::fprintf(fp, "ok,x0,x1,x2,enum_ok\n%d,%016llx,%016llx,%016llx,%d\n",
                 (ok && enum_ok) ? 1 : 0, static_cast<unsigned long long>(x0),
                 static_cast<unsigned long long>(x1), static_cast<unsigned long long>(x2), enum_ok ? 1 : 0);
    std::fclose(fp);
    if (!(ok && enum_ok)) std::exit(1);
}

static void run_s02() {
    FILE *fp = std::fopen("out/cpp/s02_parity.csv", "w");
    std::fprintf(fp, "d,walk,x0,x1,x2,fnv1a_hex\n");
    for (int d = 1; d <= 3; ++d) {
        int32_t coords[24]{}, pos_all[8][3]{};
        int ci = 0;
        for (int w = 0; w < 8; ++w) {
            Rng r(task_seed(d, 2, w));
            int32_t p[3] = {0, 0, 0};
            for (int s = 0; s < 200; ++s) { int k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
            for (int a = 0; a < 3; ++a) pos_all[w][a] = p[a];
            for (int a = 0; a < d; ++a) coords[ci++] = p[a];
        }
        uint64_t h = fnv1a_i32(coords, 8 * d);
        for (int w = 0; w < 8; ++w)
            std::fprintf(fp, "%d,%d,%d,%d,%d,%016llx\n", d, w, pos_all[w][0], pos_all[w][1], pos_all[w][2], static_cast<unsigned long long>(h));
    }
    std::fclose(fp);
}

static void run_s03() {
    FILE *fp = std::fopen("out/cpp/s03_msd.csv", "w");
    std::fprintf(fp, "d,t,msd,se,exact\n");
    for (int d = 1; d <= 3; ++d) {
        g_s03_1024[d - 1].resize(W_WALKS);
        double sum1[13]{}, sum2[13]{};
        for (int w = 0; w < W_WALKS; ++w) {
            Rng r(task_seed(d, 3, w));
            int p[3] = {0, 0, 0}, ti = 0;
            for (int t = 1; t <= 4096; ++t) {
                int k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1;
                if (t == TS_MSD[ti]) {
                    double r2 = static_cast<double>(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
                    sum1[ti] += r2; sum2[ti] += r2 * r2;
                    if (ti == 10) g_s03_1024[d - 1][w] = r2;
                    ++ti;
                }
            }
        }
        for (int ti = 0; ti < 13; ++ti) {
            double m = sum1[ti] / W_WALKS;
            double se = std::sqrt(std::max(0.0, sum2[ti] / W_WALKS - m * m) / W_WALKS);
            g_msd[d - 1][ti] = m; g_msd_se[d - 1][ti] = se;
            std::fprintf(fp, "%d,%d,%.8f,%.8f,%d\n", d, TS_MSD[ti], m, se, TS_MSD[ti]);
        }
    }
    std::fclose(fp);
}

static void run_s04() {
    FILE *fp = std::fopen("out/cpp/s04_return_sim.csv", "w");
    std::fprintf(fp, "d,T,p_hat,se\n");
    for (int d = 1; d <= 3; ++d) {
        g_s04_10k[d - 1].resize(W_WALKS);
        int hits[4]{};
        for (int w = 0; w < W_WALKS; ++w) {
            Rng r(task_seed(d, 4, w));
            int p[3] = {0, 0, 0}, ret_t = 0;
            for (int t = 1; t <= 10000; ++t) {
                int k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1;
                if ((p[0] | p[1] | p[2]) == 0) { ret_t = t; break; }
            }
            g_s04_10k[d - 1][w] = (ret_t > 0) ? 1 : 0;
            if (ret_t > 0) for (int ti = 0; ti < 4; ++ti) if (ret_t <= TS_RET[ti]) ++hits[ti];
        }
        for (int ti = 0; ti < 4; ++ti) {
            double ph = static_cast<double>(hits[ti]) / W_WALKS;
            double se = std::sqrt(std::max(0.0, ph * (1.0 - ph)) / W_WALKS);
            g_ret_sim[d - 1][ti] = ph; g_ret_se[d - 1][ti] = se;
            std::fprintf(fp, "%d,%d,%.8f,%.8f\n", d, TS_RET[ti], ph, se);
        }
    }
    std::fclose(fp);
}

static void run_s05() {
    std::vector<double> u1(T_EXACT_MAX + 1, 1.0), u2(T_EXACT_MAX + 1, 1.0), w3(T_EXACT_MAX + 1, 1.0), u3(T_EXACT_MAX + 1, 1.0), f(T_EXACT_MAX + 1, 0.0);
    w3[1] = 1.0 / 3.0;
    for (int n = 1; n < T_EXACT_MAX; ++n)
        w3[n + 1] = ((10.0*n*n + 10.0*n + 3.0)*w3[n] - static_cast<double>(n)*n*w3[n - 1]) / (9.0*(n + 1.0)*(n + 1.0));
    for (int n = 1; n <= T_EXACT_MAX; ++n) {
        u1[n] = u1[n - 1] * (2.0*n - 1.0) / (2.0*n);
        u2[n] = u1[n] * u1[n]; u3[n] = u1[n] * w3[n];
    }
    const std::vector<double> *us[3] = {&u1, &u2, &u3};
    FILE *fp = std::fopen("out/cpp/s05_return_exact.csv", "w");
    std::fprintf(fp, "d,T,F_exact\n");
    for (int d = 1; d <= 3; ++d) {
        const auto &u = *us[d - 1]; double acc = 0.0;
        for (int n = 1; n <= T_EXACT_MAX; ++n) {
            double s = u[n];
            for (int k = 1; k < n; ++k) s -= f[k] * u[n - k];
            f[n] = s; acc += s;
            for (int ti = 0; ti < 4; ++ti) if (2 * n == TS_RET[ti]) { g_ret_ex[d - 1][ti] = acc; std::fprintf(fp, "%d,%d,%.10f\n", d, TS_RET[ti], acc); }
        }
    }
    std::fclose(fp);
}

static void run_s06() {
    FILE *fp = std::fopen("out/cpp/s06_return_compare.csv", "w");
    std::fprintf(fp, "d,T,p_hat,se,F_exact,z_score,p_inf_est\n");
    double p3_inf = g_ret_sim[2][3] + (g_ret_sim[2][3] - g_ret_sim[2][2]) / (std::sqrt(10.0) - 1.0);
    for (int d = 1; d <= 3; ++d)
        for (int ti = 0; ti < 4; ++ti) {
            double z = (g_ret_se[d - 1][ti] > 1e-12) ? (g_ret_sim[d - 1][ti] - g_ret_ex[d - 1][ti]) / g_ret_se[d - 1][ti] : 0.0;
            std::fprintf(fp, "%d,%d,%.8f,%.8f,%.10f,%.4f,%.8f\n", d, TS_RET[ti], g_ret_sim[d - 1][ti], g_ret_se[d - 1][ti],
                         g_ret_ex[d - 1][ti], z, (d == 3 && ti == 3) ? p3_inf : 0.0);
        }
    std::fclose(fp);
}

static void run_s07() {
    const int Ns[6] = {10,10,10,50,50,50}, xs[6] = {2,5,8,10,25,40};
    FILE *fp = std::fopen("out/cpp/s07_ruin.csv", "w");
    std::fprintf(fp, "N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z\n");
    for (int ci = 0; ci < 6; ++ci) {
        int N = Ns[ci], x0 = xs[ci], wins = 0; double st1 = 0, st2 = 0;
        for (int w = 0; w < W_RUIN; ++w) {
            Rng r(task_seed(1, 7, ci * 100000 + w));
            int x = x0, steps = 0;
            while (x > 0 && x < N) { x += (r.rint(2) == 0) ? 1 : -1; ++steps; }
            if (x == N) ++wins;
            st1 += steps; st2 += static_cast<double>(steps) * steps;
        }
        double ph = static_cast<double>(wins) / W_RUIN, pse = std::sqrt(std::max(0.0, ph * (1.0 - ph)) / W_RUIN), pe = static_cast<double>(x0) / N;
        double th = st1 / W_RUIN, tse = std::sqrt(std::max(0.0, st2 / W_RUIN - th * th) / W_RUIN), te = static_cast<double>(x0) * (N - x0);
        g_ruin_p[ci] = ph; g_ruin_pe[ci] = pe;
        std::fprintf(fp, "%d,%d,%.8f,%.8f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                     N, x0, ph, pse, pe, (ph - pe) / pse, th, tse, te, (th - te) / tse);
    }
    std::fclose(fp);
}

static void run_s08() {
    FILE *fp = std::fopen("out/cpp/s08_endpoint.csv", "w");
    std::fprintf(fp, "d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count\n");
    int hist[40]{};
    for (int d = 1; d <= 3; ++d) {
        double m1 = 0, m2 = 0, m4 = 0; int n_ep = 20000;
        for (int w = 0; w < n_ep; ++w) {
            Rng r(task_seed(d, 8, w)); int p[3] = {0, 0, 0};
            for (int s = 0; s < 1024; ++s) { int k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
            double x = static_cast<double>(p[0]);
            m1 += x; m2 += x * x; m4 += x * x * x * x;
            if (d == 1) { int b = std::clamp(static_cast<int>(std::floor((x + 128.0) / 6.4)), 0, 39); ++hist[b]; }
        }
        m1 /= n_ep; double v = m2 / n_ep - m1 * m1;
        std::fprintf(fp, "%d,%.6f,%.6f,%.6f,%.6f,-1,0,0,%d\n", d, v, 1024.0 / d, (m4 / n_ep) / (v * v), 3.0 + (d - 3.0) / 1024.0, n_ep);
    }
    for (int b = 0; b < 40; ++b)
        std::fprintf(fp, "1,0,0,0,0,%d,%.2f,%.2f,%d\n", b, -128.0 + b * 6.4, -128.0 + (b + 1) * 6.4, hist[b]);
    std::fclose(fp);
}

static void run_s09() {
    FILE *fp = std::fopen("out/cpp/s09_ci.csv", "w");
    std::fprintf(fp, "metric,d,estimate,boot_se,ci_lo,ci_hi\n");
    std::array<double, 200> boot{};
    for (int d = 1; d <= 3; ++d) {
        Rng r1(task_seed(d, 9, 0));
        for (int b = 0; b < 200; ++b) {
            double s = 0; for (int i = 0; i < W_WALKS; ++i) s += g_s03_1024[d - 1][r1.rint(W_WALKS)];
            boot[b] = s / W_WALKS;
        }
        std::sort(boot.begin(), boot.end());
        double bm = 0, bv = 0; for (double v : boot) bm += v; bm /= 200.0;
        for (double v : boot) bv += (v - bm) * (v - bm);
        std::fprintf(fp, "msd_1024,%d,%.6f,%.6f,%.6f,%.6f\n", d, g_msd[d - 1][10], std::sqrt(bv / 199.0), boot[5], boot[194]);

        Rng r2(task_seed(d, 9, 1));
        for (int b = 0; b < 200; ++b) {
            int s = 0; for (int i = 0; i < W_WALKS; ++i) s += g_s04_10k[d - 1][r2.rint(W_WALKS)];
            boot[b] = static_cast<double>(s) / W_WALKS;
        }
        std::sort(boot.begin(), boot.end());
        bm = bv = 0; for (double v : boot) bm += v; bm /= 200.0;
        for (double v : boot) bv += (v - bm) * (v - bm);
        std::fprintf(fp, "return_10000,%d,%.6f,%.6f,%.6f,%.6f\n", d, g_ret_sim[d - 1][3], std::sqrt(bv / 199.0), boot[5], boot[194]);
    }
    std::fclose(fp);
}

static void run_s10() {
    FILE *fp = std::fopen("out/cpp/s10_fit.csv", "w");
    std::fprintf(fp, "d,slope,slope_se,intercept\n");
    for (int d = 1; d <= 3; ++d) {
        double mx = 0, my = 0;
        for (int i = 0; i < 13; ++i) { mx += std::log(TS_MSD[i]); my += std::log(g_msd[d - 1][i]); }
        mx /= 13.0; my /= 13.0;
        double num = 0, den = 0;
        for (int i = 0; i < 13; ++i) { double dx = std::log(TS_MSD[i]) - mx; num += dx * (std::log(g_msd[d - 1][i]) - my); den += dx * dx; }
        double sl = num / den, ic = my - sl * mx, rss = 0;
        for (int i = 0; i < 13; ++i) { double e = std::log(g_msd[d - 1][i]) - (sl * std::log(TS_MSD[i]) + ic); rss += e * e; }
        std::fprintf(fp, "%d,%.6f,%.6f,%.6f\n", d, sl, std::sqrt((rss / 11.0) / den), ic);
    }
    std::fclose(fp);
}

static void run_s11() {
    FILE *fp = std::fopen("out/cpp/s11_bench.csv", "w");
    std::fprintf(fp, "d,steps_per_rep,reps,ns_med,ns_mad\n");
    volatile int64_t sink = 0;
    for (int d = 1; d <= 3; ++d) {
        std::array<double, 7> vals{}, dev{};
        for (int rep = 0; rep < 7; ++rep) {
            Rng r(task_seed(d, 11, rep)); int p[3] = {0, 0, 0};
            double t0 = now_s();
            for (int s = 0; s < BENCH_STEPS; ++s) { int k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
            vals[rep] = (now_s() - t0) * 1e9 / static_cast<double>(BENCH_STEPS);
            sink += p[0] + p[1] + p[2];
        }
        std::sort(vals.begin(), vals.end()); double med = vals[3];
        for (int rep = 0; rep < 7; ++rep) dev[rep] = std::fabs(vals[rep] - med);
        std::sort(dev.begin(), dev.end());
        std::fprintf(fp, "%d,%d,7,%.4f,%.4f\n", d, BENCH_STEPS, med, dev[3]);
    }
    (void)sink;
    std::fclose(fp);
}

static void run_s12() {
    const char *cols[3] = {"#1f77b4", "#ff7f0e", "#2ca02c"};
    FILE *f1 = std::fopen("out/cpp/figures/msd.svg", "w");
    std::fprintf(f1, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for (int d = 0; d < 3; ++d) {
        std::fprintf(f1, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.6\" points=\"", cols[d]);
        for (int i = 0; i < 13; ++i)
            std::fprintf(f1, "%.1f,%.1f ", 60.0 + 500.0 * std::log10(TS_MSD[i]) / std::log10(4096.0), 310.0 - 280.0 * std::log10(std::max(1.0, g_msd[d][i])) / std::log10(4096.0));
        std::fprintf(f1, "\"/>");
    }
    std::fprintf(f1, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 MSD Log-Log (y=t)</text></svg>\n");
    std::fclose(f1);

    FILE *f2 = std::fopen("out/cpp/figures/return_prob.svg", "w");
    std::fprintf(f2, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/>");
    for (int d = 0; d < 3; ++d) {
        std::fprintf(f2, "<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"1.6\" points=\"", cols[d]);
        for (int i = 0; i < 4; ++i) std::fprintf(f2, "%.1f,%.1f ", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * g_ret_ex[d][i]);
        std::fprintf(f2, "\"/>");
        for (int i = 0; i < 4; ++i) std::fprintf(f2, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"4\" fill=\"%s\"/>", 60.0 + 500.0 * (i / 3.0), 310.0 - 280.0 * g_ret_sim[d][i], cols[d]);
    }
    std::fprintf(f2, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 Return Prob F_d(T) Exact vs Sim</text></svg>\n");
    std::fclose(f2);

    FILE *f3 = std::fopen("out/cpp/figures/ruin.svg", "w");
    std::fprintf(f3, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"30\" width=\"500\" height=\"280\" fill=\"none\" stroke=\"#ccc\"/><line x1=\"60\" y1=\"310\" x2=\"560\" y2=\"30\" stroke=\"#000\" stroke-dasharray=\"4,3\"/>");
    for (int i = 0; i < 6; ++i)
        std::fprintf(f3, "<circle cx=\"%.1f\" cy=\"%.1f\" r=\"5\" fill=\"#1f77b4\"/>", 60.0 + 500.0 * g_ruin_pe[i], 310.0 - 280.0 * g_ruin_p[i]);
    std::fprintf(f3, "<text x=\"300\" y=\"20\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">C++20 Gambler's Ruin Exact vs Sim</text></svg>\n");
    std::fclose(f3);

    double tot = 0; for (double s : g_stage_sec) tot += s;
    FILE *fj = std::fopen("out/cpp/s12_summary.json", "w");
    std::fprintf(fj, "{\"lang\":\"cpp\",\"W\":%d,\"ruin_W\":%d,\"bench_steps\":%d,\"toolchain\":\"Apple clang 16.0.0 (-std=c++20)\",\"total_sec\":%.3f}\n",
                 W_WALKS, W_RUIN, BENCH_STEPS, tot);
    std::fclose(fj);
}

int main(int argc, char **argv) {
    fs::create_directories("out/cpp/figures");
    std::string cmd = (argc > 1) ? argv[1] : "all";
    using Fn = void (*)();
    Fn fns[12] = {run_s01,run_s02,run_s03,run_s04,run_s05,run_s06,run_s07,run_s08,run_s09,run_s10,run_s11,run_s12};
    const char *names[12] = {"s01","s02","s03","s04","s05","s06","s07","s08","s09","s10","s11","s12"};
    if (cmd == "all") {
        for (int i = 0; i < 12; ++i) {
            std::printf("[cpp] tahap %d/12 (%s)\n", i + 1, names[i]);
            double t0 = now_s(); fns[i](); g_stage_sec[i] = now_s() - t0;
        }
        return 0;
    }
    for (int i = 0; i < 12; ++i) if (cmd == names[i]) { fns[i](); return 0; }
    return 1;
}
