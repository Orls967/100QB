#include "reference.hpp"
#include "variants.hpp"
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>

namespace fs = std::filesystem;
using clk = std::chrono::steady_clock;

template <class F> std::vector<double> time_runs(F&& fn, int warm = 2, int reps = 7) {
    for (int i = 0; i < warm; ++i) fn();
    std::vector<double> ts(reps);
    for (int i = 0; i < reps; ++i) { auto t0 = clk::now(); fn(); ts[i] = std::chrono::duration<double>(clk::now() - t0).count(); }
    return ts;
}
static double median_of(std::vector<double> v) { std::sort(v.begin(), v.end()); return v[v.size() / 2]; }
static double rel_err(const std::vector<double>& a, const std::vector<double>& b) {
    double mx = 0.0, df = 0.0;
    for (size_t k = 0; k < a.size(); ++k) { mx = std::max(mx, std::abs(b[k])); df = std::max(df, std::abs(a[k] - b[k])); }
    return df / mx;
}
static void run_var(const std::string& s, const Problem& prob, const Opts& opt, std::vector<double>& out) {
    if (s == "V0") solve_v0(prob, out); else if (s == "V1") solve_v1(prob, out); else if (s == "V2") solve_v2(prob, out);
    else if (s == "V3") solve_v3(prob, out, opt); else if (s == "V4") solve_v4(prob, out, opt); else solve_v5(prob, out, opt);
}
static void run_machine(const fs::path& dir) {
    FILE* f = std::fopen((dir / "machine.csv").c_str(), "w");
    std::fprintf(f, "kernel,threads,rep,seconds,metric,value\n");
    const size_t N = 1ULL << 24, iters = 20000000ULL;
    std::vector<double> a(N, 0.0), b(N, 1.1), c(N, 2.2);
    for (int T : {1, 10}) {
        auto ts = time_runs([&]() {
            std::vector<std::thread> ths;
            for (int t = 0; t < T; ++t) ths.emplace_back([&, t]() {
                for (size_t i = (t * N) / T, e = ((t + 1) * N) / T; i < e; ++i) a[i] = b[i] + 3.3 * c[i];
            });
            for (auto& th : ths) th.join();
        });
        for (int r = 0; r < 7; ++r) std::fprintf(f, "stream_triad,%d,%d,%.9f,GB_s,%.6f\n", T, r, ts[r], (24.0 * N) / (ts[r] * 1e9));
        volatile double sink = 0.0;
        auto fp = time_runs([&]() {
            std::vector<std::thread> ths; std::vector<double> out(T, 0.0);
            for (int t = 0; t < T; ++t) ths.emplace_back([&, t]() {
                double x[16]; for (int k = 0; k < 16; ++k) x[k] = 1.0 + 1e-6 * (k + t);
                for (size_t i = 0; i < iters; ++i) for (int k = 0; k < 16; ++k) x[k] = x[k] * 1.0000001 + 1e-7;
                double s = 0; for (int k = 0; k < 16; ++k) s += x[k]; out[t] = s;
            });
            for (auto& th : ths) th.join();
            double s = 0; for (double v : out) s += v; sink = s;
        });
        (void)sink;
        for (int r = 0; r < 7; ++r) std::fprintf(f, "reg_peak,%d,%d,%.9f,GFLOP_s,%.6f\n", T, r, fp[r], (32.0 * iters * T) / (fp[r] * 1e9));
    }
    std::fclose(f);
}
static void run_checks(const fs::path& dir, const char* fname, const Opts& opt) {
    FILE* f = std::fopen((dir / fname).c_str(), "w");
    std::fprintf(f, "variant,nx,nv,m,threads,price,hash,rel_err_vs_v0\n");
    const int nxs[] = {101, 201, 401, 801}, nvs[] = {51, 101, 201, 401};
    for (int g = 0; g < 4; ++g) {
        Problem prob(Params::set_a(), nxs[g], nvs[g], 200); std::vector<double> u0, u; solve_v0(prob, u0);
        for (const char* v : {"V0", "V1", "V2", "V3", "V4"}) {
            run_var(v, prob, opt, u);
            std::fprintf(f, "%s,%d,%d,200,1,%.12f,%llu,%.3e\n", v, nxs[g], nvs[g], price_at(u.data(), prob.g),
                         (unsigned long long)fnv1a_hash(u.data(), u.size()), rel_err(u, u0));
        }
        for (int t : {1, 2, 4, 10}) {
            Opts o = opt; o.threads = t; solve_v5(prob, u, o);
            std::fprintf(f, "V5_T%d,%d,%d,200,%d,%.12f,%llu,%.3e\n", t, nxs[g], nvs[g], t, price_at(u.data(), prob.g),
                         (unsigned long long)fnv1a_hash(u.data(), u.size()), rel_err(u, u0));
        }
    }
    std::fclose(f);
}
int main(int argc, char** argv) {
    std::string mode = (argc > 1) ? argv[1] : "all";
    fs::path dir = (argc > 2) ? fs::path(argv[2]) : fs::path("../../../outputs/day006d");
    fs::create_directories(dir);
    Opts opt{16, 16, 10};
    if (mode == "native") {
        run_checks(dir, "checks_native.csv", opt);
        FILE* f = std::fopen((dir / "bench_native.csv").c_str(), "w");
        std::fprintf(f, "variant,nx,nv,m,rep,seconds,B,tile,threads\n");
        Problem prob(Params::set_a(), 801, 401, 200); std::vector<double> u;
        for (const char* v : {"V0", "V1", "V2", "V3", "V4", "V5"}) {
            auto ts = time_runs([&]() { run_var(v, prob, opt, u); });
            for (int r = 0; r < 7; ++r) std::fprintf(f, "%s,801,401,200,%d,%.9f,%d,%d,%d\n", v, r, ts[r], opt.B, opt.tile, std::string(v)=="V5"?10:1);
        }
        std::fclose(f); return 0;
    }
    run_machine(dir);
    {
        FILE* f = std::fopen((dir / "scan_b.csv").c_str(), "w"); std::fprintf(f, "nx,nv,m,B,rep,seconds\n");
        double best = 1e9; int b_opt = 16;
        for (int g = 0; g < 2; ++g) {
            int nx = g ? 801 : 401, nv = g ? 401 : 201; Problem prob(Params::set_a(), nx, nv, 200); std::vector<double> u;
            for (int B : {1, 2, 4, 8, 16, 32, 64}) {
                auto ts = time_runs([&]() { solve_v3(prob, u, Opts{B, 16, 1}); });
                for (int r = 0; r < 7; ++r) std::fprintf(f, "%d,%d,200,%d,%d,%.9f\n", nx, nv, B, r, ts[r]);
                double med = median_of(ts); if (g == 1 && med < best) { best = med; b_opt = B; }
            }
        }
        opt.B = b_opt; std::fclose(f);
    }
    {
        FILE* f = std::fopen((dir / "scan_tile.csv").c_str(), "w"); std::fprintf(f, "nx,nv,m,B,tile,rep,seconds\n");
        Problem prob(Params::set_a(), 801, 401, 200); std::vector<double> u; double best = 1e9; int t_opt = 16;
        for (int tile : {1, 2, 4, 8, 16, 32, 64}) {
            auto ts = time_runs([&]() { solve_v4(prob, u, Opts{opt.B, tile, 1}); });
            for (int r = 0; r < 7; ++r) std::fprintf(f, "801,401,200,%d,%d,%d,%.9f\n", opt.B, tile, r, ts[r]);
            double med = median_of(ts); if (med < best) { best = med; t_opt = tile; }
        }
        opt.tile = t_opt; std::fclose(f);
    }
    run_checks(dir, "checks.csv", opt);
    {
        FILE* f = std::fopen((dir / "bench.csv").c_str(), "w"); std::fprintf(f, "variant,nx,nv,m,rep,seconds,B,tile,threads\n");
        const int nxs[] = {101, 201, 401, 801}, nvs[] = {51, 101, 201, 401};
        for (int g = 0; g < 4; ++g) {
            Problem prob(Params::set_a(), nxs[g], nvs[g], 200); std::vector<double> u;
            for (const char* v : {"V0", "V1", "V2", "V3", "V4", "V5"}) {
                auto ts = time_runs([&]() { run_var(v, prob, opt, u); });
                for (int r = 0; r < 7; ++r) std::fprintf(f, "%s,%d,%d,200,%d,%.9f,%d,%d,%d\n", v, nxs[g], nvs[g], r, ts[r], opt.B, opt.tile, std::string(v)=="V5"?10:1);
            }
        }
        std::fclose(f);
    }
    {
        FILE* f = std::fopen((dir / "threads.csv").c_str(), "w"); std::fprintf(f, "nx,nv,m,threads,rep,seconds,B,tile\n");
        Problem prob(Params::set_a(), 801, 401, 200); std::vector<double> u;
        for (int T : {1, 2, 3, 4, 5, 6, 8, 10}) {
            auto ts = time_runs([&]() { solve_v5(prob, u, Opts{opt.B, opt.tile, T}); });
            for (int r = 0; r < 7; ++r) std::fprintf(f, "801,401,200,%d,%d,%.9f,%d,%d\n", T, r, ts[r], opt.B, opt.tile);
        }
        std::fclose(f);
    }
    {
        FILE* f = std::fopen((dir / "conv.csv").c_str(), "w");
        std::fprintf(f, "set,nx,nv,m,dx,dv,dt,pde_price,gp_ref,abs_err,lewis_ref,gp_put_ref,parity_err\n");
        GL16 gl; const int nxs[] = {101, 201, 401, 801}, nvs[] = {51, 101, 201, 401}, ms[] = {50, 100, 200, 400};
        for (const char* sn : {"A", "B"}) {
            Params p = (sn[0] == 'A') ? Params::set_a() : Params::set_b(); HestonPsi psi{p, p.S0, p.T, p.v0};
            double ref = gp_call(psi, p.S0, p.K, p.r, p.T, gl), lew = lewis_call(psi, p.S0, p.K, p.r, p.T, gl);
            double put = gp_put(psi, p.S0, p.K, p.r, p.T, gl), par = std::abs(lew - (put + p.S0 - p.K * std::exp(-p.r * p.T)));
            for (int l = 0; l < 4; ++l) {
                Problem prob(p, nxs[l], nvs[l], ms[l]); std::vector<double> u; solve_v5(prob, u, opt);
                double pr = price_at(u.data(), prob.g);
                std::fprintf(f, "%s,%d,%d,%d,%.12f,%.12f,%.12f,%.12f,%.12f,%.12e,%.12f,%.12f,%.3e\n",
                             sn, nxs[l], nvs[l], ms[l], prob.g.dx, prob.g.dv, prob.g.dt, pr, ref, std::abs(pr - ref), lew, put, par);
            }
        }
        std::fclose(f);
    }
    {
        FILE* f = std::fopen((dir / "surface.csv").c_str(), "w"); std::fprintf(f, "S,v,pde_price,gp_ref,abs_err\n");
        Params p = Params::set_a(); Problem prob(p, 401, 201, 200); std::vector<double> u; solve_v5(prob, u, opt); GL16 gl;
        std::vector<int> js = {8}; for (int k = 1; k <= 10; ++k) if (k * 10 != 8) js.push_back(k * 10);
        std::sort(js.begin(), js.end());
        for (int j : js) for (int i = 0; i < prob.g.nx; i += 4) {
            double S = std::exp(prob.g.x_min + i * prob.g.dx), v = j * prob.g.dv;
            if (S < 50.0 || S > 200.0) continue;
            double pr = u[size_t(j) * prob.g.nx + i], ref = gp_call(HestonPsi{p, S, p.T, v}, S, p.K, p.r, p.T, gl);
            std::fprintf(f, "%.8f,%.6f,%.10f,%.10f,%.3e\n", S, v, pr, ref, std::abs(pr - ref));
        }
        std::fclose(f);
    }
    return 0;
}
