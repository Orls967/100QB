#include "ising.hpp"
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

static double now_sec() {
    using clk = std::chrono::steady_clock;
    return std::chrono::duration<double>(clk::now().time_since_epoch()).count();
}

static std::vector<double> build_temp_grid() {
    std::vector<double> raw;
    for (int i = 0; i < 25; ++i) raw.push_back(1.8 + (3.4 - 1.8) * static_cast<double>(i) / 24.0);
    for (int i = 0; i < 12; ++i) raw.push_back(2.15 + (2.40 - 2.15) * static_cast<double>(i) / 11.0);
    std::sort(raw.begin(), raw.end());
    std::vector<double> out;
    for (double t : raw) {
        if (out.empty() || std::fabs(t - out.back()) > 1e-12) out.push_back(t);
    }
    return out;
}

static void compute_stats_series(const std::vector<double> &x, double &mean_out, double &err_out, double &tau_out) {
    int n = static_cast<int>(x.size());
    double sum = 0.0;
    for (double v : x) sum += v;
    double mean = sum / static_cast<double>(n);
    double var = 0.0;
    for (double v : x) {
        double d = v - mean;
        var += d * d;
    }
    var /= static_cast<double>(n);
    double tau = 0.5;
    if (var > 1e-30) {
        int max_lag = std::min(n / 4, 2000);
        for (int t = 1; t < max_lag; ++t) {
            double c = 0.0;
            for (int i = 0; i < n - t; ++i) c += (x[i] - mean) * (x[i + t] - mean);
            double rho = (c / static_cast<double>(n - t)) / var;
            if (rho <= 0.0) break;
            tau += rho;
            if (static_cast<double>(t) >= 6.0 * tau) break;
        }
    }
    double n_eff = std::max(1.0, static_cast<double>(n) / (2.0 * tau));
    mean_out = mean;
    err_out = std::sqrt(var / n_eff);
    tau_out = tau;
}

static void run_parity(const std::string &out_dir) {
    const int L = 8;
    const double beta = 0.44;
    const int steps = 3000;
    const double p_table[5] = {1.0, 1.0, 1.0, std::exp(-4.0 * beta), std::exp(-8.0 * beta)};
    const double p_add = 1.0 - std::exp(-2.0 * beta);
    const char *names[3] = {"metropolis", "wolff", "sw"};

    FILE *fsum = std::fopen((out_dir + "/parity_cpp_summary.csv").c_str(), "w");
    std::fprintf(fsum, "algo,final_fnv1a_hex,final_E,final_M\n");

    for (int algo = 0; algo < 3; ++algo) {
        ising::Lattice lat(L);
        lat.set_cold();
        ising::Xoshiro256pp rng(ising::compute_task_seed(L, 0, algo, 0));
        FILE *fp = std::fopen((out_dir + "/parity_cpp_" + names[algo] + ".csv").c_str(), "w");
        std::fprintf(fp, "step,E_sum,M_sum\n");
        int64_t E = 0, M = 0;
        for (int s = 1; s <= steps; ++s) {
            if (algo == 0) lat.metropolis_sweep_table(rng, p_table);
            else if (algo == 1) lat.wolff_step(rng, p_add);
            else lat.sw_step(rng, p_add);
            lat.energy_mag(E, M);
            std::fprintf(fp, "%d,%lld,%lld\n", s, static_cast<long long>(E), static_cast<long long>(M));
        }
        std::fclose(fp);
        uint64_t h = ising::fnv1a_64(lat.spins.data(), static_cast<size_t>(lat.N));
        std::fprintf(fsum, "%s,%016llx,%lld,%lld\n", names[algo], static_cast<unsigned long long>(h),
                     static_cast<long long>(E), static_cast<long long>(M));
    }
    std::fclose(fsum);
}

static void run_enum(const std::string &out_dir) {
    const double betas[5] = {0.2, 0.35, 0.4406867935097715, 0.5, 0.7};
    FILE *fp = std::fopen((out_dir + "/enum_cpp.csv").c_str(), "w");
    std::fprintf(fp, "L,beta,Z,e_mean,m2_mean\n");
    for (int L = 4; L <= 5; ++L) {
        for (double b : betas) {
            auto r = ising::enumerate_exact(L, b);
            std::fprintf(fp, "%d,%.17g,%.17g,%.17g,%.17g\n", L, r.beta, r.Z, r.e_mean, r.m2_mean);
        }
    }
    std::fclose(fp);
}

struct SimTask {
    int L, t_idx, algo, n_therm, n_meas;
    double T;
    double e_mean, e_err, abs_m_mean, abs_m_err, m2_mean, m2_err, m4_mean, m4_err;
    double binder, chi, cv, tau_e, tau_m;
    uint64_t final_hash;
};

static void execute_sim_task(SimTask &tk) {
    ising::Lattice lat(tk.L);
    ising::Xoshiro256pp rng(ising::compute_task_seed(tk.L, tk.t_idx, tk.algo, 0));
    lat.set_hot(rng);

    const double beta = 1.0 / tk.T;
    const double p_table[5] = {1.0, 1.0, 1.0, std::exp(-4.0 * beta), std::exp(-8.0 * beta)};
    const double p_add = 1.0 - std::exp(-2.0 * beta);

    for (int s = 0; s < tk.n_therm; ++s) {
        if (tk.algo == 0) lat.metropolis_sweep_table(rng, p_table);
        else if (tk.algo == 1) lat.wolff_step(rng, p_add);
        else lat.sw_step(rng, p_add);
    }

    const int M = tk.n_meas;
    std::vector<double> es(M), ams(M), m2s(M), m4s(M);
    double sum_e2 = 0.0;
    const double invN = 1.0 / static_cast<double>(lat.N);

    for (int s = 0; s < M; ++s) {
        if (tk.algo == 0) lat.metropolis_sweep_table(rng, p_table);
        else if (tk.algo == 1) lat.wolff_step(rng, p_add);
        else lat.sw_step(rng, p_add);
        int64_t Ei = 0, Mi = 0;
        lat.energy_mag(Ei, Mi);
        double e = static_cast<double>(Ei) * invN;
        double m = static_cast<double>(Mi) * invN;
        double m2 = m * m;
        es[s] = e;
        ams[s] = std::fabs(m);
        m2s[s] = m2;
        m4s[s] = m2 * m2;
        sum_e2 += e * e;
    }

    double dummy_t = 0.0;
    compute_stats_series(es, tk.e_mean, tk.e_err, tk.tau_e);
    compute_stats_series(ams, tk.abs_m_mean, tk.abs_m_err, tk.tau_m);
    compute_stats_series(m2s, tk.m2_mean, tk.m2_err, dummy_t);
    compute_stats_series(m4s, tk.m4_mean, tk.m4_err, dummy_t);

    tk.binder = 1.0 - tk.m4_mean / (3.0 * tk.m2_mean * tk.m2_mean);
    tk.chi = static_cast<double>(lat.N) * (tk.m2_mean - tk.abs_m_mean * tk.abs_m_mean) / tk.T;
    tk.cv = static_cast<double>(lat.N) * (sum_e2 / static_cast<double>(M) - tk.e_mean * tk.e_mean) / (tk.T * tk.T);
    tk.final_hash = ising::fnv1a_64(lat.spins.data(), static_cast<size_t>(lat.N));
}

static void run_parallel_tasks(std::vector<SimTask> &tasks, int n_threads) {
    if (n_threads <= 1) {
        for (auto &t : tasks) execute_sim_task(t);
        return;
    }
    std::atomic<size_t> next_idx{0};
    std::vector<std::thread> pool;
    for (int t = 0; t < n_threads; ++t) {
        pool.emplace_back([&]() {
            while (true) {
                size_t idx = next_idx.fetch_add(1, std::memory_order_relaxed);
                if (idx >= tasks.size()) break;
                execute_sim_task(tasks[idx]);
            }
        });
    }
    for (auto &th : pool) th.join();
}

static void run_threads_test(const std::string &out_dir, int n_threads) {
    double temps[4] = {2.0, 2.269185314213022, 2.5, 3.0};
    std::vector<SimTask> tasks;
    for (int ti = 0; ti < 4; ++ti) {
        for (int algo = 0; algo < 3; ++algo) {
            SimTask tk{};
            tk.L = 16; tk.t_idx = ti; tk.T = temps[ti]; tk.algo = algo;
            tk.n_therm = 500; tk.n_meas = 2000;
            tasks.push_back(tk);
        }
    }
    run_parallel_tasks(tasks, n_threads);
    FILE *fp = std::fopen((out_dir + "/threads_cpp_t" + std::to_string(n_threads) + ".csv").c_str(), "w");
    std::fprintf(fp, "L,t_idx,algo,e_mean,m2_mean,final_hash\n");
    for (const auto &tk : tasks) {
        std::fprintf(fp, "%d,%d,%d,%.17g,%.17g,%016llx\n", tk.L, tk.t_idx, tk.algo,
                     tk.e_mean, tk.m2_mean, static_cast<unsigned long long>(tk.final_hash));
    }
    std::fclose(fp);
}

static void run_full(const std::string &out_dir, int n_threads) {
    const int Ls[5] = {8, 16, 32, 64, 128};
    auto temps = build_temp_grid();
    std::vector<SimTask> tasks;
    for (int L : Ls) {
        for (size_t ti = 0; ti < temps.size(); ++ti) {
            for (int algo = 0; algo < 3; ++algo) {
                SimTask tk{};
                tk.L = L; tk.t_idx = static_cast<int>(ti); tk.T = temps[ti]; tk.algo = algo;
                if (algo == 0) {
                    tk.n_therm = (L <= 64) ? 16384 : 8192;
                    tk.n_meas = (L <= 64) ? 262144 : 65536;
                } else if (algo == 1) {
                    tk.n_therm = 2048;
                    tk.n_meas = (L <= 64) ? 131072 : 32768;
                } else {
                    tk.n_therm = 2048;
                    tk.n_meas = (L <= 64) ? 131072 : 16384;
                }
                tasks.push_back(tk);
            }
        }
    }

    double t0 = now_sec();
    run_parallel_tasks(tasks, n_threads);
    double dt = now_sec() - t0;
    std::printf("C++ full run (%zu tasks, %d threads) finished in %.2f s\n", tasks.size(), n_threads, dt);

    const char *names[3] = {"metropolis", "wolff", "sw"};
    FILE *fp = std::fopen((out_dir + "/full_cpp.csv").c_str(), "w");
    std::fprintf(fp, "L,t_idx,T,algo,n_therm,n_meas,e_mean,e_err,abs_m_mean,abs_m_err,m2_mean,m2_err,m4_mean,m4_err,binder,chi,cv,tau_e,tau_m,final_hash\n");
    for (const auto &tk : tasks) {
        std::fprintf(fp, "%d,%d,%.17g,%s,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%016llx\n",
                     tk.L, tk.t_idx, tk.T, names[tk.algo], tk.n_therm, tk.n_meas,
                     tk.e_mean, tk.e_err, tk.abs_m_mean, tk.abs_m_err,
                     tk.m2_mean, tk.m2_err, tk.m4_mean, tk.m4_err,
                     tk.binder, tk.chi, tk.cv, tk.tau_e, tk.tau_m, static_cast<unsigned long long>(tk.final_hash));
    }
    std::fclose(fp);
}

static void run_bench(const std::string &out_dir) {
    const int L = 32;
    const double T = 2.269185314213022;
    const double beta = 1.0 / T;
    const double p_table[5] = {1.0, 1.0, 1.0, std::exp(-4.0 * beta), std::exp(-8.0 * beta)};
    const double p_add = 1.0 - std::exp(-2.0 * beta);
    const int reps = 7;

    FILE *fp = std::fopen((out_dir + "/bench_cpp.csv").c_str(), "w");
    std::fprintf(fp, "lang,variant,rep,ns_per_spin_update\n");

    for (int r = 0; r < reps; ++r) {
        ising::Lattice lat(L);
        ising::Xoshiro256pp rng(ising::compute_task_seed(L, 0, 0, r));
        lat.set_hot(rng);
        int sweeps = 2000;
        double t0 = now_sec();
        for (int s = 0; s < sweeps; ++s) lat.metropolis_sweep_table(rng, p_table);
        double dt = now_sec() - t0;
        std::fprintf(fp, "cpp,metropolis_table,%d,%.17g\n", r, (dt * 1e9) / (static_cast<double>(sweeps) * lat.N));

        lat.set_hot(rng);
        t0 = now_sec();
        for (int s = 0; s < sweeps; ++s) lat.metropolis_sweep_computed(rng, p_table);
        dt = now_sec() - t0;
        std::fprintf(fp, "cpp,metropolis_computed,%d,%.17g\n", r, (dt * 1e9) / (static_cast<double>(sweeps) * lat.N));

        lat.set_hot(rng);
        t0 = now_sec();
        for (int s = 0; s < sweeps; ++s) ising::metropolis_sweep_constexpr<32>(lat.spins, rng, p_table);
        dt = now_sec() - t0;
        std::fprintf(fp, "cpp,metropolis_constexpr32,%d,%.17g\n", r, (dt * 1e9) / (static_cast<double>(sweeps) * lat.N));

        lat.set_hot(rng);
        int w_steps = 2000;
        for (int s = 0; s < 200; ++s) lat.wolff_step(rng, p_add);
        t0 = now_sec();
        for (int s = 0; s < w_steps; ++s) lat.wolff_step(rng, p_add);
        dt = now_sec() - t0;
        std::fprintf(fp, "cpp,wolff,%d,%.17g\n", r, (dt * 1e9) / (static_cast<double>(w_steps) * lat.N));

        lat.set_hot(rng);
        int sw_steps = 1000;
        t0 = now_sec();
        for (int s = 0; s < sw_steps; ++s) lat.sw_step(rng, p_add);
        dt = now_sec() - t0;
        std::fprintf(fp, "cpp,sw,%d,%.17g\n", r, (dt * 1e9) / (static_cast<double>(sw_steps) * lat.N));
    }
    std::fclose(fp);
}

int main(int argc, char **argv) {
    std::string mode = "parity";
    std::string out_dir = "data";
    int n_threads = 8;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--mode") == 0 && i + 1 < argc) mode = argv[++i];
        else if (std::strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) out_dir = argv[++i];
        else if (std::strcmp(argv[i], "--threads") == 0 && i + 1 < argc) n_threads = std::atoi(argv[++i]);
    }
    if (mode == "parity") run_parity(out_dir);
    else if (mode == "enum") run_enum(out_dir);
    else if (mode == "threads") run_threads_test(out_dir, n_threads);
    else if (mode == "full") run_full(out_dir, n_threads);
    else if (mode == "bench") run_bench(out_dir);
    return 0;
}
