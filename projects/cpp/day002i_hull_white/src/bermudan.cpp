#include "hw/bermudan.hpp"

#include "hw/rng.hpp"
#include "hw/solvers.hpp"
#include "hw/trinomial_tree.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <thread>
#include <vector>

namespace hw {

namespace {

double hw_alpha_b(const NelsonSiegelCurve& curve, const HwParams& p, double t) {
    double f0t = curve.inst_forward(t);
    double b = (1.0 - std::exp(-p.a * t)) / p.a;
    return f0t + 0.5 * p.sigma * p.sigma * b * b;
}

struct ExDatePrep {
    double t0{0.0};
    int year{0};
    std::vector<double> cf{};
    std::vector<double> ln_a{};
    std::vector<double> b{};
};

std::vector<ExDatePrep> build_ex_preps(const NelsonSiegelCurve& curve,
                                       const HwParams& params,
                                       int first_yr,
                                       int last_yr,
                                       int tenor_end_yr,
                                       double fixed_rate) {
    std::vector<ExDatePrep> preps;
    for (int yr = first_yr; yr <= last_yr; ++yr) {
        ExDatePrep ep;
        ep.t0 = static_cast<double>(yr);
        ep.year = yr;
        for (int m = yr + 1; m <= tenor_end_yr; ++m) {
            double tn = static_cast<double>(m);
            double c = (m == tenor_end_yr) ? (1.0 + fixed_rate) : fixed_rate;
            ep.cf.push_back(c);
            ep.ln_a.push_back(hw_ln_a(curve, params, ep.t0, tn));
            ep.b.push_back(hw_b(params.a, ep.t0, tn));
        }
        preps.push_back(std::move(ep));
    }
    return preps;
}

double eval_exercise_value(const ExDatePrep& ep, double r_t0) {
    double pv = 0.0;
    for (std::size_t j = 0; j < ep.cf.size(); ++j) {
        pv += ep.cf[j] * std::exp(ep.ln_a[j] - ep.b[j] * r_t0);
    }
    return 1.0 - pv;
}

struct SimBatch {
    std::vector<double> x_at_ex{};      // [path * n_ex + ex_idx]
    std::vector<double> d_at_ex{};      // [path * n_ex + ex_idx]
    std::vector<double> ex_val{};       // [path * n_ex + ex_idx]
};

SimBatch simulate_bermudan_paths(const NelsonSiegelCurve& curve,
                                 const HwParams& params,
                                 const std::vector<ExDatePrep>& preps,
                                 int num_paths,
                                 int steps_per_year,
                                 std::uint64_t seed,
                                 unsigned int num_threads) {
    int n_ex = static_cast<int>(preps.size());
    int n_pairs = std::max(1, num_paths / 2);
    int actual_paths = 2 * n_pairs;

    int max_yr = preps.back().year;
    int total_steps = max_yr * steps_per_year;
    double dt = 1.0 / static_cast<double>(steps_per_year);
    double exp_adt = std::exp(-params.a * dt);
    double ou_sd = params.sigma * std::sqrt((1.0 - std::exp(-2.0 * params.a * dt)) / (2.0 * params.a));

    std::vector<double> alpha_grid(static_cast<std::size_t>(total_steps + 1));
    for (int k = 0; k <= total_steps; ++k) {
        alpha_grid[static_cast<std::size_t>(k)] = hw_alpha_b(curve, params, static_cast<double>(k) * dt);
    }

    SimBatch batch;
    batch.x_at_ex.resize(static_cast<std::size_t>(actual_paths) * static_cast<std::size_t>(n_ex));
    batch.d_at_ex.resize(static_cast<std::size_t>(actual_paths) * static_cast<std::size_t>(n_ex));
    batch.ex_val.resize(static_cast<std::size_t>(actual_paths) * static_cast<std::size_t>(n_ex));

    if (num_threads == 0) {
        num_threads = 1;
    }
    num_threads = std::min<unsigned int>(num_threads, static_cast<unsigned int>(n_pairs));

    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        workers.emplace_back([&, tid]() {
            int pair_begin = static_cast<int>(
                (static_cast<std::int64_t>(n_pairs) * tid) / num_threads);
            int pair_end = static_cast<int>(
                (static_cast<std::int64_t>(n_pairs) * (tid + 1)) / num_threads);

            for (int pair_idx = pair_begin; pair_idx < pair_end; ++pair_idx) {
                Xoshiro256StarStar rng(derive_seed(seed, static_cast<std::uint64_t>(pair_idx)));
                double x_p = 0.0;
                double x_m = 0.0;
                double r_prev_p = alpha_grid[0];
                double r_prev_m = alpha_grid[0];
                double int_p = 0.0;
                double int_m = 0.0;

                int next_ex = 0;
                int p0 = 2 * pair_idx;
                int p1 = 2 * pair_idx + 1;

                for (int k = 1; k <= total_steps; ++k) {
                    double z = rng.next_normal();
                    x_p = x_p * exp_adt + ou_sd * z;
                    x_m = x_m * exp_adt - ou_sd * z;

                    double r_curr_p = alpha_grid[static_cast<std::size_t>(k)] + x_p;
                    double r_curr_m = alpha_grid[static_cast<std::size_t>(k)] + x_m;

                    int_p += 0.5 * (r_prev_p + r_curr_p) * dt;
                    int_m += 0.5 * (r_prev_m + r_curr_m) * dt;

                    r_prev_p = r_curr_p;
                    r_prev_m = r_curr_m;

                    if (next_ex < n_ex && k == preps[static_cast<std::size_t>(next_ex)].year * steps_per_year) {
                        std::size_t idx0 = static_cast<std::size_t>(p0) * static_cast<std::size_t>(n_ex)
                                           + static_cast<std::size_t>(next_ex);
                        std::size_t idx1 = static_cast<std::size_t>(p1) * static_cast<std::size_t>(n_ex)
                                           + static_cast<std::size_t>(next_ex);

                        batch.x_at_ex[idx0] = x_p;
                        batch.d_at_ex[idx0] = std::exp(-int_p);
                        batch.ex_val[idx0] = std::max(eval_exercise_value(preps[static_cast<std::size_t>(next_ex)], r_curr_p), 0.0);

                        batch.x_at_ex[idx1] = x_m;
                        batch.d_at_ex[idx1] = std::exp(-int_m);
                        batch.ex_val[idx1] = std::max(eval_exercise_value(preps[static_cast<std::size_t>(next_ex)], r_curr_m), 0.0);
                        ++next_ex;
                    }
                }
            }
        });
    }
    for (auto& w : workers) {
        w.join();
    }
    return batch;
}

double eval_poly3(const std::array<double, 4>& beta, double x) {
    return beta[0] + x * (beta[1] + x * (beta[2] + x * beta[3]));
}

}  // namespace

BermudanLsmResult price_bermudan_lsm(const NelsonSiegelCurve& curve,
                                     const HwParams& params,
                                     int first_exercise_yr,
                                     int last_exercise_yr,
                                     int tenor_end_yr,
                                     double fixed_rate,
                                     int paths_per_pass,
                                     int steps_per_year,
                                     std::uint64_t seed,
                                     unsigned int num_threads) {
    auto preps = build_ex_preps(curve, params, first_exercise_yr, last_exercise_yr, tenor_end_yr, fixed_rate);
    int n_ex = static_cast<int>(preps.size());
    int n_pairs = std::max(1, paths_per_pass / 2);
    int actual_paths = 2 * n_pairs;

    // Pass 1: training paths
    SimBatch train = simulate_bermudan_paths(
        curve, params, preps, actual_paths, steps_per_year, seed, num_threads);

    std::vector<std::array<double, 4>> betas(static_cast<std::size_t>(n_ex));
    std::vector<double> cf_val(static_cast<std::size_t>(actual_paths));
    std::vector<int> cf_ex_idx(static_cast<std::size_t>(actual_paths), n_ex - 1);

    for (int p = 0; p < actual_paths; ++p) {
        std::size_t idx = static_cast<std::size_t>(p) * static_cast<std::size_t>(n_ex)
                          + static_cast<std::size_t>(n_ex - 1);
        cf_val[static_cast<std::size_t>(p)] = train.ex_val[idx];
        cf_ex_idx[static_cast<std::size_t>(p)] = n_ex - 1;
    }

    for (int t = n_ex - 2; t >= 0; --t) {
        // Fit degree-3 polynomial on ITM paths (or all paths if few ITM)
        std::vector<std::vector<double>> xtx(4, std::vector<double>(4, 0.0));
        std::vector<double> xty(4, 0.0);
        int itm_count = 0;

        for (int p = 0; p < actual_paths; ++p) {
            std::size_t idx = static_cast<std::size_t>(p) * static_cast<std::size_t>(n_ex)
                              + static_cast<std::size_t>(t);
            double ex = train.ex_val[idx];
            if (ex <= 0.0) {
                continue;
            }
            ++itm_count;
            double x = train.x_at_ex[idx] * 100.0;  // scale x (~0.01) for numerical conditioning
            std::size_t fut_idx = static_cast<std::size_t>(p) * static_cast<std::size_t>(n_ex)
                                  + static_cast<std::size_t>(cf_ex_idx[static_cast<std::size_t>(p)]);
            double y = cf_val[static_cast<std::size_t>(p)] * (train.d_at_ex[fut_idx] / train.d_at_ex[idx]);

            std::array<double, 4> phi{1.0, x, x * x, x * x * x};
            for (int r = 0; r < 4; ++r) {
                xty[static_cast<std::size_t>(r)] += phi[static_cast<std::size_t>(r)] * y;
                for (int c = 0; c < 4; ++c) {
                    xtx[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] +=
                        phi[static_cast<std::size_t>(r)] * phi[static_cast<std::size_t>(c)];
                }
            }
        }

        std::array<double, 4> beta{0.0, 0.0, 0.0, 0.0};
        if (itm_count >= 10) {
            for (int r = 0; r < 4; ++r) {
                xtx[static_cast<std::size_t>(r)][static_cast<std::size_t>(r)] += 1e-12;
            }
            auto solved = solvers::solve_linear(xtx, xty);
            for (int r = 0; r < 4; ++r) {
                beta[static_cast<std::size_t>(r)] = solved[static_cast<std::size_t>(r)];
            }
        }
        betas[static_cast<std::size_t>(t)] = beta;

        for (int p = 0; p < actual_paths; ++p) {
            std::size_t idx = static_cast<std::size_t>(p) * static_cast<std::size_t>(n_ex)
                              + static_cast<std::size_t>(t);
            double ex = train.ex_val[idx];
            if (ex > 0.0) {
                double x = train.x_at_ex[idx] * 100.0;
                double cont = eval_poly3(beta, x);
                if (ex >= cont) {
                    cf_val[static_cast<std::size_t>(p)] = ex;
                    cf_ex_idx[static_cast<std::size_t>(p)] = t;
                }
            }
        }
    }

    // Pass 2: independent evaluation paths (no foresight bias)
    SimBatch eval = simulate_bermudan_paths(
        curve, params, preps, actual_paths, steps_per_year, seed + 1'000'003ULL, num_threads);

    double sum_pv = 0.0;
    double sum_pv2 = 0.0;
    for (int pair_idx = 0; pair_idx < n_pairs; ++pair_idx) {
        double pair_sum = 0.0;
        for (int sub = 0; sub < 2; ++sub) {
            int p = 2 * pair_idx + sub;
            double pv = 0.0;
            for (int t = 0; t < n_ex; ++t) {
                std::size_t idx = static_cast<std::size_t>(p) * static_cast<std::size_t>(n_ex)
                                  + static_cast<std::size_t>(t);
                double ex = eval.ex_val[idx];
                if (t == n_ex - 1) {
                    pv = eval.d_at_ex[idx] * ex;
                    break;
                }
                if (ex > 0.0) {
                    double x = eval.x_at_ex[idx] * 100.0;
                    double cont = eval_poly3(betas[static_cast<std::size_t>(t)], x);
                    if (ex >= cont) {
                        pv = eval.d_at_ex[idx] * ex;
                        break;
                    }
                }
            }
            pair_sum += pv;
        }
        double pair_mean = 0.5 * pair_sum;
        sum_pv += pair_mean;
        sum_pv2 += pair_mean * pair_mean;
    }

    double n_p_d = static_cast<double>(n_pairs);
    double price = sum_pv / n_p_d;
    double var = std::max(0.0, sum_pv2 / n_p_d - price * price);
    double se = std::sqrt(var / n_p_d);

    return BermudanLsmResult{price, se, actual_paths, actual_paths};
}

Dv01Result dv01_euro_jamshidian(const NelsonSiegelCurve& curve,
                                const HwParams& params,
                                double expiry_years,
                                double tenor_end_years,
                                double fixed_rate,
                                double shift_bp) {
    double dz = shift_bp * 1e-4;
    NelsonSiegelCurve c_up = curve;
    c_up.parallel_shift += dz;
    NelsonSiegelCurve c_dn = curve;
    c_dn.parallel_shift -= dz;

    double base = hw_payer_swaption_jamshidian(curve, params, expiry_years, tenor_end_years, fixed_rate);
    double up = hw_payer_swaption_jamshidian(c_up, params, expiry_years, tenor_end_years, fixed_rate);
    double dn = hw_payer_swaption_jamshidian(c_dn, params, expiry_years, tenor_end_years, fixed_rate);
    return Dv01Result{base, up, dn, 0.5 * (up - dn)};
}

Dv01Result dv01_euro_tree(const NelsonSiegelCurve& curve,
                          const HwParams& params,
                          double expiry_years,
                          double tenor_end_years,
                          double fixed_rate,
                          int steps_per_year,
                          double shift_bp) {
    double dz = shift_bp * 1e-4;
    NelsonSiegelCurve c_up = curve;
    c_up.parallel_shift += dz;
    NelsonSiegelCurve c_dn = curve;
    c_dn.parallel_shift -= dz;

    int max_yr = static_cast<int>(std::llround(tenor_end_years));
    TrinomialTree t_base(curve, params, steps_per_year, max_yr);
    TrinomialTree t_up(c_up, params, steps_per_year, max_yr);
    TrinomialTree t_dn(c_dn, params, steps_per_year, max_yr);

    double base = t_base.price_payer_swaption(expiry_years, tenor_end_years, fixed_rate);
    double up = t_up.price_payer_swaption(expiry_years, tenor_end_years, fixed_rate);
    double dn = t_dn.price_payer_swaption(expiry_years, tenor_end_years, fixed_rate);
    return Dv01Result{base, up, dn, 0.5 * (up - dn)};
}

Dv01Result dv01_bermudan_tree(const NelsonSiegelCurve& curve,
                              const HwParams& params,
                              int first_exercise_yr,
                              int last_exercise_yr,
                              int tenor_end_yr,
                              double fixed_rate,
                              int steps_per_year,
                              double shift_bp) {
    double dz = shift_bp * 1e-4;
    NelsonSiegelCurve c_up = curve;
    c_up.parallel_shift += dz;
    NelsonSiegelCurve c_dn = curve;
    c_dn.parallel_shift -= dz;

    TrinomialTree t_base(curve, params, steps_per_year, tenor_end_yr);
    TrinomialTree t_up(c_up, params, steps_per_year, tenor_end_yr);
    TrinomialTree t_dn(c_dn, params, steps_per_year, tenor_end_yr);

    (void)first_exercise_yr;
    (void)last_exercise_yr;
    double base = t_base.price_bermudan_payer(tenor_end_yr, fixed_rate);
    double up = t_up.price_bermudan_payer(tenor_end_yr, fixed_rate);
    double dn = t_dn.price_bermudan_payer(tenor_end_yr, fixed_rate);
    return Dv01Result{base, up, dn, 0.5 * (up - dn)};
}

}  // namespace hw
