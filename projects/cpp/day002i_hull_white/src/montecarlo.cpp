#include "hw/montecarlo.hpp"

#include "hw/rng.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

namespace hw {

namespace {

// Deterministic alpha(t) = f(0, t) + (sigma^2 / (2 a^2)) * (1 - exp(-a t))^2
double hw_alpha(const NelsonSiegelCurve& curve, const HwParams& p, double t) {
    double f0t = curve.inst_forward(t);
    double b = (1.0 - std::exp(-p.a * t)) / p.a;
    return f0t + 0.5 * p.sigma * p.sigma * b * b;
}

}  // namespace

McSwaptionSuite run_mc_swaptions(const NelsonSiegelCurve& curve,
                                 const HwParams& params,
                                 const std::vector<double>& expiries,
                                 double tenor_end_years,
                                 double fixed_rate,
                                 const std::vector<double>& df_test_maturities,
                                 int num_paths,
                                 int steps_per_year,
                                 std::uint64_t seed,
                                 unsigned int num_threads) {
    // Ensure even number of paths for antithetic pairs
    int n_pairs = std::max(1, num_paths / 2);
    int actual_paths = 2 * n_pairs;

    double max_t = tenor_end_years;
    for (double t0 : expiries) {
        max_t = std::max(max_t, t0);
    }
    for (double t : df_test_maturities) {
        max_t = std::max(max_t, t);
    }
    // We only need to simulate r(t) up to max(max_expiry, max_df_test)
    double sim_end = 0.0;
    for (double t0 : expiries) {
        sim_end = std::max(sim_end, t0);
    }
    for (double t : df_test_maturities) {
        sim_end = std::max(sim_end, t);
    }

    int total_steps = static_cast<int>(std::llround(sim_end * static_cast<double>(steps_per_year)));
    double dt = 1.0 / static_cast<double>(steps_per_year);

    double exp_adt = std::exp(-params.a * dt);
    double ou_sd = params.sigma * std::sqrt((1.0 - std::exp(-2.0 * params.a * dt)) / (2.0 * params.a));

    std::vector<double> alpha_grid(static_cast<std::size_t>(total_steps + 1));
    for (int k = 0; k <= total_steps; ++k) {
        alpha_grid[static_cast<std::size_t>(k)] = hw_alpha(curve, params, static_cast<double>(k) * dt);
    }

    std::vector<int> expiry_steps(expiries.size());
    for (std::size_t i = 0; i < expiries.size(); ++i) {
        expiry_steps[i] = static_cast<int>(std::llround(expiries[i] * static_cast<double>(steps_per_year)));
    }

    std::vector<int> df_steps(df_test_maturities.size());
    for (std::size_t i = 0; i < df_test_maturities.size(); ++i) {
        df_steps[i] = static_cast<int>(std::llround(df_test_maturities[i] * static_cast<double>(steps_per_year)));
    }

    // Precompute bond coefficients A(T0, Tn), B(T0, Tn) for each expiry
    int tenor_end_int = static_cast<int>(std::llround(tenor_end_years));
    struct SwapPrep {
        double t0{0.0};
        int start_yr{0};
        std::vector<double> cf{};
        std::vector<double> ln_a{};
        std::vector<double> b{};
        std::vector<double> market_df{};
    };
    std::vector<SwapPrep> preps(expiries.size());
    for (std::size_t i = 0; i < expiries.size(); ++i) {
        preps[i].t0 = expiries[i];
        preps[i].start_yr = static_cast<int>(std::llround(expiries[i]));
        for (int yr = preps[i].start_yr + 1; yr <= tenor_end_int; ++yr) {
            double tn = static_cast<double>(yr);
            double c = (yr == tenor_end_int) ? (1.0 + fixed_rate) : fixed_rate;
            preps[i].cf.push_back(c);
            preps[i].ln_a.push_back(hw_ln_a(curve, params, preps[i].t0, tn));
            preps[i].b.push_back(hw_b(params.a, preps[i].t0, tn));
            preps[i].market_df.push_back(curve.discount(tn));
        }
    }

    if (num_threads == 0) {
        num_threads = 1;
    }
    num_threads = std::min<unsigned int>(num_threads, static_cast<unsigned int>(n_pairs));

    struct ThreadAccum {
        std::vector<double> sum_raw;
        std::vector<double> sum_raw2;
        std::vector<double> sum_cv_z;
        std::vector<double> sum_cv_z2;
        std::vector<double> sum_raw_cv_z;
        std::vector<double> sum_df;
        std::vector<double> sum_df2;
    };
    std::vector<ThreadAccum> accums(num_threads);
    for (auto& ac : accums) {
        ac.sum_raw.assign(expiries.size(), 0.0);
        ac.sum_raw2.assign(expiries.size(), 0.0);
        ac.sum_cv_z.assign(expiries.size(), 0.0);
        ac.sum_cv_z2.assign(expiries.size(), 0.0);
        ac.sum_raw_cv_z.assign(expiries.size(), 0.0);
        ac.sum_df.assign(df_test_maturities.size(), 0.0);
        ac.sum_df2.assign(df_test_maturities.size(), 0.0);
    }

    std::vector<std::thread> workers;
    workers.reserve(num_threads);
    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        workers.emplace_back([&, tid]() {
            auto& ac = accums[tid];
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

                std::size_t next_exp = 0;
                std::size_t next_df = 0;

                std::vector<double> raw_pair(expiries.size(), 0.0);
                std::vector<double> cv_pair(expiries.size(), 0.0);
                std::vector<double> df_pair(df_test_maturities.size(), 0.0);

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

                    while (next_df < df_steps.size() && df_steps[next_df] == k) {
                        double d_p = std::exp(-int_p);
                        double d_m = std::exp(-int_m);
                        df_pair[next_df] = 0.5 * (d_p + d_m);
                        ++next_df;
                    }

                    while (next_exp < expiry_steps.size() && expiry_steps[next_exp] == k) {
                        double d_p = std::exp(-int_p);
                        double d_m = std::exp(-int_m);
                        const auto& sp = preps[next_exp];

                        auto eval_path = [&](double d_t0, double r_t0, double& raw_out, double& cv_out) {
                            double pv_fixed = 0.0;
                            double cv_term = 0.0;
                            for (std::size_t j = 0; j < sp.cf.size(); ++j) {
                                double p_t0_tn = std::exp(sp.ln_a[j] - sp.b[j] * r_t0);
                                pv_fixed += sp.cf[j] * p_t0_tn;
                                cv_term += sp.cf[j] * (d_t0 * p_t0_tn - sp.market_df[j]);
                            }
                            double swap_val = 1.0 - pv_fixed;
                            double payoff = d_t0 * std::max(swap_val, 0.0);
                            raw_out = payoff;
                            cv_out = cv_term;
                        };

                        double r1 = 0.0;
                        double c1 = 0.0;
                        double r2 = 0.0;
                        double c2 = 0.0;
                        eval_path(d_p, r_curr_p, r1, c1);
                        eval_path(d_m, r_curr_m, r2, c2);

                        raw_pair[next_exp] = 0.5 * (r1 + r2);
                        cv_pair[next_exp] = 0.5 * (c1 + c2);
                        ++next_exp;
                    }
                }

                for (std::size_t i = 0; i < expiries.size(); ++i) {
                    ac.sum_raw[i] += raw_pair[i];
                    ac.sum_raw2[i] += raw_pair[i] * raw_pair[i];
                    ac.sum_cv_z[i] += cv_pair[i];
                    ac.sum_cv_z2[i] += cv_pair[i] * cv_pair[i];
                    ac.sum_raw_cv_z[i] += raw_pair[i] * cv_pair[i];
                }
                for (std::size_t i = 0; i < df_test_maturities.size(); ++i) {
                    ac.sum_df[i] += df_pair[i];
                    ac.sum_df2[i] += df_pair[i] * df_pair[i];
                }
            }
        });
    }
    for (auto& w : workers) {
        w.join();
    }

    McSwaptionSuite suite;
    double n_p_d = static_cast<double>(n_pairs);
    for (std::size_t i = 0; i < expiries.size(); ++i) {
        double s_raw = 0.0;
        double s_raw2 = 0.0;
        double s_z = 0.0;
        double s_z2 = 0.0;
        double s_rz = 0.0;
        for (const auto& ac : accums) {
            s_raw += ac.sum_raw[i];
            s_raw2 += ac.sum_raw2[i];
            s_z += ac.sum_cv_z[i];
            s_z2 += ac.sum_cv_z2[i];
            s_rz += ac.sum_raw_cv_z[i];
        }
        double mean_raw = s_raw / n_p_d;
        double var_raw = std::max(0.0, s_raw2 / n_p_d - mean_raw * mean_raw);
        double se_raw = std::sqrt(var_raw / n_p_d);

        double mean_z = s_z / n_p_d;
        double var_z = std::max(1e-30, s_z2 / n_p_d - mean_z * mean_z);
        double cov_rz = s_rz / n_p_d - mean_raw * mean_z;
        double beta = cov_rz / var_z;

        double mean_cv = mean_raw - beta * mean_z;
        double var_cv = std::max(0.0, var_raw - beta * cov_rz);
        double se_cv = std::sqrt(var_cv / n_p_d);

        suite.swaptions.push_back(McSwaptionResult{expiries[i], mean_cv, se_cv, mean_raw, se_raw});
    }

    (void)actual_paths;
    for (std::size_t i = 0; i < df_test_maturities.size(); ++i) {
        double s_df = 0.0;
        double s_df2 = 0.0;
        for (const auto& ac : accums) {
            s_df += ac.sum_df[i];
            s_df2 += ac.sum_df2[i];
        }
        double mean_df = s_df / n_p_d;
        double var_df = std::max(0.0, s_df2 / n_p_d - mean_df * mean_df);
        double se_df = std::sqrt(var_df / n_p_d);
        double mkt = curve.discount(df_test_maturities[i]);
        suite.martingale_tests.push_back(MartingaleTestResult{
            df_test_maturities[i], mkt, mean_df, se_df, std::abs(mean_df - mkt)});
    }

    return suite;
}

}  // namespace hw
