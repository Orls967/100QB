#include "hw/calibration.hpp"

#include "hw/rng.hpp"
#include "hw/solvers.hpp"

#include <algorithm>
#include <cmath>
#include <thread>
#include <vector>

namespace hw {

std::vector<SwaptionQuote> generate_coterminal_quotes(const NelsonSiegelCurve& curve,
                                                      const HwParams& true_params,
                                                      const std::vector<double>& expiries,
                                                      double tenor_end_years,
                                                      double fixed_rate) {
    std::vector<SwaptionQuote> quotes;
    quotes.reserve(expiries.size());
    for (double t0 : expiries) {
        double p = hw_payer_swaption_jamshidian(curve, true_params, t0, tenor_end_years, fixed_rate);
        quotes.push_back(SwaptionQuote{t0, tenor_end_years, fixed_rate, p});
    }
    return quotes;
}

CalibrationResult calibrate_swaptions(const NelsonSiegelCurve& curve,
                                      const std::vector<SwaptionQuote>& quotes,
                                      const HwParams& initial_guess) {
    auto objective = [&](const std::vector<double>& u) -> double {
        double a = std::exp(u[0]);
        double sigma = std::exp(u[1]);
        if (a < 1e-4 || a > 2.0 || sigma < 1e-5 || sigma > 0.25) {
            return 1e12;
        }
        HwParams p{a, sigma};
        double sse = 0.0;
        for (const auto& q : quotes) {
            double model_p = hw_payer_swaption_jamshidian(
                curve, p, q.expiry_years, q.tenor_end_years, q.fixed_rate);
            double diff = model_p - q.market_price;
            sse += diff * diff;
        }
        return sse;
    };

    std::vector<double> u0{std::log(initial_guess.a), std::log(initial_guess.sigma)};
    auto u_fit = solvers::nelder_mead(objective, u0, std::vector<double>{0.2, 0.2}, 1e-13, 800);

    CalibrationResult res;
    res.fitted = HwParams{std::exp(u_fit[0]), std::exp(u_fit[1])};
    res.iterations = 0;
    res.converged = true;

    double sse = 0.0;
    double max_rel = 0.0;
    for (const auto& q : quotes) {
        double model_p = hw_payer_swaption_jamshidian(
            curve, res.fitted, q.expiry_years, q.tenor_end_years, q.fixed_rate);
        double diff = model_p - q.market_price;
        sse += diff * diff;
        if (q.market_price > 1e-15) {
            max_rel = std::max(max_rel, std::abs(diff) / q.market_price);
        }
    }
    res.price_rmse = std::sqrt(sse / static_cast<double>(quotes.size()));
    res.max_rel_error = max_rel;
    return res;
}

double calibration_jacobian_condition(const NelsonSiegelCurve& curve,
                                      const std::vector<SwaptionQuote>& quotes,
                                      const HwParams& params) {
    const double eps_a = 1e-6 * std::max(1.0, params.a);
    const double eps_s = 1e-6 * std::max(1e-2, params.sigma);

    double jtj00 = 0.0;
    double jtj01 = 0.0;
    double jtj11 = 0.0;

    for (const auto& q : quotes) {
        double pa_p = hw_payer_swaption_jamshidian(
            curve, HwParams{params.a + eps_a, params.sigma}, q.expiry_years, q.tenor_end_years, q.fixed_rate);
        double pa_m = hw_payer_swaption_jamshidian(
            curve, HwParams{params.a - eps_a, params.sigma}, q.expiry_years, q.tenor_end_years, q.fixed_rate);
        double ps_p = hw_payer_swaption_jamshidian(
            curve, HwParams{params.a, params.sigma + eps_s}, q.expiry_years, q.tenor_end_years, q.fixed_rate);
        double ps_m = hw_payer_swaption_jamshidian(
            curve, HwParams{params.a, params.sigma - eps_s}, q.expiry_years, q.tenor_end_years, q.fixed_rate);

        double da = (pa_p - pa_m) / (2.0 * eps_a);
        double ds = (ps_p - ps_m) / (2.0 * eps_s);

        jtj00 += da * da;
        jtj01 += da * ds;
        jtj11 += ds * ds;
    }

    double tr = jtj00 + jtj11;
    double diff = jtj00 - jtj11;
    double disc = std::sqrt(diff * diff + 4.0 * jtj01 * jtj01);
    double lam_max = 0.5 * (tr + disc);
    double lam_min = 0.5 * (tr - disc);
    if (lam_min <= 0.0) {
        return 1e18;
    }
    return std::sqrt(lam_max / lam_min);
}

NoisyCalibrationSummary run_noisy_calibrations(const NelsonSiegelCurve& curve,
                                               const std::vector<SwaptionQuote>& true_quotes,
                                               const HwParams& initial_guess,
                                               double rel_noise_sd,
                                               int replications,
                                               std::uint64_t seed,
                                               unsigned int num_threads) {
    NoisyCalibrationSummary out;
    out.replications = replications;
    out.rel_noise_sd = rel_noise_sd;
    out.samples.resize(static_cast<std::size_t>(replications));

    if (num_threads == 0) {
        num_threads = 1;
    }
    num_threads = std::min<unsigned int>(num_threads, static_cast<unsigned int>(replications));

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        workers.emplace_back([&, tid]() {
            for (int r = static_cast<int>(tid); r < replications; r += static_cast<int>(num_threads)) {
                Xoshiro256StarStar rng(derive_seed(seed, static_cast<std::uint64_t>(r)));
                std::vector<SwaptionQuote> noisy = true_quotes;
                for (auto& q : noisy) {
                    double z = rng.next_normal();
                    double factor = std::max(0.01, 1.0 + rel_noise_sd * z);
                    q.market_price *= factor;
                }
                auto res = calibrate_swaptions(curve, noisy, initial_guess);
                out.samples[static_cast<std::size_t>(r)] = res.fitted;
            }
        });
    }
    for (auto& w : workers) {
        w.join();
    }

    double sum_a = 0.0;
    double sum_s = 0.0;
    for (const auto& s : out.samples) {
        sum_a += s.a;
        sum_s += s.sigma;
    }
    double n_d = static_cast<double>(replications);
    out.mean_a = sum_a / n_d;
    out.mean_sigma = sum_s / n_d;

    double var_a = 0.0;
    double var_s = 0.0;
    double cov = 0.0;
    for (const auto& s : out.samples) {
        double da = s.a - out.mean_a;
        double ds = s.sigma - out.mean_sigma;
        var_a += da * da;
        var_s += ds * ds;
        cov += da * ds;
    }
    double denom = std::max(1.0, n_d - 1.0);
    out.sd_a = std::sqrt(var_a / denom);
    out.sd_sigma = std::sqrt(var_s / denom);
    out.corr_a_sigma = (out.sd_a > 0.0 && out.sd_sigma > 0.0)
                           ? (cov / denom) / (out.sd_a * out.sd_sigma)
                           : 0.0;

    if (!true_quotes.empty()) {
        out.jacobian_cond_number = calibration_jacobian_condition(
            curve, true_quotes, HwParams{out.mean_a, out.mean_sigma});
    }
    return out;
}

}  // namespace hw
