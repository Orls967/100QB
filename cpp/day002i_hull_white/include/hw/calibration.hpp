#pragma once

#include "hw/curve.hpp"
#include "hw/hull_white.hpp"

#include <cstdint>
#include <vector>

namespace hw {

struct SwaptionQuote {
    double expiry_years{1.0};
    double tenor_end_years{10.0};
    double fixed_rate{0.045};
    double market_price{0.0};
};

struct CalibrationResult {
    HwParams fitted{};
    double price_rmse{0.0};
    double max_rel_error{0.0};
    int iterations{0};
    bool converged{false};
};

struct NoisyCalibrationSummary {
    int replications{0};
    double rel_noise_sd{0.0};
    double mean_a{0.0};
    double sd_a{0.0};
    double mean_sigma{0.0};
    double sd_sigma{0.0};
    double corr_a_sigma{0.0};
    double jacobian_cond_number{0.0};
    std::vector<HwParams> samples{};
};

/// Generates coterminal 10y payer swaption quotes under the true parameters.
std::vector<SwaptionQuote> generate_coterminal_quotes(const NelsonSiegelCurve& curve,
                                                      const HwParams& true_params,
                                                      const std::vector<double>& expiries,
                                                      double tenor_end_years,
                                                      double fixed_rate);

/// Calibrates (a, sigma) by minimizing squared price residuals via Nelder-Mead in log-space.
CalibrationResult calibrate_swaptions(const NelsonSiegelCurve& curve,
                                      const std::vector<SwaptionQuote>& quotes,
                                      const HwParams& initial_guess);

/// Computes the 2x2 condition number of J^T J (and sqrt for J) via finite differences.
double calibration_jacobian_condition(const NelsonSiegelCurve& curve,
                                      const std::vector<SwaptionQuote>& quotes,
                                      const HwParams& params);

/// Runs R noisy calibration replications in parallel using std::thread.
NoisyCalibrationSummary run_noisy_calibrations(const NelsonSiegelCurve& curve,
                                               const std::vector<SwaptionQuote>& true_quotes,
                                               const HwParams& initial_guess,
                                               double rel_noise_sd,
                                               int replications,
                                               std::uint64_t seed,
                                               unsigned int num_threads);

}  // namespace hw
