#pragma once

#include "hw/curve.hpp"
#include "hw/hull_white.hpp"

#include <cstdint>
#include <vector>

namespace hw {

struct McSwaptionResult {
    double expiry_years{0.0};
    double price{0.0};
    double std_error{0.0};
    double raw_price{0.0};
    double raw_std_error{0.0};
};

struct MartingaleTestResult {
    double maturity_years{0.0};
    double market_df{0.0};
    double mc_df{0.0};
    double std_error{0.0};
    double abs_error{0.0};
};

struct McSwaptionSuite {
    std::vector<McSwaptionResult> swaptions{};
    std::vector<MartingaleTestResult> martingale_tests{};
};

/// Runs antithetic OU Monte Carlo (dt = 1/52, trapezoidal discounting) for coterminal
/// payer swaptions with martingale control variate, and checks E[D(T)] = P(0, T).
McSwaptionSuite run_mc_swaptions(const NelsonSiegelCurve& curve,
                                 const HwParams& params,
                                 const std::vector<double>& expiries,
                                 double tenor_end_years,
                                 double fixed_rate,
                                 const std::vector<double>& df_test_maturities,
                                 int num_paths,
                                 int steps_per_year,
                                 std::uint64_t seed,
                                 unsigned int num_threads);

}  // namespace hw
