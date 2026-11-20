#pragma once

#include "hw/curve.hpp"
#include "hw/hull_white.hpp"

#include <cstdint>

namespace hw {

struct BermudanLsmResult {
    double price{0.0};
    double std_error{0.0};
    int train_paths{0};
    int eval_paths{0};
};

struct Dv01Result {
    double base_price{0.0};
    double up_price{0.0};
    double down_price{0.0};
    double dv01_per_1bp{0.0};
};

/// Two-pass Longstaff-Schwartz Monte Carlo for 10y Bermudan payer swaption
/// exercisable annually at t in {1, 2, ..., 9}.
BermudanLsmResult price_bermudan_lsm(const NelsonSiegelCurve& curve,
                                     const HwParams& params,
                                     int first_exercise_yr,
                                     int last_exercise_yr,
                                     int tenor_end_yr,
                                     double fixed_rate,
                                     int paths_per_pass,
                                     int steps_per_year,
                                     std::uint64_t seed,
                                     unsigned int num_threads);

/// Computes parallel-shift +-1bp DV01 on a 5x5 European swaption via Jamshidian.
Dv01Result dv01_euro_jamshidian(const NelsonSiegelCurve& curve,
                                const HwParams& params,
                                double expiry_years,
                                double tenor_end_years,
                                double fixed_rate,
                                double shift_bp = 1.0);

/// Computes parallel-shift +-1bp DV01 on a 5x5 European swaption via trinomial tree.
Dv01Result dv01_euro_tree(const NelsonSiegelCurve& curve,
                          const HwParams& params,
                          double expiry_years,
                          double tenor_end_years,
                          double fixed_rate,
                          int steps_per_year,
                          double shift_bp = 1.0);

/// Computes parallel-shift +-1bp DV01 on a 10y Bermudan payer swaption via trinomial tree.
Dv01Result dv01_bermudan_tree(const NelsonSiegelCurve& curve,
                              const HwParams& params,
                              int first_exercise_yr,
                              int last_exercise_yr,
                              int tenor_end_yr,
                              double fixed_rate,
                              int steps_per_year,
                              double shift_bp = 1.0);

}  // namespace hw
