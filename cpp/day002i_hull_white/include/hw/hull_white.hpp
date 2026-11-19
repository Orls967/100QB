#pragma once

#include <vector>

#include "hw/curve.hpp"

namespace hw {

struct HwParams {
  double a{0.08};
  double sigma{0.012};
};

double normal_cdf(double x) noexcept;

double hw_b(double a, double t, double big_t) noexcept;

double hw_ln_a(const NelsonSiegelCurve& curve, const HwParams& p, double t, double big_t) noexcept;

double hw_bond_price(const NelsonSiegelCurve& curve, const HwParams& p, double t, double big_t,
                     double r_t) noexcept;

double hw_zbc(const NelsonSiegelCurve& curve, const HwParams& p, double t_opt, double t_bond,
              double strike) noexcept;

double hw_zbp(const NelsonSiegelCurve& curve, const HwParams& p, double t_opt, double t_bond,
              double strike) noexcept;

double hw_cap_5y(const NelsonSiegelCurve& curve, const HwParams& p, double strike_rate) noexcept;

double hw_payer_swaption_jamshidian(const NelsonSiegelCurve& curve, const HwParams& p,
                                    double t0, int tenor_end_yr, double fixed_rate) noexcept;

} // namespace hw
