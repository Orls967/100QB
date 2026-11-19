#include "hw/hull_white.hpp"

#include <cmath>
#include <vector>

#include "hw/solvers.hpp"

namespace hw {

double normal_cdf(double x) noexcept {
  return 0.5 * std::erfc(-x * 0.70710678118654752440);
}

double hw_b(double a, double t, double big_t) noexcept {
  return (1.0 - std::exp(-a * (big_t - t))) / a;
}

double hw_ln_a(const NelsonSiegelCurve& curve, const HwParams& p, double t, double big_t) noexcept {
  const double b = hw_b(p.a, t, big_t);
  const double p0t = curve.discount(t);
  const double p0big_t = curve.discount(big_t);
  const double f0t = curve.inst_forward(t);
  const double term = (p.sigma * p.sigma / (4.0 * p.a)) * (1.0 - std::exp(-2.0 * p.a * t)) * b * b;
  return std::log(p0big_t / p0t) + b * f0t - term;
}

double hw_bond_price(const NelsonSiegelCurve& curve, const HwParams& p, double t, double big_t,
                     double r_t) noexcept {
  const double ln_a = hw_ln_a(curve, p, t, big_t);
  const double b = hw_b(p.a, t, big_t);
  return std::exp(ln_a - b * r_t);
}

double hw_zbc(const NelsonSiegelCurve& curve, const HwParams& p, double t_opt, double t_bond,
              double strike) noexcept {
  const double p0t = curve.discount(t_opt);
  const double p0s = curve.discount(t_bond);
  if (t_opt <= 0.0) {
    return std::max(p0s - strike, 0.0);
  }
  const double sig_p = (p.sigma / p.a) * (1.0 - std::exp(-p.a * (t_bond - t_opt))) *
                       std::sqrt((1.0 - std::exp(-2.0 * p.a * t_opt)) / (2.0 * p.a));
  const double h = std::log(p0s / (p0t * strike)) / sig_p + 0.5 * sig_p;
  return p0s * normal_cdf(h) - strike * p0t * normal_cdf(h - sig_p);
}

double hw_zbp(const NelsonSiegelCurve& curve, const HwParams& p, double t_opt, double t_bond,
              double strike) noexcept {
  const double p0t = curve.discount(t_opt);
  const double p0s = curve.discount(t_bond);
  if (t_opt <= 0.0) {
    return std::max(strike - p0s, 0.0);
  }
  const double sig_p = (p.sigma / p.a) * (1.0 - std::exp(-p.a * (t_bond - t_opt))) *
                       std::sqrt((1.0 - std::exp(-2.0 * p.a * t_opt)) / (2.0 * p.a));
  const double h = std::log(p0s / (p0t * strike)) / sig_p + 0.5 * sig_p;
  return strike * p0t * normal_cdf(-h + sig_p) - p0s * normal_cdf(-h);
}

double hw_cap_5y(const NelsonSiegelCurve& curve, const HwParams& p, double strike_rate) noexcept {
  const double x = 1.0 / (1.0 + strike_rate);
  double total = 0.0;
  for (int i = 1; i <= 5; ++i) {
    const double t_prev = static_cast<double>(i - 1);
    const double t_curr = static_cast<double>(i);
    total += (1.0 + strike_rate) * hw_zbp(curve, p, t_prev, t_curr, x);
  }
  return total;
}

double hw_payer_swaption_jamshidian(const NelsonSiegelCurve& curve, const HwParams& p,
                                    double t0, int tenor_end_yr, double fixed_rate) noexcept {
  const int t0_i = static_cast<int>(std::llround(t0));
  std::vector<double> pay_times;
  std::vector<double> cashflows;
  for (int yr = t0_i + 1; yr <= tenor_end_yr; ++yr) {
    pay_times.push_back(static_cast<double>(yr));
    cashflows.push_back((yr == tenor_end_yr) ? (1.0 + fixed_rate) : fixed_rate);
  }

  auto sum_bonds = [&](double r_val) {
    double s = 0.0;
    for (std::size_t i = 0; i < pay_times.size(); ++i) {
      s += cashflows[i] * hw_bond_price(curve, p, t0, pay_times[i], r_val);
    }
    return s - 1.0;
  };

  const double r_star = solvers::brent_root(sum_bonds, -0.5, 0.8);
  double price = 0.0;
  for (std::size_t i = 0; i < pay_times.size(); ++i) {
    const double x_i = hw_bond_price(curve, p, t0, pay_times[i], r_star);
    price += cashflows[i] * hw_zbp(curve, p, t0, pay_times[i], x_i);
  }
  return price;
}

} // namespace hw
