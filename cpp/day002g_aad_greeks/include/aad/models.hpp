#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "aad/cholesky.hpp"
#include "aad/dual.hpp"
#include "aad/number.hpp"

namespace aad {

struct BsAnalytic {
  double price;
  double delta;
  double d_k;
  double rho;
  double vega;
  double theta; // dC/dT (positive)
  double gamma;
};

template <typename Scalar>
Scalar bs_call(const Scalar& s, const Scalar& k, const Scalar& r, const Scalar& sigma,
               const Scalar& t) {
  const Scalar sqrt_t = math::sqrt(t);
  const Scalar vol_t = sigma * sqrt_t;
  const Scalar d1 = (math::log(s / k) + (r + Scalar(0.5) * sigma * sigma) * t) / vol_t;
  const Scalar d2 = d1 - vol_t;
  return s * math::normal_cdf(d1) - k * math::exp(-r * t) * math::normal_cdf(d2);
}

inline BsAnalytic bs_analytic_greeks(double s, double k, double r, double sigma, double t) noexcept {
  const double sqrt_t = std::sqrt(t);
  const double vol_t = sigma * sqrt_t;
  const double d1 = (std::log(s / k) + (r + 0.5 * sigma * sigma) * t) / vol_t;
  const double d2 = d1 - vol_t;
  const double nd1 = math::normal_cdf(d1);
  const double nd2 = math::normal_cdf(d2);
  const double pd1 = math::normal_pdf(d1);
  const double disc = std::exp(-r * t);

  const double price = s * nd1 - k * disc * nd2;
  const double delta = nd1;
  const double d_k = -disc * nd2;
  const double rho = k * t * disc * nd2;
  const double vega = s * pd1 * sqrt_t;
  const double theta = s * pd1 * sigma / (2.0 * sqrt_t) + r * k * disc * nd2; // dC/dT
  const double gamma = pd1 / (s * vol_t);
  return {price, delta, d_k, rho, vega, theta, gamma};
}

inline double digital_call_analytic_delta(double s, double k, double r, double sigma,
                                          double t) noexcept {
  const double sqrt_t = std::sqrt(t);
  const double vol_t = sigma * sqrt_t;
  const double d2 = (std::log(s / k) + (r - 0.5 * sigma * sigma) * t) / vol_t;
  return std::exp(-r * t) * math::normal_pdf(d2) / (s * vol_t);
}

template <typename Scalar>
Scalar geom_basket_closed_form(const std::vector<Scalar>& s0, const std::vector<Scalar>& sigma,
                               const std::vector<Scalar>& rhos, const std::vector<double>& weights,
                               const Scalar& k, const Scalar& r, const Scalar& t) {
  const std::size_t d = s0.size();
  Scalar m(0.0);
  for (std::size_t i = 0; i < d; ++i) {
    m = m + Scalar(weights[i]) * (math::log(s0[i]) + (r - Scalar(0.5) * sigma[i] * sigma[i]) * t);
  }

  // Build full correlation matrix
  std::vector<std::vector<Scalar>> corr(d, std::vector<Scalar>(d, Scalar(1.0)));
  std::size_t idx = 0;
  for (std::size_t i = 0; i < d; ++i) {
    for (std::size_t j = i + 1; j < d; ++j) {
      corr[i][j] = rhos[idx];
      corr[j][i] = rhos[idx];
      ++idx;
    }
  }

  Scalar v(0.0);
  for (std::size_t i = 0; i < d; ++i) {
    for (std::size_t j = 0; j < d; ++j) {
      v = v + Scalar(weights[i] * weights[j]) * sigma[i] * sigma[j] * corr[i][j];
    }
  }
  v = v * t;

  const Scalar sqrt_v = math::sqrt(v);
  const Scalar d1 = (m + v - math::log(k)) / sqrt_v;
  const Scalar d2 = d1 - sqrt_v;
  return math::exp(-r * t) * (math::exp(m + Scalar(0.5) * v) * math::normal_cdf(d1) -
                              k * math::normal_cdf(d2));
}

} // namespace aad
