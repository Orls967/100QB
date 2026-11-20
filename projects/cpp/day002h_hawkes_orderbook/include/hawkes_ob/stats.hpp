#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace hob::stats {

struct KsResult {
  double d_stat{0.0};
  double p_value{1.0};
};

inline double ks_p_value_asymptotic(double d, double n_eff) noexcept {
  const double sq = std::sqrt(n_eff);
  const double lam = (sq + 0.12 + 0.11 / sq) * d;
  if (lam <= 0.0) {
    return 1.0;
  }
  double sum = 0.0;
  for (int k = 1; k <= 100; ++k) {
    const double term = 2.0 * ((k % 2 == 1) ? 1.0 : -1.0) * std::exp(-2.0 * k * k * lam * lam);
    sum += term;
    if (std::abs(term) < 1e-15) {
      break;
    }
  }
  return std::clamp(sum, 0.0, 1.0);
}

inline KsResult ks_exp1(const std::vector<double>& samples) {
  if (samples.empty()) {
    return {0.0, 1.0};
  }
  std::vector<double> s = samples;
  std::sort(s.begin(), s.end());
  const double n = static_cast<double>(s.size());
  double d_max = 0.0;
  for (std::size_t i = 0; i < s.size(); ++i) {
    const double cdf = 1.0 - std::exp(-s[i]);
    const double f_hi = static_cast<double>(i + 1) / n;
    const double f_lo = static_cast<double>(i) / n;
    d_max = std::max({d_max, std::abs(f_hi - cdf), std::abs(cdf - f_lo)});
  }
  return {d_max, ks_p_value_asymptotic(d_max, n)};
}

inline KsResult ks_two_sample(const std::vector<double>& a, const std::vector<double>& b) {
  std::vector<double> sa = a;
  std::vector<double> sb = b;
  std::sort(sa.begin(), sa.end());
  std::sort(sb.begin(), sb.end());
  const std::size_t na = sa.size();
  const std::size_t nb = sb.size();
  std::size_t ia = 0, ib = 0;
  double d_max = 0.0;
  while (ia < na && ib < nb) {
    if (sa[ia] <= sb[ib]) {
      ++ia;
    } else {
      ++ib;
    }
    const double fa = static_cast<double>(ia) / static_cast<double>(na);
    const double fb = static_cast<double>(ib) / static_cast<double>(nb);
    d_max = std::max(d_max, std::abs(fa - fb));
  }
  const double n_eff =
      (static_cast<double>(na) * static_cast<double>(nb)) / static_cast<double>(na + nb);
  return {d_max, ks_p_value_asymptotic(d_max, n_eff)};
}

inline std::vector<double> acf(const std::vector<double>& xs, std::size_t max_lag) {
  const std::size_t n = xs.size();
  if (n < 2) {
    return std::vector<double>(max_lag + 1, 0.0);
  }
  double mean = 0.0;
  for (double x : xs) {
    mean += x;
  }
  mean /= static_cast<double>(n);
  double denom = 0.0;
  for (double x : xs) {
    denom += (x - mean) * (x - mean);
  }
  std::vector<double> out(max_lag + 1, 0.0);
  out[0] = 1.0;
  if (denom <= 1e-18) {
    return out;
  }
  for (std::size_t lag = 1; lag <= std::min(max_lag, n - 1); ++lag) {
    double num = 0.0;
    for (std::size_t t = 0; t + lag < n; ++t) {
      num += (xs[t] - mean) * (xs[t + lag] - mean);
    }
    out[lag] = num / denom;
  }
  return out;
}

/// Normal CDF via erfc for Wilson-Hilferty chi-square p-value approximation.
inline double normal_cdf(double x) noexcept {
  return 0.5 * std::erfc(-x * 0.70710678118654752440);
}

inline std::pair<double, double> ljung_box(const std::vector<double>& xs, std::size_t max_lag) {
  const auto r = acf(xs, max_lag);
  const double n = static_cast<double>(xs.size());
  double q = 0.0;
  for (std::size_t k = 1; k <= max_lag; ++k) {
    q += (r[k] * r[k]) / (n - static_cast<double>(k));
  }
  q *= n * (n + 2.0);
  const double df = static_cast<double>(max_lag);
  // Wilson-Hilferty cube-root approximation to chi-square tail
  const double c = 2.0 / (9.0 * df);
  const double z = (std::pow(q / df, 1.0 / 3.0) - (1.0 - c)) / std::sqrt(c);
  const double p_val = 1.0 - normal_cdf(z);
  return {q, std::clamp(p_val, 0.0, 1.0)};
}

inline double mean(const std::vector<double>& xs) noexcept {
  double s = 0.0;
  for (double x : xs) {
    s += x;
  }
  return xs.empty() ? 0.0 : s / static_cast<double>(xs.size());
}

inline double variance(const std::vector<double>& xs) noexcept {
  if (xs.size() < 2) {
    return 0.0;
  }
  const double m = mean(xs);
  double s2 = 0.0;
  for (double x : xs) {
    s2 += (x - m) * (x - m);
  }
  return s2 / static_cast<double>(xs.size() - 1);
}

inline std::vector<double> window_counts(const std::vector<double>& times, double t_end,
                                         double win_len) {
  const auto n_win = static_cast<std::size_t>(std::floor(t_end / win_len));
  std::vector<double> counts(n_win, 0.0);
  for (double t : times) {
    const auto idx = static_cast<std::size_t>(t / win_len);
    if (idx < n_win) {
      counts[idx] += 1.0;
    }
  }
  return counts;
}

inline double fano_factor(const std::vector<double>& times, double t_end, double win_len) {
  const auto c = window_counts(times, t_end, win_len);
  const double m = mean(c);
  return (m > 0.0) ? (variance(c) / m) : 0.0;
}

} // namespace hob::stats
