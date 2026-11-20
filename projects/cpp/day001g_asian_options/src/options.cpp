#include "asian/options.hpp"
#include "asian/rng.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <vector>

namespace asian {

namespace {

constexpr std::size_t kChunkSize = 4000;
constexpr double kCiZ95 = 1.96;

struct MomentAcc {
  double sum{0.0};
  double sum_sq{0.0};
};

struct CvAcc {
  double sx{0.0};
  double sy{0.0};
  double sxx{0.0};
  double syy{0.0};
  double sxy{0.0};
};

[[nodiscard]] std::size_t resolve_workers(std::size_t workers) noexcept {
  if (workers > 0) {
    return workers;
  }
  const unsigned int hc = std::thread::hardware_concurrency();
  return hc > 0 ? static_cast<std::size_t>(hc) : 1;
}

void validate_mc(const OptionParams &p, std::size_t n_paths) {
  if (p.n_steps < 1) {
    throw std::invalid_argument("n_steps must be >= 1");
  }
  if (n_paths < 2) {
    throw std::invalid_argument("n_paths must be >= 2");
  }
}

[[nodiscard]] EstimatorResult summarize_moments(std::size_t count, double sum,
                                                double sum_sq) noexcept {
  const double n = static_cast<double>(count);
  const double mean = sum / n;
  const double var = std::max((sum_sq - n * mean * mean) / (n - 1.0), 0.0);
  const double se = std::sqrt(var / n);
  return EstimatorResult{
      .price = mean,
      .se = se,
      .ci_low = mean - kCiZ95 * se,
      .ci_high = mean + kCiZ95 * se,
  };
}

template <typename Acc, typename Fn>
[[nodiscard]] std::vector<Acc> parallel_chunks(std::size_t n_items, std::size_t workers, Fn &&fn) {
  const std::size_t n_chunks = (n_items + kChunkSize - 1) / kChunkSize;
  std::vector<Acc> results(n_chunks);
  const std::size_t w = std::min(resolve_workers(workers), n_chunks);

  if (w <= 1) {
    for (std::size_t c = 0; c < n_chunks; ++c) {
      const std::size_t start = c * kChunkSize;
      const std::size_t len = std::min(n_items - start, kChunkSize);
      results[c] = fn(c, len);
    }
    return results;
  }

  std::vector<std::thread> threads;
  threads.reserve(w);
  for (std::size_t wid = 0; wid < w; ++wid) {
    threads.emplace_back([wid, w, n_chunks, n_items, &results, &fn]() {
      for (std::size_t c = wid; c < n_chunks; c += w) {
        const std::size_t start = c * kChunkSize;
        const std::size_t len = std::min(n_items - start, kChunkSize);
        results[c] = fn(c, len);
      }
    });
  }
  for (auto &th : threads) {
    th.join();
  }
  return results;
}

} // namespace

double norm_cdf(double x) noexcept {
  return 0.5 * std::erfc(-x / std::numbers::sqrt2_v<double>);
}

double bs_call(const OptionParams &p) noexcept {
  if (p.t <= 0.0 || p.sigma <= 0.0) {
    return std::max(p.s0 - p.k * std::exp(-p.r * p.t), 0.0);
  }
  const double vol_sqrt_t = p.sigma * std::sqrt(p.t);
  const double d1 = (std::log(p.s0 / p.k) + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / vol_sqrt_t;
  const double d2 = d1 - vol_sqrt_t;
  return p.s0 * norm_cdf(d1) - p.k * std::exp(-p.r * p.t) * norm_cdf(d2);
}

double bs_put(const OptionParams &p) noexcept {
  if (p.t <= 0.0 || p.sigma <= 0.0) {
    return std::max(p.k * std::exp(-p.r * p.t) - p.s0, 0.0);
  }
  const double vol_sqrt_t = p.sigma * std::sqrt(p.t);
  const double d1 = (std::log(p.s0 / p.k) + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / vol_sqrt_t;
  const double d2 = d1 - vol_sqrt_t;
  return p.k * std::exp(-p.r * p.t) * norm_cdf(-d2) - p.s0 * norm_cdf(-d1);
}

double geometric_asian_call(const OptionParams &p) {
  if (p.n_steps < 1) {
    throw std::invalid_argument("n_steps must be >= 1");
  }
  const double n = static_cast<double>(p.n_steps);
  const double mu_ln_g =
      std::log(p.s0) + (p.r - 0.5 * p.sigma * p.sigma) * p.t * (n + 1.0) / (2.0 * n);
  const double var_ln_g =
      (p.sigma * p.sigma * p.t) * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
  const double sigma_ln_g = std::sqrt(var_ln_g);
  const double expected_g = std::exp(mu_ln_g + 0.5 * var_ln_g);
  const double d1 = (mu_ln_g - std::log(p.k) + var_ln_g) / sigma_ln_g;
  const double d2 = d1 - sigma_ln_g;
  return std::exp(-p.r * p.t) * (expected_g * norm_cdf(d1) - p.k * norm_cdf(d2));
}

EstimatorResult price_european_mc(const OptionParams &p, std::size_t n_paths, std::uint64_t seed,
                                  std::size_t workers) {
  validate_mc(p, n_paths);
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);

  const auto chunks = parallel_chunks<MomentAcc>(n_paths, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    MomentAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        log_s += drift + vol * rng.next_normal();
      }
      const double payoff = std::max(std::exp(log_s) - p.k, 0.0) * disc;
      acc.sum += payoff;
      acc.sum_sq += payoff * payoff;
    }
    return acc;
  });

  double sum = 0.0;
  double sum_sq = 0.0;
  for (const auto &c : chunks) {
    sum += c.sum;
    sum_sq += c.sum_sq;
  }
  return summarize_moments(n_paths, sum, sum_sq);
}

EstimatorResult price_asian_naive(const OptionParams &p, std::size_t n_paths, std::uint64_t seed,
                                  std::size_t workers) {
  validate_mc(p, n_paths);
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  const auto chunks = parallel_chunks<MomentAcc>(n_paths, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    MomentAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      double sum_s = 0.0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        log_s += drift + vol * rng.next_normal();
        sum_s += std::exp(log_s);
      }
      const double payoff = std::max(sum_s * inv_n - p.k, 0.0) * disc;
      acc.sum += payoff;
      acc.sum_sq += payoff * payoff;
    }
    return acc;
  });

  double sum = 0.0;
  double sum_sq = 0.0;
  for (const auto &c : chunks) {
    sum += c.sum;
    sum_sq += c.sum_sq;
  }
  return summarize_moments(n_paths, sum, sum_sq);
}

EstimatorResult price_asian_antithetic(const OptionParams &p, std::size_t n_paths,
                                       std::uint64_t seed, std::size_t workers) {
  if (p.n_steps < 1) {
    throw std::invalid_argument("n_steps must be >= 1");
  }
  if (n_paths < 4 || n_paths % 2 != 0) {
    throw std::invalid_argument("n_paths must be even and >= 4");
  }
  const std::size_t n_pairs = n_paths / 2;
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  const auto chunks = parallel_chunks<MomentAcc>(n_pairs, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    MomentAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_pos = log_s0;
      double log_neg = log_s0;
      double sum_pos = 0.0;
      double sum_neg = 0.0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        const double inc = vol * rng.next_normal();
        log_pos += drift + inc;
        log_neg += drift - inc;
        sum_pos += std::exp(log_pos);
        sum_neg += std::exp(log_neg);
      }
      const double pair_avg =
          0.5 * (std::max(sum_pos * inv_n - p.k, 0.0) + std::max(sum_neg * inv_n - p.k, 0.0)) *
          disc;
      acc.sum += pair_avg;
      acc.sum_sq += pair_avg * pair_avg;
    }
    return acc;
  });

  double sum = 0.0;
  double sum_sq = 0.0;
  for (const auto &c : chunks) {
    sum += c.sum;
    sum_sq += c.sum_sq;
  }
  return summarize_moments(n_pairs, sum, sum_sq);
}

ControlVariateResult price_asian_control_variate(const OptionParams &p, std::size_t n_paths,
                                                 std::uint64_t seed, std::size_t workers) {
  validate_mc(p, n_paths);
  const double y_exact = geometric_asian_call(p);
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  const auto chunks = parallel_chunks<CvAcc>(n_paths, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    CvAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      double sum_s = 0.0;
      double sum_log_s = 0.0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        log_s += drift + vol * rng.next_normal();
        sum_s += std::exp(log_s);
        sum_log_s += log_s;
      }
      const double arith_avg = sum_s * inv_n;
      const double geom_avg = std::exp(sum_log_s * inv_n);
      const double x = std::max(arith_avg - p.k, 0.0) * disc;
      const double y = std::max(geom_avg - p.k, 0.0) * disc;
      acc.sx += x;
      acc.sy += y;
      acc.sxx += x * x;
      acc.syy += y * y;
      acc.sxy += x * y;
    }
    return acc;
  });

  double sx = 0.0;
  double sy = 0.0;
  double sxx = 0.0;
  double syy = 0.0;
  double sxy = 0.0;
  for (const auto &c : chunks) {
    sx += c.sx;
    sy += c.sy;
    sxx += c.sxx;
    syy += c.syy;
    sxy += c.sxy;
  }

  const double n = static_cast<double>(n_paths);
  const double mx = sx / n;
  const double my = sy / n;
  const double cov_xy = (sxy - n * mx * my) / (n - 1.0);
  const double var_y = (syy - n * my * my) / (n - 1.0);
  if (var_y <= 0.0) {
    throw std::runtime_error("var_y is zero; cannot compute control variate beta");
  }
  const double beta = cov_xy / var_y;
  const double z_mean = mx - beta * (my - y_exact);
  const double var_x = (sxx - n * mx * mx) / (n - 1.0);
  const double var_z = std::max(var_x + beta * beta * var_y - 2.0 * beta * cov_xy, 0.0);
  const double se = std::sqrt(var_z / n);

  return ControlVariateResult{
      .price = z_mean,
      .se = se,
      .ci_low = z_mean - kCiZ95 * se,
      .ci_high = z_mean + kCiZ95 * se,
      .beta = beta,
      .y_exact = y_exact,
  };
}

DeltaResult delta_asian_pathwise(const OptionParams &p, std::size_t n_paths, std::uint64_t seed,
                                 std::size_t workers) {
  validate_mc(p, n_paths);
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  const auto chunks = parallel_chunks<MomentAcc>(n_paths, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    MomentAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      double sum_s = 0.0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        log_s += drift + vol * rng.next_normal();
        sum_s += std::exp(log_s);
      }
      const double arith_avg = sum_s * inv_n;
      const double sample = arith_avg > p.k ? disc * (arith_avg / p.s0) : 0.0;
      acc.sum += sample;
      acc.sum_sq += sample * sample;
    }
    return acc;
  });

  double sum = 0.0;
  double sum_sq = 0.0;
  for (const auto &c : chunks) {
    sum += c.sum;
    sum_sq += c.sum_sq;
  }
  const auto res = summarize_moments(n_paths, sum, sum_sq);
  return DeltaResult{.delta = res.price, .se = res.se, .ci_low = res.ci_low, .ci_high = res.ci_high};
}

DeltaResult delta_asian_bump(const OptionParams &p, std::size_t n_paths, double h_rel,
                             std::uint64_t seed, std::size_t workers) {
  validate_mc(p, n_paths);
  if (h_rel <= 0.0) {
    throw std::invalid_argument("h_rel must be > 0");
  }
  const double h = h_rel * p.s0;
  const double scale_up = (p.s0 + h) / p.s0;
  const double scale_dn = (p.s0 - h) / p.s0;
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  const auto chunks = parallel_chunks<MomentAcc>(n_paths, workers, [&](std::size_t c, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c));
    MomentAcc acc{};
    for (std::size_t i = 0; i < len; ++i) {
      double log_s = log_s0;
      double sum_s = 0.0;
      for (std::size_t step = 0; step < p.n_steps; ++step) {
        log_s += drift + vol * rng.next_normal();
        sum_s += std::exp(log_s);
      }
      const double arith_avg = sum_s * inv_n;
      const double pay_up = std::max(arith_avg * scale_up - p.k, 0.0);
      const double pay_dn = std::max(arith_avg * scale_dn - p.k, 0.0);
      const double sample = disc * (pay_up - pay_dn) / (2.0 * h);
      acc.sum += sample;
      acc.sum_sq += sample * sample;
    }
    return acc;
  });

  double sum = 0.0;
  double sum_sq = 0.0;
  for (const auto &c : chunks) {
    sum += c.sum;
    sum_sq += c.sum_sq;
  }
  const auto res = summarize_moments(n_paths, sum, sum_sq);
  return DeltaResult{.delta = res.price, .se = res.se, .ci_low = res.ci_low, .ci_high = res.ci_high};
}

double fit_loglog_slope(std::span<const std::size_t> ns, std::span<const double> ses) {
  if (ns.size() != ses.size() || ns.size() < 2) {
    throw std::invalid_argument("ns and ses must have equal length >= 2");
  }
  const std::size_t count = ns.size();
  std::vector<double> xs(count);
  std::vector<double> ys(count);
  double sum_x = 0.0;
  double sum_y = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    xs[i] = std::log10(static_cast<double>(ns[i]));
    ys[i] = std::log10(ses[i]);
    sum_x += xs[i];
    sum_y += ys[i];
  }
  const double n = static_cast<double>(count);
  const double mx = sum_x / n;
  const double my = sum_y / n;
  double num = 0.0;
  double den = 0.0;
  for (std::size_t i = 0; i < count; ++i) {
    const double dx = xs[i] - mx;
    num += dx * (ys[i] - my);
    den += dx * dx;
  }
  if (den <= 0.0) {
    throw std::runtime_error("zero variance in x");
  }
  return num / den;
}

} // namespace asian
