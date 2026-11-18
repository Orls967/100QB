#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "asian/options.hpp"
#include "asian/rng.hpp"

#include <cmath>
#include <vector>

namespace {

asian::OptionParams default_params() {
  return asian::OptionParams{
      .s0 = 100.0,
      .k = 100.0,
      .r = 0.05,
      .sigma = 0.20,
      .t = 1.0,
      .n_steps = 252,
  };
}

} // namespace

TEST_CASE("Black-Scholes put-call parity") {
  for (const double k : {90.0, 100.0, 105.0}) {
    auto p = default_params();
    p.k = k;
    const double c = asian::bs_call(p);
    const double put = asian::bs_put(p);
    const double rhs = p.s0 - p.k * std::exp(-p.r * p.t);
    CHECK(std::abs((c - put) - rhs) < 1e-10);
  }
}

TEST_CASE("Geometric Asian closed form matches large MC") {
  const auto p = default_params();
  const double exact = asian::geometric_asian_call(p);
  constexpr std::size_t n_paths = 40000;
  const double dt = p.t / static_cast<double>(p.n_steps);
  const double drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const double vol = p.sigma * std::sqrt(dt);
  const double disc = std::exp(-p.r * p.t);
  const double log_s0 = std::log(p.s0);
  const double inv_n = 1.0 / static_cast<double>(p.n_steps);

  asian::Xoshiro256StarStar rng(42);
  double sum = 0.0;
  double sum_sq = 0.0;
  for (std::size_t i = 0; i < n_paths; ++i) {
    double log_s = log_s0;
    double sum_log_s = 0.0;
    for (std::size_t step = 0; step < p.n_steps; ++step) {
      log_s += drift + vol * rng.next_normal();
      sum_log_s += log_s;
    }
    const double geom_avg = std::exp(sum_log_s * inv_n);
    const double pay = std::max(geom_avg - p.k, 0.0) * disc;
    sum += pay;
    sum_sq += pay * pay;
  }
  const double mean = sum / static_cast<double>(n_paths);
  const double var = (sum_sq - static_cast<double>(n_paths) * mean * mean) /
                     static_cast<double>(n_paths - 1);
  const double se = std::sqrt(var / static_cast<double>(n_paths));
  CHECK(std::abs(mean - exact) < 3.0 * se);
}

TEST_CASE("Control variate gives smaller SE than naive MC") {
  const auto p = default_params();
  const auto naive = asian::price_asian_naive(p, 20000, 42);
  const auto cv = asian::price_asian_control_variate(p, 20000, 42);
  CHECK(cv.se < naive.se / 10.0);
}

TEST_CASE("Deterministic across same seed and 1 vs 4 threads") {
  auto p = default_params();
  p.n_steps = 50;
  const auto a = asian::price_asian_control_variate(p, 10000, 7, 1);
  const auto b = asian::price_asian_control_variate(p, 10000, 7, 4);
  CHECK(a.price == b.price);
  CHECK(a.se == b.se);
}

TEST_CASE("Pathwise delta matches bump-and-revalue within 3 SE") {
  const auto p = default_params();
  const auto pw = asian::delta_asian_pathwise(p, 50000, 42);
  const auto bmp = asian::delta_asian_bump(p, 50000, 0.01, 42);
  const double se_comb = std::hypot(pw.se, bmp.se);
  CHECK(std::abs(pw.delta - bmp.delta) < 3.0 * se_comb);
}

TEST_CASE("RNG standard normal mean ~ 0 and variance ~ 1") {
  asian::Xoshiro256StarStar rng(42);
  constexpr std::size_t n = 40000;
  double sum = 0.0;
  double sum_sq = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const double z = rng.next_normal();
    sum += z;
    sum_sq += z * z;
  }
  const double mean = sum / static_cast<double>(n);
  const double var = sum_sq / static_cast<double>(n) - mean * mean;
  CHECK(std::abs(mean) < 0.02);
  CHECK(std::abs(var - 1.0) < 0.03);
}
