#pragma once

#include <cstddef>
#include <cstdint>
#include <span>

namespace asian {

struct OptionParams {
  double s0{100.0};
  double k{100.0};
  double r{0.05};
  double sigma{0.20};
  double t{1.0};
  std::size_t n_steps{252};
};

struct EstimatorResult {
  double price{0.0};
  double se{0.0};
  double ci_low{0.0};
  double ci_high{0.0};
};

struct ControlVariateResult {
  double price{0.0};
  double se{0.0};
  double ci_low{0.0};
  double ci_high{0.0};
  double beta{0.0};
  double y_exact{0.0};
};

struct DeltaResult {
  double delta{0.0};
  double se{0.0};
  double ci_low{0.0};
  double ci_high{0.0};
};

[[nodiscard]] double norm_cdf(double x) noexcept;
[[nodiscard]] double bs_call(const OptionParams &p) noexcept;
[[nodiscard]] double bs_put(const OptionParams &p) noexcept;
[[nodiscard]] double geometric_asian_call(const OptionParams &p);

[[nodiscard]] EstimatorResult price_european_mc(const OptionParams &p, std::size_t n_paths,
                                                std::uint64_t seed = 42,
                                                std::size_t workers = 0);
[[nodiscard]] EstimatorResult price_asian_naive(const OptionParams &p, std::size_t n_paths,
                                                  std::uint64_t seed = 42,
                                                  std::size_t workers = 0);
[[nodiscard]] EstimatorResult price_asian_antithetic(const OptionParams &p, std::size_t n_paths,
                                                       std::uint64_t seed = 42,
                                                       std::size_t workers = 0);
[[nodiscard]] ControlVariateResult price_asian_control_variate(const OptionParams &p,
                                                               std::size_t n_paths,
                                                               std::uint64_t seed = 42,
                                                               std::size_t workers = 0);

[[nodiscard]] DeltaResult delta_asian_pathwise(const OptionParams &p, std::size_t n_paths,
                                               std::uint64_t seed = 42,
                                               std::size_t workers = 0);
[[nodiscard]] DeltaResult delta_asian_bump(const OptionParams &p, std::size_t n_paths,
                                           double h_rel = 0.01, std::uint64_t seed = 42,
                                           std::size_t workers = 0);

[[nodiscard]] double fit_loglog_slope(std::span<const std::size_t> ns,
                                      std::span<const double> ses);

} // namespace asian
