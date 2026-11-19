#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace aad {

struct McEstimate {
  double mean{0.0};
  double se{0.0};
  double ci_low{0.0};
  double ci_high{0.0};
};

struct Gbm1StepResult {
  McEstimate price;
  McEstimate delta;
  McEstimate vega;
  McEstimate rho;
  McEstimate theta;
  double time_sec{0.0};
};

struct AsianAadResult {
  McEstimate price;
  McEstimate delta;
  McEstimate delta_formula;
  McEstimate vega;
  McEstimate vega_bump;
  McEstimate rho;
  std::size_t nodes_per_path{0};
  double time_sec{0.0};
  double paths_per_sec{0.0};
};

struct DigitalEpsRow {
  double eps{0.0};
  McEstimate delta;
  double bias{0.0};
  double mse{0.0};
};

struct DigitalResult {
  double analytic_delta{0.0};
  double pathwise_raw_delta{0.0};
  std::vector<DigitalEpsRow> smoothed;
  double optimal_eps{0.0};
};

struct BasketComparison {
  std::vector<McEstimate> aad_grads; // 67 sensitivities: 10 delta, 10 vega, rho, theta, 45 rho_ij
  std::vector<McEstimate> ref_grads; // CRN FD (arith) or exact analytic AAD (geom)
  std::vector<double> diff_in_se;
  double max_diff_se{0.0};
  std::size_t violations_3se{0};
  std::size_t violations_4se{0};
  std::size_t nodes_per_path{0};
  std::size_t peak_tape_bytes{0};
  double aad_time_sec{0.0};
  double value_only_time_sec{0.0};
  double fd_time_sec{0.0};
  double speedup_vs_fd{0.0};
  double grad_to_val_ratio{0.0};
  double paths_per_sec{0.0};
};

Gbm1StepResult run_gbm_1step_aad(std::size_t n_paths, std::uint64_t seed,
                                 std::size_t n_threads = 0);

AsianAadResult run_asian_aad(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads = 0);

DigitalResult run_digital_aad(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads = 0);

BasketComparison run_basket_10_arith(std::size_t n_paths, std::uint64_t seed,
                                     std::size_t n_threads = 0, bool run_fd = true);

BasketComparison run_basket_10_geom(std::size_t n_paths, std::uint64_t seed,
                                    std::size_t n_threads = 0);

} // namespace aad
