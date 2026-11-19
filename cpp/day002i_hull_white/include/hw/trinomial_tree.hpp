#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "hw/curve.hpp"
#include "hw/hull_white.hpp"

namespace hw {

struct BranchProbs {
  int k_mid{0};
  double pd{0.0};
  double pm{0.0};
  double pu{0.0};
};

class TrinomialTree {
public:
  TrinomialTree(const NelsonSiegelCurve& curve, const HwParams& params, int steps_per_year,
                int max_years);

  [[nodiscard]] int steps_per_year() const noexcept { return m_; }
  [[nodiscard]] int total_steps() const noexcept { return n_steps_; }
  [[nodiscard]] int j_max() const noexcept { return j_max_; }
  [[nodiscard]] double dt() const noexcept { return dt_; }
  [[nodiscard]] double dx() const noexcept { return dx_; }
  [[nodiscard]] std::size_t node_count() const noexcept;
  [[nodiscard]] double max_curve_error() const noexcept { return max_curve_err_; }

  [[nodiscard]] BranchProbs branching(int j) const noexcept;
  [[nodiscard]] double rate(int step, int j) const noexcept;

  [[nodiscard]] std::vector<double> rollback_zero_bond(int mat_step, int stop_step) const;
  [[nodiscard]] double price_zbp(double t_opt, double t_bond, double strike) const;
  [[nodiscard]] double price_zbc(double t_opt, double t_bond, double strike) const;
  [[nodiscard]] double price_cap_5y(double strike_rate) const;
  [[nodiscard]] double price_payer_swaption(double t0, int tenor_end_yr, double fixed_rate) const;
  [[nodiscard]] double price_bermudan_payer(int tenor_end_yr, double fixed_rate) const;

private:
  [[nodiscard]] std::size_t idx(int j) const noexcept {
    return static_cast<std::size_t>(j + j_max_);
  }
  void build(const NelsonSiegelCurve& curve);

  HwParams params_;
  int m_{12};
  int n_steps_{0};
  int j_max_{0};
  double dt_{0.0};
  double dx_{0.0};
  double max_curve_err_{0.0};
  std::vector<double> alpha_;
  std::vector<std::vector<double>> q_;
};

} // namespace hw
