#pragma once

#include <cmath>

namespace hw {

struct NelsonSiegelCurve {
  double beta0{0.045};
  double beta1{-0.02};
  double beta2{0.015};
  double tau{2.0};
  double parallel_shift{0.0};

  [[nodiscard]] double zero_rate(double t) const noexcept {
    if (t <= 1e-12) {
      return beta0 + beta1 + parallel_shift;
    }
    const double m = t / tau;
    const double em = std::exp(-m);
    const double factor = (1.0 - em) / m;
    return beta0 + beta1 * factor + beta2 * (factor - em) + parallel_shift;
  }

  [[nodiscard]] double discount(double t) const noexcept {
    if (t <= 0.0) {
      return 1.0;
    }
    return std::exp(-zero_rate(t) * t);
  }

  [[nodiscard]] double inst_forward(double t) const noexcept {
    if (t <= 0.0) {
      return beta0 + beta1 + parallel_shift;
    }
    const double m = t / tau;
    const double em = std::exp(-m);
    return beta0 + beta1 * em + beta2 * m * em + parallel_shift;
  }
};

} // namespace hw
