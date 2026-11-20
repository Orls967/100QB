#pragma once

#include <span>
#include <string>

namespace asian {

struct ConvergencePoint {
  std::size_t n{0};
  double se_naive{0.0};
  double se_antithetic{0.0};
  double se_control_variate{0.0};
};

struct EstimatorSummary {
  std::string name;
  double price{0.0};
  double se{0.0};
  double ci_low{0.0};
  double ci_high{0.0};
  double vrf{1.0};
};

struct Slopes {
  double naive{0.0};
  double antithetic{0.0};
  double control_variate{0.0};
};

[[nodiscard]] std::string render_convergence_svg(std::span<const ConvergencePoint> points,
                                                 const Slopes &slopes);
[[nodiscard]] std::string render_estimators_svg(std::span<const EstimatorSummary> estimators);

} // namespace asian
