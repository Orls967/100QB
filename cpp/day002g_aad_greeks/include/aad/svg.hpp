#pragma once

#include <string>
#include <vector>

#include "aad/mc_engine.hpp"

namespace aad::svg {

void write_aad_vs_analytic_svg(const std::string& path, const Gbm1StepResult& gbm,
                               double an_delta, double an_vega, double an_rho, double an_theta);

void write_digital_bias_var_svg(const std::string& path, const DigitalResult& dig);

void write_basket_speedup_svg(const std::string& path, const BasketComparison& arith);

void write_tape_growth_svg(const std::string& path, std::size_t gbm_nodes,
                           std::size_t basket_nodes, std::size_t asian_nodes);

} // namespace aad::svg
