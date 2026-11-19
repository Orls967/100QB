#pragma once

#include <string>
#include <vector>

#include "hawkes_ob/hawkes.hpp"
#include "hawkes_ob/market_sim.hpp"

namespace hob::svg {

void write_fano_and_qq_svg(const std::string& path, const std::vector<double>& win_sizes,
                           const std::vector<double>& fanos, const std::vector<double>& res_exp1);

void write_mle_scatter_svg(const std::string& path, const HawkesMleMcSummary& mle);

void write_latency_svg(const std::string& path, const OrderFlowBenchmark& ob);

void write_market_stylized_svg(const std::string& path, const MarketSimResult& mkt);

} // namespace hob::svg
