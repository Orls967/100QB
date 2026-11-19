#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace hob {

struct OrderFlowBenchmark {
  std::size_t n_orders{0};
  std::size_t n_trades{0};
  std::uint64_t trade_hash{0};
  bool exact_match{false};
  std::size_t invariant_violations{0};
  std::size_t heap_allocs_after_warmup{0};
  double throughput_ops_per_sec{0.0};
  double p50_ns{0.0};
  double p90_ns{0.0};
  double p99_ns{0.0};
  double p999_ns{0.0};
  double ref_time_sec{0.0};
  double fast_time_sec{0.0};
  double speedup_vs_ref{0.0};
  std::vector<double> latencies_ns;
};

struct MarketSimResult {
  std::size_t n_market_orders{0};
  double emp_market_rate{0.0};
  double emp_market_rate_se{0.0};
  double fano_market_orders{0.0};
  double mean_spread{0.0};
  double std_spread{0.0};
  std::vector<double> depth_profile;      // avg qty at levels 1..10 from best
  std::vector<double> signed_flow_acf;    // lag 0..20
  std::vector<double> abs_ret_acf;        // lag 0..50
  double price_impact_slope{0.0};
  double price_impact_r2{0.0};
  std::size_t invariant_violations{0};
  std::int64_t net_position_sum{0};
  double time_sec{0.0};
};

OrderFlowBenchmark benchmark_orderbook(std::size_t n_orders, std::uint64_t seed,
                                       bool run_reference = true);

bool verify_equivalence_seed(std::size_t n_orders, std::uint64_t seed);

MarketSimResult run_combined_market_sim(double t_end, std::uint64_t seed);

} // namespace hob
