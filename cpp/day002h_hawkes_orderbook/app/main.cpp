#include <atomic>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <new>
#include <vector>

#include "hawkes_ob/hawkes.hpp"
#include "hawkes_ob/market_sim.hpp"
#include "hawkes_ob/stats.hpp"
#include "hawkes_ob/svg.hpp"

namespace hob {
extern std::atomic<bool> g_track_allocs;
extern std::atomic<std::size_t> g_alloc_count;
} // namespace hob

void* operator new(std::size_t sz) {
  if (hob::g_track_allocs.load(std::memory_order_relaxed)) {
    hob::g_alloc_count.fetch_add(1, std::memory_order_relaxed);
  }
  void* p = std::malloc(sz == 0 ? 1 : sz);
  if (!p) {
    throw std::bad_alloc();
  }
  return p;
}

void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace hob;

namespace {

std::filesystem::path find_repo_root() {
  auto dir = std::filesystem::current_path();
  while (!std::filesystem::is_directory(dir / "src" / "qb100")) {
    if (dir == dir.parent_path()) {
      throw std::runtime_error("Cannot find 100QB repo root (src/qb100)");
    }
    dir = dir.parent_path();
  }
  return dir;
}

} // namespace

int main() {
  const auto root = find_repo_root();
  const auto out_dir = root / "outputs" / "day002h";
  const auto fig_dir = root / "outputs" / "figures";
  std::filesystem::create_directories(out_dir);
  std::filesystem::create_directories(fig_dir);

  // Bagian A: Hawkes (T = 10,000)
  const HawkesParams p_true{0.5, 0.6, 1.2};
  const auto t_og = simulate_hawkes_ogata(p_true, 10'000.0, 42);
  const auto t_dz = simulate_hawkes_dassios_zhao(p_true, 10'000.0, 43);

  const auto win_c_dz = stats::window_counts(t_dz, 10'000.0, 50.0);
  const double rate_dz = static_cast<double>(t_dz.size()) / 10'000.0;
  const double rate_se = std::sqrt(stats::variance(win_c_dz) / static_cast<double>(win_c_dz.size())) / 50.0;
  const double fano_dz = stats::fano_factor(t_dz, 10'000.0, 50.0);
  const double fano_og = stats::fano_factor(t_og, 10'000.0, 50.0);
  const auto ks_og_dz = stats::ks_two_sample(inter_event_times(t_og), inter_event_times(t_dz));

  const auto res_good = time_rescaling_residuals(t_dz, p_true);
  const auto ks_good = stats::ks_exp1(res_good);
  const auto lb_good = stats::ljung_box(res_good, 10);

  const HawkesParams p_wrong{0.8, 0.24, 1.2};
  const auto res_wrong = time_rescaling_residuals(t_dz, p_wrong);
  const auto ks_wrong = stats::ks_exp1(res_wrong);
  const auto lb_wrong = stats::ljung_box(res_wrong, 10);

  const auto mle_mc = run_hawkes_mle_mc(p_true, 200, 5000.0, 42);

  std::vector<double> win_sizes{5, 10, 20, 50, 100, 200};
  std::vector<double> fanos;
  for (double w : win_sizes) {
    fanos.push_back(stats::fano_factor(t_dz, 10'000.0, w));
  }

  // Bagian B: Order Book 1,000,000 orders
  const auto ob_bench = benchmark_orderbook(1'000'000, 42, true);

  // Bagian C: Combined Market Simulation (T = 20,000)
  const auto mkt = run_combined_market_sim(20'000.0, 42);

  svg::write_fano_and_qq_svg((fig_dir / "day002h_hawkes-fano-qq.svg").string(), win_sizes, fanos,
                             res_good);
  svg::write_mle_scatter_svg((fig_dir / "day002h_hawkes-mle-scatter.svg").string(), mle_mc);
  svg::write_latency_svg((fig_dir / "day002h_ob-latency.svg").string(), ob_bench);
  svg::write_market_stylized_svg((fig_dir / "day002h_market-stylized-facts.svg").string(), mkt);

  std::ofstream jf(out_dir / "results.json");
  jf << std::setprecision(17);
  jf << "{\n"
     << "  \"hawkes\": {\n"
     << "    \"rate_dz\": " << rate_dz << ", \"rate_se\": " << rate_se
     << ", \"fano_dz_50\": " << fano_dz << ", \"fano_og_50\": " << fano_og
     << ", \"ks_og_dz_d\": " << ks_og_dz.d_stat << ", \"ks_og_dz_p\": " << ks_og_dz.p_value << ",\n"
     << "    \"gof_true\": {\"ks_d\": " << ks_good.d_stat << ", \"ks_p\": " << ks_good.p_value
     << ", \"lb_q10\": " << lb_good.first << ", \"lb_p\": " << lb_good.second << "},\n"
     << "    \"gof_wrong\": {\"ks_d\": " << ks_wrong.d_stat << ", \"ks_p\": " << ks_wrong.p_value
     << ", \"lb_q10\": " << lb_wrong.first << ", \"lb_p\": " << lb_wrong.second << "},\n"
     << "    \"mle_mc\": {\"r_reps\": " << mle_mc.r_reps
     << ", \"mean_mu\": " << mle_mc.mean_mu << ", \"sd_mu\": " << mle_mc.sd_mu
     << ", \"mean_se_mu\": " << mle_mc.mean_se_mu << ", \"cov_mu\": " << mle_mc.cov_mu
     << ", \"mean_alpha\": " << mle_mc.mean_alpha << ", \"sd_alpha\": " << mle_mc.sd_alpha
     << ", \"mean_se_alpha\": " << mle_mc.mean_se_alpha << ", \"cov_alpha\": " << mle_mc.cov_alpha
     << ", \"mean_beta\": " << mle_mc.mean_beta << ", \"sd_beta\": " << mle_mc.sd_beta
     << ", \"mean_se_beta\": " << mle_mc.mean_se_beta << ", \"cov_beta\": " << mle_mc.cov_beta
     << ", \"corr_alpha_beta\": " << mle_mc.corr_alpha_beta << "}\n"
     << "  },\n"
     << "  \"orderbook_1m\": {\n"
     << "    \"exact_match\": " << (ob_bench.exact_match ? "true" : "false")
     << ", \"n_trades\": " << ob_bench.n_trades
     << ", \"heap_allocs_after_warmup\": " << ob_bench.heap_allocs_after_warmup
     << ", \"throughput_ops_per_sec\": " << ob_bench.throughput_ops_per_sec
     << ", \"p50_ns\": " << ob_bench.p50_ns << ", \"p90_ns\": " << ob_bench.p90_ns
     << ", \"p99_ns\": " << ob_bench.p99_ns << ", \"p999_ns\": " << ob_bench.p999_ns
     << ", \"fast_time_sec\": " << ob_bench.fast_time_sec
     << ", \"ref_time_sec\": " << ob_bench.ref_time_sec
     << ", \"speedup_vs_ref\": " << ob_bench.speedup_vs_ref << "\n"
     << "  },\n"
     << "  \"market_sim\": {\n"
     << "    \"n_market_orders\": " << mkt.n_market_orders
     << ", \"emp_rate\": " << mkt.emp_market_rate << ", \"emp_rate_se\": " << mkt.emp_market_rate_se
     << ", \"fano\": " << mkt.fano_market_orders
     << ", \"mean_spread\": " << mkt.mean_spread << ", \"std_spread\": " << mkt.std_spread
     << ", \"flow_acf_lag1\": " << mkt.signed_flow_acf[1]
     << ", \"abs_ret_acf_lag1\": " << mkt.abs_ret_acf[1]
     << ", \"price_impact_slope\": " << mkt.price_impact_slope
     << ", \"price_impact_r2\": " << mkt.price_impact_r2
     << ", \"invariant_violations\": " << mkt.invariant_violations
     << ", \"net_position_sum\": " << mkt.net_position_sum
     << ", \"time_sec\": " << mkt.time_sec << "\n"
     << "  }\n"
     << "}\n";

  std::cout << "Hawkes rate=" << rate_dz << "±" << rate_se << " fano=" << fano_dz
            << " | KS_good p=" << ks_good.p_value << " KS_wrong p=" << ks_wrong.p_value << "\n"
            << "OrderBook 1M: match=" << ob_bench.exact_match
            << " allocs=" << ob_bench.heap_allocs_after_warmup
            << " throughput=" << ob_bench.throughput_ops_per_sec << " ops/s"
            << " speedup=" << ob_bench.speedup_vs_ref << "x\n"
            << "MarketSim: rate=" << mkt.emp_market_rate << " fano=" << mkt.fano_market_orders
            << " inv_viol=" << mkt.invariant_violations << " net_pos=" << mkt.net_position_sum << "\n";
  return 0;
}
