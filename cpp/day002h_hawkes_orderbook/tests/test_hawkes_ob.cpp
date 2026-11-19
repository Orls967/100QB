#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <cmath>
#include <vector>

#include "hawkes_ob/hawkes.hpp"
#include "hawkes_ob/market_sim.hpp"
#include "hawkes_ob/orderbook.hpp"
#include "hawkes_ob/stats.hpp"

using namespace hob;

TEST_CASE("Hawkes Ogata vs Dassios-Zhao, time-rescaling GOF and MLE") {
  const HawkesParams p_true{0.5, 0.6, 1.2};
  const auto t_og = simulate_hawkes_ogata(p_true, 2000.0, 42);
  const auto t_dz = simulate_hawkes_dassios_zhao(p_true, 2000.0, 43);

  const double rate_og = static_cast<double>(t_og.size()) / 2000.0;
  const double rate_dz = static_cast<double>(t_dz.size()) / 2000.0;
  CHECK(std::abs(rate_og - 1.0) < 0.15);
  CHECK(std::abs(rate_dz - 1.0) < 0.15);

  const auto ks2 = stats::ks_two_sample(inter_event_times(t_og), inter_event_times(t_dz));
  CHECK(ks2.p_value > 0.01);

  const auto res_good = time_rescaling_residuals(t_dz, p_true);
  const auto ks_good = stats::ks_exp1(res_good);
  CHECK(ks_good.p_value > 0.01);

  const HawkesParams p_wrong{0.8, 0.24, 1.2}; // same mean rate 1.0, wrong n=0.2
  const auto res_wrong = time_rescaling_residuals(t_dz, p_wrong);
  const auto ks_wrong = stats::ks_exp1(res_wrong);
  CHECK(ks_wrong.p_value < 0.01);

  const auto fit = fit_hawkes_mle(t_dz, 2000.0);
  CHECK(std::abs(fit.params.mu - p_true.mu) < 4.0 * fit.se_mu);
  CHECK(std::abs(fit.params.alpha - p_true.alpha) < 4.0 * fit.se_alpha);
  CHECK(std::abs(fit.params.beta - p_true.beta) < 4.0 * fit.se_beta);
}

TEST_CASE("OrderBook equivalence with ReferenceBook across 25 seeds and MODIFY priority") {
  for (std::uint64_t sd = 1; sd <= 25; ++sd) {
    CHECK(verify_equivalence_seed(1500, sd));
  }

  // Test MODIFY priority: decreasing qty preserves FIFO, increasing qty loses FIFO
  OrderBook ob(1, 1000, 100);
  std::vector<Trade> tr;
  ob.submit_limit(1, Side::Buy, 100, 10, tr);
  ob.submit_limit(2, Side::Buy, 100, 10, tr);
  ob.modify(1, 5); // keep priority ahead of id=2
  ob.submit_market(3, Side::Sell, 5, tr);
  REQUIRE(tr.size() == 1);
  CHECK(tr[0].maker_id == 1);

  ob.submit_limit(4, Side::Buy, 100, 10, tr);
  ob.modify(2, 12); // increase qty -> moves behind id=4
  ob.submit_market(5, Side::Sell, 10, tr);
  REQUIRE(tr.size() == 1);
  CHECK(tr[0].maker_id == 4);
}

TEST_CASE("Hawkes MLE Monte Carlo thread invariance and seed reproducibility") {
  const HawkesParams p_true{0.5, 0.6, 1.2};
  const auto m1 = run_hawkes_mle_mc(p_true, 8, 500.0, 42, 1);
  const auto m4 = run_hawkes_mle_mc(p_true, 8, 500.0, 42, 4);
  CHECK(m1.mean_mu == doctest::Approx(m4.mean_mu).epsilon(1e-15));
  CHECK(m1.mean_alpha == doctest::Approx(m4.mean_alpha).epsilon(1e-15));
  CHECK(m1.mean_beta == doctest::Approx(m4.mean_beta).epsilon(1e-15));
}
