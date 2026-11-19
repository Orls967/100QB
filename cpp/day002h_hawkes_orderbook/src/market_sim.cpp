#include "hawkes_ob/market_sim.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <numeric>
#include <vector>

#include "hawkes_ob/hawkes.hpp"
#include "hawkes_ob/orderbook.hpp"
#include "hawkes_ob/reference_book.hpp"
#include "hawkes_ob/rng.hpp"
#include "hawkes_ob/stats.hpp"

namespace hob {

// Atomic counter hooked by custom operator new in main.cpp to verify 0 heap allocs in hot path.
std::atomic<bool> g_track_allocs{false};
std::atomic<std::size_t> g_alloc_count{0};

namespace {

constexpr std::int64_t kMinPrice = 1;
constexpr std::int64_t kMaxPrice = 200'000;
constexpr std::int64_t kInitMid = 100'000;

struct OrderOp {
  enum Kind : std::uint8_t { Limit = 0, Cancel = 1, Market = 2, Modify = 3 };
  Kind kind{Limit};
  Side side{Side::Buy};
  std::uint64_t id{0};
  std::int64_t price{0};
  std::int64_t qty{0};
  std::uint64_t target_idx_seed{0};
};

std::uint64_t hash_trades(std::uint64_t h, const std::vector<Trade>& trades) noexcept {
  for (const auto& t : trades) {
    h ^= t.maker_id * 0x9e3779b97f4a7c15ULL;
    h = (h << 13) | (h >> 51);
    h ^= t.taker_id * 0xbf58476d1ce4e5b9ULL;
    h ^= static_cast<std::uint64_t>(t.price) * 0x94d049bb133111ebULL;
    h ^= static_cast<std::uint64_t>(t.qty) * 0x517cc1b727220a95ULL;
  }
  return h;
}

std::vector<OrderOp> generate_random_flow(std::size_t n_orders, std::uint64_t seed) {
  Xoshiro256StarStar rng(seed);
  std::vector<OrderOp> ops;
  ops.reserve(n_orders);
  std::int64_t mid_est = kInitMid;

  for (std::size_t i = 0; i < n_orders; ++i) {
    const double u = rng.next_uniform();
    OrderOp op;
    op.id = i;
    op.side = (rng.next_uniform() < 0.5) ? Side::Buy : Side::Sell;
    op.target_idx_seed = rng.next_u64();

    if (u < 0.45) {
      op.kind = OrderOp::Limit;
      const std::int64_t depth = rng.next_geometric(0.25);
      const std::int64_t offset = (rng.next_uniform() < 0.10) ? -1 : depth;
      op.price = (op.side == Side::Buy) ? (mid_est - offset) : (mid_est + offset);
      op.price = std::clamp(op.price, kMinPrice + 10, kMaxPrice - 10);
      op.qty = rng.next_geometric(0.2);
    } else if (u < 0.70) {
      op.kind = OrderOp::Cancel;
    } else if (u < 0.90) {
      op.kind = OrderOp::Market;
      op.qty = rng.next_geometric(0.2);
    } else {
      op.kind = OrderOp::Modify;
      op.qty = rng.next_geometric(0.2);
    }
    ops.push_back(op);
  }
  return ops;
}

double now_sec() noexcept {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

} // namespace

bool verify_equivalence_seed(std::size_t n_orders, std::uint64_t seed) {
  const auto ops = generate_random_flow(n_orders, seed);
  OrderBook fast(kMinPrice, kMaxPrice, n_orders + 16);
  ReferenceBook ref;

  std::vector<Trade> fast_trades, ref_trades;
  fast_trades.reserve(64);
  ref_trades.reserve(64);

  for (const auto& op : ops) {
    if (op.kind == OrderOp::Limit) {
      fast.submit_limit(op.id, op.side, op.price, op.qty, fast_trades);
      ref.submit_limit(op.id, op.side, op.price, op.qty, ref_trades);
    } else if (op.kind == OrderOp::Market) {
      fast.submit_market(op.id, op.side, op.qty, fast_trades);
      ref.submit_market(op.id, op.side, op.qty, ref_trades);
    } else if (op.kind == OrderOp::Cancel) {
      if (!fast.active_ids().empty()) {
        const auto& act = fast.active_ids();
        const std::uint64_t tid = act[op.target_idx_seed % act.size()];
        const bool ok1 = fast.cancel(tid);
        const bool ok2 = ref.cancel(tid);
        if (ok1 != ok2) {
          return false;
        }
      }
      fast_trades.clear();
      ref_trades.clear();
    } else {
      if (!fast.active_ids().empty()) {
        const auto& act = fast.active_ids();
        const std::uint64_t tid = act[op.target_idx_seed % act.size()];
        const bool ok1 = fast.modify(tid, op.qty);
        const bool ok2 = ref.modify(tid, op.qty);
        if (ok1 != ok2) {
          return false;
        }
      }
      fast_trades.clear();
      ref_trades.clear();
    }

    if (fast_trades != ref_trades || !fast.check_invariants()) {
      return false;
    }
  }
  return fast.snapshot() == ref.snapshot();
}

OrderFlowBenchmark benchmark_orderbook(std::size_t n_orders, std::uint64_t seed,
                                       bool run_reference) {
  const auto ops = generate_random_flow(n_orders, seed);
  OrderBook fast(kMinPrice, kMaxPrice, n_orders + 16);
  std::vector<Trade> trades;
  trades.reserve(256);

  // Warmup 1000 ops then reset
  for (std::size_t i = 0; i < std::min<std::size_t>(1000, n_orders); ++i) {
    fast.submit_limit(i, ops[i].side, kInitMid + (ops[i].side == Side::Buy ? -1 : 1), 10, trades);
  }
  fast.reset();

  std::vector<double> lats_ns(n_orders, 0.0);
  std::size_t n_trades = 0;
  std::uint64_t t_hash = 0xcbf29ce484222325ULL;

  g_alloc_count.store(0, std::memory_order_relaxed);
  g_track_allocs.store(true, std::memory_order_seq_cst);
  const double t0 = now_sec();

  for (std::size_t i = 0; i < n_orders; ++i) {
    const auto& op = ops[i];
    const auto ts0 = std::chrono::steady_clock::now();
    if (op.kind == OrderOp::Limit) {
      fast.submit_limit(op.id, op.side, op.price, op.qty, trades);
    } else if (op.kind == OrderOp::Market) {
      fast.submit_market(op.id, op.side, op.qty, trades);
    } else if (op.kind == OrderOp::Cancel) {
      if (!fast.active_ids().empty()) {
        const auto& act = fast.active_ids();
        fast.cancel(act[op.target_idx_seed % act.size()]);
      }
      trades.clear();
    } else {
      if (!fast.active_ids().empty()) {
        const auto& act = fast.active_ids();
        fast.modify(act[op.target_idx_seed % act.size()], op.qty);
      }
      trades.clear();
    }
    const auto ts1 = std::chrono::steady_clock::now();
    lats_ns[i] = static_cast<double>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(ts1 - ts0).count());
    n_trades += trades.size();
    t_hash = hash_trades(t_hash, trades);
  }

  const double fast_time = now_sec() - t0;
  g_track_allocs.store(false, std::memory_order_seq_cst);
  const std::size_t allocs = g_alloc_count.load(std::memory_order_relaxed);

  std::vector<double> sorted_lats = lats_ns;
  std::sort(sorted_lats.begin(), sorted_lats.end());
  auto pct = [&](double p) {
    return sorted_lats[static_cast<std::size_t>(p * static_cast<double>(sorted_lats.size() - 1))];
  };

  OrderFlowBenchmark out;
  out.n_orders = n_orders;
  out.n_trades = n_trades;
  out.trade_hash = t_hash;
  out.heap_allocs_after_warmup = allocs;
  out.throughput_ops_per_sec = static_cast<double>(n_orders) / fast_time;
  out.p50_ns = pct(0.50);
  out.p90_ns = pct(0.90);
  out.p99_ns = pct(0.99);
  out.p999_ns = pct(0.999);
  out.fast_time_sec = fast_time;
  out.latencies_ns = std::move(sorted_lats);

  if (run_reference) {
    // Run exact comparison with target IDs recorded from a deterministic OrderBook pass
    OrderBook fast_check(kMinPrice, kMaxPrice, n_orders + 16);
    ReferenceBook ref;
    std::vector<Trade> f_tr, r_tr;
    f_tr.reserve(256);
    r_tr.reserve(256);

    // Pre-resolve target IDs so reference book executes exact same cancel/modify IDs
    std::vector<std::uint64_t> resolved_targets(n_orders, UINT64_MAX);
    for (std::size_t i = 0; i < n_orders; ++i) {
      const auto& op = ops[i];
      if (op.kind == OrderOp::Limit) {
        fast_check.submit_limit(op.id, op.side, op.price, op.qty, f_tr);
      } else if (op.kind == OrderOp::Market) {
        fast_check.submit_market(op.id, op.side, op.qty, f_tr);
      } else if ((op.kind == OrderOp::Cancel || op.kind == OrderOp::Modify) &&
                 !fast_check.active_ids().empty()) {
        const auto& act = fast_check.active_ids();
        const std::uint64_t tid = act[op.target_idx_seed % act.size()];
        resolved_targets[i] = tid;
        if (op.kind == OrderOp::Cancel) {
          fast_check.cancel(tid);
        } else {
          fast_check.modify(tid, op.qty);
        }
      }
    }

    const double t_ref0 = now_sec();
    for (std::size_t i = 0; i < n_orders; ++i) {
      const auto& op = ops[i];
      if (op.kind == OrderOp::Limit) {
        ref.submit_limit(op.id, op.side, op.price, op.qty, r_tr);
      } else if (op.kind == OrderOp::Market) {
        ref.submit_market(op.id, op.side, op.qty, r_tr);
      } else if (op.kind == OrderOp::Cancel && resolved_targets[i] != UINT64_MAX) {
        ref.cancel(resolved_targets[i]);
      } else if (op.kind == OrderOp::Modify && resolved_targets[i] != UINT64_MAX) {
        ref.modify(resolved_targets[i], op.qty);
      }
    }
    out.ref_time_sec = now_sec() - t_ref0;
    out.speedup_vs_ref = out.ref_time_sec / std::max(fast_time, 1e-9);
    out.exact_match = (fast_check.snapshot() == ref.snapshot());
  }
  return out;
}

MarketSimResult run_combined_market_sim(double t_end, std::uint64_t seed) {
  const double t0 = now_sec();
  const HawkesParams hp{0.5, 0.6, 1.2};
  const auto mkt_times = simulate_hawkes_dassios_zhao(hp, t_end, seed);

  const std::size_t max_orders = 500'000;
  OrderBook ob(kMinPrice, kMaxPrice, max_orders);
  std::vector<std::uint8_t> order_trader(max_orders, 0); // 0=MM, 1=Informed
  std::int64_t net_pos[2] = {0, 0};

  std::uint64_t next_id = 0;
  std::vector<Trade> trades;
  trades.reserve(128);

  // Initial book: 20 levels per side with qty 50
  for (std::int64_t lvl = 1; lvl <= 20; ++lvl) {
    order_trader[next_id] = 0;
    ob.submit_limit(next_id++, Side::Buy, kInitMid - lvl, 50, trades);
    order_trader[next_id] = 0;
    ob.submit_limit(next_id++, Side::Sell, kInitMid + lvl, 50, trades);
  }

  Xoshiro256StarStar rng(derive_seed(seed, 1));
  Side last_mkt_side = Side::Buy;
  std::size_t mkt_idx = 0;
  double t = 0.0;

  // Sampling grid for stylized facts (every 1.0 time unit)
  const std::size_t n_bins = static_cast<std::size_t>(t_end);
  std::vector<double> spreads;
  std::vector<double> mid_series;
  std::vector<double> signed_vol_win(n_bins, 0.0);
  std::vector<double> signed_flow;
  std::vector<double> depth_sum(10, 0.0);
  std::size_t depth_samples = 0;
  std::size_t inv_viol = 0;
  std::size_t next_sample_t = 1;

  while (t < t_end && next_id + 10 < max_orders) {
    const double next_mkt_t =
        (mkt_idx < mkt_times.size()) ? mkt_times[mkt_idx] : (t_end + 1.0);
    const double rate_lim = 6.0; // 3.0 per side
    const double rate_can = 0.05 * static_cast<double>(ob.active_order_count());
    const double rate_pois = rate_lim + rate_can;
    const double next_pois_t = t + rng.next_exp(rate_pois);

    if (next_mkt_t <= next_pois_t && next_mkt_t < t_end) {
      t = next_mkt_t;
      ++mkt_idx;
      if (rng.next_uniform() >= 0.7) {
        last_mkt_side = (last_mkt_side == Side::Buy) ? Side::Sell : Side::Buy;
      }
      const std::int64_t qty = rng.next_geometric(0.2); // mean 5
      order_trader[next_id] = 1;
      ob.submit_market(next_id++, last_mkt_side, qty, trades);
      const double sign = (last_mkt_side == Side::Buy) ? 1.0 : -1.0;
      signed_flow.push_back(sign * static_cast<double>(qty));
      const auto bin = std::min(static_cast<std::size_t>(t), n_bins - 1);
      signed_vol_win[bin] += sign * static_cast<double>(qty);

      for (const auto& tr : trades) {
        const auto mt = order_trader[tr.maker_id];
        const auto tt = order_trader[tr.taker_id];
        if (last_mkt_side == Side::Buy) {
          net_pos[tt] += tr.qty;
          net_pos[mt] -= tr.qty;
        } else {
          net_pos[tt] -= tr.qty;
          net_pos[mt] += tr.qty;
        }
      }
    } else if (next_pois_t < t_end) {
      t = next_pois_t;
      if (rng.next_uniform() * rate_pois < rate_lim) {
        const Side s = (rng.next_uniform() < 0.5) ? Side::Buy : Side::Sell;
        const std::int64_t ref_p =
            (s == Side::Buy)
                ? (ob.best_ask() <= kMaxPrice ? ob.best_ask() : kInitMid)
                : (ob.best_bid() >= kMinPrice ? ob.best_bid() : kInitMid);
        const std::int64_t depth = rng.next_geometric(0.2);
        const std::int64_t price =
            std::clamp((s == Side::Buy) ? (ref_p - depth) : (ref_p + depth), kMinPrice + 10, kMaxPrice - 10);
        const std::int64_t qty = rng.next_geometric(0.1);
        order_trader[next_id] = 0;
        ob.submit_limit(next_id++, s, price, qty, trades);
        for (const auto& tr : trades) {
          const auto mt = order_trader[tr.maker_id];
          const auto tt = order_trader[tr.taker_id];
          if (s == Side::Buy) {
            net_pos[tt] += tr.qty;
            net_pos[mt] -= tr.qty;
          } else {
            net_pos[tt] -= tr.qty;
            net_pos[mt] += tr.qty;
          }
        }
      } else if (!ob.active_ids().empty()) {
        const auto& act = ob.active_ids();
        ob.cancel(act[rng.next_u64() % act.size()]);
      }
    } else {
      break;
    }

    while (t >= static_cast<double>(next_sample_t) && next_sample_t <= n_bins) {
      if (ob.best_bid() >= kMinPrice && ob.best_ask() <= kMaxPrice) {
        spreads.push_back(static_cast<double>(ob.best_ask() - ob.best_bid()));
        mid_series.push_back(0.5 * static_cast<double>(ob.best_ask() + ob.best_bid()));
      }
      if (next_sample_t % 100 == 0) {
        if (!ob.check_invariants() || (net_pos[0] + net_pos[1] != 0)) {
          ++inv_viol;
        }
        const auto snap = ob.snapshot();
        for (const auto& lv : snap) {
          const std::int64_t dist =
              (lv.side == Side::Buy) ? (ob.best_bid() - lv.price) : (lv.price - ob.best_ask());
          if (dist >= 0 && dist < 10) {
            depth_sum[static_cast<std::size_t>(dist)] += 0.5 * static_cast<double>(lv.qty);
          }
        }
        ++depth_samples;
      }
      ++next_sample_t;
    }
  }

  std::vector<double> abs_ret;
  std::vector<double> delta_mid;
  for (std::size_t i = 1; i < mid_series.size(); ++i) {
    const double dm = mid_series[i] - mid_series[i - 1];
    delta_mid.push_back(dm);
    abs_ret.push_back(std::abs(dm));
  }

  // Price impact regression: delta_mid vs signed_vol_win[1..]
  const std::size_t m_reg = std::min(delta_mid.size(), signed_vol_win.size() - 1);
  double mx = 0.0, my = 0.0;
  for (std::size_t i = 0; i < m_reg; ++i) {
    mx += signed_vol_win[i + 1];
    my += delta_mid[i];
  }
  mx /= static_cast<double>(m_reg);
  my /= static_cast<double>(m_reg);
  double sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (std::size_t i = 0; i < m_reg; ++i) {
    const double dx = signed_vol_win[i + 1] - mx;
    const double dy = delta_mid[i] - my;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  const double slope = (sxx > 0.0) ? (sxy / sxx) : 0.0;
  const double r2 = (sxx > 0.0 && syy > 0.0) ? ((sxy * sxy) / (sxx * syy)) : 0.0;

  const auto win_c = stats::window_counts(mkt_times, t_end, 50.0);
  const double m_rate = static_cast<double>(mkt_times.size()) / t_end;
  const double m_se = std::sqrt(stats::variance(win_c) / static_cast<double>(win_c.size())) / 50.0;

  for (double& d : depth_sum) {
    if (depth_samples > 0) {
      d /= static_cast<double>(depth_samples);
    }
  }

  MarketSimResult res;
  res.n_market_orders = mkt_times.size();
  res.emp_market_rate = m_rate;
  res.emp_market_rate_se = m_se;
  res.fano_market_orders = stats::fano_factor(mkt_times, t_end, 50.0);
  res.mean_spread = stats::mean(spreads);
  res.std_spread = std::sqrt(stats::variance(spreads));
  res.depth_profile = std::move(depth_sum);
  res.signed_flow_acf = stats::acf(signed_flow, 20);
  res.abs_ret_acf = stats::acf(abs_ret, 50);
  res.price_impact_slope = slope;
  res.price_impact_r2 = r2;
  res.invariant_violations = inv_viol;
  res.net_position_sum = net_pos[0] + net_pos[1];
  res.time_sec = now_sec() - t0;
  return res;
}

} // namespace hob
