#include "aad/mc_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <thread>
#include <vector>

#include "aad/cholesky.hpp"
#include "aad/models.hpp"
#include "aad/number.hpp"
#include "aad/rng.hpp"
#include "aad/tape.hpp"

namespace aad {

namespace {

constexpr std::size_t kChunkSize = 4000;

struct WelfordAcc {
  std::size_t count{0};
  double mean{0.0};
  double m2{0.0};

  void push(double x) noexcept {
    ++count;
    const double d = x - mean;
    mean += d / static_cast<double>(count);
    m2 += d * (x - mean);
  }

  void merge(const WelfordAcc& other) noexcept {
    if (other.count == 0) {
      return;
    }
    if (count == 0) {
      *this = other;
      return;
    }
    const std::size_t total = count + other.count;
    const double delta = other.mean - mean;
    const double scale =
        (static_cast<double>(count) * static_cast<double>(other.count)) / static_cast<double>(total);
    mean += delta * (static_cast<double>(other.count) / static_cast<double>(total));
    m2 += other.m2 + delta * delta * scale;
    count = total;
  }

  [[nodiscard]] McEstimate to_estimate() const noexcept {
    const double n = static_cast<double>(count);
    const double var = (count > 1) ? std::max(m2 / (n - 1.0), 0.0) : 0.0;
    const double se = std::sqrt(var / n);
    return {mean, se, mean - 1.96 * se, mean + 1.96 * se};
  }
};

std::size_t resolve_threads(std::size_t n_threads) noexcept {
  if (n_threads > 0) {
    return n_threads;
  }
  const unsigned hw = std::thread::hardware_concurrency();
  return hw > 0 ? static_cast<std::size_t>(hw) : 1;
}

template <typename ChunkFn>
auto run_chunked(std::size_t n_paths, std::size_t n_threads, ChunkFn&& fn) {
  const std::size_t n_chunks = (n_paths + kChunkSize - 1) / kChunkSize;
  const std::size_t eff_threads = std::min(resolve_threads(n_threads), n_chunks);
  using ResultT = decltype(fn(std::size_t{0}, std::size_t{0}));
  std::vector<ResultT> results(n_chunks);

  if (eff_threads <= 1) {
    for (std::size_t c = 0; c < n_chunks; ++c) {
      const std::size_t len = std::min(kChunkSize, n_paths - c * kChunkSize);
      results[c] = fn(c, len);
    }
    return results;
  }

  std::vector<std::thread> workers;
  workers.reserve(eff_threads);
  for (std::size_t wid = 0; wid < eff_threads; ++wid) {
    workers.emplace_back([&, wid]() {
      for (std::size_t c = wid; c < n_chunks; c += eff_threads) {
        const std::size_t len = std::min(kChunkSize, n_paths - c * kChunkSize);
        results[c] = fn(c, len);
      }
    });
  }
  for (auto& w : workers) {
    w.join();
  }
  return results;
}

double now_sec() noexcept {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double>(clock::now().time_since_epoch()).count();
}

} // namespace

Gbm1StepResult run_gbm_1step_aad(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads) {
  const double t0 = now_sec();
  struct ChunkOut {
    WelfordAcc acc[5];
  };

  auto chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
    Tape tape;
    TapeGuard guard(tape);
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));

    const Number s0(100.0);
    const Number k(100.0);
    const Number r(0.03);
    const Number sigma(0.20);
    const Number t(1.0);
    const Number disc = math::exp(-r * t);
    const Number drift = (r - Number(0.5) * sigma * sigma) * t;
    const Number vol = sigma * math::sqrt(t);
    tape.mark();

    ChunkOut out;
    for (std::size_t p = 0; p < len; ++p) {
      const double z = rng.next_normal();
      const Number s_t = s0 * math::exp(drift + vol * Number(z));
      const Number payoff = disc * math::max(s_t - k, Number(0.0));
      tape.backprop_to_mark(payoff.index());
      tape.backprop_from_mark();

      out.acc[0].push(payoff.value());
      out.acc[1].push(s0.adjoint());
      out.acc[2].push(sigma.adjoint());
      out.acc[3].push(r.adjoint());
      out.acc[4].push(t.adjoint());
      tape.zero_all_adjoints();
      tape.rewind_to_mark();
    }
    return out;
  });

  WelfordAcc merged[5];
  for (const auto& ch : chunks) {
    for (int i = 0; i < 5; ++i) {
      merged[i].merge(ch.acc[i]);
    }
  }
  return {merged[0].to_estimate(), merged[1].to_estimate(), merged[2].to_estimate(),
          merged[3].to_estimate(), merged[4].to_estimate(), now_sec() - t0};
}

AsianAadResult run_asian_aad(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads) {
  const double t0 = now_sec();
  constexpr std::size_t n_steps = 252;
  struct ChunkOut {
    WelfordAcc acc[6]; // price, delta, delta_formula, vega, vega_bump, rho
    std::size_t nodes_per_path{0};
  };

  auto chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
    Tape tape;
    tape.reserve(8192);
    TapeGuard guard(tape);
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));

    const Number s0(100.0);
    const Number k(100.0);
    const Number r(0.05);
    const Number sigma(0.20);
    const Number t(1.0);
    const Number dt = t / Number(static_cast<double>(n_steps));
    const Number disc = math::exp(-r * t);
    const Number drift = (r - Number(0.5) * sigma * sigma) * dt;
    const Number vol = sigma * math::sqrt(dt);
    tape.mark();

    constexpr double h_sig = 1e-4;
    const double dt_d = 1.0 / static_cast<double>(n_steps);
    const double disc_d = std::exp(-0.05);
    const double drift_p = (0.05 - 0.5 * (0.20 + h_sig) * (0.20 + h_sig)) * dt_d;
    const double drift_m = (0.05 - 0.5 * (0.20 - h_sig) * (0.20 - h_sig)) * dt_d;
    const double vol_p = (0.20 + h_sig) * std::sqrt(dt_d);
    const double vol_m = (0.20 - h_sig) * std::sqrt(dt_d);

    std::vector<double> zs(n_steps);
    ChunkOut out;
    for (std::size_t p = 0; p < len; ++p) {
      Number s = s0;
      Number sum_s(0.0);
      for (std::size_t step = 0; step < n_steps; ++step) {
        const double z = rng.next_normal();
        zs[step] = z;
        s = s * math::exp(drift + vol * Number(z));
        sum_s = sum_s + s;
      }
      const Number avg = sum_s / Number(static_cast<double>(n_steps));
      const Number payoff = disc * math::max(avg - k, Number(0.0));
      if (p == 0) {
        out.nodes_per_path = tape.size() - tape.mark_index();
      }
      tape.backprop_to_mark(payoff.index());
      tape.backprop_from_mark();

      const double df = (avg.value() > 100.0) ? (disc_d * avg.value() / 100.0) : 0.0;

      // CRN central bump for vega
      double sp = 100.0, sm = 100.0, sump = 0.0, summ = 0.0;
      for (std::size_t step = 0; step < n_steps; ++step) {
        sp *= std::exp(drift_p + vol_p * zs[step]);
        sm *= std::exp(drift_m + vol_m * zs[step]);
        sump += sp;
        summ += sm;
      }
      const double vp = disc_d * std::max(sump / static_cast<double>(n_steps) - 100.0, 0.0);
      const double vm = disc_d * std::max(summ / static_cast<double>(n_steps) - 100.0, 0.0);

      out.acc[0].push(payoff.value());
      out.acc[1].push(s0.adjoint());
      out.acc[2].push(df);
      out.acc[3].push(sigma.adjoint());
      out.acc[4].push((vp - vm) / (2.0 * h_sig));
      out.acc[5].push(r.adjoint());
      tape.zero_all_adjoints();
      tape.rewind_to_mark();
    }
    return out;
  });

  WelfordAcc merged[6];
  std::size_t nodes = 0;
  for (const auto& ch : chunks) {
    nodes = std::max(nodes, ch.nodes_per_path);
    for (int i = 0; i < 6; ++i) {
      merged[i].merge(ch.acc[i]);
    }
  }
  const double dt_s = now_sec() - t0;
  return {merged[0].to_estimate(),
          merged[1].to_estimate(),
          merged[2].to_estimate(),
          merged[3].to_estimate(),
          merged[4].to_estimate(),
          merged[5].to_estimate(),
          nodes,
          dt_s,
          static_cast<double>(n_paths) / dt_s};
}

DigitalResult run_digital_aad(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads) {
  const double an_delta = digital_call_analytic_delta(100.0, 100.0, 0.03, 0.20, 1.0);
  const std::vector<double> eps_list{0.1, 0.5, 1.0, 2.0, 5.0};

  struct ChunkOut {
    WelfordAcc raw_acc;
    WelfordAcc eps_acc[5];
  };

  auto chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
    Tape tape;
    TapeGuard guard(tape);
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));

    const Number s0(100.0);
    const Number k(100.0);
    const Number r(0.03);
    const Number sigma(0.20);
    const Number t(1.0);
    const Number disc = math::exp(-r * t);
    const Number drift = (r - Number(0.5) * sigma * sigma) * t;
    const Number vol = sigma * math::sqrt(t);
    tape.mark();

    ChunkOut out;
    for (std::size_t p = 0; p < len; ++p) {
      const double z = rng.next_normal();
      const Number s_t = s0 * math::exp(drift + vol * Number(z));
      // Discontinuous payoff indicator has zero pathwise derivative everywhere
      out.raw_acc.push(0.0);
      for (std::size_t e = 0; e < eps_list.size(); ++e) {
        const Number pay = disc * math::normal_cdf((s_t - k) / Number(eps_list[e]));
        tape.backprop(pay.index());
        out.eps_acc[e].push(s0.adjoint());
      }
      tape.zero_all_adjoints();
      tape.rewind_to_mark();
    }
    return out;
  });

  WelfordAcc raw_m;
  WelfordAcc eps_m[5];
  for (const auto& ch : chunks) {
    raw_m.merge(ch.raw_acc);
    for (std::size_t e = 0; e < 5; ++e) {
      eps_m[e].merge(ch.eps_acc[e]);
    }
  }

  DigitalResult res;
  res.analytic_delta = an_delta;
  res.pathwise_raw_delta = raw_m.mean;
  double best_mse = 1e100;
  for (std::size_t e = 0; e < 5; ++e) {
    const McEstimate est = eps_m[e].to_estimate();
    const double bias = est.mean - an_delta;
    const double mse = bias * bias + est.se * est.se;
    res.smoothed.push_back({eps_list[e], est, bias, mse});
    if (mse < best_mse) {
      best_mse = mse;
      res.optimal_eps = eps_list[e];
    }
  }
  return res;
}

namespace {

constexpr std::size_t kDim = 10;
constexpr std::size_t kNRhos = 45;
constexpr std::size_t kNParams = 67; // 10 S0 + 10 sigma + r + T + 45 rho

template <bool Geometric>
double eval_basket_double(const double* s0, const double* sig, double r, double t,
                          const std::vector<std::vector<double>>& chol, const double* z) noexcept {
  double acc = Geometric ? 1.0 : 0.0;
  for (std::size_t i = 0; i < kDim; ++i) {
    double corr_z = 0.0;
    for (std::size_t j = 0; j <= i; ++j) {
      corr_z += chol[i][j] * z[j];
    }
    const double s_i =
        s0[i] * std::exp((r - 0.5 * sig[i] * sig[i]) * t + sig[i] * std::sqrt(t) * corr_z);
    if constexpr (Geometric) {
      acc *= std::pow(s_i, 0.1);
    } else {
      acc += 0.1 * s_i;
    }
  }
  return std::exp(-r * t) * std::max(acc - 100.0, 0.0);
}

template <bool Geometric>
BasketComparison run_basket_impl(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads,
                                 bool run_fd) {
  std::vector<double> sig_base(kDim);
  for (std::size_t i = 0; i < kDim; ++i) {
    sig_base[i] = 0.15 + 0.20 * static_cast<double>(i) / static_cast<double>(kDim - 1);
  }
  const std::vector<double> rho_base(kNRhos, 0.3);

  // 1. AAD run
  const double t_aad0 = now_sec();
  struct ChunkOut {
    WelfordAcc grads[kNParams];
    std::size_t nodes_per_path{0};
    std::size_t peak_bytes{0};
  };

  auto chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
    Tape tape;
    tape.reserve(2048);
    TapeGuard guard(tape);
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));

    std::vector<Number> s0(kDim), sig(kDim), rhos(kNRhos);
    for (std::size_t i = 0; i < kDim; ++i) {
      s0[i] = Number(100.0);
    }
    for (std::size_t i = 0; i < kDim; ++i) {
      sig[i] = Number(sig_base[i]);
    }
    const Number r(0.03);
    const Number t(1.0);
    for (std::size_t k = 0; k < kNRhos; ++k) {
      rhos[k] = Number(rho_base[k]);
    }

    const auto chol = cholesky_from_rhos(rhos, kDim);
    const Number sqrt_t = math::sqrt(t);
    const Number disc = math::exp(-r * t);
    std::vector<Number> drifts(kDim), vols(kDim);
    for (std::size_t i = 0; i < kDim; ++i) {
      drifts[i] = (r - Number(0.5) * sig[i] * sig[i]) * t;
      vols[i] = sig[i] * sqrt_t;
    }
    tape.mark();

    ChunkOut out;
    double z[kDim]{};
    for (std::size_t p = 0; p < len; ++p) {
      for (std::size_t i = 0; i < kDim; ++i) {
        z[i] = rng.next_normal();
      }
      Number basket = Geometric ? Number(1.0) : Number(0.0);
      for (std::size_t i = 0; i < kDim; ++i) {
        Number corr_z(0.0);
        for (std::size_t j = 0; j <= i; ++j) {
          corr_z = corr_z + chol[i][j] * Number(z[j]);
        }
        const Number s_i = s0[i] * math::exp(drifts[i] + vols[i] * corr_z);
        if constexpr (Geometric) {
          basket = basket * math::pow(s_i, 0.1);
        } else {
          basket = basket + Number(0.1) * s_i;
        }
      }
      const Number payoff = disc * math::max(basket - Number(100.0), Number(0.0));
      if (p == 0) {
        out.nodes_per_path = tape.size() - tape.mark_index();
        out.peak_bytes = tape.capacity_bytes();
      }
      tape.backprop_to_mark(payoff.index());
      tape.backprop_from_mark();

      for (std::size_t i = 0; i < kDim; ++i) {
        out.grads[i].push(s0[i].adjoint());
      }
      for (std::size_t i = 0; i < kDim; ++i) {
        out.grads[kDim + i].push(sig[i].adjoint());
      }
      out.grads[20].push(r.adjoint());
      out.grads[21].push(t.adjoint());
      for (std::size_t k = 0; k < kNRhos; ++k) {
        out.grads[22 + k].push(rhos[k].adjoint());
      }
      tape.zero_all_adjoints();
      tape.rewind_to_mark();
    }
    return out;
  });
  const double aad_time = now_sec() - t_aad0;

  WelfordAcc merged[kNParams];
  std::size_t nodes_per_path = 0;
  std::size_t peak_bytes = 0;
  for (const auto& ch : chunks) {
    nodes_per_path = std::max(nodes_per_path, ch.nodes_per_path);
    peak_bytes = std::max(peak_bytes, ch.peak_bytes);
    for (std::size_t k = 0; k < kNParams; ++k) {
      merged[k].merge(ch.grads[k]);
    }
  }

  // 2. Value-only timing (for grad/value ratio)
  const double t_val0 = now_sec();
  const auto chol_d = cholesky_from_rhos(rho_base, kDim);
  double s0_d[kDim];
  std::fill(std::begin(s0_d), std::end(s0_d), 100.0);
  auto val_chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
    Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));
    double sum = 0.0;
    double z[kDim]{};
    for (std::size_t p = 0; p < len; ++p) {
      for (std::size_t i = 0; i < kDim; ++i) {
        z[i] = rng.next_normal();
      }
      sum += eval_basket_double<Geometric>(s0_d, sig_base.data(), 0.03, 1.0, chol_d, z);
    }
    return sum;
  });
  volatile double sink = std::accumulate(val_chunks.begin(), val_chunks.end(), 0.0);
  (void)sink;
  const double val_time = now_sec() - t_val0;

  BasketComparison cmp;
  cmp.nodes_per_path = nodes_per_path;
  cmp.peak_tape_bytes = peak_bytes;
  cmp.aad_time_sec = aad_time;
  cmp.value_only_time_sec = val_time;
  cmp.grad_to_val_ratio = aad_time / std::max(val_time, 1e-9);
  cmp.paths_per_sec = static_cast<double>(n_paths) / aad_time;
  for (std::size_t k = 0; k < kNParams; ++k) {
    cmp.aad_grads.push_back(merged[k].to_estimate());
  }

  if constexpr (Geometric) {
    // Exact AAD on closed-form geometric basket
    Tape tape;
    TapeGuard guard(tape);
    std::vector<Number> s0(kDim), sig(kDim), rhos(kNRhos);
    for (std::size_t i = 0; i < kDim; ++i) {
      s0[i] = Number(100.0);
      sig[i] = Number(sig_base[i]);
    }
    const Number r(0.03), t(1.0), k_num(100.0);
    for (std::size_t idx = 0; idx < kNRhos; ++idx) {
      rhos[idx] = Number(rho_base[idx]);
    }
    const std::vector<double> w(kDim, 0.1);
    const Number c_val = geom_basket_closed_form(s0, sig, rhos, w, k_num, r, t);
    tape.backprop(c_val.index());

    auto push_exact = [&](double v) { cmp.ref_grads.push_back({v, 0.0, v, v}); };
    for (std::size_t i = 0; i < kDim; ++i) {
      push_exact(s0[i].adjoint());
    }
    for (std::size_t i = 0; i < kDim; ++i) {
      push_exact(sig[i].adjoint());
    }
    push_exact(r.adjoint());
    push_exact(t.adjoint());
    for (std::size_t idx = 0; idx < kNRhos; ++idx) {
      push_exact(rhos[idx].adjoint());
    }
  } else if (run_fd) {
    // Central FD + CRN (2 * 67 re-evaluations)
    const double t_fd0 = now_sec();
    constexpr double h_s = 1e-2;
    constexpr double h_sig = 1e-4;
    constexpr double h_r = 1e-4;
    constexpr double h_t = 1e-4;
    constexpr double h_rho = 1e-4;

    std::vector<std::vector<std::vector<double>>> chols_p(kNRhos), chols_m(kNRhos);
    for (std::size_t k = 0; k < kNRhos; ++k) {
      auto rp = rho_base;
      auto rm = rho_base;
      rp[k] += h_rho;
      rm[k] -= h_rho;
      chols_p[k] = cholesky_from_rhos(rp, kDim);
      chols_m[k] = cholesky_from_rhos(rm, kDim);
    }

    struct FdChunk {
      WelfordAcc grads[kNParams];
    };
    auto fd_chunks = run_chunked(n_paths, n_threads, [&](std::size_t c_idx, std::size_t len) {
      Xoshiro256StarStar rng(derive_chunk_seed(seed, c_idx));
      FdChunk out;
      double z[kDim]{};
      double s0_tmp[kDim]{};
      double sig_tmp[kDim]{};
      for (std::size_t p = 0; p < len; ++p) {
        for (std::size_t i = 0; i < kDim; ++i) {
          z[i] = rng.next_normal();
        }
        std::copy(std::begin(s0_d), std::end(s0_d), std::begin(s0_tmp));
        std::copy(sig_base.begin(), sig_base.end(), std::begin(sig_tmp));

        for (std::size_t i = 0; i < kDim; ++i) {
          s0_tmp[i] = 100.0 + h_s;
          const double vp = eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chol_d, z);
          s0_tmp[i] = 100.0 - h_s;
          const double vm = eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chol_d, z);
          s0_tmp[i] = 100.0;
          out.grads[i].push((vp - vm) / (2.0 * h_s));
        }
        for (std::size_t i = 0; i < kDim; ++i) {
          sig_tmp[i] = sig_base[i] + h_sig;
          const double vp = eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chol_d, z);
          sig_tmp[i] = sig_base[i] - h_sig;
          const double vm = eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chol_d, z);
          sig_tmp[i] = sig_base[i];
          out.grads[kDim + i].push((vp - vm) / (2.0 * h_sig));
        }
        const double rp_v =
            eval_basket_double<false>(s0_tmp, sig_tmp, 0.03 + h_r, 1.0, chol_d, z);
        const double rm_v =
            eval_basket_double<false>(s0_tmp, sig_tmp, 0.03 - h_r, 1.0, chol_d, z);
        out.grads[20].push((rp_v - rm_v) / (2.0 * h_r));

        const double tp_v =
            eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0 + h_t, chol_d, z);
        const double tm_v =
            eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0 - h_t, chol_d, z);
        out.grads[21].push((tp_v - tm_v) / (2.0 * h_t));

        for (std::size_t k = 0; k < kNRhos; ++k) {
          const double vp =
              eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chols_p[k], z);
          const double vm =
              eval_basket_double<false>(s0_tmp, sig_tmp, 0.03, 1.0, chols_m[k], z);
          out.grads[22 + k].push((vp - vm) / (2.0 * h_rho));
        }
      }
      return out;
    });
    cmp.fd_time_sec = now_sec() - t_fd0;
    cmp.speedup_vs_fd = cmp.fd_time_sec / std::max(cmp.aad_time_sec, 1e-9);

    WelfordAcc fd_merged[kNParams];
    for (const auto& ch : fd_chunks) {
      for (std::size_t k = 0; k < kNParams; ++k) {
        fd_merged[k].merge(ch.grads[k]);
      }
    }
    for (std::size_t k = 0; k < kNParams; ++k) {
      cmp.ref_grads.push_back(fd_merged[k].to_estimate());
    }
  }

  if (!cmp.ref_grads.empty()) {
    for (std::size_t k = 0; k < kNParams; ++k) {
      const double se_comb = std::hypot(cmp.aad_grads[k].se, cmp.ref_grads[k].se);
      const double z_val =
          (se_comb > 0.0) ? std::abs(cmp.aad_grads[k].mean - cmp.ref_grads[k].mean) / se_comb : 0.0;
      cmp.diff_in_se.push_back(z_val);
      cmp.max_diff_se = std::max(cmp.max_diff_se, z_val);
      if (z_val > 3.0) {
        ++cmp.violations_3se;
      }
      if (z_val > 4.0) {
        ++cmp.violations_4se;
      }
    }
  }
  return cmp;
}

} // namespace

BasketComparison run_basket_10_arith(std::size_t n_paths, std::uint64_t seed, std::size_t n_threads,
                                     bool run_fd) {
  return run_basket_impl<false>(n_paths, seed, n_threads, run_fd);
}

BasketComparison run_basket_10_geom(std::size_t n_paths, std::uint64_t seed,
                                    std::size_t n_threads) {
  return run_basket_impl<true>(n_paths, seed, n_threads, false);
}

} // namespace aad
