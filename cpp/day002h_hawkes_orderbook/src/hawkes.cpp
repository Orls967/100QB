#include "hawkes_ob/hawkes.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <thread>

#include "hawkes_ob/rng.hpp"
#include "hawkes_ob/stats.hpp"

namespace hob {

std::vector<double> simulate_hawkes_ogata(const HawkesParams& p, double t_end, std::uint64_t seed) {
  Xoshiro256StarStar rng(seed);
  std::vector<double> times;
  times.reserve(static_cast<std::size_t>(t_end * (p.mu / std::max(1.0 - p.alpha / p.beta, 0.1)) * 1.2));
  double t = 0.0;
  double ex = 0.0;
  while (t < t_end) {
    const double lam_bar = p.mu + ex;
    const double dt = rng.next_exp(lam_bar);
    t += dt;
    if (t >= t_end) {
      break;
    }
    ex *= std::exp(-p.beta * dt);
    const double lam_t = p.mu + ex;
    if (rng.next_uniform() * lam_bar <= lam_t) {
      times.push_back(t);
      ex += p.alpha;
    }
  }
  return times;
}

std::vector<double> simulate_hawkes_dassios_zhao(const HawkesParams& p, double t_end,
                                                 std::uint64_t seed) {
  Xoshiro256StarStar rng(seed);
  std::vector<double> times;
  times.reserve(static_cast<std::size_t>(t_end * (p.mu / std::max(1.0 - p.alpha / p.beta, 0.1)) * 1.2));
  double t = 0.0;
  double ex = 0.0;
  while (t < t_end) {
    const double s1 = rng.next_exp(p.mu);
    double s2 = std::numeric_limits<double>::infinity();
    if (ex > 0.0) {
      const double u = std::max(rng.next_uniform(), 1e-16);
      const double d = 1.0 + p.beta * std::log(u) / ex;
      if (d > 0.0) {
        s2 = -std::log(d) / p.beta;
      }
    }
    const double dt = std::min(s1, s2);
    t += dt;
    if (t >= t_end) {
      break;
    }
    times.push_back(t);
    ex = ex * std::exp(-p.beta * dt) + p.alpha;
  }
  return times;
}

std::vector<double> inter_event_times(const std::vector<double>& times) {
  std::vector<double> dts;
  if (times.empty()) {
    return dts;
  }
  dts.reserve(times.size());
  double prev = 0.0;
  for (double t : times) {
    dts.push_back(t - prev);
    prev = t;
  }
  return dts;
}

std::vector<double> time_rescaling_residuals(const std::vector<double>& times,
                                             const HawkesParams& p) {
  std::vector<double> res;
  res.reserve(times.size());
  double t_prev = 0.0;
  double ex_prev = 0.0;
  for (double t : times) {
    const double dt = t - t_prev;
    const double decay = std::exp(-p.beta * dt);
    const double lam = p.mu * dt + (ex_prev / p.beta) * (1.0 - decay);
    res.push_back(lam);
    ex_prev = ex_prev * decay + p.alpha;
    t_prev = t;
  }
  return res;
}

namespace {

void ll_and_grad_log_params(const std::vector<double>& times, double t_end, double u_mu,
                            double u_alpha, double u_beta, double& out_nll,
                            std::array<double, 3>& out_grad) {
  const double mu = std::exp(u_mu);
  const double alpha = std::exp(u_alpha);
  const double beta = std::exp(u_beta);
  const std::size_t n = times.size();

  double r_k = 0.0;
  double dr_dbeta = 0.0;
  double sum_log_lam = 0.0;
  double d_ll_dmu = -t_end;
  double d_ll_dalpha = 0.0;
  double d_ll_dbeta = 0.0;

  double sum_comp = 0.0;
  double d_comp_dbeta = 0.0;

  for (std::size_t i = 0; i < n; ++i) {
    if (i > 0) {
      const double dt = times[i] - times[i - 1];
      const double e_dt = std::exp(-beta * dt);
      const double prev_1pr = 1.0 + r_k;
      r_k = e_dt * prev_1pr;
      dr_dbeta = e_dt * (dr_dbeta - dt * prev_1pr);
    }
    const double lam = mu + alpha * r_k;
    sum_log_lam += std::log(lam);
    d_ll_dmu += 1.0 / lam;
    d_ll_dalpha += r_k / lam;
    d_ll_dbeta += (alpha * dr_dbeta) / lam;

    const double dt_end = t_end - times[i];
    const double e_end = std::exp(-beta * dt_end);
    sum_comp += 1.0 - e_end;
    d_comp_dbeta += dt_end * e_end;
  }

  const double ll = sum_log_lam - mu * t_end - (alpha / beta) * sum_comp;
  d_ll_dalpha -= sum_comp / beta;
  d_ll_dbeta -= alpha * (d_comp_dbeta * beta - sum_comp) / (beta * beta);

  out_nll = -ll;
  out_grad[0] = -d_ll_dmu * mu;
  out_grad[1] = -d_ll_dalpha * alpha;
  out_grad[2] = -d_ll_dbeta * beta;
}

bool invert_3x3(const std::array<std::array<double, 3>, 3>& a,
                std::array<std::array<double, 3>, 3>& inv) noexcept {
  const double c00 = a[1][1] * a[2][2] - a[1][2] * a[2][1];
  const double c01 = -(a[1][0] * a[2][2] - a[1][2] * a[2][0]);
  const double c02 = a[1][0] * a[2][1] - a[1][1] * a[2][0];
  const double det = a[0][0] * c00 + a[0][1] * c01 + a[0][2] * c02;
  if (std::abs(det) < 1e-18) {
    return false;
  }
  const double idet = 1.0 / det;
  inv[0][0] = c00 * idet;
  inv[1][0] = c01 * idet;
  inv[2][0] = c02 * idet;
  inv[0][1] = -(a[0][1] * a[2][2] - a[0][2] * a[2][1]) * idet;
  inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * idet;
  inv[2][1] = -(a[0][0] * a[2][1] - a[0][1] * a[2][0]) * idet;
  inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * idet;
  inv[1][2] = -(a[0][0] * a[1][2] - a[0][2] * a[1][0]) * idet;
  inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * idet;
  return true;
}

} // namespace

double hawkes_log_likelihood(const std::vector<double>& times, double t_end,
                             const HawkesParams& p) {
  double nll = 0.0;
  std::array<double, 3> g{};
  ll_and_grad_log_params(times, t_end, std::log(p.mu), std::log(p.alpha), std::log(p.beta), nll, g);
  return -nll;
}

HawkesMleFit fit_hawkes_mle(const std::vector<double>& times, double t_end) {
  const double rate = static_cast<double>(times.size()) / t_end;
  std::array<double, 3> x{std::log(0.5 * rate), std::log(0.5), std::log(1.0)};
  std::array<std::array<double, 3>, 3> h_inv{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};

  double f = 0.0;
  std::array<double, 3> g{};
  ll_and_grad_log_params(times, t_end, x[0], x[1], x[2], f, g);

  for (int it = 0; it < 100; ++it) {
    const double gnorm = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    if (gnorm < 1e-8) {
      break;
    }
    std::array<double, 3> p_dir{};
    for (int i = 0; i < 3; ++i) {
      p_dir[i] = -(h_inv[i][0] * g[0] + h_inv[i][1] * g[1] + h_inv[i][2] * g[2]);
    }
    double slope = p_dir[0] * g[0] + p_dir[1] * g[1] + p_dir[2] * g[2];
    if (slope >= 0.0) {
      p_dir = {-g[0], -g[1], -g[2]};
      slope = -(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    }

    double step = 1.0;
    std::array<double, 3> x_next{};
    double f_next = 0.0;
    std::array<double, 3> g_next{};
    for (int ls = 0; ls < 30; ++ls) {
      for (int i = 0; i < 3; ++i) {
        x_next[i] = std::clamp(x[i] + step * p_dir[i], -8.0, 6.0);
      }
      ll_and_grad_log_params(times, t_end, x_next[0], x_next[1], x_next[2], f_next, g_next);
      if (f_next <= f + 1e-4 * step * slope) {
        break;
      }
      step *= 0.5;
    }

    std::array<double, 3> s_vec{x_next[0] - x[0], x_next[1] - x[1], x_next[2] - x[2]};
    std::array<double, 3> y_vec{g_next[0] - g[0], g_next[1] - g[1], g_next[2] - g[2]};
    const double sy = s_vec[0] * y_vec[0] + s_vec[1] * y_vec[1] + s_vec[2] * y_vec[2];
    if (sy > 1e-12) {
      const double rho = 1.0 / sy;
      std::array<std::array<double, 3>, 3> i_sy{};
      std::array<std::array<double, 3>, 3> i_ys{};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          i_sy[i][j] = (i == j ? 1.0 : 0.0) - rho * s_vec[i] * y_vec[j];
          i_ys[i][j] = (i == j ? 1.0 : 0.0) - rho * y_vec[i] * s_vec[j];
        }
      }
      std::array<std::array<double, 3>, 3> tmp{};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          tmp[i][j] = i_sy[i][0] * h_inv[0][j] + i_sy[i][1] * h_inv[1][j] + i_sy[i][2] * h_inv[2][j];
        }
      }
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          h_inv[i][j] = tmp[i][0] * i_ys[0][j] + tmp[i][1] * i_ys[1][j] + tmp[i][2] * i_ys[2][j] +
                        rho * s_vec[i] * s_vec[j];
        }
      }
    }
    x = x_next;
    f = f_next;
    g = g_next;
  }

  // Numerical Hessian of NLL in (mu, alpha, beta) for asymptotic SE
  const HawkesParams p_hat{std::exp(x[0]), std::exp(x[1]), std::exp(x[2])};
  const double theta[3] = {p_hat.mu, p_hat.alpha, p_hat.beta};
  const double eps[3] = {1e-4 * p_hat.mu, 1e-4 * p_hat.alpha, 1e-4 * p_hat.beta};
  auto nll_orig = [&](double m, double a, double b) {
    return -hawkes_log_likelihood(times, t_end, {m, a, b});
  };
  std::array<std::array<double, 3>, 3> hess{};
  const double f0 = nll_orig(theta[0], theta[1], theta[2]);
  for (int i = 0; i < 3; ++i) {
    double tp[3] = {theta[0], theta[1], theta[2]};
    double tm[3] = {theta[0], theta[1], theta[2]};
    tp[i] += eps[i];
    tm[i] -= eps[i];
    hess[i][i] = (nll_orig(tp[0], tp[1], tp[2]) - 2.0 * f0 + nll_orig(tm[0], tm[1], tm[2])) /
                 (eps[i] * eps[i]);
    for (int j = i + 1; j < 3; ++j) {
      double tpp[3] = {theta[0], theta[1], theta[2]};
      double tpm[3] = {theta[0], theta[1], theta[2]};
      double tmp[3] = {theta[0], theta[1], theta[2]};
      double tmm[3] = {theta[0], theta[1], theta[2]};
      tpp[i] += eps[i];
      tpp[j] += eps[j];
      tpm[i] += eps[i];
      tpm[j] -= eps[j];
      tmp[i] -= eps[i];
      tmp[j] += eps[j];
      tmm[i] -= eps[i];
      tmm[j] -= eps[j];
      const double hij = (nll_orig(tpp[0], tpp[1], tpp[2]) - nll_orig(tpm[0], tpm[1], tpm[2]) -
                          nll_orig(tmp[0], tmp[1], tmp[2]) + nll_orig(tmm[0], tmm[1], tmm[2])) /
                         (4.0 * eps[i] * eps[j]);
      hess[i][j] = hij;
      hess[j][i] = hij;
    }
  }

  std::array<std::array<double, 3>, 3> cov{};
  double se_mu = 0.0, se_alpha = 0.0, se_beta = 0.0, corr_ab = 0.0;
  if (invert_3x3(hess, cov)) {
    se_mu = std::sqrt(std::max(cov[0][0], 0.0));
    se_alpha = std::sqrt(std::max(cov[1][1], 0.0));
    se_beta = std::sqrt(std::max(cov[2][2], 0.0));
    if (se_alpha > 0.0 && se_beta > 0.0) {
      corr_ab = cov[1][2] / (se_alpha * se_beta);
    }
  }
  return {p_hat, -f, se_mu, se_alpha, se_beta, corr_ab};
}

HawkesMleMcSummary run_hawkes_mle_mc(const HawkesParams& true_p, std::size_t r_reps, double t_end,
                                     std::uint64_t master_seed, std::size_t n_threads) {
  const unsigned hw = std::thread::hardware_concurrency();
  const std::size_t eff_t =
      std::min((n_threads > 0) ? n_threads : (hw > 0 ? static_cast<std::size_t>(hw) : 1), r_reps);
  std::vector<HawkesMleFit> fits(r_reps);

  if (eff_t <= 1) {
    for (std::size_t r = 0; r < r_reps; ++r) {
      const auto ts = simulate_hawkes_dassios_zhao(true_p, t_end, derive_seed(master_seed, r));
      fits[r] = fit_hawkes_mle(ts, t_end);
    }
  } else {
    std::vector<std::thread> workers;
    workers.reserve(eff_t);
    for (std::size_t wid = 0; wid < eff_t; ++wid) {
      workers.emplace_back([&, wid]() {
        for (std::size_t r = wid; r < r_reps; r += eff_t) {
          const auto ts = simulate_hawkes_dassios_zhao(true_p, t_end, derive_seed(master_seed, r));
          fits[r] = fit_hawkes_mle(ts, t_end);
        }
      });
    }
    for (auto& w : workers) {
      w.join();
    }
  }

  std::vector<double> mus, alphas, betas, se_m, se_a, se_b;
  std::vector<HawkesParams> ests;
  std::size_t cov_m = 0, cov_a = 0, cov_b = 0;
  for (const auto& f : fits) {
    mus.push_back(f.params.mu);
    alphas.push_back(f.params.alpha);
    betas.push_back(f.params.beta);
    se_m.push_back(f.se_mu);
    se_a.push_back(f.se_alpha);
    se_b.push_back(f.se_beta);
    ests.push_back(f.params);
    if (std::abs(f.params.mu - true_p.mu) <= 1.96 * f.se_mu) {
      ++cov_m;
    }
    if (std::abs(f.params.alpha - true_p.alpha) <= 1.96 * f.se_alpha) {
      ++cov_a;
    }
    if (std::abs(f.params.beta - true_p.beta) <= 1.96 * f.se_beta) {
      ++cov_b;
    }
  }

  const double ma = stats::mean(alphas);
  const double mb = stats::mean(betas);
  const double sda = std::sqrt(stats::variance(alphas));
  const double sdb = std::sqrt(stats::variance(betas));
  double cov_ab = 0.0;
  for (std::size_t i = 0; i < r_reps; ++i) {
    cov_ab += (alphas[i] - ma) * (betas[i] - mb);
  }
  cov_ab /= static_cast<double>(r_reps - 1);

  HawkesMleMcSummary out;
  out.r_reps = r_reps;
  out.mean_mu = stats::mean(mus);
  out.sd_mu = std::sqrt(stats::variance(mus));
  out.bias_mu = out.mean_mu - true_p.mu;
  out.mean_se_mu = stats::mean(se_m);
  out.cov_mu = static_cast<double>(cov_m) / static_cast<double>(r_reps);
  out.mean_alpha = ma;
  out.sd_alpha = sda;
  out.bias_alpha = ma - true_p.alpha;
  out.mean_se_alpha = stats::mean(se_a);
  out.cov_alpha = static_cast<double>(cov_a) / static_cast<double>(r_reps);
  out.mean_beta = mb;
  out.sd_beta = sdb;
  out.bias_beta = mb - true_p.beta;
  out.mean_se_beta = stats::mean(se_b);
  out.cov_beta = static_cast<double>(cov_b) / static_cast<double>(r_reps);
  out.corr_alpha_beta = cov_ab / (sda * sdb);
  out.estimates = std::move(ests);
  return out;
}

} // namespace hob
