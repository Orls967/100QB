#include "hw/trinomial_tree.hpp"

#include <algorithm>
#include <cmath>

namespace hw {

TrinomialTree::TrinomialTree(const NelsonSiegelCurve& curve, const HwParams& params,
                             int steps_per_year, int max_years)
    : params_(params),
      m_(steps_per_year),
      n_steps_(steps_per_year * max_years),
      dt_(1.0 / static_cast<double>(steps_per_year)),
      dx_(params.sigma * std::sqrt(3.0 * dt_)) {
  const double raw = 0.184 / (params_.a * dt_);
  j_max_ = static_cast<int>(std::floor(raw)) + 1;
  alpha_.resize(static_cast<std::size_t>(n_steps_ + 1), 0.0);
  q_.assign(static_cast<std::size_t>(n_steps_ + 1),
            std::vector<double>(static_cast<std::size_t>(2 * j_max_ + 1), 0.0));
  build(curve);
}

BranchProbs TrinomialTree::branching(int j) const noexcept {
  const double mu = -params_.a * static_cast<double>(j) * dt_;
  if (j <= -j_max_) {
    return {j + 1, 7.0 / 6.0 + 0.5 * (mu * mu - 3.0 * mu), -1.0 / 3.0 - mu * mu + 2.0 * mu,
            1.0 / 6.0 + 0.5 * (mu * mu - mu)};
  }
  if (j >= j_max_) {
    return {j - 1, 1.0 / 6.0 + 0.5 * (mu * mu + mu), -1.0 / 3.0 - mu * mu - 2.0 * mu,
            7.0 / 6.0 + 0.5 * (mu * mu + 3.0 * mu)};
  }
  return {j, 1.0 / 6.0 + 0.5 * (mu * mu - mu), 2.0 / 3.0 - mu * mu,
          1.0 / 6.0 + 0.5 * (mu * mu + mu)};
}

void TrinomialTree::build(const NelsonSiegelCurve& curve) {
  q_[0][idx(0)] = 1.0;
  max_curve_err_ = 0.0;

  for (int i = 0; i < n_steps_; ++i) {
    const int j_lo = std::max(-j_max_, -i);
    const int j_hi = std::min(j_max_, i);
    double sum_qe = 0.0;
    for (int j = j_lo; j <= j_hi; ++j) {
      sum_qe += q_[i][idx(j)] * std::exp(-static_cast<double>(j) * dx_ * dt_);
    }
    const double p_next = curve.discount(static_cast<double>(i + 1) * dt_);
    alpha_[i] = (std::log(sum_qe) - std::log(p_next)) / dt_;

    for (int j = j_lo; j <= j_hi; ++j) {
      const double q_disc =
          q_[i][idx(j)] * std::exp(-(alpha_[i] + static_cast<double>(j) * dx_) * dt_);
      const BranchProbs br = branching(j);
      q_[i + 1][idx(br.k_mid - 1)] += q_disc * br.pd;
      q_[i + 1][idx(br.k_mid)] += q_disc * br.pm;
      q_[i + 1][idx(br.k_mid + 1)] += q_disc * br.pu;
    }

    double sum_q_next = 0.0;
    for (double qv : q_[i + 1]) {
      sum_q_next += qv;
    }
    max_curve_err_ = std::max(max_curve_err_, std::abs(sum_q_next - p_next));
  }
  alpha_[n_steps_] = alpha_[n_steps_ - 1];
}

double TrinomialTree::rate(int step, int j) const noexcept {
  return alpha_[step] + static_cast<double>(j) * dx_;
}

std::size_t TrinomialTree::node_count() const noexcept {
  std::size_t cnt = 0;
  for (int i = 0; i <= n_steps_; ++i) {
    const int span = std::min(i, j_max_);
    cnt += static_cast<std::size_t>(2 * span + 1);
  }
  return cnt;
}

std::vector<double> TrinomialTree::rollback_zero_bond(int mat_step, int stop_step) const {
  const std::size_t width = static_cast<std::size_t>(2 * j_max_ + 1);
  std::vector<double> v(width, 1.0);
  std::vector<double> v_prev(width, 0.0);
  for (int i = mat_step - 1; i >= stop_step; --i) {
    const int j_lo = std::max(-j_max_, -i);
    const int j_hi = std::min(j_max_, i);
    for (int j = j_lo; j <= j_hi; ++j) {
      const BranchProbs br = branching(j);
      const double cont =
          br.pd * v[idx(br.k_mid - 1)] + br.pm * v[idx(br.k_mid)] + br.pu * v[idx(br.k_mid + 1)];
      v_prev[idx(j)] = cont * std::exp(-rate(i, j) * dt_);
    }
    v.swap(v_prev);
  }
  return v;
}

double TrinomialTree::price_zbp(double t_opt, double t_bond, double strike) const {
  const int opt_step = static_cast<int>(std::llround(t_opt * m_));
  const int mat_step = static_cast<int>(std::llround(t_bond * m_));
  std::vector<double> bond = rollback_zero_bond(mat_step, opt_step);
  const std::size_t width = static_cast<std::size_t>(2 * j_max_ + 1);
  std::vector<double> v(width, 0.0);
  for (int j = std::max(-j_max_, -opt_step); j <= std::min(j_max_, opt_step); ++j) {
    v[idx(j)] = std::max(strike - bond[idx(j)], 0.0);
  }
  std::vector<double> v_prev(width, 0.0);
  for (int i = opt_step - 1; i >= 0; --i) {
    for (int j = std::max(-j_max_, -i); j <= std::min(j_max_, i); ++j) {
      const BranchProbs br = branching(j);
      const double cont =
          br.pd * v[idx(br.k_mid - 1)] + br.pm * v[idx(br.k_mid)] + br.pu * v[idx(br.k_mid + 1)];
      v_prev[idx(j)] = cont * std::exp(-rate(i, j) * dt_);
    }
    v.swap(v_prev);
  }
  return v[idx(0)];
}

double TrinomialTree::price_zbc(double t_opt, double t_bond, double strike) const {
  const int opt_step = static_cast<int>(std::llround(t_opt * m_));
  const int mat_step = static_cast<int>(std::llround(t_bond * m_));
  std::vector<double> bond = rollback_zero_bond(mat_step, opt_step);
  const std::size_t width = static_cast<std::size_t>(2 * j_max_ + 1);
  std::vector<double> v(width, 0.0);
  for (int j = std::max(-j_max_, -opt_step); j <= std::min(j_max_, opt_step); ++j) {
    v[idx(j)] = std::max(bond[idx(j)] - strike, 0.0);
  }
  std::vector<double> v_prev(width, 0.0);
  for (int i = opt_step - 1; i >= 0; --i) {
    for (int j = std::max(-j_max_, -i); j <= std::min(j_max_, i); ++j) {
      const BranchProbs br = branching(j);
      const double cont =
          br.pd * v[idx(br.k_mid - 1)] + br.pm * v[idx(br.k_mid)] + br.pu * v[idx(br.k_mid + 1)];
      v_prev[idx(j)] = cont * std::exp(-rate(i, j) * dt_);
    }
    v.swap(v_prev);
  }
  return v[idx(0)];
}

double TrinomialTree::price_cap_5y(double strike_rate) const {
  const double x = 1.0 / (1.0 + strike_rate);
  double total = 0.0;
  for (int yr = 1; yr <= 5; ++yr) {
    total += (1.0 + strike_rate) * price_zbp(static_cast<double>(yr - 1), static_cast<double>(yr), x);
  }
  return total;
}

double TrinomialTree::price_payer_swaption(double t0, int tenor_end_yr, double fixed_rate) const {
  const int opt_step = static_cast<int>(std::llround(t0 * m_));
  const int end_step = tenor_end_yr * m_;
  const std::size_t width = static_cast<std::size_t>(2 * j_max_ + 1);

  std::vector<double> c_bond(width, 1.0 + fixed_rate);
  std::vector<double> v_prev(width, 0.0);
  for (int i = end_step - 1; i >= opt_step; --i) {
    const int j_lo = std::max(-j_max_, -i);
    const int j_hi = std::min(j_max_, i);
    for (int j = j_lo; j <= j_hi; ++j) {
      const BranchProbs br = branching(j);
      const double cont = br.pd * c_bond[idx(br.k_mid - 1)] + br.pm * c_bond[idx(br.k_mid)] +
                          br.pu * c_bond[idx(br.k_mid + 1)];
      v_prev[idx(j)] = cont * std::exp(-rate(i, j) * dt_);
    }
    c_bond.swap(v_prev);
    if (i > opt_step && (i % m_ == 0)) {
      for (int j = j_lo; j <= j_hi; ++j) {
        c_bond[idx(j)] += fixed_rate;
      }
    }
  }

  std::vector<double> v(width, 0.0);
  for (int j = std::max(-j_max_, -opt_step); j <= std::min(j_max_, opt_step); ++j) {
    v[idx(j)] = std::max(1.0 - c_bond[idx(j)], 0.0);
  }
  for (int i = opt_step - 1; i >= 0; --i) {
    for (int j = std::max(-j_max_, -i); j <= std::min(j_max_, i); ++j) {
      const BranchProbs br = branching(j);
      const double cont =
          br.pd * v[idx(br.k_mid - 1)] + br.pm * v[idx(br.k_mid)] + br.pu * v[idx(br.k_mid + 1)];
      v_prev[idx(j)] = cont * std::exp(-rate(i, j) * dt_);
    }
    v.swap(v_prev);
  }
  return v[idx(0)];
}

double TrinomialTree::price_bermudan_payer(int tenor_end_yr, double fixed_rate) const {
  const int end_step = tenor_end_yr * m_;
  const std::size_t width = static_cast<std::size_t>(2 * j_max_ + 1);

  std::vector<double> c_bond(width, 1.0 + fixed_rate);
  std::vector<double> v_opt(width, 0.0);
  std::vector<double> c_prev(width, 0.0);
  std::vector<double> v_prev(width, 0.0);

  for (int i = end_step - 1; i >= 0; --i) {
    const int j_lo = std::max(-j_max_, -i);
    const int j_hi = std::min(j_max_, i);
    for (int j = j_lo; j <= j_hi; ++j) {
      const BranchProbs br = branching(j);
      const double disc = std::exp(-rate(i, j) * dt_);
      const double cb_cont = br.pd * c_bond[idx(br.k_mid - 1)] + br.pm * c_bond[idx(br.k_mid)] +
                             br.pu * c_bond[idx(br.k_mid + 1)];
      c_prev[idx(j)] = cb_cont * disc;

      const double opt_cont = br.pd * v_opt[idx(br.k_mid - 1)] + br.pm * v_opt[idx(br.k_mid)] +
                              br.pu * v_opt[idx(br.k_mid + 1)];
      v_prev[idx(j)] = opt_cont * disc;
    }
    c_bond.swap(c_prev);
    v_opt.swap(v_prev);

    if (i > 0 && (i % m_ == 0)) {
      for (int j = j_lo; j <= j_hi; ++j) {
        const double ex_val = std::max(1.0 - c_bond[idx(j)], 0.0);
        v_opt[idx(j)] = std::max(v_opt[idx(j)], ex_val);
        c_bond[idx(j)] += fixed_rate;
      }
    }
  }
  return v_opt[idx(0)];
}

} // namespace hw
