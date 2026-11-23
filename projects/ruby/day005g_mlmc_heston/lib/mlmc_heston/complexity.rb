# frozen_string_literal: true

module MlmcHeston
  module Complexity
    module_function

    def giles_asymptotic_slope(beta_var, alpha_bias, gamma: 1.0)
      return 2.0 if beta_var > gamma

      eff_alpha = [alpha_bias, 0.5].max
      2.0 + ((gamma - beta_var) / eff_alpha)
    end

    def compare_mlmc_vs_standard(model:, params:, eps_grid: [0.08, 0.04, 0.02, 0.01],
                                 seed: MAIN_SEED, workers: nil, n0: 4, rate_est: nil)
      probe = Mlmc.level_stats(model: model, params: params, level: 0, n_paths: 8_000,
                               seed: seed, workers: workers, n0: n0)
      v0 = probe[:var_y]

      rows = eps_grid.map do |eps|
        res = Mlmc.run_adaptive(model: model, params: params, eps: eps,
                                seed: seed, workers: workers, n0: n0)
        steps_std = (n0 * 2.0 * (0.08 / eps)).ceil
        n_std = (2.0 * v0 / (eps * eps)).ceil
        std_cost = steps_std.to_f * n_std
        res.merge(std_mc_cost: std_cost, speedup: std_cost / res[:total_cost],
                  giles_formula_cost: giles_theorem_cost(res[:vars], n0, eps))
      end

      summarize_complexity(rows, rate_est)
    end

    def giles_theorem_cost(vars, n0, eps)
      sum_vc = vars.each_with_index.sum(0.0) do |v, l|
        c_l = (l.zero? ? 1.0 : 1.5) * n0 * (2**l)
        Math.sqrt(v * c_l)
      end
      2.0 * sum_vc * sum_vc / (eps * eps)
    end

    def summarize_complexity(rows, rate_est)
      log_inv_eps = rows.map { |r| Math.log(1.0 / r[:eps]) }
      mlmc_fit = ols_fit(log_inv_eps, rows.map { |r| Math.log(r[:total_cost]) })
      std_fit = ols_fit(log_inv_eps, rows.map { |r| Math.log(r[:std_mc_cost]) })
      giles_pre_slope = Mlmc.ols_slope(log_inv_eps, rows.map { |r| Math.log(r[:giles_formula_cost]) })

      finest_vars = rows.last[:vars]
      beta_from_vars = estimate_beta_from_vars(finest_vars)
      beta_used = rate_est ? rate_est[:beta_strong] : beta_from_vars
      alpha_used = rate_est ? rate_est[:alpha_weak] : 1.0

      {
        rows: rows, mlmc_cost_slope: mlmc_fit[:slope], std_mc_cost_slope: std_fit[:slope],
        beta_var_levels: beta_from_vars, beta_rate_part2: beta_used, alpha_bias_part2: alpha_used,
        giles_asymptotic_slope: giles_asymptotic_slope(beta_used, alpha_used),
        giles_preasymptotic_slope: giles_pre_slope,
        v1_to_v0_ratio: finest_vars[1] / [finest_vars[0], 1.0e-18].max,
        crossover_eps_extrapolated: extrapolated_crossover_eps(mlmc_fit, std_fit)
      }
    end
    private_class_method :summarize_complexity

    def estimate_beta_from_vars(vars)
      return 1.0 if vars.length < 3

      ls = (1...(vars.length)).map(&:to_f)
      log2_vs = (1...(vars.length)).map { |l| Math.log2([vars[l], 1.0e-18].max) }
      -Mlmc.ols_slope(ls, log2_vs)
    end
    private_class_method :estimate_beta_from_vars

    def ols_fit(xs, ys)
      slope = Mlmc.ols_slope(xs, ys)
      mx = xs.sum(0.0) / xs.length
      my = ys.sum(0.0) / ys.length
      { slope: slope, intercept: my - (slope * mx) }
    end
    private_class_method :ols_fit

    def extrapolated_crossover_eps(mlmc_fit, std_fit)
      dp = std_fit[:slope] - mlmc_fit[:slope]
      return nil if dp <= 0.05

      log_inv_eps = (mlmc_fit[:intercept] - std_fit[:intercept]) / dp
      Math.exp(-log_inv_eps)
    end
    private_class_method :extrapolated_crossover_eps
  end
end
