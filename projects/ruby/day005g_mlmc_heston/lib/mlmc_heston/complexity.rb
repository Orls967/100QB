# frozen_string_literal: true

module MlmcHeston
  module Complexity
    module_function

    def compare_mlmc_vs_standard(model:, params:, eps_grid: [0.08, 0.04, 0.02, 0.01],
                                 seed: MAIN_SEED, workers: nil, n0: 4)
      probe = Mlmc.level_stats(model: model, params: params, level: 0, n_paths: 8_000,
                               seed: seed, workers: workers, n0: n0)
      v0 = probe[:var_y]

      rows = eps_grid.map do |eps|
        res = Mlmc.run_adaptive(model: model, params: params, eps: eps,
                                seed: seed, workers: workers, n0: n0)
        steps_std = (n0 * 2.0 * (0.08 / eps)).ceil
        n_std = (2.0 * v0 / (eps * eps)).ceil
        std_cost = steps_std.to_f * n_std
        res.merge(std_mc_cost: std_cost, speedup: std_cost / res[:total_cost])
      end

      log_inv_eps = rows.map { |r| Math.log(1.0 / r[:eps]) }
      mlmc_slope = Mlmc.ols_slope(log_inv_eps, rows.map { |r| Math.log(r[:total_cost]) })
      std_slope = Mlmc.ols_slope(log_inv_eps, rows.map { |r| Math.log(r[:std_mc_cost]) })
      { rows: rows, mlmc_cost_slope: mlmc_slope, std_mc_cost_slope: std_slope }
    end
  end
end
