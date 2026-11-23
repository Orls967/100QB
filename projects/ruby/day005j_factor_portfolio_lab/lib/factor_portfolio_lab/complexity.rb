# frozen_string_literal: true

module FactorPortfolioLab
  # Computational complexity scaling benchmark across N in {10, 25, 50, 100, 200}.
  module Complexity
    module_function

    def run(master_seed:, sizes: [10, 25, 50, 100, 200])
      rows = sizes.each_with_index.map do |n_assets, idx|
        world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 1300, idx))
        sim = FactorWorld.simulate(world, t_len: [n_assets * 2, 60].max, seed: Rng.derive_seed(master_seed, 1301, idx))
        sigma = world[:sigma]

        t_inv = measure_ms { LinAlg.invert(sigma) }
        t_rp = measure_ms { RiskParity.solve_newton(sigma) }
        t_qp = measure_ms { Markowitz.solve_qp(sigma, mu: world[:mu], kappa: 0.3, max_iter: 400) }
        t_lw = measure_ms { Shrinkage.shrink_cov_lw(sim[:returns]) }

        { n_assets: n_assets, invert_ms: t_inv, rp_newton_ms: t_rp, qp_ms: t_qp, lw_ms: t_lw }
      end

      {
        rows: rows,
        slope_invert: log_log_slope(rows, :invert_ms),
        slope_rp_newton: log_log_slope(rows, :rp_newton_ms),
        slope_qp: log_log_slope(rows, :qp_ms),
        slope_lw: log_log_slope(rows, :lw_ms)
      }
    end

    def measure_ms
      t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      yield
      (Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0) * 1000.0
    end

    def log_log_slope(rows, key)
      x_mat = rows.map { |r| [1.0, Math.log(r[:n_assets].to_f)] }
      y_vec = rows.map { |r| Math.log([r[key], 1.0e-4].max) }
      Ols.fit(x_mat, y_vec)[:beta][1]
    end
  end
end
