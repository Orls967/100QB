# frozen_string_literal: true

module FactorPortfolioLab
  # Comparison of Risk Parity solvers (Newton Spinu vs CCD) and concentration vs Markowitz/1/N.
  module RpStudy
    module_function

    def run(master_seed:, sizes: [10, 25, 50])
      comparisons = sizes.each_with_index.map do |n_assets, idx|
        world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 500, idx))
        sigma = world[:sigma]

        t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
        newton_res = RiskParity.solve_newton(sigma)
        t_newton = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

        t1 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
        ccd_res = RiskParity.solve_ccd(sigma)
        t_ccd = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t1

        max_w_diff = newton_res[:weights].each_with_index.map { |w, i| (w - ccd_res[:weights][i]).abs }.max
        eq_w = Array.new(n_assets, 1.0 / n_assets)
        eq_rc = RiskParity.summarize_solution(sigma, eq_w, 0, true)

        {
          n_assets: n_assets,
          newton_iters: newton_res[:iterations],
          newton_time_s: t_newton,
          newton_rc_max_dev: newton_res[:rc_max_dev],
          ccd_iters: ccd_res[:iterations],
          ccd_time_s: t_ccd,
          ccd_rc_max_dev: ccd_res[:rc_max_dev],
          max_weight_diff: max_w_diff,
          rp_vol: newton_res[:portfolio_vol],
          eq_vol: eq_rc[:portfolio_vol],
          rp_hhi: newton_res[:weights].sum { |w| w * w },
          eq_rc_max_dev: eq_rc[:rc_max_dev]
        }
      end

      { comparisons: comparisons }
    end
  end
end
