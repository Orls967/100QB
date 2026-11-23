# frozen_string_literal: true

module FactorPortfolioLab
  # Factor stress testing (linear + quadratic gamma convexity) and closed-form Mahalanobis reverse stress test.
  module Stress
    module_function

    def scenario_pnl(weights:, betas:, delta_f:, gamma_f: nil)
      g = LinAlg.mat_vec(LinAlg.transpose(betas), weights)
      linear = LinAlg.dot(g, delta_f)
      quad = gamma_f ? 0.5 * LinAlg.quad_form(delta_f, gamma_f) : 0.0
      {
        linear_pnl: linear,
        quadratic_pnl: quad,
        total_pnl: linear + quad,
        factor_exposure: g
      }
    end

    # Closed-form minimum-Mahalanobis-norm shock df* minimizing df' Sigma_f^-1 df s.t. g' df = -target_loss.
    def reverse_stress_test(weights:, betas:, sigma_f:, target_loss:)
      g = LinAlg.mat_vec(LinAlg.transpose(betas), weights)
      sf_g = LinAlg.mat_vec(sigma_f, g)
      denom = LinAlg.dot(g, sf_g)
      raise ArgumentError, 'zero factor exposure' if denom <= 1.0e-18

      scale = -target_loss.to_f / denom
      df_star = sf_g.map { |v| scale * v }
      sf_inv = LinAlg.invert(sigma_f)
      mahal_sq = LinAlg.quad_form(df_star, sf_inv)

      {
        shock_df: df_star,
        achieved_pnl: LinAlg.dot(g, df_star),
        mahalanobis_sq: mahal_sq,
        mahalanobis_dist: Math.sqrt([mahal_sq, 0.0].max),
        theory_mahalanobis_sq: (target_loss * target_loss) / denom
      }
    end

    def run_study(master_seed:, n_assets: 20)
      world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 1100, 0))
      w_rp = RiskParity.solve_newton(world[:sigma])[:weights]
      w_mv = Markowitz.solve_qp(world[:sigma], mu: world[:mu], kappa: 0.6)[:weights]
      vols = world[:sigma_f].each_with_index.map { |r, i| Math.sqrt(r[i]) }
      gamma_f = [[-2.5, 0.0, 0.0], [0.0, -1.5, 0.0], [0.0, 0.0, -1.0]]

      scenarios = {
        'MarketCrash3Sigma' => [-3.0 * vols[0], -1.5 * vols[1], -1.0 * vols[2]],
        'Stagflation' => [-2.0 * vols[0], 2.5 * vols[1], -1.5 * vols[2]],
        'FactorReversal' => [-1.5 * vols[0], -2.5 * vols[1], 2.0 * vols[2]]
      }

      scen_out = scenarios.transform_values do |df|
        {
          rp: scenario_pnl(weights: w_rp, betas: world[:betas], delta_f: df, gamma_f: gamma_f),
          mv: scenario_pnl(weights: w_mv, betas: world[:betas], delta_f: df, gamma_f: gamma_f)
        }
      end

      rev_rp = reverse_stress_test(weights: w_rp, betas: world[:betas], sigma_f: world[:sigma_f], target_loss: 0.15)
      rev_mv = reverse_stress_test(weights: w_mv, betas: world[:betas], sigma_f: world[:sigma_f], target_loss: 0.15)

      { scenarios: scen_out, reverse_rp: rev_rp, reverse_mv: rev_mv }
    end
  end
end
