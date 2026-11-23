# frozen_string_literal: true

module FactorPortfolioLab
  # Efficient frontier study comparing unconstrained, long-only, and weight-capped Markowitz QP.
  module MvStudy
    module_function

    def run(master_seed:, n_assets: 20, n_points: 45)
      world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 700, 0))
      sigma = world[:sigma]
      mu = world[:mu]

      kappas = Array.new(n_points) { |i| (i.to_f / (n_points - 1)) * 1.5 }
      long_only = build_frontier(sigma, mu, kappas, Array.new(n_assets, 0.0), Array.new(n_assets, 1.0))
      capped = build_frontier(sigma, mu, kappas, Array.new(n_assets, 0.0), Array.new(n_assets, 0.12))

      gmv_unc = Markowitz.unconstrained_gmv(sigma)
      tan_unc = Markowitz.unconstrained_tangency(sigma, mu)
      unc_targets = Array.new(n_points) do |i|
        mu_g = LinAlg.dot(gmv_unc, mu)
        mu_t = LinAlg.dot(tan_unc, mu)
        mu_g + ((i.to_f / (n_points - 1)) * [(mu_t - mu_g).abs * 1.3, 0.006].max)
      end
      unconstrained = unc_targets.map do |tgt|
        w = Markowitz.unconstrained_frontier_point(sigma, mu, tgt)
        { expected_return: LinAlg.dot(w, mu), volatility: Math.sqrt([LinAlg.quad_form(w, sigma), 0.0].max), weights: w }
      end

      {
        n_assets: n_assets,
        n_points: n_points,
        unconstrained: unconstrained,
        long_only: long_only,
        capped: capped,
        max_kkt_long_only: long_only.map { |p| p[:max_kkt] }.max,
        max_kkt_capped: capped.map { |p| p[:max_kkt] }.max
      }
    end

    def build_frontier(sigma, mu, kappas, lower, upper)
      kappas.map do |kappa|
        sol = Markowitz.solve_qp(sigma, mu: mu, kappa: kappa, lower: lower, upper: upper)
        {
          kappa: kappa,
          expected_return: sol[:expected_return],
          volatility: sol[:volatility],
          max_kkt: sol[:kkt][:max_kkt],
          n_active: sol[:weights].count { |w| w > 1.0e-4 },
          max_weight: sol[:weights].max,
          weights: sol[:weights]
        }
      end
    end
  end
end
