# frozen_string_literal: true

module FactorPortfolioLab
  # Exact Euler risk attribution and factor/idiosyncratic decomposition (including misspecified K).
  module Attribution
    module_function

    def euler_decomposition(weights, sigma)
      sw = LinAlg.mat_vec(sigma, weights)
      port_var = LinAlg.dot(weights, sw)
      port_vol = Math.sqrt([port_var, 0.0].max)
      rc = weights.each_with_index.map { |w, i| (w * sw[i]) / [port_vol, 1.0e-18].max }
      sum_rc = rc.sum
      {
        portfolio_vol: port_vol,
        risk_contributions: rc,
        sum_rc: sum_rc,
        euler_error: (sum_rc - port_vol).abs,
        rel_rc: rc.map { |r| r / [port_vol, 1.0e-18].max }
      }
    end

    def factor_risk_decomposition(weights, betas, sigma_f, idio_var)
      port_exposure = LinAlg.mat_vec(LinAlg.transpose(betas), weights)
      sys_var = LinAlg.quad_form(port_exposure, sigma_f)
      idio = weights.each_with_index.sum { |w, i| w * w * idio_var[i] }
      tot_var = sys_var + idio
      {
        factor_exposures: port_exposure,
        systematic_var: sys_var,
        idiosyncratic_var: idio,
        total_var: tot_var,
        systematic_share: tot_var <= 1e-18 ? 0.0 : sys_var / tot_var,
        idiosyncratic_share: tot_var <= 1e-18 ? 0.0 : idio / tot_var
      }
    end

    def misspecified_k_study(master_seed:, n_assets: 20, t_len: 360)
      world = FactorWorld.build_world(n_assets: n_assets, n_factors: 3, seed: Rng.derive_seed(master_seed, 1000, 0))
      sim = FactorWorld.simulate(world, t_len: t_len, seed: Rng.derive_seed(master_seed, 1001, 0))
      w = RiskParity.solve_newton(world[:sigma])[:weights]

      euler_true = euler_decomposition(w, world[:sigma])
      exact_fac = factor_risk_decomposition(w, world[:betas], world[:sigma_f], world[:idio_var])

      k_models = [1, 2, 3].map do |k_fit|
        f_sub = sim[:factors].map { |row| row[0, k_fit] }
        ts = TimeSeriesRegression.fit(sim[:returns], f_sub)
        f_bar = Pca.column_means(f_sub, k_fit)
        sf_sub = Grs.sample_cov_mle(f_sub, f_bar, k_fit)
        idio_est = ts[:sigma_eps_unbiased].each_with_index.map { |row, i| row[i] }
        dec = factor_risk_decomposition(w, ts[:betas], sf_sub, idio_est)
        dec.merge(k_fit: k_fit)
      end

      {
        euler_error: euler_true[:euler_error],
        true_systematic_share: exact_fac[:systematic_share],
        k_models: k_models
      }
    end
  end
end
