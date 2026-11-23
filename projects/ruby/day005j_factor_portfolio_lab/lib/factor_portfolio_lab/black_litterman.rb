# frozen_string_literal: true

module FactorPortfolioLab
  # Black-Litterman Bayesian asset allocation (He-Litterman & Idzorek confidence calibration).
  module BlackLitterman
    module_function

    def implied_equilibrium(sigma, w_mkt, risk_aversion: 2.5)
      sw = LinAlg.mat_vec(sigma, w_mkt)
      sw.map { |v| risk_aversion * v }
    end

    def he_litterman_omega(sigma, p_mat, tau: 0.05)
      k_views = p_mat.length
      tau_sigma = scale_mat(sigma, tau)
      p_ts_pt = LinAlg.mat_mul(LinAlg.mat_mul(p_mat, tau_sigma), LinAlg.transpose(p_mat))
      Array.new(k_views) do |i|
        Array.new(k_views) { |j| i == j ? [p_ts_pt[i][i], 1.0e-12].max : 0.0 }
      end
    end

    def idzorek_omega(sigma, p_mat, confidences, tau: 0.05)
      k_views = p_mat.length
      tau_sigma = scale_mat(sigma, tau)
      Array.new(k_views) do |i|
        c = confidences[i].clamp(1.0e-4, 1.0 - 1.0e-4)
        p_row = p_mat[i]
        v_i = LinAlg.quad_form(p_row, tau_sigma)
        omega_ii = ((1.0 - c) / c) * [v_i, 1.0e-12].max
        Array.new(k_views) { |j| i == j ? omega_ii : 0.0 }
      end
    end

    def posterior(sigma:, pi:, p_mat:, q_vec:, omega:, tau: 0.05)
      tau_sigma = scale_mat(sigma, tau)
      pt = LinAlg.transpose(p_mat)
      ts_pt = LinAlg.mat_mul(tau_sigma, pt)
      inner = LinAlg.mat_mul(p_mat, ts_pt)
      inner.each_with_index { |row, i| row.each_index { |j| row[j] += omega[i][j] } }
      inner_inv = LinAlg.invert(inner)

      p_pi = LinAlg.mat_vec(p_mat, pi)
      diff = q_vec.each_with_index.map { |q, i| q - p_pi[i] }
      gain = LinAlg.mat_mul(ts_pt, inner_inv)
      adj = LinAlg.mat_vec(gain, diff)
      mu_post = pi.each_with_index.map { |p_i, i| p_i + adj[i] }

      cov_adj = LinAlg.mat_mul(gain, LinAlg.mat_mul(p_mat, tau_sigma))
      n = sigma.length
      sigma_mu_post = Array.new(n) do |i|
        Array.new(n) { |j| tau_sigma[i][j] - cov_adj[i][j] }
      end
      sigma_bl = Array.new(n) do |i|
        Array.new(n) { |j| sigma[i][j] + sigma_mu_post[i][j] }
      end

      {
        mu_post: mu_post,
        sigma_mu_post: sigma_mu_post,
        sigma_bl: sigma_bl,
        pi: pi
      }
    end

    def scale_mat(mat, scalar)
      mat.map { |row| row.map { |v| v * scalar } }
    end
  end
end
