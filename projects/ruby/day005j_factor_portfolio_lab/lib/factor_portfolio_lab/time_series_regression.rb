# frozen_string_literal: true

module FactorPortfolioLab
  # Simultaneous multi-asset time-series factor regression R_{t,i} = alpha_i + beta_i' f_t + eps_{t,i}.
  module TimeSeriesRegression
    module_function

    def fit(returns, factors)
      t_len = returns.length
      n_assets = returns[0].length
      k_factors = factors[0].length
      x_mat = factors.map { |f| [1.0] + f }
      p_dim = k_factors + 1

      xtx = Ols.gram_matrix(x_mat, p_dim)
      xtx_inv = LinAlg.invert(xtx)
      coefs = solve_multi_rhs(x_mat, returns, xtx_inv, p_dim, n_assets)
      alphas = coefs[0]
      betas = Array.new(n_assets) { |i| Array.new(k_factors) { |k| coefs[k + 1][i] } }

      residuals = compute_multi_residuals(x_mat, returns, coefs, t_len, n_assets, p_dim)
      dof = [t_len - p_dim, 1].max
      sigma_eps_mle = residual_cov(residuals, t_len, n_assets, t_len.to_f)
      sigma_eps_unbiased = residual_cov(residuals, t_len, n_assets, dof.to_f)
      se_alpha, t_alpha, se_alpha_hc1, t_alpha_hc1 = alpha_inference(
        x_mat, residuals, xtx_inv, sigma_eps_unbiased, alphas, t_len, n_assets, p_dim
      )

      {
        alphas: alphas,
        betas: betas,
        residuals: residuals,
        sigma_eps_mle: sigma_eps_mle,
        sigma_eps_unbiased: sigma_eps_unbiased,
        se_alpha: se_alpha,
        t_alpha: t_alpha,
        se_alpha_hc1: se_alpha_hc1,
        t_alpha_hc1: t_alpha_hc1,
        r2: compute_r2(returns, residuals, t_len, n_assets),
        n_obs: t_len,
        n_assets: n_assets,
        n_factors: k_factors
      }
    end

    def solve_multi_rhs(x_mat, returns, xtx_inv, p_dim, n_assets)
      xty = Array.new(p_dim) { Array.new(n_assets, 0.0) }
      x_mat.each_with_index do |x_t, t|
        r_t = returns[t]
        p_dim.times do |j|
          xj = x_t[j]
          n_assets.times { |i| xty[j][i] += xj * r_t[i] }
        end
      end
      LinAlg.mat_mul(xtx_inv, xty)
    end

    def compute_multi_residuals(x_mat, returns, coefs, t_len, n_assets, p_dim)
      Array.new(t_len) do |t|
        x_t = x_mat[t]
        r_t = returns[t]
        Array.new(n_assets) do |i|
          pred = 0.0
          p_dim.times { |j| pred += x_t[j] * coefs[j][i] }
          r_t[i] - pred
        end
      end
    end

    def residual_cov(residuals, t_len, n_assets, denom)
      cov = Array.new(n_assets) { Array.new(n_assets, 0.0) }
      t_len.times do |t|
        e_t = residuals[t]
        n_assets.times do |i|
          ei = e_t[i]
          j = i
          while j < n_assets
            cov[i][j] += ei * e_t[j]
            j += 1
          end
        end
      end
      inv_d = 1.0 / denom
      n_assets.times do |i|
        (i...n_assets).each do |j|
          v = cov[i][j] * inv_d
          cov[i][j] = v
          cov[j][i] = v
        end
      end
      cov
    end

    def alpha_inference(x_mat, residuals, xtx_inv, sigma_eps_unbiased, alphas, t_len, n_assets, p_dim)
      h00 = xtx_inv[0][0]
      se_ols = Array.new(n_assets) { |i| Math.sqrt([h00 * sigma_eps_unbiased[i][i], 0.0].max) }
      t_ols = Array.new(n_assets) { |i| se_ols[i] <= 1e-15 ? 0.0 : alphas[i] / se_ols[i] }

      b0 = xtx_inv[0]
      weights = x_mat.map { |x_t| LinAlg.dot(b0, x_t) }
      scale = t_len.to_f / [t_len - p_dim, 1].max
      se_hc1 = Array.new(n_assets) do |i|
        var = 0.0
        t_len.times do |t|
          wt = weights[t] * residuals[t][i]
          var += wt * wt
        end
        Math.sqrt([var * scale, 0.0].max)
      end
      t_hc1 = Array.new(n_assets) { |i| se_hc1[i] <= 1e-15 ? 0.0 : alphas[i] / se_hc1[i] }
      [se_ols, t_ols, se_hc1, t_hc1]
    end

    def compute_r2(returns, residuals, t_len, n_assets)
      Array.new(n_assets) do |i|
        mean_r = returns.sum { |r| r[i] } / t_len
        tss = returns.sum { |r| (r[i] - mean_r)**2 }
        rss = residuals.sum { |e| e[i]**2 }
        tss <= 1e-18 ? 0.0 : 1.0 - (rss / tss)
      end
    end
  end
end
