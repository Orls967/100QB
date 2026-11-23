# frozen_string_literal: true

module FactorPortfolioLab
  # Two-pass Fama-MacBeth (1973) cross-sectional asset pricing estimator.
  module FamaMacbeth
    module_function

    def fit(returns, factors, hac_lags: 0, add_intercept: true)
      ts = TimeSeriesRegression.fit(returns, factors)
      betas = ts[:betas]
      t_len = returns.length
      n_assets = returns[0].length
      k_factors = factors[0].length

      x_cs = betas.map { |b| add_intercept ? [1.0] + b : b.dup }
      p_dim = x_cs[0].length
      xtx_inv = LinAlg.invert(Ols.gram_matrix(x_cs, p_dim))
      proj = LinAlg.mat_mul(xtx_inv, LinAlg.transpose(x_cs))

      gammas = Array.new(t_len) { |t| LinAlg.mat_vec(proj, returns[t]) }
      lambda_hat = Pca.column_means(gammas, p_dim)
      vcov_fm = fm_vcov(gammas, lambda_hat, t_len, p_dim)
      vcov_hac = hac_lags.positive? ? gammas_hac_vcov(gammas, lambda_hat, t_len, p_dim, hac_lags) : vcov_fm

      se_fm = Ols.diag_sqrt(vcov_fm)
      se_hac = Ols.diag_sqrt(vcov_hac)
      t_fm = compute_t_stats(lambda_hat, se_fm)
      t_hac = compute_t_stats(lambda_hat, se_hac)

      mean_r = Pca.column_means(returns, n_assets)
      cs_r2 = cross_sectional_r2(mean_r, x_cs, lambda_hat)

      {
        lambda_hat: lambda_hat,
        gammas: gammas,
        betas: betas,
        vcov_fm: vcov_fm,
        vcov_hac: vcov_hac,
        se_fm: se_fm,
        se_hac: se_hac,
        t_fm: t_fm,
        t_hac: t_hac,
        cs_r2: cs_r2,
        n_obs: t_len,
        n_assets: n_assets,
        n_factors: k_factors,
        add_intercept: add_intercept
      }
    end

    def fm_vcov(gammas, lambda_hat, t_len, p_dim)
      cov = Array.new(p_dim) { Array.new(p_dim, 0.0) }
      gammas.each do |g|
        p_dim.times do |i|
          di = g[i] - lambda_hat[i]
          j = i
          while j < p_dim
            cov[i][j] += di * (g[j] - lambda_hat[j])
            j += 1
          end
        end
      end
      scale = 1.0 / (t_len * [t_len - 1, 1].max)
      p_dim.times do |i|
        (i...p_dim).each do |j|
          v = cov[i][j] * scale
          cov[i][j] = v
          cov[j][i] = v
        end
      end
      cov
    end

    def gammas_hac_vcov(gammas, lambda_hat, t_len, p_dim, lags)
      centered = gammas.map { |g| Array.new(p_dim) { |i| g[i] - lambda_hat[i] } }
      vcov = Array.new(p_dim) { Array.new(p_dim, 0.0) }
      t_len.times { |t| Ols.add_outer!(vcov, centered[t], centered[t], 1.0, p_dim) }
      (1..lags).each do |lag|
        w = 1.0 - (lag.to_f / (lags + 1))
        (lag...t_len).each { |t| Ols.add_sym_outer!(vcov, centered[t], centered[t - lag], w, p_dim) }
      end
      scale = 1.0 / (t_len * t_len)
      vcov.map { |row| row.map { |v| v * scale } }
    end

    def compute_t_stats(est, ses)
      est.each_with_index.map { |val, i| ses[i] <= 1e-15 ? 0.0 : val / ses[i] }
    end

    def cross_sectional_r2(mean_r, x_cs, lambda_hat)
      preds = x_cs.map { |row| LinAlg.dot(row, lambda_hat) }
      avg_r = mean_r.sum / mean_r.length
      tss = mean_r.sum { |r| (r - avg_r)**2 }
      rss = mean_r.each_with_index.sum { |r, i| (r - preds[i])**2 }
      tss <= 1e-18 ? 0.0 : 1.0 - (rss / tss)
    end
  end
end
