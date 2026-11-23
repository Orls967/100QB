# frozen_string_literal: true

module FactorPortfolioLab
  # Gibbons-Ross-Shanken (1989) exact finite-sample F test and asymptotic Wald chi-squared test.
  module Grs
    module_function

    def test(returns, factors)
      fit = TimeSeriesRegression.fit(returns, factors)
      t_len = fit[:n_obs]
      n_assets = fit[:n_assets]
      k_factors = fit[:n_factors]
      df2 = t_len - n_assets - k_factors
      raise ArgumentError, "insufficient degrees of freedom: T=#{t_len}, N=#{n_assets}, K=#{k_factors}" if df2 <= 0

      f_bar = Pca.column_means(factors, k_factors)
      sigma_f_mle = sample_cov_mle(factors, f_bar, k_factors)
      sharpe_f_sq = LinAlg.quad_form(f_bar, LinAlg.invert(sigma_f_mle))

      alphas = fit[:alphas]
      alpha_quad = LinAlg.quad_form(alphas, LinAlg.invert(fit[:sigma_eps_mle]))
      denom = 1.0 + sharpe_f_sq

      grs_f = (df2.to_f / n_assets) * (alpha_quad / denom)
      wald_chi2 = t_len.to_f * (alpha_quad / denom)

      {
        grs_f: grs_f,
        df1: n_assets,
        df2: df2,
        p_value_f: Distributions.f_sf(grs_f, n_assets, df2),
        wald_chi2: wald_chi2,
        p_value_wald: Distributions.chi2_sf(wald_chi2, n_assets),
        sharpe_f_sq: sharpe_f_sq,
        mean_abs_alpha: alphas.sum(&:abs) / n_assets,
        alphas: alphas
      }
    end

    def sample_cov_mle(rows, means, dim)
      t_len = rows.length
      cov = Array.new(dim) { Array.new(dim, 0.0) }
      rows.each do |r|
        dim.times do |i|
          di = r[i] - means[i]
          j = i
          while j < dim
            cov[i][j] += di * (r[j] - means[j])
            j += 1
          end
        end
      end
      inv_t = 1.0 / t_len
      dim.times do |i|
        (i...dim).each do |j|
          v = cov[i][j] * inv_t
          cov[i][j] = v
          cov[j][i] = v
        end
      end
      cov
    end
  end
end
