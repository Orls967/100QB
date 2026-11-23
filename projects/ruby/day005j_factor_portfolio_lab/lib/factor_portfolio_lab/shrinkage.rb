# frozen_string_literal: true

module FactorPortfolioLab
  # Covariance and expected-return shrinkage estimators (Ledoit-Wolf 2004, diagonal, James-Stein).
  module Shrinkage
    module_function

    # Ledoit-Wolf (2004) analytical shrinkage toward scaled identity mu * I.
    def shrink_cov_lw(rows)
      t_len = rows.length
      n = rows[0].length
      means = Pca.column_means(rows, n)
      centered = rows.map { |r| Array.new(n) { |j| r[j] - means[j] } }
      s_mat = Grs.sample_cov_mle(centered, Array.new(n, 0.0), n)

      mu = LinAlg.trace(s_mat) / n
      d2 = frobenius_sq_diff_diag(s_mat, mu, n)
      b2_raw = estimate_b2_bar(centered, s_mat, t_len, n)
      b2 = [b2_raw, d2].min
      delta = d2 <= 1.0e-18 ? 0.0 : (b2 / d2).clamp(0.0, 1.0)

      shrunk = Array.new(n) do |i|
        Array.new(n) { |j| ((1.0 - delta) * s_mat[i][j]) + (i == j ? delta * mu : 0.0) }
      end
      { sigma: shrunk, delta: delta, target_mu: mu, sample_cov: s_mat }
    end

    def frobenius_sq_diff_diag(s_mat, mu, n)
      total = 0.0
      n.times do |i|
        n.times do |j|
          diff = s_mat[i][j] - (i == j ? mu : 0.0)
          total += diff * diff
        end
      end
      total
    end

    def estimate_b2_bar(centered, s_mat, t_len, n)
      sum_sq = 0.0
      centered.each do |row|
        n.times do |i|
          ri = row[i]
          n.times do |j|
            diff = (ri * row[j]) - s_mat[i][j]
            sum_sq += diff * diff
          end
        end
      end
      sum_sq / (t_len.to_f * t_len)
    end

    def shrink_cov_diag(sample_cov, delta: 0.3)
      n = sample_cov.length
      Array.new(n) do |i|
        Array.new(n) { |j| i == j ? sample_cov[i][i] : (1.0 - delta) * sample_cov[i][j] }
      end
    end

    # Jorion / James-Stein shrinkage of sample mean vector toward grand mean.
    def shrink_mu_js(mu_hat, sigma, t_len)
      n = mu_hat.length
      return mu_hat.dup if n <= 2

      ones = Array.new(n, 1.0)
      inv_s_mu = LinAlg.solve_spd(sigma, mu_hat)
      inv_s_1 = LinAlg.solve_spd(sigma, ones)
      mu_gmv = inv_s_1.each_with_index.sum { |v, i| v * mu_hat[i] } / [inv_s_1.sum, 1.0e-15].max
      diff = mu_hat.map { |m| m - mu_gmv }
      inv_s_diff = LinAlg.solve_spd(sigma, diff)
      quad = LinAlg.dot(diff, inv_s_diff)
      phi = ((n - 2).to_f / [t_len * quad, 1.0e-12].max).clamp(0.0, 1.0)
      _ = inv_s_mu
      mu_hat.each_with_index.map { |m, _i| ((1.0 - phi) * m) + (phi * mu_gmv) }
    end
  end
end
