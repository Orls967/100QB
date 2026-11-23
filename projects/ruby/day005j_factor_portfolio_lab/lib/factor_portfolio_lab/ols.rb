# frozen_string_literal: true

module FactorPortfolioLab
  # Ordinary Least Squares regression with classical, White HC0/HC1, and Newey-West HAC standard errors.
  module Ols
    module_function

    # Fits y = X * beta + eps. X is T x P, y is length T.
    def fit(x_mat, y_vec, hac_lags: 0)
      t_len = y_vec.length
      p_dim = x_mat[0].length
      xtx = gram_matrix(x_mat, p_dim)
      xty = cross_vec(x_mat, y_vec, p_dim)
      xtx_inv = LinAlg.invert(xtx)
      beta = LinAlg.mat_vec(xtx_inv, xty)
      residuals = compute_residuals(x_mat, y_vec, beta)
      dof = [t_len - p_dim, 1].max
      sigma2 = residuals.sum { |e| e * e } / dof
      vcov_ols = scale_matrix(xtx_inv, sigma2)
      vcov_hc1 = white_vcov(x_mat, residuals, xtx_inv, t_len, p_dim)
      vcov_hac = hac_lags.positive? ? newey_west_vcov(x_mat, residuals, xtx_inv, hac_lags, p_dim) : vcov_hc1

      {
        beta: beta,
        residuals: residuals,
        sigma2: sigma2,
        r2: r_squared(y_vec, residuals),
        xtx_inv: xtx_inv,
        vcov_ols: vcov_ols,
        vcov_hc1: vcov_hc1,
        vcov_hac: vcov_hac,
        se_ols: diag_sqrt(vcov_ols),
        se_hc1: diag_sqrt(vcov_hc1),
        se_hac: diag_sqrt(vcov_hac)
      }
    end

    def gram_matrix(x_mat, p_dim)
      gram = Array.new(p_dim) { Array.new(p_dim, 0.0) }
      x_mat.each do |row|
        p_dim.times do |i|
          ri = row[i]
          j = i
          while j < p_dim
            gram[i][j] += ri * row[j]
            j += 1
          end
        end
      end
      p_dim.times do |i|
        ((i + 1)...p_dim).each { |j| gram[j][i] = gram[i][j] }
      end
      gram
    end

    def cross_vec(x_mat, y_vec, p_dim)
      vec = Array.new(p_dim, 0.0)
      x_mat.each_with_index do |row, t|
        yt = y_vec[t]
        p_dim.times { |j| vec[j] += row[j] * yt }
      end
      vec
    end

    def compute_residuals(x_mat, y_vec, beta)
      y_vec.each_with_index.map do |yt, t|
        row = x_mat[t]
        pred = 0.0
        row.each_with_index { |xij, j| pred += xij * beta[j] }
        yt - pred
      end
    end

    def r_squared(y_vec, residuals)
      mean_y = y_vec.sum / y_vec.length
      tss = y_vec.sum { |y| (y - mean_y)**2 }
      rss = residuals.sum { |e| e * e }
      tss <= 1e-18 ? 0.0 : 1.0 - (rss / tss)
    end

    def white_vcov(x_mat, residuals, xtx_inv, t_len, p_dim)
      meat = Array.new(p_dim) { Array.new(p_dim, 0.0) }
      scale = t_len.to_f / [t_len - p_dim, 1].max
      x_mat.each_with_index do |row, t|
        e2 = residuals[t] * residuals[t] * scale
        p_dim.times do |i|
          ri_e2 = row[i] * e2
          j = i
          while j < p_dim
            meat[i][j] += ri_e2 * row[j]
            j += 1
          end
        end
      end
      symmetrize!(meat, p_dim)
      sandwich(xtx_inv, meat)
    end

    def newey_west_vcov(x_mat, residuals, xtx_inv, lags, p_dim)
      t_len = x_mat.length
      meat = Array.new(p_dim) { Array.new(p_dim, 0.0) }
      t_len.times do |t|
        add_outer!(meat, x_mat[t], x_mat[t], residuals[t] * residuals[t], p_dim)
      end
      (1..lags).each do |lag|
        weight = 1.0 - (lag.to_f / (lags + 1))
        (lag...t_len).each do |t|
          coef = weight * residuals[t] * residuals[t - lag]
          add_sym_outer!(meat, x_mat[t], x_mat[t - lag], coef, p_dim)
        end
      end
      sandwich(xtx_inv, meat)
    end

    def add_outer!(mat, u, v, weight, p_dim)
      p_dim.times do |i|
        ui = u[i] * weight
        p_dim.times { |j| mat[i][j] += ui * v[j] }
      end
    end

    def add_sym_outer!(mat, u, v, weight, p_dim)
      p_dim.times do |i|
        ui = u[i] * weight
        vi = v[i] * weight
        p_dim.times { |j| mat[i][j] += (ui * v[j]) + (vi * u[j]) }
      end
    end

    def symmetrize!(mat, p_dim)
      p_dim.times do |i|
        ((i + 1)...p_dim).each { |j| mat[j][i] = mat[i][j] }
      end
    end

    def sandwich(bread, meat)
      LinAlg.mat_mul(LinAlg.mat_mul(bread, meat), bread)
    end

    def scale_matrix(mat, scalar)
      mat.map { |row| row.map { |v| v * scalar } }
    end

    def diag_sqrt(mat)
      mat.each_with_index.map { |row, i| Math.sqrt([row[i], 0.0].max) }
    end
  end
end
