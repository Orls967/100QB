# frozen_string_literal: true

module FactorPortfolioLab
  # Shanken (1992) errors-in-variables correction for Fama-MacBeth two-pass standard errors.
  module Shanken
    module_function

    def correct(fm_fit, factors)
      t_len = factors.length
      k_factors = factors[0].length
      add_int = fm_fit[:add_intercept]
      offset = add_int ? 1 : 0

      lambda_f = fm_fit[:lambda_hat][offset, k_factors]
      f_bar = Pca.column_means(factors, k_factors)
      sigma_f = Grs.sample_cov_mle(factors, f_bar, k_factors)
      c_sh = LinAlg.quad_form(lambda_f, LinAlg.invert(sigma_f))
      mult = 1.0 + c_sh

      v_fm = fm_fit[:vcov_fm]
      p_dim = v_fm.length
      v_sh = Array.new(p_dim) do |i|
        Array.new(p_dim) do |j|
          sf_ij = i >= offset && j >= offset ? sigma_f[i - offset][j - offset] / t_len : 0.0
          (mult * (v_fm[i][j] - sf_ij)) + sf_ij
        end
      end

      se_sh = Ols.diag_sqrt(v_sh)
      t_sh = FamaMacbeth.compute_t_stats(fm_fit[:lambda_hat], se_sh)

      {
        c_shanken: c_sh,
        multiplier: mult,
        vcov_shanken: v_sh,
        se_shanken: se_sh,
        t_shanken: t_sh,
        se_ratio: se_sh.each_with_index.map { |s, i| fm_fit[:se_fm][i] <= 1e-15 ? 1.0 : s / fm_fit[:se_fm][i] }
      }
    end
  end
end
