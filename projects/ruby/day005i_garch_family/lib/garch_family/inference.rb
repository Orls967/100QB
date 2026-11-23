# frozen_string_literal: true

module GarchFamily
  module Inference
    module_function

    def neg_hessian_from_score(model, params, returns, dist)
      keys = model.param_names(dist)
      d = keys.length
      raw = Array.new(d) { Array.new(d, 0.0) }
      keys.each_with_index do |k, j|
        h = param_step(k, params[k])
        gp = model.score(params.merge(k => params[k] + h), returns, dist)
        gm = model.score(params.merge(k => params[k] - h), returns, dist)
        d.times { |i| raw[i][j] = -(gp[i] - gm[i]) / (2.0 * h) }
      end
      Array.new(d) { |i| Array.new(d) { |j| 0.5 * (raw[i][j] + raw[j][i]) } }
    end

    def opg_matrix(model, params, returns, dist)
      rows = model.score_contributions(params, returns, dist)
      d = rows.first.length
      b = Array.new(d) { Array.new(d, 0.0) }
      rows.each do |g|
        d.times do |i|
          gi = g[i]
          d.times { |j| b[i][j] += gi * g[j] }
        end
      end
      b
    end

    def param_step(key, val)
      return [val.abs * 1.0e-4, 1.0e-9].max if key == :omega

      [val.abs * 1.0e-5, 1.0e-6].max
    end
    private_class_method :param_step

    def covariance_report(model, params, returns, dist)
      a = neg_hessian_from_score(model, params, returns, dist)
      b = opg_matrix(model, params, returns, dist)
      a_inv = safe_invert(a)
      return nil if a_inv.nil?

      cov_sw = Linalg.mat_mul(Linalg.mat_mul(a_inv, b), a_inv)
      se_naive = Array.new(a.length) { |i| Math.sqrt([a_inv[i][i], 0.0].max) }
      se_sw = Array.new(a.length) { |i| Math.sqrt([cov_sw[i][i], 0.0].max) }
      { a_mat: a, b_mat: b, cov_naive: a_inv, cov_sandwich: cov_sw,
        se_naive: se_naive, se_sandwich: se_sw }
    end

    def wald_test(est_vec, cov, indices, null_vec = nil)
      null_v = null_vec || Array.new(indices.length, 0.0)
      diff = Array.new(indices.length) { |i| est_vec[indices[i]] - null_v[i] }
      sub = indices.map { |r| indices.map { |c| cov[r][c] } }
      inv_sub = safe_invert(sub)
      return { stat: Float::NAN, p_value: Float::NAN } if inv_sub.nil?

      stat = Linalg.quad_form(diff, inv_sub)
      { stat: stat, p_value: Special.chi2_sf(stat, indices.length) }
    end

    def lr_test(ll_full, ll_null, df: 1, mixture_zero_one: false)
      stat = [2.0 * (ll_full - ll_null), 0.0].max
      pval = mixture_zero_one ? Special.chi2_mixture_01_sf(stat) : Special.chi2_sf(stat, df)
      { stat: stat, p_value: pval }
    end

    def safe_invert(mat)
      Linalg.invert(mat)
    rescue StandardError
      nil
    end
    private_class_method :safe_invert
  end
end
