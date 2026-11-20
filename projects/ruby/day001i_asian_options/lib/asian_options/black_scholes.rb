# frozen_string_literal: true

module AsianOptions
  module BlackScholes
    module_function

    def bs_call(params)
      params.validate!
      return intrinsic_call(params) if params.t <= 0.0 || params.sigma <= 0.0

      d1, d2 = d1_d2(params)
      (params.s0 * Stats.norm_cdf(d1)) - (params.k * Math.exp(-params.r * params.t) * Stats.norm_cdf(d2))
    end

    def bs_put(params)
      params.validate!
      return intrinsic_put(params) if params.t <= 0.0 || params.sigma <= 0.0

      d1, d2 = d1_d2(params)
      (params.k * Math.exp(-params.r * params.t) * Stats.norm_cdf(-d2)) - (params.s0 * Stats.norm_cdf(-d1))
    end

    # Rumus tertutup opsi call Asia geometrik monitoring diskret (Kemna-Vorst).
    def geometric_asian_call(params)
      params.validate!
      mu_ln_g, var_ln_g = geometric_log_moments(params)
      sigma_ln_g = Math.sqrt(var_ln_g)
      expected_g = Math.exp(mu_ln_g + (0.5 * var_ln_g))
      d1 = (mu_ln_g - Math.log(params.k) + var_ln_g) / sigma_ln_g
      d2 = d1 - sigma_ln_g
      Math.exp(-params.r * params.t) * ((expected_g * Stats.norm_cdf(d1)) - (params.k * Stats.norm_cdf(d2)))
    end

    def geometric_log_moments(params)
      n = params.n_steps.to_f
      sigma_sq = params.sigma * params.sigma
      mu_ln_g = Math.log(params.s0) + ((params.r - (0.5 * sigma_sq)) * params.t * (n + 1.0) / (2.0 * n))
      var_ln_g = sigma_sq * params.t * (n + 1.0) * ((2.0 * n) + 1.0) / (6.0 * n * n)
      [mu_ln_g, var_ln_g]
    end
    private_class_method :geometric_log_moments

    def intrinsic_call(params)
      [params.s0 - (params.k * Math.exp(-params.r * params.t)), 0.0].max
    end
    private_class_method :intrinsic_call

    def intrinsic_put(params)
      [(params.k * Math.exp(-params.r * params.t)) - params.s0, 0.0].max
    end
    private_class_method :intrinsic_put

    def d1_d2(params)
      vol_sqrt_t = params.sigma * Math.sqrt(params.t)
      drift_term = (params.r + (0.5 * params.sigma * params.sigma)) * params.t
      d1 = (Math.log(params.s0 / params.k) + drift_term) / vol_sqrt_t
      [d1, d1 - vol_sqrt_t]
    end
    private_class_method :d1_d2
  end
end
