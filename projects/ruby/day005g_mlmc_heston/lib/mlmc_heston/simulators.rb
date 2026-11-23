# frozen_string_literal: true

module MlmcHeston
  module Simulators
    module_function

    def cir_exact_mean(v0, kappa, theta, t)
      theta + ((v0 - theta) * Math.exp(-kappa * t))
    end

    def cir_exact_var(v0, kappa, theta, xi, t)
      ek = Math.exp(-kappa * t)
      e2k = ek * ek
      (v0 * (xi * xi / kappa) * (ek - e2k)) + ((theta * xi * xi / (2.0 * kappa)) * ((1.0 - ek)**2))
    end

    # Pasangan jalur terkopel (fine, coarse) untuk GBM (skema :euler atau :milstein).
    def gbm_coupled_pair(params, strike, maturity, level, scheme:, rng:, n0: 4)
      spot = params.fetch(:spot)
      rate = params.fetch(:rate)
      sigma = params.fetch(:sigma)
      n_f = n0 * (2**level)
      dt_f = maturity / n_f
      sqrt_dt_f = Math.sqrt(dt_f)
      disc = Math.exp(-rate * maturity)

      s_f = spot
      s_c = spot
      dw_c = 0.0
      dt_c = 2.0 * dt_f

      i = 0
      while i < n_f
        dw = sqrt_dt_f * rng.next_normal
        s_f = step_gbm(s_f, rate, sigma, dt_f, dw, scheme)
        dw_c += dw
        if i.odd?
          s_c = step_gbm(s_c, rate, sigma, dt_c, dw_c, scheme) if level.positive?
          dw_c = 0.0
        end
        i += 1
      end

      pf = disc * [s_f - strike, 0.0].max
      pc = level.positive? ? disc * [s_c - strike, 0.0].max : 0.0
      [pf, pf - pc]
    end

    def step_gbm(spot, rate, sigma, dt, dw, scheme)
      mil = scheme == :milstein ? 0.5 * sigma * sigma * spot * ((dw * dw) - dt) : 0.0
      spot + (rate * spot * dt) + (sigma * spot * dw) + mil
    end
    private_class_method :step_gbm

    # Pasangan jalur terkopel Heston Full-Truncation Euler (atau Milstein pada S) untuk European & Asian Call.
    def heston_coupled_pair(params, strike, maturity, level, scheme:, payoff_type:, rng:, n0: 4)
      n_f = n0 * (2**level)
      dt_f = maturity / n_f
      sqrt_dt_f = Math.sqrt(dt_f)
      rho = params.fetch(:rho)
      s1r = Math.sqrt(1.0 - (rho * rho))

      state_f = [Math.log(params.fetch(:spot)), params.fetch(:v0), 0.0]
      state_c = [state_f[0], state_f[1], 0.0]
      dw1_c = 0.0
      dw2_c = 0.0

      i = 0
      while i < n_f
        z1 = rng.next_normal
        z2 = (rho * z1) + (s1r * rng.next_normal)
        dw1 = sqrt_dt_f * z1
        dw2 = sqrt_dt_f * z2
        step_heston_ft!(state_f, params, dt_f, dw1, dw2, scheme)
        dw1_c += dw1
        dw2_c += dw2
        if i.odd?
          step_heston_ft!(state_c, params, 2.0 * dt_f, dw1_c, dw2_c, scheme) if level.positive?
          dw1_c = 0.0
          dw2_c = 0.0
        end
        i += 1
      end

      evaluate_heston_payoffs(state_f, state_c, params, strike, maturity, level, n_f, payoff_type)
    end

    def step_heston_ft!(st, p, dt, dw1, dw2, scheme)
      v_pos = st[1].positive? ? st[1] : 0.0
      sqrt_v = Math.sqrt(v_pos)
      # Log-spot Euler dengan volatilitas full-truncation (+ opsional koreksi Milstein pada S dalam koordinat S)
      st[0] += ((p.fetch(:rate) - (0.5 * v_pos)) * dt) + (sqrt_v * dw1)
      mil_v = scheme == :milstein ? 0.25 * p.fetch(:xi) * p.fetch(:xi) * ((dw2 * dw2) - dt) : 0.0
      st[1] += (p.fetch(:kappa) * (p.fetch(:theta) - v_pos) * dt) + (p.fetch(:xi) * sqrt_v * dw2) + mil_v
      st[2] += Math.exp(st[0])
    end
    private_class_method :step_heston_ft!

    def evaluate_heston_payoffs(sf, sc, params, strike, maturity, level, n_f, payoff_type)
      disc = Math.exp(-params.fetch(:rate) * maturity)
      s_f = payoff_type == :asian ? sf[2] / n_f : Math.exp(sf[0])
      pf = disc * [s_f - strike, 0.0].max
      return [pf, pf] if level.zero?

      s_c = payoff_type == :asian ? sc[2] / (n_f / 2) : Math.exp(sc[0])
      pc = disc * [s_c - strike, 0.0].max
      [pf, pf - pc]
    end
    private_class_method :evaluate_heston_payoffs
  end
end
