# frozen_string_literal: true

module GarchFamily
  class Garch11 < VolModel::Base
    DELTA = 1.0e-4

    def name = 'garch11'

    def param_names(dist)
      base = %i[mu omega alpha beta]
      dist.is_a?(Innovations::StudentT) ? base + [:nu] : base
    end

    def unconditional_variance(params)
      den = [1.0 - persistence(params), 1.0e-12].max
      params[:omega] / den
    end

    def persistence(params) = params[:alpha] + params[:beta]

    def fourth_moment_slack(params, kappa_z: 3.0)
      a = params[:alpha]
      b = params[:beta]
      1.0 - (kappa_z * a * a) - (2.0 * a * b) - (b * b)
    end

    def theoretical_kurtosis(params)
      a = params[:alpha]
      b = params[:beta]
      p = a + b
      den = 1.0 - (p * p) - (2.0 * a * a)
      return Float::INFINITY if den <= 0.0

      3.0 * (1.0 - (p * p)) / den
    end

    def theoretical_acf(params, lag)
      return Float::NAN if fourth_moment_slack(params) <= 0.0

      a = params[:alpha]
      b = params[:beta]
      p = a + b
      rho1 = a * (1.0 - (b * p)) / (1.0 - (2.0 * b * p) + (b * b))
      lag == 1 ? rho1 : rho1 * (p**(lag - 1))
    end

    def forecast_variance(params, next_var, horizon)
      v_bar = unconditional_variance(params)
      p = persistence(params)
      v_bar + ((p**(horizon - 1)) * (next_var - v_bar))
    end

    def term_structure(params, cur_var, horizons)
      horizons.map do |h|
        v_h = forecast_variance(params, cur_var, h)
        Math.sqrt(v_h * 252.0)
      end
    end

    def news_impact(params, eps_grid, base_var: nil)
      h_base = base_var || unconditional_variance(params)
      eps_grid.map { |e| params[:omega] + (params[:alpha] * e * e) + (params[:beta] * h_base) }
    end

    def simulate(params, n_obs, dist, rng, burn: 500)
      mu = params.fetch(:mu, 0.0)
      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      h = unconditional_variance(params)
      rets = Array.new(n_obs)
      vars = Array.new(n_obs)
      eps = 0.0
      (n_obs + burn).times do |t|
        h = params[:omega] + (params[:alpha] * eps * eps) + (params[:beta] * h) if t.positive?
        eps = Math.sqrt(h) * dist.sample(rng, nu)
        if t >= burn
          rets[t - burn] = mu + eps
          vars[t - burn] = h
        end
      end
      { returns: rets, true_vars: vars }
    end

    def filter(params, returns)
      mu = params.fetch(:mu, 0.0)
      h = unconditional_variance(params)
      vars = Array.new(returns.length, h)
      (1...returns.length).each do |t|
        ep = returns[t - 1] - mu
        h = params[:omega] + (params[:alpha] * ep * ep) + (params[:beta] * h)
        vars[t] = h
      end
      vars
    end

    def score_contributions(params, returns, dist)
      mu = params.fetch(:mu, 0.0)
      w = params[:omega]
      a = params[:alpha]
      b = params[:beta]
      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      den = [1.0 - a - b, 1.0e-12].max
      h = w / den
      dh = [0.0, 1.0 / den, h / den, h / den]

      returns.each_with_index.map do |r, t|
        h, dh = step_garch_derivs(t, returns, mu, w, a, b, h, dh)
        obs_grad_from_dh(r - mu, h, dh, dist, nu)
      end
    end

    def transform(params, dist)
      lim = 1.0 - DELTA
      s = [1.0 - ((params[:alpha] + params[:beta]) / lim), 1.0e-6].max
      u = [params.fetch(:mu, 0.0), Math.log([params[:omega], 1.0e-12].max),
           Math.log([params[:alpha] / (lim * s), 1.0e-6].max),
           Math.log([params[:beta] / (lim * s), 1.0e-6].max)]
      u << Math.log([params.fetch(:nu, dist.nu) - 2.05, 1.0e-4].max) if dist.is_a?(Innovations::StudentT)
      u
    end

    def untransform(theta, dist)
      lim = 1.0 - DELTA
      ea = Math.exp(theta[2].clamp(-25.0, 25.0))
      eb = Math.exp(theta[3].clamp(-25.0, 25.0))
      norm = 1.0 + ea + eb
      out = { mu: theta[0], omega: Math.exp(theta[1].clamp(-30.0, 5.0)),
              alpha: lim * ea / norm, beta: lim * eb / norm }
      out[:nu] = 2.05 + Math.exp(theta[4].clamp(-12.0, 8.0)) if dist.is_a?(Innovations::StudentT)
      out
    end

    private

    def step_garch_derivs(t, returns, mu, w, a, b, h, dh)
      return [h, dh] if t.zero?

      ep = returns[t - 1] - mu
      [w + (a * ep * ep) + (b * h),
       [(-2.0 * a * ep) + (b * dh[0]), 1.0 + (b * dh[1]), (ep * ep) + (b * dh[2]), h + (b * dh[3])]]
    end
  end
end
