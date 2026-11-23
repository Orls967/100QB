# frozen_string_literal: true

module GarchFamily
  class GjrGarch11 < VolModel::Base
    DELTA = 1.0e-4

    def name = 'gjr_garch11'

    def param_names(dist)
      base = %i[mu omega alpha gamma beta]
      dist.is_a?(Innovations::StudentT) ? base + [:nu] : base
    end

    def persistence(params) = params[:alpha] + (0.5 * params[:gamma]) + params[:beta]

    def unconditional_variance(params)
      den = [1.0 - persistence(params), 1.0e-12].max
      params[:omega] / den
    end

    def forecast_variance(params, next_var, horizon)
      v_bar = unconditional_variance(params)
      p = persistence(params)
      v_bar + ((p**(horizon - 1)) * (next_var - v_bar))
    end

    def news_impact(params, eps_grid, base_var: nil)
      h_base = base_var || unconditional_variance(params)
      eps_grid.map do |e|
        ind = e.negative? ? 1.0 : 0.0
        params[:omega] + ((params[:alpha] + (params[:gamma] * ind)) * e * e) + (params[:beta] * h_base)
      end
    end

    def simulate(params, n_obs, dist, rng, burn: 500)
      mu = params.fetch(:mu, 0.0)
      w = params[:omega]
      a = params[:alpha]
      gm = params[:gamma]
      b = params[:beta]
      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      h = unconditional_variance(params)
      rets = Array.new(n_obs)
      vars = Array.new(n_obs)
      eps = 0.0

      (n_obs + burn).times do |t|
        if t.positive?
          ind = eps.negative? ? 1.0 : 0.0
          h = w + ((a + (gm * ind)) * eps * eps) + (b * h)
        end
        z = dist.sample(rng, nu)
        eps = Math.sqrt(h) * z
        if t >= burn
          rets[t - burn] = mu + eps
          vars[t - burn] = h
        end
      end
      { returns: rets, true_vars: vars }
    end

    def filter(params, returns)
      mu = params.fetch(:mu, 0.0)
      w = params[:omega]
      a = params[:alpha]
      gm = params[:gamma]
      b = params[:beta]
      h = unconditional_variance(params)
      vars = Array.new(returns.length, h)
      returns.each_with_index do |_r, t|
        next if t.zero?

        eps = returns[t - 1] - mu
        ind = eps.negative? ? 1.0 : 0.0
        h = w + ((a + (gm * ind)) * eps * eps) + (b * h)
        vars[t] = h
      end
      vars
    end

    def score_contributions(params, returns, dist)
      mu = params.fetch(:mu, 0.0)
      den = [1.0 - persistence(params), 1.0e-12].max
      h = params[:omega] / den
      dh = [0.0, 1.0 / den, h / den, 0.5 * h / den, h / den]

      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      returns.each_with_index.map do |r, t|
        h, dh = step_gjr_derivs(params, t, returns, mu, h, dh)
        obs_grad_from_dh(r - mu, h, dh, dist, nu)
      end
    end

    def transform(params, dist)
      lim = 1.0 - DELTA
      s = [1.0 - (persistence(params) / lim), 1.0e-6].max
      u = [params.fetch(:mu, 0.0), Math.log([params[:omega], 1.0e-12].max),
           Math.log([params[:alpha] / (lim * s), 1.0e-6].max),
           Math.log([(0.5 * params[:gamma]) / (lim * s), 1.0e-6].max),
           Math.log([params[:beta] / (lim * s), 1.0e-6].max)]
      u << Math.log([params.fetch(:nu, dist.nu) - 2.05, 1.0e-4].max) if dist.is_a?(Innovations::StudentT)
      u
    end

    def untransform(theta, dist)
      lim = 1.0 - DELTA
      ea = Math.exp(theta[2].clamp(-25.0, 25.0))
      eg = Math.exp(theta[3].clamp(-25.0, 25.0))
      eb = Math.exp(theta[4].clamp(-25.0, 25.0))
      norm = 1.0 + ea + eg + eb
      out = { mu: theta[0], omega: Math.exp(theta[1].clamp(-30.0, 5.0)),
              alpha: lim * ea / norm, gamma: 2.0 * lim * eg / norm, beta: lim * eb / norm }
      out[:nu] = 2.05 + Math.exp(theta[5].clamp(-12.0, 8.0)) if dist.is_a?(Innovations::StudentT)
      out
    end

    private

    def step_gjr_derivs(p, t, returns, mu, h, dh)
      return [h, dh] if t.zero?

      ep = returns[t - 1] - mu
      ind = ep.negative? ? 1.0 : 0.0
      coef = p[:alpha] + (p[:gamma] * ind)
      [p[:omega] + (coef * ep * ep) + (p[:beta] * h),
       [(-2.0 * coef * ep) + (p[:beta] * dh[0]), 1.0 + (p[:beta] * dh[1]),
        (ep * ep) + (p[:beta] * dh[2]), (ind * ep * ep) + (p[:beta] * dh[3]), h + (p[:beta] * dh[4])]]
    end
  end
end
