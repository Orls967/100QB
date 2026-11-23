# frozen_string_literal: true

module GarchFamily
  class Egarch11 < VolModel::Base
    def name = 'egarch11'

    def param_names(dist)
      base = %i[mu omega alpha gamma beta]
      dist.is_a?(Innovations::StudentT) ? base + [:nu] : base
    end

    def persistence(params) = params[:beta]

    def unconditional_variance(params)
      w, a, gm, b = params.values_at(:omega, :alpha, :gamma, :beta)
      prod = 1.0
      bk = 1.0
      200.times do
        prod *= normal_abs_mgf(bk * a, bk * gm) * Math.exp(-bk * a * Special::SQRT_2_OVER_PI)
        bk *= b
        break if bk.abs < 1.0e-12
      end
      Math.exp(w / [1.0 - b, 1.0e-12].max) * prod
    end

    def forecast_variance(params, next_var, horizon)
      lv = Math.log(next_var)
      w = params[:omega]
      b = params[:beta]
      b_pow = b**(horizon - 1)
      log_h = ((w / (1.0 - b)) * (1.0 - b_pow)) + (b_pow * lv)
      Math.exp(log_h.clamp(-35.0, 20.0))
    end

    def news_impact(params, eps_grid, base_var: nil)
      h_base = base_var || unconditional_variance(params)
      sqrt_h = Math.sqrt(h_base)
      lh_base = Math.log(h_base)
      ez = Special::SQRT_2_OVER_PI
      eps_grid.map do |e|
        z = e / sqrt_h
        lh = params[:omega] + (params[:alpha] * (z.abs - ez)) + (params[:gamma] * z) + (params[:beta] * lh_base)
        Math.exp(lh.clamp(-35.0, 20.0))
      end
    end

    def simulate(params, n_obs, dist, rng, burn: 500)
      mu = params.fetch(:mu, 0.0)
      ez = Special::SQRT_2_OVER_PI
      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      lh = params[:omega] / [1.0 - params[:beta], 1.0e-12].max
      rets = Array.new(n_obs)
      vars = Array.new(n_obs)
      z = 0.0
      (n_obs + burn).times do |t|
        lh = next_log_h(params, z, ez, lh) if t.positive?
        h = Math.exp(lh)
        z = dist.sample(rng, nu)
        next if t < burn

        rets[t - burn] = mu + (Math.sqrt(h) * z)
        vars[t - burn] = h
      end
      { returns: rets, true_vars: vars }
    end

    def filter(params, returns)
      mu = params.fetch(:mu, 0.0)
      ez = Special::SQRT_2_OVER_PI
      lh = (params[:omega] / [1.0 - params[:beta], 1.0e-12].max).clamp(-35.0, 20.0)
      vars = Array.new(returns.length, Math.exp(lh))
      (1...returns.length).each do |t|
        lh = next_log_h(params, (returns[t - 1] - mu) / Math.sqrt(vars[t - 1]), ez, lh)
        vars[t] = Math.exp(lh)
      end
      vars
    end

    def score_contributions(params, returns, dist)
      vars = filter(params, returns)
      mu = params.fetch(:mu, 0.0)
      nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
      ez = Special::SQRT_2_OVER_PI
      dlh = init_egarch_derivs(params, returns.length)
      1.upto(returns.length - 1) do |t|
        dlh[t] = step_egarch_deriv(params, returns[t - 1] - mu, vars[t - 1], ez, dlh[t - 1])
      end
      returns.each_with_index.map { |r, t| egarch_obs_grad(r - mu, vars[t], dlh[t], dist, nu) }
    end

    def transform(params, dist)
      u = [params.fetch(:mu, 0.0), params[:omega], params[:alpha], params[:gamma],
           Math.atanh(params[:beta].clamp(-0.998, 0.998))]
      dist.is_a?(Innovations::StudentT) ? u + [Math.log([params.fetch(:nu, dist.nu) - 2.05, 1.0e-4].max)] : u
    end

    def untransform(theta, dist)
      out = { mu: theta[0], omega: theta[1], alpha: theta[2], gamma: theta[3],
              beta: 0.999 * Math.tanh(theta[4].clamp(-10.0, 10.0)) }
      out[:nu] = 2.05 + Math.exp(theta[5].clamp(-12.0, 8.0)) if dist.is_a?(Innovations::StudentT)
      out
    end

    private

    def next_log_h(p, z, ez, lh)
      (p[:omega] + (p[:alpha] * (z.abs - ez)) + (p[:gamma] * z) + (p[:beta] * lh)).clamp(-35.0, 20.0)
    end

    def normal_abs_mgf(c1, c2)
      (Math.exp(0.5 * ((c1 + c2)**2)) * Special.norm_cdf(c1 + c2)) +
        (Math.exp(0.5 * ((c1 - c2)**2)) * Special.norm_cdf(c1 - c2))
    end

    def init_egarch_derivs(p, n_obs)
      den = [1.0 - p[:beta], 1.0e-12].max
      out = Array.new(n_obs)
      out[0] = [0.0, 1.0 / den, 0.0, 0.0, p[:omega] / (den * den)]
      out
    end

    def step_egarch_deriv(p, ep, prev_var, ez, prev)
      sq = Math.sqrt(prev_var)
      z = ep / sq
      sgn = z >= 0.0 ? 1.0 : -1.0
      b_eff = p[:beta] - (0.5 * ((p[:alpha] * sgn) + p[:gamma]) * z)
      [-(((p[:alpha] * sgn) + p[:gamma]) / sq) + (b_eff * prev[0]), 1.0 + (b_eff * prev[1]),
       (z.abs - ez) + (b_eff * prev[2]), z + (b_eff * prev[3]), Math.log(prev_var) + (b_eff * prev[4])]
    end

    def egarch_obs_grad(eps, h, dlh_t, dist, nu)
      dl_dlh = h * dist.score_wrt_var(eps, h, nu)
      dl_dmu = dist.score_wrt_mu_direct(eps, h, nu)
      g = [dl_dmu + (dl_dlh * dlh_t[0]), dl_dlh * dlh_t[1], dl_dlh * dlh_t[2], dl_dlh * dlh_t[3], dl_dlh * dlh_t[4]]
      g << dist.score_wrt_nu(eps, h, nu) if dist.is_a?(Innovations::StudentT)
      g
    end
  end
end
