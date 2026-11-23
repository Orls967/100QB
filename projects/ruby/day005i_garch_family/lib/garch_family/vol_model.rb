# frozen_string_literal: true

module GarchFamily
  module VolModel
    class Base
      def simulate(_params, _n_obs, _dist, _rng) = raise(NotImplementedError)

      def filter(_params, _returns) = raise(NotImplementedError)

      def loglik(params, returns, dist)
        vars = filter(params, returns)
        mu = params.fetch(:mu, 0.0)
        nu = dist.is_a?(Innovations::StudentT) ? params.fetch(:nu, dist.nu) : nil
        returns.each_with_index.sum(0.0) { |r, t| dist.obs_loglik(r - mu, vars[t], nu) }
      end

      def score(params, returns, dist)
        rows = score_contributions(params, returns, dist)
        Array.new(rows.first.length) { |j| rows.sum(0.0) { |row| row[j] } }
      end

      def obs_grad_from_dh(eps, h, dh, dist, nu)
        dl_dh = dist.score_wrt_var(eps, h, nu)
        g = [dist.score_wrt_mu_direct(eps, h, nu) + (dl_dh * dh[0])]
        (1...dh.length).each { |j| g << (dl_dh * dh[j]) }
        g << dist.score_wrt_nu(eps, h, nu) if dist.is_a?(Innovations::StudentT)
        g
      end

      def unconditional_variance(_params) = raise(NotImplementedError)

      def persistence(_params) = raise(NotImplementedError)

      def forecast_variance(_params, _last_var, _horizon) = raise(NotImplementedError)

      def transform(_params, _dist) = raise(NotImplementedError)

      def untransform(_theta, _dist) = raise(NotImplementedError)
    end
  end
end
