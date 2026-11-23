# frozen_string_literal: true

module BacktestOverfitting
  module DataGenerator
    module_function

    def iid_gaussian(n_steps, sigma:, rng:, mu: 0.0)
      Array.new(n_steps) { mu + (sigma * rng.next_normal) }
    end

    # AR(1): r_t = mu + phi * (r_{t-1} - mu) + eps_t.
    # Bila marginal_sigma true, simpangan baku stasioner r_t sama persis dengan sigma.
    def ar1_gaussian(n_steps, phi:, sigma:, rng:, mu: 0.0, marginal_sigma: true)
      raise ArgumentError, '|phi| harus < 1' if phi.abs >= 1.0

      eps_sd = marginal_sigma ? sigma * Math.sqrt(1.0 - (phi * phi)) : sigma
      stat_sd = marginal_sigma ? sigma : sigma / Math.sqrt(1.0 - (phi * phi))
      r = mu + (stat_sd * rng.next_normal)
      out = Array.new(n_steps)
      n_steps.times do |t|
        r = mu + (phi * (r - mu)) + (eps_sd * rng.next_normal) if t.positive?
        out[t] = r
      end
      out
    end

    def two_state_regime(n_steps, p00:, p11:, sigma0:, sigma1:, rng:, mu0: 0.0, mu1: 0.0)
      state = rng.next_float < 0.5 ? 0 : 1
      out = Array.new(n_steps)
      n_steps.times do |t|
        u = rng.next_float
        stay_prob = state.zero? ? p00 : p11
        state = u < stay_prob ? state : (1 - state)
        sd = state.zero? ? sigma0 : sigma1
        mu = state.zero? ? mu0 : mu1
        out[t] = mu + (sd * rng.next_normal)
      end
      out
    end

    def student_t_returns(n_steps, df:, sigma:, rng:, mu: 0.0)
      scale = sigma * Math.sqrt((df - 2.0) / df)
      Array.new(n_steps) { mu + (scale * rng.next_student_t(df)) }
    end

    def lognormal_returns(n_steps, mu_log:, sigma_log:, rng:)
      Array.new(n_steps) { Math.exp(mu_log + (sigma_log * rng.next_normal)) - 1.0 }
    end

    def gbm_simple_returns(n_steps, mu:, sigma:, dt:, rng:)
      drift = (mu - (0.5 * sigma * sigma)) * dt
      vol = sigma * Math.sqrt(dt)
      Array.new(n_steps) { Math.expm1(drift + (vol * rng.next_normal)) }
    end
  end
end
