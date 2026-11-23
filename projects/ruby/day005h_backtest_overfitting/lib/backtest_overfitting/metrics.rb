# frozen_string_literal: true

module BacktestOverfitting
  module Metrics
    module_function

    def basic_moments(returns)
      n = returns.length
      raise ArgumentError, 'returns kosong' if n < 2

      m = returns.sum(0.0) / n
      m2 = 0.0
      m3 = 0.0
      m4 = 0.0
      returns.each do |r|
        d = r - m
        d2 = d * d
        m2 += d2
        m3 += d2 * d
        m4 += d2 * d2
      end
      var = m2 / (n - 1)
      sd = Math.sqrt([var, 0.0].max)
      [m, sd, n, m2 / n, m3 / n, m4 / n]
    end

    def sharpe_per_period(returns)
      m, sd, = basic_moments(returns)
      sd.positive? ? m / sd : 0.0
    end

    def annualized_sharpe(returns, periods_per_year: 252)
      sharpe_per_period(returns) * Math.sqrt(periods_per_year)
    end

    def max_drawdown(returns, init_equity: 1.0)
      eq = init_equity
      peak = init_equity
      mdd = 0.0
      returns.each do |r|
        eq *= (1.0 + r)
        peak = eq if eq > peak
        dd = peak.positive? ? (peak - eq) / peak : 0.0
        mdd = dd if dd > mdd
      end
      mdd
    end

    def summary(returns, periods_per_year: 252)
      m, sd, n, pop_var, m3, m4 = basic_moments(returns)
      sr = sd.positive? ? m / sd : 0.0
      skew = pop_var.positive? ? m3 / (pop_var**1.5) : 0.0
      kurt = pop_var.positive? ? m4 / (pop_var * pop_var) : 3.0
      {
        n_obs: n, mean: m, std: sd, sharpe_per_period: sr,
        annualized_sharpe: sr * Math.sqrt(periods_per_year),
        skewness: skew, kurtosis: kurt, max_drawdown: max_drawdown(returns)
      }
    end
  end
end
