# frozen_string_literal: true

module FactorPortfolioLab
  # Portfolio risk and performance metrics (annualized assuming monthly periods by default).
  module Metrics
    module_function

    def summarize_returns(port_returns, turnovers: [], hhis: [], periods_per_year: 12.0)
      t_len = port_returns.length
      mean_r = port_returns.sum / t_len
      var_r = port_returns.sum { |r| (r - mean_r)**2 } / [t_len - 1, 1].max
      sd_r = Math.sqrt([var_r, 0.0].max)
      ann_ret = mean_r * periods_per_year
      ann_vol = sd_r * Math.sqrt(periods_per_year)
      sharpe = ann_vol <= 1e-12 ? 0.0 : ann_ret / ann_vol

      downside = port_returns.select(&:negative?)
      down_sd = downside.empty? ? 1.0e-12 : Math.sqrt(downside.sum { |r| r * r } / t_len)
      sortino = (mean_r * periods_per_year) / (down_sd * Math.sqrt(periods_per_year))

      mdd = max_drawdown(port_returns)
      var_95, cvar_95 = var_cvar(port_returns, alpha: 0.05)

      {
        ann_return: ann_ret,
        ann_vol: ann_vol,
        sharpe: sharpe,
        sortino: sortino,
        max_drawdown: mdd,
        calmar: mdd <= 1e-12 ? 0.0 : ann_ret / mdd,
        var_95: var_95,
        cvar_95: cvar_95,
        mean_turnover: turnovers.empty? ? 0.0 : turnovers.sum / turnovers.length,
        mean_hhi: hhis.empty? ? 0.0 : hhis.sum / hhis.length
      }
    end

    def max_drawdown(returns)
      wealth = 1.0
      peak = 1.0
      mdd = 0.0
      returns.each do |r|
        wealth *= (1.0 + r)
        peak = wealth if wealth > peak
        dd = (peak - wealth) / peak
        mdd = dd if dd > mdd
      end
      mdd
    end

    def var_cvar(returns, alpha: 0.05)
      sorted = returns.sort
      idx = [(sorted.length * alpha).floor, 1].max
      var = -sorted[idx - 1]
      cvar = -sorted[0...idx].sum / idx
      [var, cvar]
    end

    def turnover(w_prev, w_curr)
      w_prev.each_with_index.sum { |wp, i| (w_curr[i] - wp).abs }
    end

    def diversification_ratio(weights, sigma)
      port_vol = Math.sqrt([LinAlg.quad_form(weights, sigma), 1.0e-18].max)
      weighted_vol = weights.each_with_index.sum { |w, i| w * Math.sqrt([sigma[i][i], 0.0].max) }
      weighted_vol / port_vol
    end
  end
end
