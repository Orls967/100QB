# frozen_string_literal: true

module FactorPortfolioLab
  # Strictly causal rolling-window portfolio backtest engine with transaction costs.
  module Backtest
    module_function

    def run(returns:, strategy:, window: 60, rebalance_every: 3, tc_bps: 5.0, factors: nil)
      t_len = returns.length
      n_assets = returns[0].length
      raise ArgumentError, 't_len <= window' if t_len <= window

      tc_rate = tc_bps * 1.0e-4
      w_curr = Array.new(n_assets, 1.0 / n_assets)
      port_returns = []
      turnovers = []
      hhis = []
      weight_history = []

      (window...t_len).each do |t|
        if ((t - window) % rebalance_every).zero?
          r_win = returns[(t - window)...t]
          f_win = factors ? factors[(t - window)...t] : nil
          w_new = strategy.weights(r_win, factors_window: f_win)
          to = t == window ? 0.0 : Metrics.turnover(w_curr, w_new)
          cost = tc_rate * to
          w_curr = w_new
        else
          to = 0.0
          cost = 0.0
        end

        gross = LinAlg.dot(w_curr, returns[t])
        port_returns << (gross - cost)
        turnovers << to
        hhis << w_curr.sum { |w| w * w }
        weight_history << w_curr if ((t - window) % rebalance_every).zero?
        w_curr = drift_weights(w_curr, returns[t])
      end

      {
        strategy: strategy.name,
        tc_bps: tc_bps,
        port_returns: port_returns,
        turnovers: turnovers,
        weight_history: weight_history,
        metrics: Metrics.summarize_returns(port_returns, turnovers: turnovers, hhis: hhis)
      }
    end

    def drift_weights(weights, r_t)
      grown = weights.each_with_index.map { |w, i| w * (1.0 + r_t[i]) }
      tot = grown.sum
      tot <= 1e-12 ? weights.dup : grown.map { |g| g / tot }
    end
  end
end
