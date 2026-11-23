# frozen_string_literal: true

module BacktestOverfitting
  module BacktestEngine
    module_function

    def run(returns, positions, cost_bps: 0.0, init_equity: 1.0)
      n = returns.length
      raise ArgumentError, 'Panjang returns dan positions tidak sama' unless positions.length == n

      rate = cost_bps * 1.0e-4
      gross = Array.new(n, 0.0)
      net = Array.new(n, 0.0)
      turnovers = Array.new(n, 0.0)
      equity = Array.new(n + 1, init_equity)
      prev_pos = 0.0

      n.times do |t|
        pos = positions[t].to_f
        raise ArgumentError, "Posisi tidak valid: #{pos}" unless [-1.0, 0.0, 1.0].include?(pos)

        to = (pos - prev_pos).abs
        g = pos * returns[t]
        nt = g - (to * rate)
        gross[t] = g
        net[t] = nt
        turnovers[t] = to
        equity[t + 1] = equity[t] * (1.0 + nt)
        prev_pos = pos
      end

      { gross_returns: gross, net_returns: net, turnovers: turnovers,
        total_turnover: turnovers.sum(0.0), equity_curve: equity, final_equity: equity.last }
    end

    def evaluate_library(strategies, returns, cost_bps: 0.0)
      strategies.map do |strat|
        pos = strat.generate_positions(returns)
        res = run(returns, pos, cost_bps: cost_bps)
        net = res[:net_returns]
        { name: strat.name, net_returns: net, sharpe: Metrics.sharpe_per_period(net),
          total_turnover: res[:total_turnover], final_equity: res[:final_equity] }
      end
    end

    def verify_no_lookahead(strategy, returns, split_idx, rng)
      pos_orig = strategy.generate_positions(returns)
      mutated = returns.dup
      (split_idx...returns.length).each { |i| mutated[i] = 0.05 * rng.next_normal }
      pos_mut = strategy.generate_positions(mutated)
      (0..split_idx).all? { |t| pos_orig[t] == pos_mut[t] }
    end

    def verify_equity_identity(run_result, init_equity: 1.0)
      prod = run_result[:net_returns].reduce(init_equity) { |acc, r| acc * (1.0 + r) }
      (run_result[:final_equity] - prod).abs
    end
  end
end
