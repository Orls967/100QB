# frozen_string_literal: true

module FactorPortfolioLab
  # Large-scale Monte Carlo rolling two-regime backtest across 8 strategies and 3 cost tiers.
  module BacktestStudy
    module_function

    def run(master_seed:, n_worlds: 300, n_assets: 15, t_len: 180, window: 60, tc_levels: [0.0, 5.0, 20.0])
      world_base = FactorWorld.build_world(
        n_assets: n_assets,
        n_factors: 3,
        seed: Rng.derive_seed(master_seed, 900, 0)
      )
      strategies = Strategy.all_strategies

      world_results = Parallel.map_items(n_items: n_worlds, chunk_size: 10) do |w_idx|
        seed = Rng.derive_seed(master_seed, 901, w_idx)
        sim = FactorWorld.simulate_two_regime(world_base, t_len: t_len, seed: seed)
        run_single_world(sim, strategies, window, tc_levels)
      end

      aggregate(strategies, tc_levels, world_results, n_worlds)
    end

    def run_single_world(sim, strategies, window, tc_levels)
      strategies.to_h do |strat|
        base = Backtest.run(
          returns: sim[:returns],
          strategy: strat,
          window: window,
          rebalance_every: 3,
          tc_bps: 0.0,
          factors: sim[:factors]
        )
        by_tc = tc_levels.to_h do |tc|
          if tc.zero?
            [tc, base[:metrics]]
          else
            net_rets = apply_tc(base[:port_returns], base[:turnovers], tc)
            [tc, Metrics.summarize_returns(net_rets, turnovers: base[:turnovers], hhis: [base[:metrics][:mean_hhi]])]
          end
        end
        [strat.name, by_tc]
      end
    end

    def apply_tc(gross_rets, turnovers, tc_bps)
      rate = tc_bps * 1.0e-4
      gross_rets.each_with_index.map { |r, i| r - (rate * turnovers[i]) }
    end

    def aggregate(strategies, tc_levels, world_results, n_worlds)
      summary = {}
      strategies.each do |strat|
        summary[strat.name] = {}
        tc_levels.each do |tc|
          rows = world_results.map { |wr| wr[strat.name][tc] }
          summary[strat.name][tc] = average_metrics(rows, n_worlds)
        end
      end
      { n_worlds: n_worlds, tc_levels: tc_levels, summary: summary }
    end

    def average_metrics(rows, n_worlds)
      keys = rows[0].keys
      keys.to_h do |k|
        vals = rows.map { |r| r[k] }
        mean = vals.sum / n_worlds
        Math.sqrt(vals.sum { |v| (v - mean)**2 } / [n_worlds - 1, 1].max / n_worlds)
        [k, mean]
      end.merge(sharpe_se: compute_se(rows.map { |r| r[:sharpe] }, n_worlds))
    end

    def compute_se(vals, n_worlds)
      mean = vals.sum / n_worlds
      Math.sqrt(vals.sum { |v| (v - mean)**2 } / [n_worlds - 1, 1].max / n_worlds)
    end
  end
end
