# frozen_string_literal: true

module FactorPortfolioLab
  # Heavy-tail Student-t(4) robustness study across GRS size, Shanken coverage, and portfolio CVaR.
  module Robustness
    module_function

    def run(master_seed:, n_grs_reps: 800, n_port_worlds: 120)
      world = FactorWorld.build_world(n_assets: 20, n_factors: 3, seed: Rng.derive_seed(master_seed, 1200, 0))

      gauss_grs = grs_size(world, master_seed, 1201, n_grs_reps, t_len: 120, df: nil)
      t4_grs = grs_size(world, master_seed, 1202, n_grs_reps, t_len: 120, df: 4.0)

      gauss_port = portfolio_tail(world, master_seed, 1203, n_port_worlds, df: nil)
      t4_port = portfolio_tail(world, master_seed, 1204, n_port_worlds, df: 4.0)

      {
        grs_gaussian_rej: gauss_grs,
        grs_t4_rej: t4_grs,
        portfolio_gaussian: gauss_port,
        portfolio_t4: t4_port
      }
    end

    def grs_size(world, master_seed, stream_id, n_reps, t_len:, df:)
      reps = Parallel.map_items(n_items: n_reps, chunk_size: 40) do |rep|
        seed = Rng.derive_seed(master_seed, stream_id, rep)
        sim = FactorWorld.simulate(world, t_len: t_len, seed: seed, df: df)
        res = Grs.test(sim[:returns], sim[:factors])
        [res[:p_value_f] < 0.05 ? 1.0 : 0.0, res[:p_value_wald] < 0.05 ? 1.0 : 0.0]
      end
      f_rate = reps.sum { |r| r[0] } / n_reps
      w_rate = reps.sum { |r| r[1] } / n_reps
      { grs_rej_05: f_rate, wald_rej_05: w_rate, n_reps: n_reps }
    end

    def portfolio_tail(world, master_seed, stream_id, n_worlds, df:)
      strats = [Strategy::EqualWeight.new, Strategy::MinVarShrink.new, Strategy::RiskParitySample.new]
      rows = Parallel.map_items(n_items: n_worlds, chunk_size: 10) do |w_idx|
        seed = Rng.derive_seed(master_seed, stream_id, w_idx)
        sim = FactorWorld.simulate(world, t_len: 150, seed: seed, df: df)
        strats.to_h do |s|
          bt = Backtest.run(returns: sim[:returns], strategy: s, window: 60, rebalance_every: 3, tc_bps: 5.0)
          [s.name, bt[:metrics]]
        end
      end

      strats.to_h do |s|
        cvar = rows.sum { |r| r[s.name][:cvar_95] } / n_worlds
        sharpe = rows.sum { |r| r[s.name][:sharpe] } / n_worlds
        mdd = rows.sum { |r| r[s.name][:max_drawdown] } / n_worlds
        [s.name, { cvar_95: cvar, sharpe: sharpe, max_drawdown: mdd }]
      end
    end
  end
end
