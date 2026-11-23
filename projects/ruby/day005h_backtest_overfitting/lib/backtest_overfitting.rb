# frozen_string_literal: true

require_relative 'backtest_overfitting/rng'
require_relative 'backtest_overfitting/parallel'
require_relative 'backtest_overfitting/special'
require_relative 'backtest_overfitting/data_generator'
require_relative 'backtest_overfitting/strategy'
require_relative 'backtest_overfitting/backtest_engine'
require_relative 'backtest_overfitting/metrics'
require_relative 'backtest_overfitting/sharpe_stats'
require_relative 'backtest_overfitting/selection_bias'
require_relative 'backtest_overfitting/multiple_testing'
require_relative 'backtest_overfitting/cscv_pbo'
require_relative 'backtest_overfitting/purged_cv'
require_relative 'backtest_overfitting/json_out'
require_relative 'backtest_overfitting/svg'
require_relative 'backtest_overfitting/criteria'
require_relative 'backtest_overfitting/figures'
require_relative 'backtest_overfitting/runner'

module BacktestOverfitting
  MAIN_SEED = 20_261_123
end
