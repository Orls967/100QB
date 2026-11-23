# frozen_string_literal: true

require_relative 'garch_family/rng'
require_relative 'garch_family/parallel'
require_relative 'garch_family/special'
require_relative 'garch_family/linalg'
require_relative 'garch_family/innovations'
require_relative 'garch_family/vol_model'
require_relative 'garch_family/garch11'
require_relative 'garch_family/gjr_garch11'
require_relative 'garch_family/egarch11'
require_relative 'garch_family/ewma'
require_relative 'garch_family/moments_study'
require_relative 'garch_family/optim'
require_relative 'garch_family/inference'
require_relative 'garch_family/mle'
require_relative 'garch_family/innovation_study'
require_relative 'garch_family/risk_backtest'
require_relative 'garch_family/model_comparison'
require_relative 'garch_family/dcc'
require_relative 'garch_family/dcc_studies'
require_relative 'garch_family/har_rv'
require_relative 'garch_family/misspec_portfolio'
require_relative 'garch_family/json_out'
require_relative 'garch_family/svg'
require_relative 'garch_family/figures'
require_relative 'garch_family/criteria'
require_relative 'garch_family/runner'

module GarchFamily
  MAIN_SEED = 20_261_123
end
