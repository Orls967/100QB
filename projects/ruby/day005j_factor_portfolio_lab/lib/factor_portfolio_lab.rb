# frozen_string_literal: true

require_relative 'factor_portfolio_lab/rng'
require_relative 'factor_portfolio_lab/parallel'
require_relative 'factor_portfolio_lab/distributions'
require_relative 'factor_portfolio_lab/root_finder'
require_relative 'factor_portfolio_lab/lin_alg'
require_relative 'factor_portfolio_lab/factor_world'
require_relative 'factor_portfolio_lab/pca'
require_relative 'factor_portfolio_lab/spiked_study'
require_relative 'factor_portfolio_lab/ols'
require_relative 'factor_portfolio_lab/time_series_regression'
require_relative 'factor_portfolio_lab/grs'
require_relative 'factor_portfolio_lab/grs_study'
require_relative 'factor_portfolio_lab/fama_macbeth'
require_relative 'factor_portfolio_lab/shanken'
require_relative 'factor_portfolio_lab/fm_study'
require_relative 'factor_portfolio_lab/risk_parity'
require_relative 'factor_portfolio_lab/rp_study'
require_relative 'factor_portfolio_lab/black_litterman'
require_relative 'factor_portfolio_lab/bl_study'
require_relative 'factor_portfolio_lab/constraints'
require_relative 'factor_portfolio_lab/markowitz'
require_relative 'factor_portfolio_lab/mv_study'

module FactorPortfolioLab
  MAIN_SEED = 20_261_123
end
