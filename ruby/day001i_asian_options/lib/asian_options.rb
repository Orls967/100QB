# frozen_string_literal: true

require_relative 'asian_options/types'
require_relative 'asian_options/rng'
require_relative 'asian_options/stats'
require_relative 'asian_options/black_scholes'
require_relative 'asian_options/parallel'
require_relative 'asian_options/path_kernels'
require_relative 'asian_options/chunk_simulator'
require_relative 'asian_options/estimators'
require_relative 'asian_options/greeks'
require_relative 'asian_options/convergence_svg'
require_relative 'asian_options/svg'
require_relative 'asian_options/runner'

module AsianOptions
  DEFAULT_PARAMS = Params.new(
    s0: 100.0,
    k: 100.0,
    r: 0.05,
    sigma: 0.20,
    t: 1.0,
    n_steps: 252
  )

  DEFAULT_N_PATHS = 200_000
  DEFAULT_SEED = 42
  CONVERGENCE_GRID = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000].freeze
end
