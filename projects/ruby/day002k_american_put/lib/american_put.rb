# frozen_string_literal: true

require_relative 'american_put/types'
require_relative 'american_put/rng'
require_relative 'american_put/stats'
require_relative 'american_put/linear_solver'
require_relative 'american_put/parallel'
require_relative 'american_put/crr'
require_relative 'american_put/lsm_simulator'
require_relative 'american_put/lsm'
require_relative 'american_put/psor_grid'
require_relative 'american_put/psor'
require_relative 'american_put/svg'
require_relative 'american_put/runner'

module AmericanPut
  DEFAULT_PARAMS = Params.new(
    s0: 100.0,
    k: 100.0,
    r: 0.05,
    sigma: 0.20,
    t: 1.0
  )

  DEFAULT_SEED = 42
  DEFAULT_CRR_STEPS = 4_000
  DEFAULT_LSM_PATHS = 100_000
  DEFAULT_LSM_STEPS = 50
  DEFAULT_PSOR_NX = 400
  DEFAULT_PSOR_NT = 400

  CRR_STEPS_GRID = [50, 100, 200, 500, 1_000, 2_000, 4_000].freeze
  LSM_PATHS_GRID = [5_000, 10_000, 25_000, 50_000, 100_000].freeze
  PSOR_GRID = [[100, 100], [200, 200], [400, 400], [800, 800]].freeze
  SPOT_GRID = [80.0, 90.0, 100.0, 110.0, 120.0].freeze
end
