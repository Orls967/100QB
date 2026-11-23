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

module GarchFamily
  MAIN_SEED = 20_261_123
end
