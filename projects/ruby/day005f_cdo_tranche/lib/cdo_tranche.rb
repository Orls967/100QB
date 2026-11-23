# frozen_string_literal: true

require_relative 'cdo_tranche/rng'
require_relative 'cdo_tranche/parallel'
require_relative 'cdo_tranche/special'
require_relative 'cdo_tranche/models'
require_relative 'cdo_tranche/asb'
require_relative 'cdo_tranche/lhp'
require_relative 'cdo_tranche/pricing'
require_relative 'cdo_tranche/base_corr'
require_relative 'cdo_tranche/importance'
require_relative 'cdo_tranche/model_risk'
require_relative 'cdo_tranche/json_out'
require_relative 'cdo_tranche/svg'
require_relative 'cdo_tranche/criteria'
require_relative 'cdo_tranche/figures'
require_relative 'cdo_tranche/runner'

module CdoTranche
  MAIN_SEED = 20_261_123
end
