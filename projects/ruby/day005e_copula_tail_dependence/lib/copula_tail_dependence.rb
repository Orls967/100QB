# frozen_string_literal: true

require_relative 'copula_tail_dependence/rng'
require_relative 'copula_tail_dependence/parallel'
require_relative 'copula_tail_dependence/special'
require_relative 'copula_tail_dependence/linalg'
require_relative 'copula_tail_dependence/stats'
require_relative 'copula_tail_dependence/copulas'
require_relative 'copula_tail_dependence/dependence'
require_relative 'copula_tail_dependence/euler'
require_relative 'copula_tail_dependence/fit'
require_relative 'copula_tail_dependence/importance'

module CopulaTailDependence
  MAIN_SEED = 20_261_123
end
