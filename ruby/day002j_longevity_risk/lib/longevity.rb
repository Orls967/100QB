# frozen_string_literal: true

RubyVM::YJIT.enable if defined?(RubyVM::YJIT) && RubyVM::YJIT.respond_to?(:enable)

require_relative 'longevity/rng'
require_relative 'longevity/parallel'
require_relative 'longevity/stats'
require_relative 'longevity/linalg'
require_relative 'longevity/optimizers'
require_relative 'longevity/lee_carter'
require_relative 'longevity/gompertz_makeham'
require_relative 'longevity/forecast'
require_relative 'longevity/annuity'
require_relative 'longevity/portfolio'
require_relative 'longevity/svg'
