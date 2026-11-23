# frozen_string_literal: true

module GarchFamily
  class Ewma < VolModel::Base
    def initialize(lambda_decay: 0.94)
      super()
      @lambda = lambda_decay
    end

    def name = 'ewma'

    def persistence(_params = {}) = @lambda

    def unconditional_variance(_params = {}) = Float::INFINITY

    def forecast_variance(_params, next_var, _horizon) = next_var

    def news_impact(_params, eps_grid, base_var:)
      eps_grid.map { |e| ((1.0 - @lambda) * e * e) + (@lambda * base_var) }
    end

    def filter(params, returns)
      mu = params.fetch(:mu, 0.0)
      init_n = [20, returns.length].min
      h = returns[0, init_n].sum(0.0) { |r| (r - mu)**2 } / init_n
      vars = Array.new(returns.length, h)
      returns.each_with_index do |_r, t|
        next if t.zero?

        eps = returns[t - 1] - mu
        h = (@lambda * h) + ((1.0 - @lambda) * eps * eps)
        vars[t] = h
      end
      vars
    end
  end
end
