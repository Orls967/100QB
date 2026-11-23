# frozen_string_literal: true

module CdoTranche
  module Models
    TRANCHES = [
      { name: '0-3%', attach: 0.00, detach: 0.03 },
      { name: '3-7%', attach: 0.03, detach: 0.07 },
      { name: '7-10%', attach: 0.07, detach: 0.10 },
      { name: '10-15%', attach: 0.10, detach: 0.15 },
      { name: '15-30%', attach: 0.15, detach: 0.30 },
      { name: '30-100%', attach: 0.30, detach: 1.00 }
    ].freeze

    PAYMENT_TIMES = (1..20).map { |q| 0.25 * q }.freeze
    RISK_FREE_RATE = 0.02

    module_function

    def homogeneous_portfolio(n_names: 125, hazard: 0.015, recovery: 0.40)
      notional = 1.0 / n_names
      lgd = 1.0 - recovery
      Array.new(n_names) do |i|
        { id: i, notional: notional, hazard: hazard, recovery: recovery, lgd: notional * lgd }
      end
    end

    # Portofolio heterogen 125 nama (5 sektor dengan hazard & recovery berbeda).
    def heterogeneous_portfolio(n_names: 125)
      hazards = [0.006, 0.010, 0.015, 0.022, 0.030]
      recoveries = [0.50, 0.45, 0.40, 0.35, 0.30]
      notional = 1.0 / n_names
      Array.new(n_names) do |i|
        s = i % 5
        { id: i, notional: notional, hazard: hazards[s], recovery: recoveries[s],
          lgd: notional * (1.0 - recoveries[s]) }
      end
    end

    def default_prob(hazard, t)
      1.0 - Math.exp(-hazard * t)
    end

    def tranche_payoff(port_loss, attach, detach)
      (port_loss - attach).clamp(0.0, detach - attach)
    end
  end
end
