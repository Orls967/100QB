# frozen_string_literal: true

module AsianOptions
  Params = Data.define(:s0, :k, :r, :sigma, :t, :n_steps) do
    def validate!
      raise ArgumentError, 's0 harus > 0' unless s0.positive?
      raise ArgumentError, 'k harus > 0' unless k.positive?
      raise ArgumentError, 'n_steps harus >= 1' unless n_steps >= 1
    end

    def to_h
      {
        s0: s0.to_f,
        k: k.to_f,
        r: r.to_f,
        sigma: sigma.to_f,
        t: t.to_f,
        n_steps: n_steps.to_i
      }
    end
  end

  EstimatorResult = Data.define(:price, :se, :ci_low, :ci_high) do
    def to_h
      { price: price, se: se, ci_low: ci_low, ci_high: ci_high }
    end
  end

  ControlVariateResult = Data.define(:price, :se, :ci_low, :ci_high, :beta, :y_exact) do
    def to_h
      {
        price: price,
        se: se,
        ci_low: ci_low,
        ci_high: ci_high,
        beta: beta,
        y_exact: y_exact
      }
    end
  end

  DeltaResult = Data.define(:delta, :se, :ci_low, :ci_high) do
    def to_h
      { delta: delta, se: se, ci_low: ci_low, ci_high: ci_high }
    end
  end

  ConvergencePoint = Data.define(:n, :se_naive, :se_antithetic, :se_control_variate) do
    def to_h
      {
        n: n,
        se_naive: se_naive,
        se_antithetic: se_antithetic,
        se_control_variate: se_control_variate
      }
    end
  end

  Slopes = Data.define(:naive, :antithetic, :control_variate) do
    def to_h
      { naive: naive, antithetic: antithetic, control_variate: control_variate }
    end
  end

  EstimatorSummary = Data.define(:name, :price, :se, :ci_low, :ci_high, :vrf, :beta, :y_exact) do
    def to_h
      base = {
        name: name,
        price: price,
        se: se,
        ci_low: ci_low,
        ci_high: ci_high,
        vrf: vrf
      }
      return base if beta.nil? && y_exact.nil?

      base.merge(beta: beta, y_exact: y_exact)
    end
  end
end
