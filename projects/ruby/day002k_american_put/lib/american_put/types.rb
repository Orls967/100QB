# frozen_string_literal: true

module AmericanPut
  Params = Data.define(:s0, :k, :r, :sigma, :t) do
    def validate!
      raise ArgumentError, 's0 harus > 0' unless s0.positive?
      raise ArgumentError, 'k harus > 0' unless k.positive?
      raise ArgumentError, 'sigma harus >= 0' if sigma.negative?
      raise ArgumentError, 't harus >= 0' if t.negative?
    end

    def with_s0(new_s0)
      Params.new(s0: new_s0.to_f, k: k, r: r, sigma: sigma, t: t)
    end
  end

  CrrResult = Data.define(
    :n_steps,
    :american_price,
    :european_price,
    :early_exercise_premium,
    :delta,
    :gamma,
    :theta
  )

  LsmResult = Data.define(
    :n_paths,
    :n_steps,
    :degree,
    :antithetic,
    :price,
    :se,
    :ci_low,
    :ci_high,
    :early_exercise_premium
  )

  PsorResult = Data.define(
    :n_x,
    :n_t,
    :american_price,
    :european_price,
    :early_exercise_premium,
    :delta,
    :gamma,
    :avg_sweeps
  )

  SpotComparisonPoint = Data.define(
    :s0,
    :intrinsic,
    :bs_european_put,
    :crr_american,
    :lsm_american,
    :lsm_se,
    :psor_american,
    :early_exercise_premium
  )
end
