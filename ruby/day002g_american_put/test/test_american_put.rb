# frozen_string_literal: true

require_relative 'test_helper'

class TestAmericanPut < Minitest::Test
  def setup
    @params = AmericanPut::DEFAULT_PARAMS
  end

  def test_crr_american_vs_european_and_greeks
    res = AmericanPut::Crr.price(@params, n_steps: 500)
    bs_eu = AmericanPut::Stats.bs_put(@params)

    assert_in_delta bs_eu, res.european_price, 0.01
    assert_operator res.american_price, :>, res.european_price
    assert_in_delta 6.090, res.american_price, 0.01
    assert_operator res.delta, :<, 0.0
    assert_operator res.delta, :>, -1.0
    assert_operator res.gamma, :>, 0.0
  end

  def test_psor_matches_crr_and_bs
    crr = AmericanPut::Crr.price(@params, n_steps: 1_000)
    psor = AmericanPut::Psor.price(@params, n_x: 400, n_t: 400)
    bs_eu = AmericanPut::Stats.bs_put(@params)

    assert_in_delta bs_eu, psor.european_price, 0.005
    assert_in_delta crr.american_price, psor.american_price, 0.005
    assert_in_delta crr.delta, psor.delta, 0.01
    assert_in_delta crr.gamma, psor.gamma, 0.005
  end

  def test_lsm_deterministic_across_workers_and_consistent
    res1 = AmericanPut::Lsm.price(@params, n_paths: 8_000, n_steps: 25, seed: 42, workers: 1)
    res4 = AmericanPut::Lsm.price(@params, n_paths: 8_000, n_steps: 25, seed: 42, workers: 4)

    assert_in_delta res1.price, res4.price, 1e-15
    assert_in_delta res1.se, res4.se, 1e-15

    crr = AmericanPut::Crr.price(@params, n_steps: 500)
    z = AmericanPut::Stats.z_score(res4.price - crr.american_price, res4.se)
    assert_operator z, :<, 3.0
  end
end
