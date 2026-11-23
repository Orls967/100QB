# frozen_string_literal: true

require 'minitest/autorun'
require 'mlmc_heston'

class TestCore < Minitest::Test
  HESTON_BASE = {
    spot: 100.0, rate: 0.02, v0: 0.04, kappa: 1.5, theta: 0.04, xi: 0.5, rho: -0.7
  }.freeze

  def test_heston_cf_reduces_to_bs_when_xi_zero
    p = HESTON_BASE.merge(xi: 0.0)
    c_h = MlmcHeston::HestonCf.call_price(p, 100.0, 1.0)
    c_b = MlmcHeston::Special.bs_call(100.0, 100.0, 0.02, 0.2, 1.0)
    assert_in_delta c_b, c_h, 1.0e-10
  end

  def test_heston_put_call_parity
    [85.0, 100.0, 115.0].each do |k|
      c = MlmcHeston::HestonCf.call_price(HESTON_BASE, k, 1.0)
      pt = MlmcHeston::HestonCf.put_price(HESTON_BASE, k, 1.0)
      rhs = 100.0 - (k * Math.exp(-0.02))
      assert_in_delta rhs, c - pt, 1.0e-12
    end
  end

  def test_heston_cf_quadrature_convergence
    c96 = MlmcHeston::HestonCf.call_price(HESTON_BASE, 100.0, 1.0, u_max: 140.0, n_gl: 96)
    c128 = MlmcHeston::HestonCf.call_price(HESTON_BASE, 100.0, 1.0, u_max: 180.0, n_gl: 128)
    assert_in_delta c128, c96, 1.0e-9
  end

  def test_cir_moments_match_full_truncation_mc
    rng = MlmcHeston::Rng::Sampler.new(42)
    n_paths = 15_000
    n_steps = 100
    dt = 1.0 / n_steps
    sqrt_dt = Math.sqrt(dt)
    v_sum = 0.0
    v_sum2 = 0.0
    n_paths.times do
      v = 0.04
      n_steps.times do
        vp = v.positive? ? v : 0.0
        v += (1.5 * (0.04 - vp) * dt) + (0.5 * Math.sqrt(vp) * sqrt_dt * rng.next_normal)
      end
      vp = v.positive? ? v : 0.0
      v_sum += vp
      v_sum2 += vp * vp
    end
    m = v_sum / n_paths
    var = (v_sum2 / n_paths) - (m * m)
    assert_in_delta MlmcHeston::Simulators.cir_exact_mean(0.04, 1.5, 0.04, 1.0), m, 0.002
    assert_in_delta MlmcHeston::Simulators.cir_exact_var(0.04, 1.5, 0.04, 0.5, 1.0), var, 0.0005
  end

  def test_parallel_worker_invariance_1_vs_4
    run_one = lambda do |w|
      MlmcHeston::Parallel.map_chunks(n_items: 400, workers: w, chunk_size: 100) do |cidx, len|
        rng = MlmcHeston::Rng::Sampler.new(MlmcHeston::Rng.derive_seed(99, 1, cidx))
        s = 0.0
        len.times { s += rng.next_normal }
        s
      end.sum(0.0)
    end
    assert_equal [run_one.call(1)].pack('D'), [run_one.call(4)].pack('D')
  end
end
