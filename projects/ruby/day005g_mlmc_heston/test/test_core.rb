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

  def test_splitmix64_seed_determinism_and_avalanche
    s1 = MlmcHeston::Rng.derive_seed(20_261_123, 100, 0)
    s2 = MlmcHeston::Rng.derive_seed(20_261_123, 100, 1)
    assert_equal s1, MlmcHeston::Rng.derive_seed(20_261_123, 100, 0)
    refute_equal s1, s2
    assert_operator (s1 ^ s2).to_s(2).count('1'), :>=, 15
  end

  def test_box_muller_sampler_moments
    rng = MlmcHeston::Rng::Sampler.new(12_345)
    draws = Array.new(20_000) { rng.next_normal }
    m = draws.sum(0.0) / draws.length
    var = (draws.sum(0.0) { |x| (x - m)**2 } / draws.length)
    assert_in_delta 0.0, m, 0.03
    assert_in_delta 1.0, var, 0.04
  end

  def test_gauss_legendre_integrates_polynomial_degree_7_exactly
    nodes, weights = MlmcHeston::Special.gauss_legendre(8)
    integral = nodes.zip(weights).sum(0.0) { |x, w| w * ((x**6) + (2.0 * (x**4)) + (3.0 * (x**3)) + 1.0) }
    exact = (2.0 / 7.0) + (4.0 / 5.0) + 2.0
    assert_in_delta exact, integral, 1.0e-13
  end

  def test_complex_helpers_identities
    z = Complex(-0.7, 1.3)
    ez = MlmcHeston::Special.cexp(z)
    lz = MlmcHeston::Special.clog(ez)
    assert_operator (lz - z).abs, :<, 1.0e-14
    sz = MlmcHeston::Special.csqrt(z)
    assert_operator ((sz * sz) - z).abs, :<, 1.0e-14
  end

  def test_heston_cf_monotone_in_strike_and_vol
    c90 = MlmcHeston::HestonCf.call_price(HESTON_BASE, 90.0, 1.0)
    c100 = MlmcHeston::HestonCf.call_price(HESTON_BASE, 100.0, 1.0)
    c110 = MlmcHeston::HestonCf.call_price(HESTON_BASE, 110.0, 1.0)
    assert_operator c90, :>, c100
    assert_operator c100, :>, c110
    c_hi_v = MlmcHeston::HestonCf.call_price(HESTON_BASE.merge(v0: 0.09, theta: 0.09), 100.0, 1.0)
    assert_operator c_hi_v, :>, c100
  end

  def test_coarse_fine_coupling_uses_exact_summed_brownian_increments
    zs = Array.new(16) { |i| ((i % 5) - 2) * 0.45 }
    rec = Struct.new(:arr, :idx) do
      def next_normal
        v = arr[idx]
        self.idx += 1
        v
      end
    end
    params = { spot: 100.0, rate: 0.05, sigma: 0.25 }
    pf, y = MlmcHeston::Simulators.gbm_coupled_pair(params, 100.0, 1.0, 2, scheme: :milstein,
                                                                           rng: rec.new(zs, 0), n0: 4)
    pc_from_pair = pf - y
    coarse_zs = Array.new(8) { |j| (zs[2 * j] + zs[(2 * j) + 1]) / Math.sqrt(2.0) }
    pc_direct, = MlmcHeston::Simulators.gbm_coupled_pair(params, 100.0, 1.0, 1, scheme: :milstein,
                                                                                rng: rec.new(coarse_zs, 0), n0: 4)
    assert_in_delta pc_direct, pc_from_pair, 1.0e-12
  end
end
