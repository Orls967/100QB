# frozen_string_literal: true

require 'minitest/autorun'
require 'cdo_tranche'

class TestCore < Minitest::Test
  def test_asb_matches_exact_homogeneous_binomial
    n = 125
    hazard = 0.015
    recovery = 0.40
    t = 5.0
    rho = 0.25
    obs = CdoTranche::Models.homogeneous_portfolio(n_names: n, hazard: hazard, recovery: recovery)

    CdoTranche::Models::TRANCHES.each do |tr|
      el_asb = CdoTranche::Asb.tranche_el_asb(obs, t, rho, tr[:attach], tr[:detach], loss_unit: 0.0048)
      el_bin = CdoTranche::Asb.homogeneous_tranche_el_exact(n, hazard, recovery, t, rho, tr[:attach], tr[:detach])
      assert_in_delta el_bin, el_asb, 1.0e-12, tr[:name]
    end
  end

  def test_asb_matches_bruteforce_small_heterogeneous
    obs = CdoTranche::Models.heterogeneous_portfolio(n_names: 12)
    t = 5.0
    rho = 0.30
    probs = obs.map { |ob| CdoTranche::Asb.conditional_prob_1f(CdoTranche::Models.default_prob(ob[:hazard], t), rho, 0.5) }
    dist, unit = CdoTranche::Asb.conditional_loss_distribution(obs, probs, loss_unit: 0.0005)

    # Brute-force 2^12 enumerasi
    bf_el = 0.0
    (1 << 12).times do |mask|
      pr = 1.0
      loss = 0.0
      12.times do |i|
        if mask.nobits?(1 << i)
          pr *= 1.0 - probs[i]
        else
          pr *= probs[i]
          loss += obs[i][:lgd]
        end
      end
      bf_el += pr * CdoTranche::Models.tranche_payoff(loss, 0.05, 0.15)
    end

    asb_el = 0.0
    dist.each_with_index do |pr, k|
      asb_el += pr * CdoTranche::Models.tranche_payoff(k * unit, 0.05, 0.15)
    end
    assert_in_delta bf_el, asb_el, 1.0e-4
  end

  def test_lhp_quadrature_matches_bivariate_closed_form
    p = CdoTranche::Models.default_prob(0.015, 5.0)
    [0.03, 0.07, 0.10, 0.15, 0.30].each do |d|
      el_q = CdoTranche::Lhp.equity_el(d, p, 0.40, 0.25)
      el_b = CdoTranche::Lhp.equity_el_bivariate(d, p, 0.40, 0.25)
      assert_in_delta el_b, el_q, 1.0e-9, "detach=#{d}"
    end
  end

  def test_tranche_sum_conservation
    obs = CdoTranche::Models.heterogeneous_portfolio(n_names: 125)
    t = 5.0
    rho = 0.25
    sum_tr = CdoTranche::Models::TRANCHES.sum(0.0) do |tr|
      CdoTranche::Asb.tranche_el_asb(obs, t, rho, tr[:attach], tr[:detach], loss_unit: 0.0004)
    end
    port_el = obs.sum(0.0) { |ob| CdoTranche::Models.default_prob(ob[:hazard], t) * ob[:lgd] }
    assert_in_delta port_el, sum_tr, 1.0e-4
  end

  def test_parallel_worker_invariance_1_vs_4
    run_one = lambda do |w|
      CdoTranche::Parallel.map_chunks(n_items: 400, workers: w, chunk_size: 100) do |cidx, len|
        rng = CdoTranche::Rng::Sampler.new(CdoTranche::Rng.derive_seed(123, 1, cidx))
        s = 0.0
        len.times { s += rng.next_normal }
        s
      end.sum(0.0)
    end
    assert_equal [run_one.call(1)].pack('D'), [run_one.call(4)].pack('D')
  end

  def test_gauss_legendre_matches_analytic_integrals
    nodes, weights = CdoTranche::Special.gauss_legendre(16)
    poly = nodes.each_with_index.sum(0.0) { |x, i| weights[i] * (x**4) }
    assert_in_delta 0.4, poly, 1.0e-14
    exp_int = nodes.each_with_index.sum(0.0) { |x, i| 0.5 * weights[i] * Math.exp(0.5 * (x + 1.0)) }
    assert_in_delta Math.exp(1.0) - 1.0, exp_int, 1.0e-14
  end

  def test_asb_n1_and_n2_exact_enumeration
    obs1 = [{ id: 0, notional: 1.0, hazard: 0.02, recovery: 0.4, lgd: 0.6 }]
    dist1, u1 = CdoTranche::Asb.conditional_loss_distribution(obs1, [0.1], loss_unit: 0.1)
    assert_in_delta 0.9, dist1[0], 1.0e-14
    assert_in_delta 0.1, dist1[(0.6 / u1).round], 1.0e-14

    obs2 = [
      { id: 0, notional: 0.5, hazard: 0.02, recovery: 0.4, lgd: 0.3 },
      { id: 1, notional: 0.5, hazard: 0.04, recovery: 0.5, lgd: 0.2 }
    ]
    dist2, u2 = CdoTranche::Asb.conditional_loss_distribution(obs2, [0.1, 0.2], loss_unit: 0.1)
    assert_in_delta 0.72, dist2[0], 1.0e-14
    assert_in_delta 0.18, dist2[(0.2 / u2).round], 1.0e-14
    assert_in_delta 0.08, dist2[(0.3 / u2).round], 1.0e-14
    assert_in_delta 0.02, dist2[(0.5 / u2).round], 1.0e-14
  end

  def test_lhp_convergence_as_n_grows
    gran = CdoTranche::Pricing.granularity_study(ns: [25, 125, 500])
    assert_operator gran[2][:abs_err], :<, gran[1][:abs_err]
    assert_operator gran[1][:abs_err], :<, gran[0][:abs_err]
  end
end
