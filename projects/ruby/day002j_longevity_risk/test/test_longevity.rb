# frozen_string_literal: true

require 'minitest/autorun'
require_relative '../lib/longevity'

class TestLongevity < Minitest::Test
  def test_svd_reconstruction_and_orthogonality
    [[6, 4, 11], [10, 7, 29], [4, 7, 53]].each do |m, n, seed|
      rng = Longevity::Rng.new(seed)
      mat = Array.new(m) { Array.new(n) { rng.next_normal } }
      res = Longevity::Linalg.svd(mat)
      u = res[:u]
      s = res[:s]
      v = res[:v]

      s.each_cons(2) { |a, b| assert_operator a, :>=, b }
      s.each { |sv| assert_operator sv, :>=, 0.0 }

      rec = Longevity::Linalg.reconstruct_from_svd(u, s, v)
      diff = Array.new(m) { |i| Array.new(n) { |j| mat[i][j] - rec[i][j] } }
      assert_operator Longevity::Linalg.frobenius_norm(diff), :<, 1e-12

      utu = Longevity::Linalg.matmul(Longevity::Linalg.transpose(u), u)
      vtv = Longevity::Linalg.matmul(Longevity::Linalg.transpose(v), v)
      k = s.size
      eye_k = Longevity::Linalg.eye(k)
      err_u = Array.new(k) { |i| Array.new(k) { |j| utu[i][j] - eye_k[i][j] } }
      err_v = Array.new(k) { |i| Array.new(k) { |j| vtv[i][j] - eye_k[i][j] } }
      assert_operator Longevity::Linalg.frobenius_norm(err_u), :<, 1e-12
      assert_operator Longevity::Linalg.frobenius_norm(err_v), :<, 1e-12
    end
  end

  def test_lee_carter_noiseless
    truth = Longevity::LeeCarter.true_parameters(seed: 42)
    ln_m = Longevity::LeeCarter.simulate_log_mortality(truth, noise_sd: 0.0, seed: 42)
    fit = Longevity::LeeCarter.fit(ln_m)

    max_a = truth[:a].each_index.map { |i| (fit[:a][i] - truth[:a][i]).abs }.max
    max_b = truth[:b].each_index.map { |i| (fit[:b][i] - truth[:b][i]).abs }.max
    max_k = truth[:k].each_index.map { |t| (fit[:k][t] - truth[:k][t]).abs }.max

    assert_operator max_a, :<, 1e-9
    assert_operator max_b, :<, 1e-9
    assert_operator max_k, :<, 1e-9
  end

  def test_gompertz_makeham_noiseless
    truth = Longevity::LeeCarter.true_parameters(seed: 42)
    ln_m = Longevity::LeeCarter.simulate_log_mortality(truth, noise_sd: 0.0, seed: 42)
    fit = Longevity::LeeCarter.fit(ln_m)
    a_unshifted = fit[:a].each_with_index.map { |ax, i| ax - (fit[:b][i] * truth[:k_raw_mean]) }
    gm = Longevity::GompertzMakeham.fit(truth[:ages], a_unshifted)

    assert_operator ((gm[:a] - Longevity::LeeCarter::TRUE_GM_A) / Longevity::LeeCarter::TRUE_GM_A).abs, :<, 1e-6
    assert_operator ((gm[:b] - Longevity::LeeCarter::TRUE_GM_B) / Longevity::LeeCarter::TRUE_GM_B).abs, :<, 1e-6
    assert_operator ((gm[:c] - Longevity::LeeCarter::TRUE_GM_C) / Longevity::LeeCarter::TRUE_GM_C).abs, :<, 1e-6
  end

  def test_forecast_moments_and_worker_invariance
    truth = Longevity::LeeCarter.true_parameters(seed: 42)
    rw = Longevity::Forecast.estimate_rw_drift(truth[:k])
    p1 = Longevity::Forecast.simulate_k_paths(rw, horizon: 40, num_paths: 4000, base_seed: 777, workers: 1)
    p4 = Longevity::Forecast.simulate_k_paths(rw, horizon: 40, num_paths: 4000, base_seed: 777, workers: 4)
    assert_equal p1, p4

    moms = Longevity::Forecast.check_moments(p1, rw, horizons: [1, 10, 40], param_uncertainty: false)
    moms.each do |row|
      assert_operator row[:mean_z].abs, :<, 3.0
      assert_operator row[:var_z].abs, :<, 3.0
    end
  end

  def test_annuity_exact_and_simulation
    truth = Longevity::LeeCarter.true_parameters(seed: 42)
    rw = Longevity::Forecast.estimate_rw_drift(truth[:k])
    k_med = Array.new(35) { |j| rw[:k_last] + ((j + 1) * rw[:d_hat]) }
    qx = Longevity::Annuity.cohort_qx(truth[:a], truth[:b], k_med)
    exact = Longevity::Annuity.deterministic_moments(qx)

    assert_operator (exact[:mean_backward] - exact[:mean_forward]).abs, :<, 1e-12

    sim1 = Longevity::Annuity.simulate_individual_lives(qx, num_lives: 20_000, base_seed: 99, workers: 1)
    sim3 = Longevity::Annuity.simulate_individual_lives(qx, num_lives: 20_000, base_seed: 99, workers: 3)
    assert_equal sim1, sim3
    assert_operator (sim1[:mean] - exact[:mean_backward]).abs, :<, 3.0 * sim1[:mean_se]
    assert_operator (sim1[:variance] - exact[:variance]).abs, :<, 3.0 * sim1[:var_se]
  end

  def test_portfolio_total_variance_law
    truth = Longevity::LeeCarter.true_parameters(seed: 42)
    rw = Longevity::Forecast.estimate_rw_drift(truth[:k])
    k_paths = Longevity::Forecast.simulate_k_paths(rw, horizon: 35, num_paths: 800, base_seed: 555, workers: 2)
    scen = Longevity::Portfolio.scenario_moments(truth[:a], truth[:b], k_paths, workers: 2)
    ltv = Longevity::Portfolio.total_variance_summary(scen, [50])
    sim = Longevity::Portfolio.simulate_portfolio_full(
      truth[:a], truth[:b], k_paths, num_policies: 50, base_seed: 888, workers: 2
    )

    row50 = ltv[:by_n][0]
    assert_operator (sim[:mean] - row50[:mean]).abs, :<, 3.0 * sim[:mean_se]
    assert_operator (sim[:variance] - row50[:var]).abs, :<, 3.0 * sim[:var_se]
  end
end
