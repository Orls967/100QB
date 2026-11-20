# frozen_string_literal: true

require_relative 'test_helper'

class TestAsianOptions < Minitest::Test
  BS_ANCHOR = 10.450583572185565
  GEOM_ANCHOR = 5.5655088313493559

  def setup
    @params = AsianOptions::DEFAULT_PARAMS
  end

  def test_bs_put_call_parity_and_anchor
    call = AsianOptions::BlackScholes.bs_call(@params)
    put = AsianOptions::BlackScholes.bs_put(@params)
    rhs = @params.s0 - (@params.k * Math.exp(-@params.r * @params.t))

    assert_in_delta rhs, call - put, 1e-12
    rel_err = (call - BS_ANCHOR).abs / BS_ANCHOR
    assert_operator rel_err, :<, 1e-12, "Rel error BS terlalu besar: #{rel_err}"
  end

  def test_geometric_asian_closed_form_vs_anchor_and_mc
    exact = AsianOptions::BlackScholes.geometric_asian_call(@params)
    rel_err = (exact - GEOM_ANCHOR).abs / GEOM_ANCHOR
    assert_operator rel_err, :<, 1e-12, "Rel error Kemna-Vorst terlalu besar: #{rel_err}"

    spec = AsianOptions::Estimators.build_gbm_spec(@params)
    chunks = AsianOptions::Parallel.map_chunks(n_items: 40_000, workers: 2) do |chunk_idx, len|
      simulate_geometric_chunk(spec, chunk_idx, len, 123)
    end
    mc = AsianOptions::Stats.summarize_univariate(AsianOptions::Stats.merge_univariate(chunks))
    z = AsianOptions::Stats.z_score(mc.price - exact, mc.se)
    assert_operator z, :<, 3.0, "Geometric Asian MC #{mc.price} melenceng #{z} SE dari #{exact}"
  end

  def test_control_variate_beats_naive_variance
    n_paths = 20_000
    naive = AsianOptions::Estimators.price_asian_naive(@params, n_paths: n_paths, seed: 42, workers: 2)
    cv = AsianOptions::Estimators.price_asian_control_variate(@params, n_paths: n_paths, seed: 42, workers: 2)

    assert_operator cv.se, :<, naive.se * 0.1
    vrf = AsianOptions::Stats.variance_reduction_factor(naive.se, cv.se)
    assert_operator vrf, :>, 500.0
    assert_in_delta 1.035, cv.beta, 0.05
  end

  def test_deterministic_seed_and_single_vs_multiprocess_identical
    n_paths = 12_000
    assert_equal(
      AsianOptions::Estimators.price_asian_naive(@params, n_paths: n_paths, seed: 42, workers: 1),
      AsianOptions::Estimators.price_asian_naive(@params, n_paths: n_paths, seed: 42, workers: 4)
    )
    assert_equal(
      AsianOptions::Estimators.price_asian_antithetic(@params, n_paths: n_paths, seed: 42, workers: 1),
      AsianOptions::Estimators.price_asian_antithetic(@params, n_paths: n_paths, seed: 42, workers: 4)
    )
    assert_equal(
      AsianOptions::Estimators.price_asian_control_variate(@params, n_paths: n_paths, seed: 42, workers: 1),
      AsianOptions::Estimators.price_asian_control_variate(@params, n_paths: n_paths, seed: 42, workers: 4)
    )
    assert_equal(
      AsianOptions::Greeks.delta_pathwise(@params, n_paths: n_paths, seed: 42, workers: 1),
      AsianOptions::Greeks.delta_pathwise(@params, n_paths: n_paths, seed: 42, workers: 4)
    )
  end

  def test_delta_pathwise_matches_bump_within_se
    n_paths = 20_000
    pw = AsianOptions::Greeks.delta_pathwise(@params, n_paths: n_paths, seed: 42, workers: 2)
    bmp = AsianOptions::Greeks.delta_bump(@params, n_paths: n_paths, h_rel: 0.01, seed: 42, workers: 2)

    z = AsianOptions::Stats.z_score(pw.delta - bmp.delta, pw.se, bmp.se)
    assert_operator z, :<, 3.0, "Selisih delta pathwise vs bump terlalu besar (#{z} SE)"
    assert_in_delta 0.59, pw.delta, 0.02
  end

  private

  def simulate_geometric_chunk(spec, chunk_idx, len, seed)
    rng = AsianOptions::Rng::NormalGenerator.new(AsianOptions::Rng.derive_chunk_seed(seed, chunk_idx))
    normals = Array.new(spec.n_steps, 0.0)
    mean = 0.0
    m2 = 0.0

    len.times do |idx|
      rng.fill_normals!(normals)
      val = geometric_payoff(normals, spec)
      delta = val - mean
      mean += delta / (idx + 1)
      m2 += delta * (val - mean)
    end
    AsianOptions::Stats::UnivariateChunk.new(sample_count: len, mean: mean, m2: m2)
  end

  def geometric_payoff(normals, spec)
    log_s = spec.log_s0
    sum_log_s = 0.0
    spec.n_steps.times do |step|
      log_s += spec.drift + (spec.vol * normals[step])
      sum_log_s += log_s
    end
    geom_avg = Math.exp(sum_log_s * spec.inv_n)
    geom_avg > spec.k ? (geom_avg - spec.k) * spec.disc : 0.0
  end
end
