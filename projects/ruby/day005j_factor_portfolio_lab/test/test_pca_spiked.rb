# frozen_string_literal: true

require_relative 'test_helper'

class TestPcaSpiked < Minitest::Test
  include FactorPortfolioLab

  def test_factor_world_exact_covariance_structure
    world = FactorWorld.build_world(n_assets: 8, n_factors: 2, seed: 101)
    assert_equal 8, world[:sigma].length
    8.times do |i|
      8.times { |j| assert_in_delta world[:sigma][i][j], world[:sigma][j][i], 1e-12 }
      assert_operator world[:sigma][i][i], :>, world[:idio_var][i]
    end
  end

  def test_factor_world_student_t_variance_preserved
    world = FactorWorld.build_world(n_assets: 5, n_factors: 2, seed: 202)
    sim = FactorWorld.simulate(world, t_len: 5000, seed: 203, df: 5.0)
    cov_hat = LinAlg.cov_matrix(sim[:returns])
    5.times { |i| assert_in_delta world[:sigma][i][i], cov_hat[i][i], 0.0015 }
  end

  def test_pca_trace_and_effective_rank
    world = FactorWorld.build_world(n_assets: 12, n_factors: 3, seed: 303)
    sim = FactorWorld.simulate(world, t_len: 300, seed: 304)
    fit = Pca.fit(sim[:returns])
    assert_in_delta LinAlg.trace(fit[:cov]), fit[:eigenvalues].sum, 1e-9
    eff = Pca.effective_rank(fit[:eigenvalues])
    assert_operator eff, :>=, 1.0
    assert_operator eff, :<=, 12.0
  end

  def test_marchenko_pastur_bounds_and_bbp_threshold
    mp = Pca.marchenko_pastur_bounds(0.25, sigma2: 1.0)
    assert_in_delta 0.25, mp[:lambda_minus], 1e-12
    assert_in_delta 2.25, mp[:lambda_plus], 1e-12
    assert_in_delta 0.5, mp[:bbp_threshold], 1e-12

    sub = Pca.bbp_prediction(0.3, 0.25)
    sup = Pca.bbp_prediction(2.0, 0.25)
    refute sub[:above_threshold]
    assert_equal 0.0, sub[:alignment_sq]
    assert sup[:above_threshold]
    assert_in_delta 3.375, sup[:sample_eigenvalue], 1e-12
    assert_in_delta(3.75 / 4.5, sup[:alignment_sq], 1e-12)
  end

  def test_bbp_super_threshold_spike_accuracy
    res = SpikedStudy.run(master_seed: MAIN_SEED, n_assets: 30, t_len: 120, n_reps: 80, spikes: [2.5])
    lvl = res[:levels][0]
    assert_in_delta lvl[:theory_lambda1], lvl[:emp_lambda1_mean], 0.15
    assert_in_delta lvl[:theory_alignment_sq], lvl[:emp_alignment_sq_mean], 0.05
  end

  def test_bbp_sub_threshold_spike_stays_at_edge
    res = SpikedStudy.run(master_seed: MAIN_SEED, n_assets: 30, t_len: 120, n_reps: 80, spikes: [0.2])
    lvl = res[:levels][0]
    assert_in_delta lvl[:theory_lambda1], lvl[:emp_lambda1_mean], 0.15
    assert_operator lvl[:emp_alignment_sq_mean], :<, 0.15
  end

  def test_spiked_study_multiple_levels
    res = SpikedStudy.run(master_seed: MAIN_SEED, n_assets: 20, t_len: 80, n_reps: 20, spikes: [0.0, 1.5])
    assert_equal 2, res[:levels].length
    assert_operator res[:levels][1][:emp_lambda1_mean], :>, res[:levels][0][:emp_lambda1_mean]
  end

  def test_pca_centering_invariant_to_constant_shift
    world = FactorWorld.build_world(n_assets: 6, n_factors: 2, seed: 305)
    sim = FactorWorld.simulate(world, t_len: 80, seed: 306)
    shifted = sim[:returns].map { |r| r.each_with_index.map { |v, j| v + (0.1 * (j + 1)) } }
    fit1 = Pca.fit(sim[:returns], center: true)
    fit2 = Pca.fit(shifted, center: true)
    6.times { |i| assert_in_delta fit1[:eigenvalues][i], fit2[:eigenvalues][i], 1e-12 }
  end
end
