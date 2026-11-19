# frozen_string_literal: true

require_relative 'test_helper'
require 'tmpdir'

class TestStatsAndRng < Minitest::Test
  def test_normal_rng_statistics_and_determinism
    gen_a = AsianOptions::Rng::NormalGenerator.new(99)
    gen_b = AsianOptions::Rng::NormalGenerator.new(99)
    n = 100_000
    buf_a = Array.new(n, 0.0)
    gen_a.fill_normals!(buf_a)

    20.times do |idx|
      assert_equal buf_a[idx], gen_b.next_normal
    end

    mean = buf_a.sum / n
    var = buf_a.sum { |x| (x - mean)**2 } / (n - 1)
    assert_in_delta 0.0, mean, 0.02
    assert_in_delta 1.0, var, 0.02
  end

  def test_loglog_slope_and_svg_generation
    ns = [1_000, 10_000, 100_000]
    ses = ns.map { |n| 10.0 / Math.sqrt(n) }
    slope = AsianOptions::Stats.fit_loglog_slope(ns, ses)
    assert_in_delta(-0.5, slope, 1e-12)

    pts = [
      AsianOptions::ConvergencePoint.new(n: 1_000, se_naive: 0.25, se_antithetic: 0.18, se_control_variate: 0.006),
      AsianOptions::ConvergencePoint.new(n: 10_000, se_naive: 0.08, se_antithetic: 0.056, se_control_variate: 0.002)
    ]
    slopes = AsianOptions::Slopes.new(naive: -0.5, antithetic: -0.5, control_variate: -0.49)
    ests = [
      AsianOptions::EstimatorSummary.new(
        name: 'Naive MC (200k)', price: 5.78, se: 0.018,
        ci_low: 5.74, ci_high: 5.82, vrf: 1.0, beta: nil, y_exact: nil
      )
    ]

    Dir.mktmpdir do |dir|
      p1 = File.join(dir, 'conv.svg')
      p2 = File.join(dir, 'est.svg')
      AsianOptions::Svg.write_convergence_svg(p1, pts, slopes)
      AsianOptions::Svg.write_estimators_svg(p2, ests)
      assert_includes File.read(p1), '<svg'
      assert_includes File.read(p2), 'VRF 1.0x'
    end
  end
end
