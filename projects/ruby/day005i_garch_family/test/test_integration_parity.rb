# frozen_string_literal: true

require 'minitest/autorun'
require 'tmpdir'
require 'garch_family'

class TestIntegrationAndParity < Minitest::Test
  def test_structural_break_causes_spurious_igarch_that_splits_resolve
    sim = GarchFamily::MisspecPortfolio.simulate_break_series(n_obs: 3000, seed: 3001)
    g = GarchFamily::Garch11.new
    nd = GarchFamily::Innovations::Normal.new
    f_full = GarchFamily::Optim.fit(g, sim[:returns], dist: nd)
    f_pre = GarchFamily::Optim.fit(g, sim[:returns][0...1500], dist: nd)
    assert_operator g.persistence(f_full[:params]), :>, g.persistence(f_pre[:params]) + 0.04
  end

  def test_student_t_lr_mixture_critical_value_beats_naive_chi2
    mix_p = GarchFamily::Special.chi2_mixture_01_sf(2.705543)
    chi_p = GarchFamily::Special.chi2_sf(2.705543, 1)
    assert_in_delta 0.05, mix_p, 1.0e-5
    assert_operator chi_p, :>, 0.09
  end

  def test_parallel_map_items_returns_exact_order_and_matches_single_worker
    r1 = GarchFamily::Parallel.map_items(n_items: 45, workers: 1, chunk_size: 7) { |i| (i * i) + 0.5 }
    r4 = GarchFamily::Parallel.map_items(n_items: 45, workers: 4, chunk_size: 7) { |i| (i * i) + 0.5 }
    assert_equal 45, r1.length
    assert_equal r1, r4
  end

  def test_json_out_sanitizes_nan_and_inf
    san = GarchFamily::JsonOut.sanitize({ a: Float::NAN, b: Float::INFINITY, c: [1.0, -Float::INFINITY] })
    assert_equal 'NaN', san[:a]
    assert_equal 'Infinity', san[:b]
    assert_equal [1.0, '-Infinity'], san[:c]
  end

  def test_svg_writes_valid_line_and_bar_charts_including_negative_bars
    Dir.mktmpdir do |dir|
      lp = File.join(dir, 'line.svg')
      bp = File.join(dir, 'bar.svg')
      GarchFamily::Svg.write_line_chart(lp, title: 'T', x_label: 'X', y_label: 'Y',
                                            series: [{ label: 'S', x: [1.0, 2.0], y: [3.0, 4.0], color: '#000' }])
      GarchFamily::Svg.write_bar_chart(bp, title: 'B', y_label: 'V',
                                           groups: [{ label: 'G', bars: [{ value: -0.5, color: '#f00' }] }])
      assert_includes File.read(lp), '</svg>'
      assert_includes File.read(bp), '<rect'
    end
  end

  def test_nic_skew_ratio_exceeds_one_only_when_gamma_nonzero
    gjr = GarchFamily::GjrGarch11.new
    p_sym = { mu: 0.0, omega: 5.0e-6, alpha: 0.08, gamma: 0.0, beta: 0.85 }
    p_asym = { mu: 0.0, omega: 5.0e-6, alpha: 0.03, gamma: 0.10, beta: 0.85 }
    s_nic = gjr.news_impact(p_sym, [-0.03, 0.03])
    a_nic = gjr.news_impact(p_asym, [-0.03, 0.03])
    assert_in_delta 1.0, s_nic[0] / s_nic[1], 1.0e-12
    assert_operator a_nic[0] / a_nic[1], :>, 1.9
  end

  def test_find_repo_root_locates_project_root
    root = GarchFamily::Runner.find_repo_root(__dir__)
    assert File.exist?(File.join(root, 'docs', 'ROADMAP.md'))
  end
end
