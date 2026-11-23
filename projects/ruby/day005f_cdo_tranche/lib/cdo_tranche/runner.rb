# frozen_string_literal: true

module CdoTranche
  module RunnerParts
    module_function

    def run_part1
      obs = Models.homogeneous_portfolio(n_names: 125)
      max_bin_err = Models::TRANCHES.map do |tr|
        el_a = Asb.tranche_el_asb(obs, 5.0, 0.25, tr[:attach], tr[:detach], loss_unit: 0.0048)
        el_b = Asb.homogeneous_tranche_el_exact(125, 0.015, 0.40, 5.0, 0.25, tr[:attach], tr[:detach])
        (el_a - el_b).abs
      end.max
      p5 = Models.default_prob(0.015, 5.0)
      max_lhp_err = [0.03, 0.07, 0.10, 0.15, 0.30].map do |d|
        (Lhp.equity_el(d, p5, 0.40, 0.25) - Lhp.equity_el_bivariate(d, p5, 0.40, 0.25)).abs
      end.max
      { max_asb_vs_binomial_err: max_bin_err, max_lhp_bivariate_err: max_lhp_err,
        granularity: Pricing.granularity_study }
    end

    def run_part2
      homog = Models::TRANCHES.map do |tr|
        Pricing.homogeneous_exact_pv(125, 0.015, 0.40, 0.25, tr[:attach], tr[:detach]).merge(name: tr[:name])
      end
      homog_full = Pricing.homogeneous_exact_pv(125, 0.015, 0.40, 0.25, 0.0, 1.0)[:pv_prot]
      obs_het = Models.heterogeneous_portfolio(n_names: 125)
      hetero = Models::TRANCHES.map do |tr|
        Pricing.asb_tranche_pv(obs_het, 0.25, tr[:attach], tr[:detach], loss_unit: 0.0004).merge(name: tr[:name])
      end
      het_full = Pricing.asb_tranche_pv(obs_het, 0.25, 0.0, 1.0, loss_unit: 0.0004)[:pv_prot]
      {
        homogeneous: homog, heterogeneous: hetero,
        homog_pv_sum_err: (homog.sum(0.0) { |r| r[:pv_prot] } - homog_full).abs,
        hetero_pv_sum_err: (hetero.sum(0.0) { |r| r[:pv_prot] } - het_full).abs
      }
    end

    def run_part3
      detaches = [0.03, 0.07, 0.10, 0.15, 0.30]
      market_els = detaches.map { |d| BaseCorr.t_copula_equity_el(125, 0.015, 0.40, 5.0, 0.25, 4.0, d) }
      el_fine = BaseCorr.t_copula_equity_el(125, 0.015, 0.40, 5.0, 0.25, 4.0, 0.10, m_gl: 80, w_gl: 40)
      quad_diff = (market_els[2] - el_fine).abs
      base_corrs = BaseCorr.calibrate_base_correlation(market_els, detaches)
      prev_d = 0.0
      norm_els = base_corrs.map do |row|
        w = row[:detach] - prev_d
        prev_d = row[:detach]
        row[:implied_tranche_el] / w
      end
      scans = build_compound_scans(market_els, detaches)
      {
        market_t4_equity_els: market_els, quadrature_conv_diff: quad_diff,
        base_correlations: base_corrs, normalized_tranche_els: norm_els,
        skew_sensitivity: BaseCorr.t_copula_sensitivity_grid(detaches: detaches),
        compound_scans: scans,
        mezz_compound_roots: scans[1][:roots], senior_compound_roots: scans[2][:roots]
      }
    end

    def build_compound_scans(market_els, detaches)
      prev_el = 0.0
      prev_d = 0.0
      detaches.each_with_index.map do |d, idx|
        tr_el = market_els[idx] - prev_el
        roots = BaseCorr.compound_correlation_scan(tr_el, prev_d, d)
        res = { tranche: "#{(prev_d * 100).round}-#{(d * 100).round}%", n_roots: roots.length, roots: roots }
        prev_el = market_els[idx]
        prev_d = d
        res
      end
    end
  end

  class Runner
    def initialize(seed: MAIN_SEED, workers: nil, r_is: 100_000, r_mr: 30_000)
      @seed = seed
      @workers = Parallel.resolve_workers(workers)
      @r_is = r_is
      @r_mr = r_mr
    end

    def self.find_repo_root(start_dir)
      dir = File.expand_path(start_dir)
      loop do
        return dir if File.exist?(File.join(dir, 'docs', 'ROADMAP.md'))

        parent = File.dirname(dir)
        return File.expand_path('../../..', start_dir) if parent == dir

        dir = parent
      end
    end

    def execute
      t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      p1, t1 = timed { RunnerParts.run_part1 }
      p2, t2 = timed { RunnerParts.run_part2 }
      p3, t3 = timed { RunnerParts.run_part3 }
      p4, t4 = timed { run_part4 }
      p5, t5 = timed { ModelRisk.compare_models(n_paths: @r_mr, seed: @seed, workers: @workers) }
      bit_diffs = check_worker_invariance
      tot = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      {
        metadata: build_metadata([t1, t2, t3, t4, t5], tot),
        criteria: Criteria.evaluate(p1, p2, p3, p4, p5, bit_diffs, tot),
        part1_anchors: p1, part2_valuation: p2, part3_base_corr: p3,
        part4_importance: p4, part5_model_risk: p5
      }
    end

    def write_outputs(payload, repo_root:)
      json_path = File.join(repo_root, 'outputs', 'day005f', 'results.json')
      JsonOut.write(json_path, payload)
      { json: json_path, svgs: Figures.write_all(payload, File.join(repo_root, 'outputs', 'figures')) }
    end

    private

    def timed
      ts = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      [yield, Process.clock_gettime(Process::CLOCK_MONOTONIC) - ts]
    end

    def build_metadata(ts, tot)
      {
        day: 5, session: 'f', topic_source: 'backup_cdo_tranche_lhp_asb_exponential_twisting',
        language: 'ruby', ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?,
        main_seed: @seed, workers_used: @workers, run_completed_without_errors: true,
        elapsed_seconds: { part1_anchors: ts[0], part2_valuation: ts[1], part3_base_corr: ts[2],
                           part4_importance: ts[3], part5_model_risk: ts[4], total: tot }
      }
    end

    def run_part4
      [[0.10, 0.15], [0.15, 0.30], [0.30, 1.00]].map do |a, d|
        Importance.compare_naive_vs_is(attach: a, detach: d, n_paths: @r_is, seed: @seed, workers: @workers)
      end
    end

    def check_worker_invariance
      a = Importance.compare_naive_vs_is(attach: 0.15, detach: 0.30, n_paths: 4_000, seed: @seed, workers: 1)
      b = Importance.compare_naive_vs_is(attach: 0.15, detach: 0.30, n_paths: 4_000, seed: @seed, workers: @workers)
      [a[:is_el]].pack('D') == [b[:is_el]].pack('D') ? 0 : 1
    end
  end
end
