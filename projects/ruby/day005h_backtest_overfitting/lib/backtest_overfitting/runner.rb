# frozen_string_literal: true

module BacktestOverfitting
  class Runner
    def initialize(seed: MAIN_SEED, workers: nil)
      @seed = seed
      @workers = Parallel.resolve_workers(workers)
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
      p1, t1 = timed { run_part1_anchors }
      p2, t2 = timed { SharpeStats.run_mc_study(seed: @seed, workers: @workers) }
      p3, t3 = timed { SelectionBias.run_study(seed: @seed, workers: @workers) }
      p4, t4 = timed { MultipleTesting.run_study(seed: @seed, workers: @workers) }
      p5, t5 = timed { CscvPbo.run_study(seed: @seed, workers: @workers) }
      p6, t6 = timed { PurgedCv.run_study(seed: @seed, workers: @workers) }
      bit_diffs = check_worker_invariance
      tot = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      parts = {
        part1_anchors: p1, part2_sharpe_stats: p2, part3_selection_bias: p3,
        part4_multiple_testing: p4, part5_cscv_pbo: p5, part6_purged_cv: p6
      }
      { metadata: build_metadata([t1, t2, t3, t4, t5, t6], tot),
        criteria: Criteria.evaluate(parts, bit_diffs, tot) }.merge(parts)
    end

    def write_outputs(payload, repo_root:)
      json_path = File.join(repo_root, 'outputs', 'day005h', 'results.json')
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
        day: 5, session: 'h', topic_source: 'roadmap_091_backtest_overfitting',
        language: 'ruby', ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?,
        main_seed: @seed, workers_used: @workers, run_completed_without_errors: true,
        elapsed_seconds: { part1_anchors: ts[0], part2_sharpe_stats: ts[1], part3_selection_bias: ts[2],
                           part4_multiple_testing: ts[3], part5_cscv_pbo: ts[4], part6_purged_cv: ts[5],
                           total: tot }
      }
    end

    def run_part1_anchors
      rng = Rng::Sampler.new(Rng.derive_seed(@seed, 10, 0))
      rets = DataGenerator.ar1_gaussian(200_000, phi: 0.1, sigma: 0.01, rng: rng)
      strat = Strategy::Momentum.new(lookback: 1)
      pos = strat.generate_positions(rets)
      r0 = BacktestEngine.run(rets, pos, cost_bps: 0.0)
      emp_ar1 = r0[:net_returns][1..].sum(0.0) / (rets.length - 1)
      theo_ar1 = 0.1 * 0.01 * Math.sqrt(2.0 / Math::PI)
      gbm_anchor = check_gbm_bh_anchor(rng)
      {
        equity_identity_err: BacktestEngine.verify_equity_identity(r0),
        zero_cost_err: r0[:net_returns].each_with_index.map { |v, i| (v - (pos[i] * rets[i])).abs }.max,
        no_lookahead_passed: BacktestEngine.verify_no_lookahead(Strategy.build_library.first, rets[0, 300], 150, rng),
        ar1_sign_emp: emp_ar1, ar1_sign_theo: theo_ar1,
        ar1_sign_z: (emp_ar1 - theo_ar1).abs / (0.01 / Math.sqrt(rets.length - 1)),
        gbm_bh_log_emp: gbm_anchor[:emp], gbm_bh_log_theo: gbm_anchor[:theo], gbm_bh_log_z: gbm_anchor[:z]
      }
    end

    def check_gbm_bh_anchor(rng)
      n_reps = 4_000
      n_steps = 252
      dt = 1.0 / n_steps
      pos = Array.new(n_steps, 1.0)
      logs = Array.new(n_reps) do
        r = DataGenerator.gbm_simple_returns(n_steps, mu: 0.08, sigma: 0.20, dt: dt, rng: rng)
        Math.log(BacktestEngine.run(r, pos, cost_bps: 0.0)[:final_equity])
      end
      emp = logs.sum(0.0) / n_reps
      theo = 0.08 - (0.5 * 0.04)
      { emp: emp, theo: theo, z: (emp - theo).abs / (0.20 / Math.sqrt(n_reps)) }
    end

    def check_worker_invariance
      a = SharpeStats.gaussian_coverage_mc(n_obs: 250, sr_true: 0.1, n_reps: 1_000, seed: @seed, workers: 1)
      b = SharpeStats.gaussian_coverage_mc(n_obs: 250, sr_true: 0.1, n_reps: 1_000, seed: @seed, workers: @workers)
      [a[:emp_sd]].pack('D') == [b[:emp_sd]].pack('D') ? 0 : 1
    end
  end
end
