# frozen_string_literal: true

require 'pathname'

module GarchFamily
  class Runner
    def self.find_repo_root(start_dir = Dir.pwd)
      cur = Pathname.new(File.expand_path(start_dir))
      loop do
        return cur.to_s if (cur / 'docs' / 'ROADMAP.md').exist?
        break if cur.parent == cur

        cur = cur.parent
      end
      File.expand_path('../../..', __dir__)
    end

    def initialize(seed: MAIN_SEED, workers: nil)
      @seed = seed
      @workers = Parallel.resolve_workers(workers)
    end

    def execute
      t0 = now_sec
      parts, timings = run_all_parts
      parity_diffs = check_worker_parity
      timings[:total] = now_sec - t0
      payload = parts.merge(metadata: build_meta(timings, parity_diffs))
      payload[:criteria] = Criteria.evaluate_all(payload)
      payload
    end

    def write_outputs(payload, repo_root: self.class.find_repo_root)
      out_json = File.join(repo_root, 'outputs', 'day005i', 'results.json')
      fig_dir = File.join(repo_root, 'outputs', 'figures')
      JsonOut.write(out_json, payload)
      Figures.generate_all(payload, out_dir: fig_dir)
      { json: out_json, figures_dir: fig_dir }
    end

    private

    def run_all_parts
      t = {}
      parts = {
        part1_moments: timed(t, :part1_moments) { MomentsStudy.run_study(seed: @seed, workers: @workers) },
        part2_mle: timed(t, :part2_mle) { Mle.run_study(seed: @seed, workers: @workers) },
        part3_innovations: timed(t, :part3_innovations) { InnovationStudy.run_study(seed: @seed, workers: @workers) },
        part4_backtest: timed(t, :part4_backtest) { RiskBacktest.run_study(seed: @seed, workers: @workers) },
        part5_comparison: timed(t, :part5_comparison) { ModelComparison.run_study(seed: @seed, workers: @workers) },
        part6_dcc: timed(t, :part6_dcc) { Dcc.run_study(seed: @seed, workers: @workers) },
        part7_har_rv: timed(t, :part7_har_rv) { HarRv.run_study(seed: @seed) },
        part8_structural_break: timed(t, :part8_break) do
          MisspecPortfolio.run_structural_break_study(seed: @seed, workers: @workers)
        end,
        part9_portfolio_stress: timed(t, :part9_port) { MisspecPortfolio.run_five_asset_portfolio(seed: @seed) }
      }
      [parts, t]
    end

    def timed(timings, key)
      t0 = now_sec
      res = yield
      timings[key] = now_sec - t0
      res
    end

    def check_worker_parity
      seq = Mle.run_qmle_cell(n_obs: 500, n_reps: 8, true_dist: Innovations::Normal.new,
                              stream_id: 999, seed: @seed, workers: 1)
      par = Mle.run_qmle_cell(n_obs: 500, n_reps: 8, true_dist: Innovations::Normal.new,
                              stream_id: 999, seed: @seed, workers: @workers)
      s_vec = [seq[:params][:alpha][:mean], seq[:params][:beta][:rmse]]
      p_vec = [par[:params][:alpha][:mean], par[:params][:beta][:rmse]]
      s_vec.pack('D*') == p_vec.pack('D*') ? 0 : 1
    end

    def build_meta(timings, parity_diffs)
      {
        day: 5, session: 'i', topic_source: 'garch_family_complete_pipeline',
        language: 'ruby', ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) ? RubyVM::YJIT.enabled? : false,
        main_seed: @seed, workers_used: @workers,
        worker_parity_bit_diffs: parity_diffs,
        run_completed_without_errors: true,
        elapsed_seconds: timings
      }
    end

    def now_sec = Process.clock_gettime(Process::CLOCK_MONOTONIC)
  end
end
