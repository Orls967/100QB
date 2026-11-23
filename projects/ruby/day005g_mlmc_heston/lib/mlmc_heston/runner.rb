# frozen_string_literal: true

module MlmcHeston
  class Runner
    GBM_BASE = { spot: 100.0, rate: 0.05, sigma: 0.20, strike: 100.0, maturity: 1.0 }.freeze
    HESTON_BASE = {
      spot: 100.0, rate: 0.02, v0: 0.04, kappa: 1.5, theta: 0.04, xi: 0.50, rho: -0.7,
      strike: 100.0, maturity: 1.0
    }.freeze

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
      p2, t2 = timed { run_part2_rates }
      p3, t3 = timed { run_part3_complexity(p2) }
      p4, t4 = timed { run_part4_antithetic }
      bit_diffs = check_worker_invariance
      tot = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      {
        metadata: build_metadata([t1, t2, t3, t4], tot),
        criteria: Criteria.evaluate(p1, p2, p3, p4, bit_diffs, tot),
        part1_anchors: p1, part2_rates: p2, part3_complexity: p3, part4_antithetic: p4
      }
    end

    def write_outputs(payload, repo_root:)
      json_path = File.join(repo_root, 'outputs', 'day005g', 'results.json')
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
        day: 5, session: 'g', topic_source: 'backup_mlmc_giles_cir_heston',
        language: 'ruby', ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?,
        main_seed: @seed, workers_used: @workers, run_completed_without_errors: true,
        elapsed_seconds: { part1_anchors: ts[0], part2_rates: ts[1], part3_complexity: ts[2],
                           part4_antithetic: ts[3], total: tot }
      }
    end

    def run_part1_anchors
      bs_ex = Special.bs_call(100.0, 100.0, 0.05, 0.20, 1.0)
      hc_ex = HestonCf.call_price(HESTON_BASE, 100.0, 1.0)
      hc_xi0 = HestonCf.call_price(HESTON_BASE.merge(xi: 0.0), 100.0, 1.0)
      bs_xi0 = Special.bs_call(100.0, 100.0, 0.02, 0.20, 1.0)
      hc_put = HestonCf.put_price(HESTON_BASE, 100.0, 1.0)
      hc_128 = HestonCf.call_price(HESTON_BASE, 100.0, 1.0, u_max: 180.0, n_gl: 128)
      {
        anchors: {
          bs_exact: bs_ex, heston_cf_exact: hc_ex,
          xi_zero_bs_err: (hc_xi0 - bs_xi0).abs,
          parity_err: ((hc_ex - hc_put) - (100.0 - (100.0 * Math.exp(-0.02)))).abs,
          cf_quad_err: (hc_ex - hc_128).abs
        }
      }
    end

    def run_part2_rates
      {
        gbm_euler: Mlmc.estimate_rates(model: :gbm, params: GBM_BASE.merge(scheme: :euler),
                                       seed: @seed, workers: @workers),
        gbm_milstein: Mlmc.estimate_rates(model: :gbm, params: GBM_BASE.merge(scheme: :milstein),
                                          seed: @seed, workers: @workers),
        heston_feller_ok: Mlmc.estimate_rates(model: :heston, params: HESTON_BASE.merge(xi: 0.25, scheme: :euler),
                                              seed: @seed, workers: @workers),
        heston_feller_violated: Mlmc.estimate_rates(model: :heston, params: HESTON_BASE.merge(scheme: :euler),
                                                    seed: @seed, workers: @workers)
      }
    end

    def run_part3_complexity(p2)
      {
        gbm_milstein: Complexity.compare_mlmc_vs_standard(
          model: :gbm, params: GBM_BASE.merge(scheme: :milstein), seed: @seed, workers: @workers,
          rate_est: p2[:gbm_milstein]
        ),
        heston_european: Complexity.compare_mlmc_vs_standard(
          model: :heston, params: HESTON_BASE.merge(scheme: :euler), seed: @seed, workers: @workers,
          rate_est: p2[:heston_feller_violated]
        ),
        heston_asian: Complexity.compare_mlmc_vs_standard(
          model: :heston, params: HESTON_BASE.merge(scheme: :euler, payoff_type: :asian),
          seed: @seed, workers: @workers
        )
      }
    end

    def run_part4_antithetic
      [1, 2, 3, 4].map do |lvl|
        Antithetic.heston_antithetic_var(HESTON_BASE, 100.0, 1.0, level: lvl, n_paths: 12_000,
                                                                  seed: @seed, workers: @workers)
      end
    end

    def check_worker_invariance
      a = Mlmc.run_adaptive(model: :heston, params: HESTON_BASE, eps: 0.05, seed: @seed, workers: 1)
      b = Mlmc.run_adaptive(model: :heston, params: HESTON_BASE, eps: 0.05, seed: @seed, workers: @workers)
      [a[:price]].pack('D') == [b[:price]].pack('D') ? 0 : 1
    end
  end
end
