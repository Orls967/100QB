# frozen_string_literal: true

require 'fileutils'
require 'json'

module AsianOptions
  class Runner
    def initialize(params: DEFAULT_PARAMS, n_paths: DEFAULT_N_PATHS, seed: DEFAULT_SEED, workers: nil)
      @params = params
      @n_paths = n_paths
      @seed = seed
      @workers = workers
    end

    def execute
      core = run_core_benchmark
      estimators = build_estimator_summaries(core[:naive], core[:anti], core[:cv])
      convergence, t_conv = timed { run_convergence }
      slopes = fit_all_slopes(convergence)
      payload = {
        params: @params.to_h, n_paths: @n_paths, seed: @seed, sanity_check: core[:sanity],
        estimators: estimators.map(&:to_h), delta: build_delta_section(core[:pw], core[:bmp]),
        convergence: convergence.map(&:to_h), slopes: slopes.to_h,
        timings_ms: core[:timings].merge(convergence: t_conv.round(3))
      }
      { payload: payload, estimators: estimators, convergence: convergence, slopes: slopes }
    end

    def write_outputs(result, repo_root:)
      out_dir = File.join(repo_root, 'outputs', 'day001i')
      fig_dir = File.join(repo_root, 'outputs', 'figures')
      FileUtils.mkdir_p([out_dir, fig_dir])
      paths = {
        json: File.join(out_dir, 'results.json'),
        svg_convergence: File.join(fig_dir, 'day001i_asian-convergence-loglog.svg'),
        svg_estimators: File.join(fig_dir, 'day001i_asian-estimators-ci.svg')
      }
      File.write(paths[:json], "#{JSON.pretty_generate(result[:payload])}\n")
      Svg.write_convergence_svg(paths[:svg_convergence], result[:convergence], result[:slopes])
      Svg.write_estimators_svg(paths[:svg_estimators], result[:estimators])
      paths
    end

    def self.find_repo_root(start_dir = Dir.pwd)
      dir = File.expand_path(start_dir)
      until File.directory?(File.join(dir, 'src', 'qb100'))
        raise 'Tidak menemukan root repository 100QB (src/qb100)' if dir == File.dirname(dir)

        dir = File.dirname(dir)
      end
      dir
    end

    private

    def run_core_benchmark
      opts = { n_paths: @n_paths, seed: @seed, workers: @workers }
      core_start = monotonic_ms
      sanity, t_sanity = timed { run_sanity_check }
      naive, t_naive = timed { Estimators.price_asian_naive(@params, **opts) }
      anti, t_anti = timed { Estimators.price_asian_antithetic(@params, **opts) }
      cv_res, t_cv = timed { Estimators.price_asian_control_variate(@params, **opts) }
      pw, t_pw = timed { Greeks.delta_pathwise(@params, **opts) }
      bmp, t_bmp = timed { Greeks.delta_bump(@params, h_rel: 0.01, **opts) }
      timings = {
        sanity_check: t_sanity, naive: t_naive, antithetic: t_anti, control_variate: t_cv,
        delta_pathwise: t_pw, delta_bump: t_bmp, total_200k_core: monotonic_ms - core_start
      }.transform_values { |ms| ms.round(3) }
      { sanity: sanity, naive: naive, anti: anti, cv: cv_res, pw: pw, bmp: bmp, timings: timings }
    end

    def monotonic_ms = Process.clock_gettime(Process::CLOCK_MONOTONIC) * 1000.0

    def timed
      t0 = monotonic_ms
      [yield, monotonic_ms - t0]
    end

    def run_sanity_check
      bs_euro = BlackScholes.bs_call(@params)
      mc_euro = Estimators.price_european_mc(@params, n_paths: @n_paths, seed: @seed, workers: @workers)
      diff = mc_euro.price - bs_euro
      { bs_call: bs_euro, mc_call: mc_euro.to_h, diff: diff, z_score: Stats.z_score(diff, mc_euro.se) }
    end

    def build_estimator_summaries(naive, anti, cv_res)
      [
        to_summary('Naive MC (200k)', naive, 1.0),
        to_summary('Antithetic (100k pasang)', anti, Stats.variance_reduction_factor(naive.se, anti.se)),
        to_summary('Control Variate (Geom. Asia)', cv_res, Stats.variance_reduction_factor(naive.se, cv_res.se),
                   beta: cv_res.beta, y_exact: cv_res.y_exact)
      ]
    end

    def to_summary(name, res, vrf, beta: nil, y_exact: nil)
      EstimatorSummary.new(name: name, price: res.price, se: res.se, ci_low: res.ci_low,
                           ci_high: res.ci_high, vrf: vrf, beta: beta, y_exact: y_exact)
    end

    def build_delta_section(pathwise, bump)
      diff = pathwise.delta - bump.delta
      { pathwise: pathwise.to_h, bump_1pct: bump.to_h, diff: diff,
        z_score: Stats.z_score(diff, pathwise.se, bump.se) }
    end

    def run_convergence
      CONVERGENCE_GRID.map do |n_val|
        opts = { n_paths: n_val, seed: @seed, workers: @workers }
        ConvergencePoint.new(n: n_val, se_naive: Estimators.price_asian_naive(@params, **opts).se,
                             se_antithetic: Estimators.price_asian_antithetic(@params, **opts).se,
                             se_control_variate: Estimators.price_asian_control_variate(@params, **opts).se)
      end
    end

    def fit_all_slopes(convergence)
      counts = convergence.map(&:n)
      Slopes.new(naive: Stats.fit_loglog_slope(counts, convergence.map(&:se_naive)),
                 antithetic: Stats.fit_loglog_slope(counts, convergence.map(&:se_antithetic)),
                 control_variate: Stats.fit_loglog_slope(counts, convergence.map(&:se_control_variate)))
    end
  end
end
