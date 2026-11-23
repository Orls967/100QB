# frozen_string_literal: true

module FactorPortfolioLab
  # Main experiment runner executing Parts 1-12, writing results.json and 26 SVG figures.
  module Runner
    module_function

    def run_all(root_dir:)
      t0 = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      data = {
        spiked: SpikedStudy.run(master_seed: MAIN_SEED, n_reps: 500),
        grs: GrsStudy.run(master_seed: MAIN_SEED, n_reps: 5000),
        fama_macbeth: FmStudy.run(master_seed: MAIN_SEED, n_reps: 1000),
        risk_parity: RpStudy.run(master_seed: MAIN_SEED, sizes: [10, 25, 50]),
        black_litterman: BlStudy.run(master_seed: MAIN_SEED, n_mc: 500_000),
        markowitz: MvStudy.run(master_seed: MAIN_SEED, n_assets: 20, n_points: 45),
        resampling: Resampling.run_michaud(master_seed: MAIN_SEED, n_boot: 500),
        backtest: BacktestStudy.run(master_seed: MAIN_SEED, n_worlds: 300),
        attribution: Attribution.misspecified_k_study(master_seed: MAIN_SEED),
        stress: Stress.run_study(master_seed: MAIN_SEED),
        robustness: Robustness.run(master_seed: MAIN_SEED, n_grs_reps: 800, n_port_worlds: 120),
        complexity: Complexity.run(master_seed: MAIN_SEED, sizes: [10, 25, 50, 100, 200])
      }

      criteria = Criteria.evaluate(data)
      fig_dir = File.join(root_dir, 'outputs', 'figures')
      figures = Figures.generate_all(data, fig_dir)
      elapsed = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      payload = {
        project: 'day005j_factor_portfolio_lab',
        ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?,
        main_seed: MAIN_SEED,
        elapsed_seconds: elapsed,
        n_figures: figures.length,
        criteria: criteria,
        studies: strip_heavy_arrays(data)
      }

      out_json = File.join(root_dir, 'outputs', 'day005j', 'results.json')
      JsonOut.write(out_json, payload)
      print_summary(criteria, figures.length, elapsed, out_json)
      payload
    end

    def strip_heavy_arrays(data)
      out = data.dup
      mv = data[:markowitz].dup
      mv[:unconstrained] = mv[:unconstrained].map { |p| p.except(:weights) }
      mv[:long_only] = mv[:long_only].map { |p| p.except(:weights) }
      mv[:capped] = mv[:capped].map { |p| p.except(:weights) }
      out[:markowitz] = mv
      out
    end

    def print_summary(criteria, n_figs, elapsed, out_json)
      passed = criteria.values.count { |v| v[:status] == 'tercapai' }
      puts "day005j complete in #{elapsed.round(2)}s | #{passed}/#{criteria.size} criteria tercapai | " \
           "#{n_figs} SVGs -> #{out_json}"
      criteria.each do |k, v|
        next if v[:status] == 'tercapai'

        puts "  [BELUM TERCAPAI] #{k}: observed=#{v[:observed]} threshold=#{v[:threshold]} (#{v[:note]})"
      end
    end
  end
end
