# frozen_string_literal: true

module CopulaTailDependence
  module RunnerParts
    module_function

    def run_part1(seed, workers, r_part1, copula_specs)
      u_grid = [0.9, 0.99, 0.999]
      {
        copulas: copula_specs.map { |spec| eval_copula(spec, u_grid, seed, workers, r_part1) },
        gaussian_chi_grid: Copulas.gaussian_chi_grid(0.6),
        max_orthant_err: max_orthant_err
      }
    end

    def eval_copula(spec, u_grid, seed, workers, r_part1)
      key = spec[:key]
      params = spec[:params]
      id = { gaussian: 1, t: 2, clayton: 3, gumbel: 4 }.fetch(key)
      chunks = Parallel.map_chunks(n_items: r_part1, workers: workers, chunk_size: 1_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(seed, 100 + id, cidx))
        u1 = Array.new(len, 0.0)
        u2 = Array.new(len, 0.0)
        len.times { |i| u1[i], u2[i] = Copulas.sample_pair(key, params, rng) }
        [Dependence.kendall_tau(u1, u2), Dependence.spearman_rho(u1, u2),
         u_grid.map { |u| Dependence.empirical_tail(u1, u2, u, side: spec[:tail_side]) }]
      end
      summarize_copula(spec, u_grid, chunks, r_part1)
    end

    def summarize_copula(spec, u_grid, chunks, r_part1)
      key = spec[:key]
      params = spec[:params]
      taus = chunks.map { |c| c[0] }
      tau_emp = Stats.mean(taus)
      tau_se = Stats.jackknife_se(taus)
      tau_ex = Copulas.exact_kendall_tau(key, params)
      tail_ex = if spec[:tail_side] == :lower
                  Copulas.exact_lower_tail(key, params)
                else
                  Copulas.exact_upper_tail(key, params)
                end
      {
        copula: key.to_s, label: spec[:label], tau_exact: tau_ex, tau_emp: tau_emp,
        tau_se: tau_se, tau_z: (tau_emp - tau_ex) / tau_se,
        spearman_exact: Copulas.exact_spearman_rho(key, params),
        spearman_emp: Stats.mean(chunks.map { |c| c[1] }),
        tail_side: spec[:tail_side].to_s, tail_exact: tail_ex,
        tail_by_u: build_tail_rows(key, params, u_grid, chunks, tail_ex, r_part1)
      }
    end

    def build_tail_rows(key, params, u_grid, chunks, tail_ex, r_part1)
      u_grid.each_with_index.map do |u, ui|
        emp = Stats.mean(chunks.map { |c| c[2][ui] })
        p_hat = emp * (1.0 - u)
        se_chi = Math.sqrt([p_hat * (1.0 - p_hat), 0.0].max / r_part1) / (1.0 - u)
        ex_u = key == :gaussian ? Copulas.exact_gaussian_chi(u, params[:rho]) : nil
        { u: u, empirical: emp, se: se_chi, exact_finite_u: ex_u,
          z_finite_u: ex_u ? (emp - ex_u) / [se_chi, 1.0e-15].max : nil,
          z_asymp: (emp - tail_ex) / [se_chi, 1.0e-15].max }
      end
    end

    def max_orthant_err
      sp = Special
      [
        (sp.bivariate_norm_cdf(-1.2, 0.7, 0.0) - (sp.norm_cdf(-1.2) * sp.norm_cdf(0.7))).abs,
        (sp.bivariate_norm_cdf(-1.5, 0.3, 1.0 - 1.0e-15) - sp.norm_cdf(-1.5)).abs,
        (sp.bivariate_norm_cdf(0.5, 0.5, -1.0 + 1.0e-15) - (sp.norm_cdf(0.5) + sp.norm_cdf(0.5) - 1.0)).abs,
        (sp.bivariate_norm_cdf(-1.0, 0.5, 0.6) - sp.bivariate_norm_cdf(0.5, -1.0, 0.6)).abs,
        (sp.bivariate_norm_cdf(0.0, 0.0, 0.5) - (0.25 + (Math.asin(0.5) / (2.0 * Math::PI)))).abs
      ].max
    end
  end

  class Runner
    COPULA_SPECS = [
      { key: :gaussian, label: 'Gaussian(rho=0.6)', params: { rho: 0.6 }, tail_side: :upper },
      { key: :t, label: 'Student-t(df=4, rho=0.6)', params: { rho: 0.6, df: 4.0 }, tail_side: :upper },
      { key: :clayton, label: 'Clayton(theta=2.0)', params: { theta: 2.0 }, tail_side: :lower },
      { key: :gumbel, label: 'Gumbel(theta=2.0)', params: { theta: 2.0 }, tail_side: :upper }
    ].freeze

    def initialize(seed: MAIN_SEED, workers: nil, r_part1: 200_000, r_euler: 200_000)
      @seed = seed
      @workers = Parallel.resolve_workers(workers)
      @r_part1 = r_part1
      @r_euler = r_euler
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
      p1, t_p1 = timed { RunnerParts.run_part1(@seed, @workers, @r_part1, COPULA_SPECS) }
      p2, t_p2 = timed { run_part2_euler }
      p3, t_p3 = timed { run_part3_estimation }
      p4, t_p4 = timed { run_part4_importance }
      bit_diffs = check_worker_invariance
      total_sec = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      {
        metadata: build_metadata(t_p1, t_p2, t_p3, t_p4, total_sec),
        criteria: Criteria.evaluate(p1, p2, p3, p4, bit_diffs, total_sec),
        part1_dependence: p1, part2_euler: p2, part3_estimation: p3, part4_importance: p4
      }
    end

    def write_outputs(payload, repo_root:)
      json_path = File.join(repo_root, 'outputs', 'day005e', 'results.json')
      JsonOut.write(json_path, payload)
      fig_dir = File.join(repo_root, 'outputs', 'figures')
      { json: json_path, svgs: Figures.write_all(payload, fig_dir, seed: @seed, copula_specs: COPULA_SPECS) }
    end

    private

    def timed
      t_start = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      [yield, Process.clock_gettime(Process::CLOCK_MONOTONIC) - t_start]
    end

    def build_metadata(t1, t2, t3, t4, total_sec)
      {
        day: 5, session: 'e', topic_source: 'kopula_ketergantungan_ekor_dan_alokasi_euler',
        language: 'ruby', ruby_version: RUBY_VERSION,
        yjit_enabled: defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?,
        main_seed: @seed, workers_used: @workers, run_completed_without_errors: true,
        elapsed_seconds: { part1_dependence: t1, part2_euler: t2, part3_estimation: t3,
                           part4_importance: t4, total: total_sec }
      }
    end

    def run_part2_euler
      ex = Euler.gaussian_exact
      fd = Euler.gaussian_es_finite_diff
      mc_g = Euler.mc_euler_es(copula: :gaussian, n_paths: @r_euler, seed: @seed, workers: @workers)
      max_z = ex[:es_contribs].each_with_index.map { |c, i| (mc_g[:contribs][i] - c).abs / mc_g[:contrib_ses][i] }.max
      {
        gaussian_exact: ex, gaussian_finite_diff: fd,
        max_fd_err: ex[:es_contribs].each_with_index.map { |c, i| (fd[i] - c).abs }.max,
        mc_gaussian: mc_g, max_mc_gaussian_z: max_z,
        mc_t4: Euler.mc_euler_es(copula: :t, df: 4.0, n_paths: @r_euler, seed: @seed, workers: @workers),
        mc_clayton: Euler.mc_euler_es(copula: :clayton, theta: 2.0, n_paths: @r_euler, seed: @seed, workers: @workers)
      }
    end

    def run_part3_estimation
      {
        tau_vs_mle: Fit.compare_tau_vs_mle(r_reps: 500, ns: [250, 1000, 5000], seed: @seed, workers: @workers),
        misspecification: Fit.misspecification_study(r_reps: 200, n_obs: 2000, seed: @seed, workers: @workers),
        marginal_vs_copula_es99: {
          t4_marginal_gaussian_copula: Fit.equicorr_portfolio_es99(
            copula: :gaussian, rho: 0.5, df: 4.0, marginal: :t, n_paths: 60_000, seed: @seed + 2
          ),
          normal_marginal_t4_copula: Fit.equicorr_portfolio_es99(
            copula: :t, rho: 0.5, df: 4.0, marginal: :normal, n_paths: 60_000, seed: @seed + 3
          )
        }
      }
    end

    def run_part4_importance
      [[2, 0.5], [2, 0.8], [5, 0.5], [5, 0.8]].map do |d, rho|
        Importance.compare_naive_vs_is(d: d, a: -3.0, rho: rho, n_paths: 100_000, seed: @seed, workers: @workers)
      end
    end

    def check_worker_invariance
      a = Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: @seed, workers: 1)
      b = Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: @seed, workers: @workers)
      [a[:is_mean]].pack('D') == [b[:is_mean]].pack('D') ? 0 : 1
    end
  end
end
