# frozen_string_literal: true

module CopulaTailDependence
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
      p1, t_p1 = timed { run_part1_dependence }
      p2, t_p2 = timed { run_part2_euler }
      p3, t_p3 = timed { run_part3_estimation }
      p4, t_p4 = timed { run_part4_importance }
      bit_diffs = check_worker_invariance
      total_sec = Process.clock_gettime(Process::CLOCK_MONOTONIC) - t0

      criteria = evaluate_criteria(p1, p2, p3, p4, bit_diffs, total_sec)
      {
        metadata: build_metadata(t_p1, t_p2, t_p3, t_p4, total_sec),
        criteria: criteria,
        part1_dependence: p1,
        part2_euler: p2,
        part3_estimation: p3,
        part4_importance: p4
      }
    end

    def write_outputs(payload, repo_root:)
      json_path = File.join(repo_root, 'outputs', 'day005e', 'results.json')
      JsonOut.write(json_path, payload)
      fig_dir = File.join(repo_root, 'outputs', 'figures')
      svg_paths = write_figures(payload, fig_dir)
      { json: json_path, svgs: svg_paths }
    end

    private

    def timed
      t_start = Process.clock_gettime(Process::CLOCK_MONOTONIC)
      res = yield
      [res, Process.clock_gettime(Process::CLOCK_MONOTONIC) - t_start]
    end

    def build_metadata(t1, t2, t3, t4, total_sec)
      yjit = defined?(RubyVM::YJIT) && RubyVM::YJIT.enabled?
      {
        day: 5,
        session: 'e',
        topic_source: 'kopula_ketergantungan_ekor_dan_alokasi_euler',
        language: 'ruby',
        ruby_version: RUBY_VERSION,
        yjit_enabled: yjit,
        main_seed: @seed,
        workers_used: @workers,
        run_completed_without_errors: true,
        elapsed_seconds: {
          part1_dependence: t1,
          part2_euler: t2,
          part3_estimation: t3,
          part4_importance: t4,
          total: total_sec
        }
      }
    end

    def run_part1_dependence
      u_grid = [0.9, 0.99, 0.999]
      copula_rows = COPULA_SPECS.map { |spec| evaluate_copula_mc(spec, u_grid) }
      orthant_err = max_orthant_anchor_error
      { copulas: copula_rows, max_orthant_err: orthant_err }
    end

    def evaluate_copula_mc(spec, u_grid)
      key = spec[:key]
      params = spec[:params]
      chunks = Parallel.map_chunks(n_items: @r_part1, workers: @workers, chunk_size: 1_000) do |cidx, len|
        rng = Rng::Sampler.new(Rng.derive_seed(@seed, 100 + copula_id(key), cidx))
        u1 = Array.new(len, 0.0)
        u2 = Array.new(len, 0.0)
        len.times do |i|
          a, b = Copulas.sample_pair(key, params, rng)
          u1[i] = a
          u2[i] = b
        end
        tau_b = Dependence.kendall_tau(u1, u2)
        rho_b = Dependence.spearman_rho(u1, u2)
        tails = u_grid.map { |u| Dependence.empirical_tail(u1, u2, u, side: spec[:tail_side]) }
        [tau_b, rho_b, tails]
      end

      taus = chunks.map { |c| c[0] }
      rhos = chunks.map { |c| c[1] }
      tau_emp = Stats.mean(taus)
      tau_se = Stats.jackknife_se(taus)
      tau_ex = Copulas.exact_kendall_tau(key, params)
      rho_emp = Stats.mean(rhos)
      rho_ex = Copulas.exact_spearman_rho(key, params)
      tail_ex = if spec[:tail_side] == :lower
                  Copulas.exact_lower_tail(key, params)
                else
                  Copulas.exact_upper_tail(key, params)
                end
      tail_emp = u_grid.each_with_index.map do |u, ui|
        { u: u, empirical: Stats.mean(chunks.map { |c| c[2][ui] }) }
      end

      {
        copula: key.to_s, label: spec[:label],
        tau_exact: tau_ex, tau_emp: tau_emp, tau_se: tau_se, tau_z: (tau_emp - tau_ex) / tau_se,
        spearman_exact: rho_ex, spearman_emp: rho_emp,
        tail_side: spec[:tail_side].to_s, tail_exact: tail_ex, tail_by_u: tail_emp
      }
    end

    def copula_id(key)
      { gaussian: 1, t: 2, clayton: 3, gumbel: 4 }.fetch(key)
    end

    def max_orthant_anchor_error
      sp = Special
      errs = [
        (sp.bivariate_norm_cdf(-1.2, 0.7, 0.0) - (sp.norm_cdf(-1.2) * sp.norm_cdf(0.7))).abs,
        (sp.bivariate_norm_cdf(-1.5, 0.3, 1.0 - 1.0e-15) - sp.norm_cdf(-1.5)).abs,
        (sp.bivariate_norm_cdf(0.5, 0.5, -1.0 + 1.0e-15) - (sp.norm_cdf(0.5) + sp.norm_cdf(0.5) - 1.0)).abs,
        (sp.bivariate_norm_cdf(-1.0, 0.5, 0.6) - sp.bivariate_norm_cdf(0.5, -1.0, 0.6)).abs,
        (sp.bivariate_norm_cdf(0.0, 0.0, 0.5) - (0.25 + (Math.asin(0.5) / (2.0 * Math::PI)))).abs
      ]
      errs.max
    end

    def run_part2_euler
      ex = Euler.gaussian_exact
      fd = Euler.gaussian_es_finite_diff
      mc_g = Euler.mc_euler_es(copula: :gaussian, n_paths: @r_euler, seed: @seed, workers: @workers)
      mc_t = Euler.mc_euler_es(copula: :t, df: 4.0, n_paths: @r_euler, seed: @seed, workers: @workers)
      mc_c = Euler.mc_euler_es(copula: :clayton, theta: 2.0, n_paths: @r_euler, seed: @seed, workers: @workers)
      max_z = ex[:es_contribs].each_with_index.map do |c, i|
        (mc_g[:contribs][i] - c).abs / mc_g[:contrib_ses][i]
      end.max
      max_fd_err = ex[:es_contribs].each_with_index.map { |c, i| (fd[i] - c).abs }.max
      {
        gaussian_exact: ex,
        gaussian_finite_diff: fd,
        max_fd_err: max_fd_err,
        mc_gaussian: mc_g,
        max_mc_gaussian_z: max_z,
        mc_t4: mc_t,
        mc_clayton: mc_c
      }
    end

    def run_part3_estimation
      tau_vs_mle = Fit.compare_tau_vs_mle(r_reps: 500, ns: [250, 1000, 5000], seed: @seed, workers: @workers)
      misspec = Fit.misspecification_study(r_reps: 200, n_obs: 2000, seed: @seed, workers: @workers)
      es_t_marg_g_cop = Fit.equicorr_portfolio_es99(copula: :gaussian, rho: 0.5, df: 4.0, marginal: :t,
                                                    n_paths: 60_000, seed: @seed + 2)
      es_n_marg_t_cop = Fit.equicorr_portfolio_es99(copula: :t, rho: 0.5, df: 4.0, marginal: :normal,
                                                    n_paths: 60_000, seed: @seed + 3)
      {
        tau_vs_mle: tau_vs_mle,
        misspecification: misspec,
        marginal_vs_copula_es99: {
          t4_marginal_gaussian_copula: es_t_marg_g_cop,
          normal_marginal_t4_copula: es_n_marg_t_cop
        }
      }
    end

    def run_part4_importance
      [
        Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.5, n_paths: 100_000, seed: @seed, workers: @workers),
        Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.8, n_paths: 100_000, seed: @seed, workers: @workers),
        Importance.compare_naive_vs_is(d: 5, a: -3.0, rho: 0.5, n_paths: 100_000, seed: @seed, workers: @workers),
        Importance.compare_naive_vs_is(d: 5, a: -3.0, rho: 0.8, n_paths: 100_000, seed: @seed, workers: @workers)
      ]
    end

    def check_worker_invariance
      a = Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: @seed, workers: 1)
      b = Importance.compare_naive_vs_is(d: 2, a: -3.0, rho: 0.5, n_paths: 4_000, seed: @seed, workers: @workers)
      [a[:is_mean]].pack('D') == [b[:is_mean]].pack('D') ? 0 : 1
    end

    def evaluate_criteria(p1, p2, p3, p4, bit_diffs, total_sec)
      evaluate_part1_2_criteria(p1, p2).merge(evaluate_part3_4_criteria(p3, p4, bit_diffs, total_sec))
    end

    def evaluate_part1_2_criteria(p1, p2)
      max_tau_z = p1[:copulas].map { |c| c[:tau_z].abs }.max
      gauss_row = p1[:copulas].find { |c| c[:copula] == 'gaussian' }
      sp_err = (gauss_row[:spearman_emp] - gauss_row[:spearman_exact]).abs
      max_tail_err = p1[:copulas].map { |c| (c[:tail_by_u].last[:empirical] - c[:tail_exact]).abs }.max
      g_chi = gauss_row[:tail_by_u].map { |r| r[:empirical] }
      g_mono = g_chi[0] > g_chi[1] && g_chi[1] > g_chi[2]
      ex_sum_err = (p2[:gaussian_exact][:es] - p2[:gaussian_exact][:es_contribs].sum(0.0)).abs

      {
        c1_tau_and_spearman_vs_closed_form: crit(
          max_tau_z <= 3.0 && sp_err <= 0.015,
          { max_tau_z: max_tau_z, gaussian_spearman_err: sp_err },
          'max |tau_z| <= 3.0 dan |spearman_err| <= 0.015 (R=200k)',
          'Kendall tau (SE jackknife) dan Spearman cocok dengan rumus tertutup'
        ),
        c2_tail_dependence_convergence: crit(
          max_tail_err < 0.08 && g_mono,
          { max_tail_err_u999: max_tail_err, gaussian_chi_decreasing: g_mono },
          'max |lambda_hat(0.999) - lambda_exact| < 0.08 dan chi_Gauss(u) turun monoton',
          'Konvergensi koefisien ketergantungan ekor ke nilai analitik dan bias berhingga Gaussian'
        ),
        c3_bivariate_orthant_accuracy: crit(
          p1[:max_orthant_err] < 1.0e-10, { max_orthant_err: p1[:max_orthant_err] },
          'max_orthant_err < 1e-10',
          'Probabilitas ortan bivariat Gaussian via Gauss-Legendre akurat terhadap kasus batas'
        ),
        c4_euler_gaussian_exact_and_mc: crit(
          ex_sum_err < 1.0e-10 && p2[:max_fd_err] < 1.0e-6 && p2[:max_mc_gaussian_z] <= 3.0,
          { exact_sum_err: ex_sum_err, max_fd_err: p2[:max_fd_err], max_mc_z: p2[:max_mc_gaussian_z] },
          'exact_sum_err < 1e-10, max_fd_err < 1e-6, max_mc_z <= 3.0',
          'Alokasi Euler Gaussian memenuhi teorema Euler, beda hingga, dan MC'
        )
      }
    end

    def evaluate_part3_4_criteria(p3, p4, bit_diffs, total_sec)
      n5000 = p3[:tau_vs_mle].last
      ms = p3[:misspecification]
      max_is_z = p4.map { |r| r[:is_z].abs }.max
      min_vr = p4.map { |r| r[:variance_reduction] }.min

      {
        c5_cml_vs_tau_inversion: crit(
          n5000[:rho_mle_bias].abs < 0.01 && n5000[:rho_tau_bias].abs < 0.01,
          { n5000_rho_mle_bias: n5000[:rho_mle_bias], n5000_rho_tau_bias: n5000[:rho_tau_bias] },
          '|rho_mle_bias| < 0.01 dan |rho_tau_bias| < 0.01 pada n=5000 (R=500)',
          'CML Gaussian dan inversi tau konsisten serta tak-bias'
        ),
        c6_misspecification_aic_lrt_es99: crit(
          ms[:single_sample][:gauss_rejected_at_5pct] && ms[:es99_fitted_gauss_mean] < ms[:es99_true_t4],
          {
            lrt_t_vs_gauss: ms[:single_sample][:lrt_t_vs_gauss], es99_true_t4: ms[:es99_true_t4],
            es99_fitted_gauss_mean: ms[:es99_fitted_gauss_mean], es99_rel_bias: ms[:es99_rel_bias]
          },
          'LRT > 3.841 (Gaussian ditolak) dan ES99_gauss < ES99_true_t4',
          'Kopula Gaussian ditolak oleh LRT dan meremehkan ES99 portofolio di bawah kebenaran t(4)'
        ),
        c7_importance_sampling_vs_anchor: crit(
          max_is_z <= 3.0 && min_vr > 50.0,
          { max_is_z: max_is_z, min_variance_reduction: min_vr },
          'max |is_z| <= 3.0 dan min_variance_reduction > 50.0 (R=100k)',
          'Importance sampling pergeseran-rerata + conditional MC tak-bias dan hemat variansi besar'
        ),
        c8_worker_invariance_1_vs_n: crit(
          bit_diffs.zero?, { bit_diffs: bit_diffs }, 'bit_diffs == 0',
          'Determinisme paralel Process.fork bit-identik 1 vs N worker'
        ),
        c9_runtime_under_300s: crit(
          total_sec < 300.0, { total_seconds: total_sec }, '< 300.0s',
          'Anggaran waktu eksekusi runner Ruby penuh'
        )
      }
    end

    def crit(ok, measured, threshold, note)
      {
        status: ok ? 'tercapai' : 'belum tercapai',
        measured: measured,
        threshold: threshold,
        note: note
      }
    end

    def write_figures(payload, fig_dir)
      paths = {}
      paths[:scatter] = write_scatter_fig(fig_dir)
      paths[:chi_u] = write_chi_fig(payload[:part1_dependence], fig_dir)
      paths[:lambda_param] = write_lambda_param_fig(fig_dir)
      paths[:euler] = write_euler_fig(payload[:part2_euler], fig_dir)
      paths[:nu_profile] = write_nu_profile_fig(payload[:part3_estimation][:misspecification], fig_dir)
      paths[:es99_spec] = write_es99_spec_fig(payload[:part3_estimation], fig_dir)
      paths[:vr_is] = write_vr_fig(payload[:part4_importance], fig_dir)
      paths
    end

    def write_scatter_fig(fig_dir)
      rng = Rng::Sampler.new(@seed + 900)
      colors = { gaussian: '#2563eb', t: '#dc2626', clayton: '#059669', gumbel: '#d97706' }
      panels = COPULA_SPECS.map do |spec|
        u1 = Array.new(600, 0.0)
        u2 = Array.new(600, 0.0)
        600.times do |i|
          a, b = Copulas.sample_pair(spec[:key], spec[:params], rng)
          u1[i] = a
          u2[i] = b
        end
        { label: spec[:label], u1: u1, u2: u2, color: colors[spec[:key]] }
      end
      out = File.join(fig_dir, 'day005e_copula_scatter.svg')
      Svg.write_scatter_panels(out, title: 'Sampel Kopula Bivariat (N=600 per panel)', panels: panels)
      out
    end

    def write_chi_fig(p1, fig_dir)
      colors = ['#2563eb', '#dc2626', '#059669', '#d97706']
      us = [0.9, 0.99, 0.999]
      series = p1[:copulas].each_with_index.map do |row, idx|
        {
          label: row[:copula],
          x: us,
          y: row[:tail_by_u].map { |r| r[:empirical] },
          color: colors[idx],
          dashed: false
        }
      end
      out = File.join(fig_dir, 'day005e_tail_chi_u.svg')
      Svg.write_line_chart(out, title: 'Konvergensi Empiris chi(u) vs Ambang u', x_label: 'Ambang u',
                                y_label: 'chi(u)', series: series)
      out
    end

    def write_lambda_param_fig(fig_dir)
      rhos = (1..18).map { |i| i * 0.05 }
      l_g = rhos.map { 0.0 }
      l_t4 = rhos.map { |r| Copulas.exact_upper_tail(:t, rho: r, df: 4.0) }
      l_c = rhos.map { |r| Copulas.exact_lower_tail(:clayton, theta: (2.0 * r) / (1.0 - r + 0.1)) }
      out = File.join(fig_dir, 'day005e_lambda_vs_param.svg')
      Svg.write_line_chart(out, title: 'Koefisien Ketergantungan Ekor vs Parameter', x_label: 'Parameter (rho / skala)',
                                y_label: 'lambda', series: [
                                  { label: 'Gaussian', x: rhos, y: l_g, color: '#2563eb', dashed: true },
                                  { label: 'Student-t (df=4)', x: rhos, y: l_t4, color: '#dc2626', dashed: false },
                                  { label: 'Clayton (lambda_L)', x: rhos, y: l_c, color: '#059669', dashed: false }
                                ])
      out
    end

    def write_euler_fig(p2, fig_dir)
      assets = [1.0, 2.0, 3.0, 4.0, 5.0]
      out = File.join(fig_dir, 'day005e_euler_contribs.svg')
      Svg.write_line_chart(out, title: 'Alokasi Modal Euler ES99 per Aset (d=5)', x_label: 'Indeks Aset',
                                y_label: 'Kontribusi ES99', series: [
                                  { label: 'Gaussian Eksak', x: assets, y: p2[:gaussian_exact][:es_contribs],
                                    color: '#111827', dashed: true },
                                  { label: 'Gaussian MC', x: assets, y: p2[:mc_gaussian][:contribs],
                                    color: '#2563eb', dashed: false },
                                  { label: 'Student-t (df=4)', x: assets, y: p2[:mc_t4][:contribs],
                                    color: '#dc2626', dashed: false },
                                  { label: 'Clayton', x: assets, y: p2[:mc_clayton][:contribs],
                                    color: '#059669', dashed: false }
                                ])
      out
    end

    def write_nu_profile_fig(ms, fig_dir)
      profs = ms[:single_sample][:t_profiles].reject { |r| r[:df].infinite? }
      out = File.join(fig_dir, 'day005e_nu_profile.svg')
      Svg.write_line_chart(out, title: 'Profil Log-Likelihood Kopula Student-t terhadap nu',
                                x_label: 'Derajat Kebebasan nu', y_label: 'Log-Likelihood', series: [
                                  { label: 'Profile LogLik', x: profs.map { |r| r[:df] },
                                    y: profs.map { |r| r[:loglik] }, color: '#2563eb', dashed: false }
                                ])
      out
    end

    def write_es99_spec_fig(p3, fig_dir)
      ms = p3[:misspecification]
      mv = p3[:marginal_vs_copula_es99]
      xs = [1.0, 2.0, 3.0, 4.0]
      ys = [ms[:es99_fitted_gauss_mean], ms[:es99_true_t4], mv[:normal_marginal_t4_copula],
            mv[:t4_marginal_gaussian_copula]]
      out = File.join(fig_dir, 'day005e_es99_misspec.svg')
      Svg.write_line_chart(out, title: 'ES99 Portofolio (1:Gauss, 2:t4, 3:Norm+t4, 4:t4+Gauss)',
                                x_label: 'Spesifikasi Model (1..4)', y_label: 'ES99', series: [
                                  { label: 'ES99 Portofolio', x: xs, y: ys, color: '#dc2626', dashed: false }
                                ])
      out
    end

    def write_vr_fig(p4, fig_dir)
      xs = [1.0, 2.0, 3.0, 4.0]
      ys = p4.map { |r| Math.log10(r[:variance_reduction]) }
      out = File.join(fig_dir, 'day005e_is_variance_reduction.svg')
      Svg.write_line_chart(out, title: 'Log10 Faktor Penghematan Variansi IS (1..4: d2r5, d2r8, d5r5, d5r8)',
                                x_label: 'Skenario', y_label: 'log10(VR Factor)', series: [
                                  { label: 'log10(VR)', x: xs, y: ys, color: '#059669', dashed: false }
                                ])
      out
    end
  end
end
