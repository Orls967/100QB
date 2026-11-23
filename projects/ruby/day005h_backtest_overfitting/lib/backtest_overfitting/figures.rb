# frozen_string_literal: true

module BacktestOverfitting
  module Figures
    module_function

    def write_all(payload, fig_dir)
      {
        best_sharpe: write_best_sharpe(payload[:part3_selection_bias], fig_dir),
        coverage: write_coverage(payload[:part2_sharpe_stats], fig_dir),
        multiple_testing: write_multiple_testing(payload[:part4_multiple_testing], fig_dir),
        cscv_lambda: write_cscv_lambda(payload[:part5_cscv_pbo], fig_dir),
        is_oos_cost: write_is_oos_cost(payload[:part5_cscv_pbo], fig_dir),
        purged_cv: write_purged_cv(payload[:part6_purged_cv], fig_dir)
      }
    end

    def write_best_sharpe(p3, fig_dir)
      tbl = p3[:expected_max_table]
      xs = tbl.map { |r| Math.log10(r[:n_trials]) }
      out = File.join(fig_dir, 'day005h_best_sharpe_vs_emax.svg')
      Svg.write_line_chart(out, title: 'E[max Z] Eksak vs Aproksimasi Euler-Mascheroni (log10 N)',
                                x_label: 'log10(N strategi)', y_label: 'E[max Z]', series: [
                                  { label: 'Gauss-Legendre Eksak', x: xs, y: tbl.map { |r| r[:exact] },
                                    color: '#111827', dashed: false },
                                  { label: 'Aproksimasi Bailey-LdP', x: xs, y: tbl.map { |r| r[:approx] },
                                    color: '#dc2626', dashed: true }
                                ])
      out
    end

    def write_coverage(p2, fig_dir)
      out = File.join(fig_dir, 'day005h_sharpe_coverage.svg')
      g = p2[:gaussian_table]
      ng = p2[:non_gaussian_table]
      ac = p2[:autocorr_lo]
      groups = [
        { label: 'Gauss 250', bars: [{ value: g[0][:coverage_95], color: '#2563eb' }] },
        { label: 'Gauss 1k', bars: [{ value: g[1][:coverage_95], color: '#2563eb' }] },
        { label: 't(5) Naif/Adj', bars: [{ value: ng[0][:coverage_95_naive_gauss], color: '#dc2626' },
                                         { value: ng[0][:coverage_95_moment_adj], color: '#059669' }] },
        { label: 'AR(1) Naif/Boot', bars: [{ value: ac[:naive_coverage_95], color: '#dc2626' },
                                           { value: ac[:stationary_bootstrap_coverage_95], color: '#059669' }] }
      ]
      Svg.write_bar_chart(out, title: 'Cakupan Interval 95% Sharpe: Naif vs Moment/Stationary Bootstrap',
                               y_label: 'Cakupan Aktual', groups: groups)
      out
    end

    def write_multiple_testing(p4, fig_dir)
      out = File.join(fig_dir, 'day005h_multiple_testing_fdr.svg')
      ind = p4[:independent][:methods]
      cor = p4[:correlated_rho04][:methods]
      groups = %i[unadjusted bonferroni holm bh].map do |m|
        { label: m.to_s.upcase, bars: [
          { value: ind[m][:fdr], color: '#2563eb' },
          { value: ind[m][:power], color: '#059669' },
          { value: cor[m][:fdr], color: '#dc2626' }
        ] }
      end
      Svg.write_bar_chart(out, title: 'Pengujian Ganda (alpha=0.05): FDR & Power (Biru/Hijau=Indep, Merah=Corr)',
                               y_label: 'Tingkat (0..1)', groups: groups)
      out
    end

    def write_cscv_lambda(p5, fig_dir)
      out = File.join(fig_dir, 'day005h_cscv_logit_lambda.svg')
      Svg.write_histogram(out, title: 'Distribusi Logit OOS lambda CSCV (12.870 Pembagian C(16,8))',
                               x_label: 'Logit lambda = ln(omega / (1 - omega))', y_label: 'Proporsi',
                               series_vals: [
                                 { label: 'Noise Murni (PBO~0.5)',
                                   values: p5[:noise_anchor][:first_mat][:lambdas_sample],
                                   color: '#dc2626', dashed: false },
                                 { label: 'Sinyal Ditanam (PBO<<0.5)',
                                   values: p5[:planted_signal][:lambdas_sample],
                                   color: '#059669', dashed: false }
                               ])
      out
    end

    def write_is_oos_cost(p5, fig_dir)
      out = File.join(fig_dir, 'day005h_is_vs_oos_degradation.svg')
      rows = p5[:cost_study][:rows]
      xs = rows.map { |r| r[:cost_bps] }
      Svg.write_line_chart(out, title: 'Degradasi Biaya Transaksi: Turnover Tinggi vs Rendah (0..20 bps)',
                                x_label: 'Biaya Transaksi c (bps)', y_label: 'Sharpe Harian Rata-rata', series: [
                                  { label: 'High Turnover (>=0.40)',
                                    x: xs, y: rows.map { |r| r[:hi_turnover_mean_sharpe] },
                                    color: '#dc2626', dashed: false },
                                  { label: 'Low Turnover (<0.08)',
                                    x: xs, y: rows.map { |r| r[:lo_turnover_mean_sharpe] },
                                    color: '#059669', dashed: false },
                                  { label: 'Strategi Terbaik',
                                    x: xs, y: rows.map { |r| r[:best_sharpe] },
                                    color: '#2563eb', dashed: true }
                                ])
      out
    end

    def write_purged_cv(p6, fig_dir)
      out = File.join(fig_dir, 'day005h_purged_vs_naive_cv.svg')
      Svg.write_histogram(out, title: 'Korelasi OOS Label Tumpang-Tindih (h=10): Naive CV vs Purged+Embargo',
                               x_label: 'Korelasi OOS pada Fitur Noise Independen', y_label: 'Proporsi',
                               series_vals: [
                                 { label: 'Naive Shuffled K-Fold (Bocor)', values: p6[:naive_sample],
                                   color: '#dc2626', dashed: false },
                                 { label: 'Purged + Embargo K-Fold (Bersih)', values: p6[:purged_sample],
                                   color: '#059669', dashed: false }
                               ])
      out
    end
  end
end
