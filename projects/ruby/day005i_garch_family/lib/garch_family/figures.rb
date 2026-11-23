# frozen_string_literal: true

module GarchFamily
  module Figures
    module_function

    def generate_all(payload, out_dir:)
      FileUtils.mkdir_p(out_dir)
      emit_part1_part2(payload, out_dir)
      emit_part3_part4_part5(payload, out_dir)
      emit_part6_to_part9(payload, out_dir)
    end

    def emit_part1_part2(payload, out_dir)
      p1 = payload[:part1_moments]
      p2 = payload[:part2_mle]
      draw_moments_chart(p1, "#{out_dir}/day005i_moments_var_kurt.svg")
      draw_acf_term_chart(p1, "#{out_dir}/day005i_garch_acf_term_structure.svg")
      draw_qmle_rmse_chart(p2, "#{out_dir}/day005i_qmle_rmse_convergence.svg")
      draw_sandwich_coverage_chart(p2, "#{out_dir}/day005i_sandwich_vs_naive_coverage.svg")
    end
    private_class_method :emit_part1_part2

    def emit_part3_part4_part5(payload, out_dir)
      p3 = payload[:part3_innovations]
      p4 = payload[:part4_backtest]
      p5 = payload[:part5_comparison]
      draw_student_t_chart(p3, "#{out_dir}/day005i_student_t_efficiency_lr.svg")
      draw_var_es_chart(p4, "#{out_dir}/day005i_var_es_backtest_comparison.svg")
      draw_nic_chart(p5, "#{out_dir}/day005i_news_impact_curves.svg")
      draw_model_sel_chart(p5, "#{out_dir}/day005i_qlike_dm_model_selection.svg")
    end
    private_class_method :emit_part3_part4_part5

    def emit_part6_to_part9(payload, out_dir)
      draw_dcc_chart(payload[:part6_dcc], "#{out_dir}/day005i_dcc_dynamic_correlation.svg")
      draw_har_chart(payload[:part7_har_rv], "#{out_dir}/day005i_har_rv_infill_and_forecast.svg")
      draw_break_chart(payload[:part8_structural_break], "#{out_dir}/day005i_structural_break_igarch.svg")
      draw_port_comp_chart(payload[:part9_portfolio_stress], "#{out_dir}/day005i_portfolio_component_risk.svg")
      draw_port_stress_chart(payload[:part9_portfolio_stress], "#{out_dir}/day005i_portfolio_stress_var_es.svg")
    end
    private_class_method :emit_part6_to_part9

    def draw_moments_chart(p1, path)
      groups = p1[:garch_grid].each_with_index.map do |r, i|
        { label: "G#{i + 1}", bars: [
          { value: r[:kurt_theo], color: '#2563eb' }, { value: r[:kurt_emp], color: '#16a34a' }
        ] }
      end
      Svg.write_bar_chart(path, title: 'GARCH(1,1) Kurtosis: Teoretis (Biru) vs MC 400k (Hijau)',
                                y_label: 'Kurtosis', groups: groups)
    end
    private_class_method :draw_moments_chart

    def draw_acf_term_chart(p1, path)
      lags = p1[:acf_sq_lags]
      series = [
        { label: 'ACF Kuadrat Teori', x: lags.map { |r| r[:lag] }, y: lags.map { |r| r[:theo] },
          color: '#1d4ed8', dashed: false },
        { label: 'ACF Kuadrat MC', x: lags.map { |r| r[:lag] }, y: lags.map { |r| r[:emp] },
          color: '#dc2626', dashed: true }
      ]
      Svg.write_line_chart(path, title: 'Fungsi Autokorelasi Return Kuadrat GARCH(1,1): Teori vs MC',
                                 x_label: 'Lag k', y_label: 'ACF(r^2)', series: series)
    end
    private_class_method :draw_acf_term_chart

    def draw_qmle_rmse_chart(p2, path)
      ts = p2[:normal_qmle].map { |r| r[:n_obs].to_f }
      series = [
        { label: 'Norm RMSE(alpha)', x: ts, y: p2[:normal_qmle].map { |r| r[:params][:alpha][:rmse] },
          color: '#2563eb', dashed: false },
        { label: 't(5) QMLE RMSE(alpha)', x: ts, y: p2[:student_t5_qmle].map { |r| r[:params][:alpha][:rmse] },
          color: '#dc2626', dashed: false },
        { label: 'Norm RMSE(beta)', x: ts, y: p2[:normal_qmle].map { |r| r[:params][:beta][:rmse] },
          color: '#16a34a', dashed: true }
      ]
      Svg.write_line_chart(path, title: 'Konvergensi Root-T QMLE (T = 500, 1000, 2500)',
                                 x_label: 'Panjang Sampel T', y_label: 'RMSE Parameter', series: series)
    end
    private_class_method :draw_qmle_rmse_chart

    def draw_sandwich_coverage_chart(p2, path)
      groups = %i[omega alpha beta].map do |k|
        nr = p2[:normal_qmle].last[:params][k]
        tr = p2[:student_t5_qmle].last[:params][k]
        { label: k.to_s, bars: [
          { value: nr[:coverage_naive_95], color: '#60a5fa' },
          { value: tr[:coverage_naive_95], color: '#ef4444' },
          { value: tr[:coverage_sandwich_95], color: '#16a34a' }
        ] }
      end
      Svg.write_bar_chart(path, title: 'Cakupan 95% (T=2500): Norm Naif vs t(5) Naif (Merah) & Sandwich (Hijau)',
                                y_label: 'Cakupan Empiris', groups: groups)
    end
    private_class_method :draw_sandwich_coverage_chart

    def draw_student_t_chart(p3, path)
      eff = p3[:efficiency_comparison].last[:params]
      groups = %i[omega alpha beta].map do |k|
        { label: "Efisiensi #{k}", bars: [{ value: eff[k][:mse_eff_ratio], color: '#2563eb' }] }
      end
      groups << { label: 'Null LR Mix 5%',
                  bars: [{ value: p3[:lr_boundary_null_normal][:rej_rate_mixture_5pct], color: '#16a34a' }] }
      Svg.write_bar_chart(path, title: 'Efisiensi Student-t MLE vs QMLE & Ukuran LR Campuran (T=2500)',
                                y_label: 'Rasio MSE / Rejection Rate', groups: groups)
    end
    private_class_method :draw_student_t_chart

    def draw_var_es_chart(p4, path)
      ro = p4[:rolling_oos]
      mods = %w[const ewma garch_norm garch_t garch_fhs]
      keys = %w[const_norm_1pct ewma_norm_1pct garch_norm_1pct garch_t_1pct garch_fhs_1pct]
      groups = mods.zip(keys).map do |lbl, key|
        { label: lbl, bars: [
          { value: ro[key][:hit_rate] * 100.0, color: '#2563eb' },
          { value: ro[key][:acerbi_z2], color: '#dc2626' }
        ] }
      end
      Svg.write_bar_chart(path, title: 'Rolling OOS 1%: Hit Rate % (Biru, Target 1.0) & Acerbi Z2 (Merah)',
                                y_label: 'Hit Rate % / Z2 ES', groups: groups)
    end
    private_class_method :draw_var_es_chart

    def draw_nic_chart(p5, path)
      x = p5[:asymmetric_dgp][:nic_grid]
      curves = p5[:asymmetric_dgp][:nic_curves]
      series = [
        { label: 'GARCH(1,1)', x: x, y: curves[:garch11][:values], color: '#6b7280', dashed: true },
        { label: 'GJR-GARCH(1,1)', x: x, y: curves[:gjr][:values], color: '#dc2626', dashed: false },
        { label: 'EGARCH(1,1)', x: x, y: curves[:egarch][:values], color: '#2563eb', dashed: false }
      ]
      Svg.write_line_chart(path, title: 'News Impact Curves pada DGP Asimetris (Leverage Negatif)',
                                 x_label: 'Kejutan eps_{t-1}', y_label: 'Variansi Berikut h_t', series: series)
    end
    private_class_method :draw_nic_chart

    def draw_model_sel_chart(p5, path)
      q = p5[:asymmetric_dgp][:oos_qlike]
      groups = %i[garch11 gjr egarch ewma].map do |k|
        { label: k.to_s, bars: [{ value: q[k], color: k == :gjr ? '#16a34a' : '#2563eb' }] }
      end
      Svg.write_bar_chart(path, title: 'OOS QLIKE Loss pada DGP Asimetris (Lebih Rendah = Lebih Baik)',
                                y_label: 'Mean OOS QLIKE', groups: groups)
    end
    private_class_method :draw_model_sel_chart

    def draw_dcc_chart(p6, path)
      anch = p6[:ccc_pd_anchor]
      t_idx = (1..anch[:sample_true_rhos].length).map(&:to_f)
      series = [
        { label: 'Korelasi DCC Sejati', x: t_idx, y: anch[:sample_true_rhos], color: '#111827', dashed: false },
        { label: 'Korelasi DCC Terfilter', x: t_idx, y: anch[:sample_filt_rhos], color: '#2563eb', dashed: true }
      ]
      Svg.write_line_chart(path, title: 'Korelasi Dinamis DCC(1,1): Sejati vs Filter Stage-2',
                                 x_label: 'Hari t', y_label: 'rho_{12,t}', series: series)
    end
    private_class_method :draw_dcc_chart

    def draw_har_chart(p7, path)
      grid = p7[:infill_anchor][:grid]
      series = [
        { label: 'MSE(RV - IV)', x: grid.map { |r| Math.log10(r[:m_bars]) },
          y: grid.map { |r| Math.log10(r[:mse_rv_iv]) }, color: '#16a34a', dashed: false }
      ]
      Svg.write_line_chart(path, title: 'Penskalaan Infill Intraday log10 MSE(RV - IV) vs log10(M)',
                                 x_label: 'log10(Jumlah Bar Intraday M)', y_label: 'log10 MSE', series: series)
    end
    private_class_method :draw_har_chart

    def draw_break_chart(p8, path)
      groups = [
        { label: 'Sejati', bars: [{ value: p8[:true_persistence], color: '#16a34a' }] },
        { label: 'Full Pooled (Spurious)', bars: [{ value: p8[:full_sample_persistence_mean], color: '#dc2626' }] },
        { label: 'Pre/Post Split', bars: [{ value: p8[:subsample_persistence_mean], color: '#2563eb' }] }
      ]
      Svg.write_bar_chart(path, title: 'Break Struktural Omega: Persistensi Semu IGARCH vs Pemulihan Split',
                                y_label: 'Estimasi Persistensi (a+b)', groups: groups)
    end
    private_class_method :draw_break_chart

    def draw_port_comp_chart(p9, path)
      st = p9[:stress_test]
      groups = 5.times.map do |i|
        { label: "Aset #{i + 1}", bars: [
          { value: st[:calm][:rc_pct][i] * 100.0, color: '#2563eb' },
          { value: st[:stress][:rc_pct][i] * 100.0, color: '#dc2626' }
        ] }
      end
      Svg.write_bar_chart(path, title: 'Dekomposisi Risiko Euler (% Kontribusi): Tenang vs Krisis',
                                y_label: 'Persentase Risiko %', groups: groups)
    end
    private_class_method :draw_port_comp_chart

    def draw_port_stress_chart(p9, path)
      st = p9[:stress_test]
      groups = [
        { label: 'Tenang (GJR-DCC)', bars: [
          { value: st[:calm][:port_var1] * 100.0, color: '#2563eb' },
          { value: st[:calm][:port_es1] * 100.0, color: '#60a5fa' }
        ] },
        { label: 'Krisis (rho=0.85, 2.5x)', bars: [
          { value: st[:stress][:port_var1] * 100.0, color: '#dc2626' },
          { value: st[:stress][:port_es1] * 100.0, color: '#f97316' }
        ] }
      ]
      Svg.write_bar_chart(path, title: 'Stres Portofolio 5-Aset: VaR 1% & ES 1% Harian (%)',
                                y_label: 'Persen Kerugian (%)', groups: groups)
    end
    private_class_method :draw_port_stress_chart
  end
end
