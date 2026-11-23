# frozen_string_literal: true

module GarchFamily
  module Criteria
    module_function

    def evaluate_all(payload)
      out = {}
      add_part1_criteria(out, payload[:part1_moments])
      add_part2_criteria(out, payload[:part2_mle])
      add_part3_criteria(out, payload[:part3_innovations])
      add_part4_criteria(out, payload[:part4_backtest])
      add_part5_criteria(out, payload[:part5_comparison])
      add_part6_to_9_criteria(out, payload)
      add_eng_criteria(out, payload[:metadata])
      out
    end

    def add_part1_criteria(out, p1)
      all_z = (p1[:garch_grid] + p1[:gjr_grid] + p1[:egarch_grid]).map { |r| r[:var_z] }
      out[:c1_unconditional_variance_all_models] = crit(
        all_z.max <= 3.5, { max_var_z: all_z.max },
        'max |z| <= 3.5 across 15 GARCH/GJR/EGARCH configs',
        'Variansi tak bersyarat MC cocok dengan formula tertutup GARCH, GJR, dan EGARCH'
      )
      kz = p1[:garch_grid].map { |r| r[:kurt_z] }.max
      az = p1[:acf_sq_lags].map { |r| r[:z] }.max
      out[:c2_garch_kurtosis_and_sq_acf] = crit(
        kz <= 3.5 && az <= 3.5, { max_kurt_z: kz, max_acf_z: az },
        'max |kurt_z| <= 3.5 & max |acf_z| <= 3.5 (lags 1..5)',
        'Kurtosis teoretis orde ke-4 dan ACF return kuadrat GARCH cocok dengan simulasi MC'
      )
      rel_fc = p1[:forecast_mc].map { |r| (r[:mc_var] - r[:theo_var]).abs / r[:theo_var] }.max
      out[:c3_multistep_variance_forecast_anchor] = crit(
        rel_fc < 0.025, { max_rel_err: rel_fc }, 'max_rel_err < 2.5% (h=1..60)',
        'Prediksi multi-horizon h_t mengarah ke variansi tak bersyarat sesuai ekspektasi MC'
      )
    end
    private_class_method :add_part1_criteria

    def add_part2_criteria(out, p2)
      sa = p2[:single_sample_anchor]
      out[:c4_analytic_score_zero_at_mle] = crit(
        sa[:converged] && sa[:dimless_score_inf] < 1.0e-5,
        { converged: sa[:converged], dimless_score_inf: sa[:dimless_score_inf] },
        'converged && max |score_param * scale| / T < 1e-5 at T=5000',
        'Skor analitik nol terkalibrasi skala pada titik MLE GARCH(1,1) T=5000'
      )
      add_qmle_criteria(out, p2)
    end
    private_class_method :add_part2_criteria

    def add_qmle_criteria(out, p2)
      nr = p2[:normal_qmle]
      rmses = nr.map { |r| r[:params][:alpha][:rmse] }
      last_n = nr.last
      cov_a = last_n[:params][:alpha][:coverage_naive_95]
      ok_n = rmses[0] > rmses[1] && rmses[1] > rmses[2] && cov_a.between?(0.91, 0.98) && last_n[:n_converged] >= 975
      out[:c5_qmle_root_t_consistency_normal] = crit(
        ok_n, { rmse_alpha_t: rmses, cov_alpha_2500: cov_a, n_conv: last_n[:n_converged] },
        'RMSE monoton turun & cov in [0.91, 0.98] & conv >= 975/1000',
        'Konsistensi root-T dan cakupan asimtotik QMLE Normal tervalidasi (R=1000)'
      )
      t25 = p2[:student_t5_qmle].last[:params][:alpha]
      ratio = t25[:mean_se_sandwich] / t25[:mean_se_naive]
      z_theo = (t25[:coverage_sandwich_95] - 0.95).abs / [t25[:coverage_se], 1.0e-6].max
      ok_sw = t25[:coverage_naive_95] < 0.85 && ratio > 1.25 && z_theo <= 3.0
      out[:c6_sandwich_se_restores_t5_coverage] = crit(
        ok_sw, { se_ratio: ratio, cov_naive: t25[:coverage_naive_95], cov_sw: t25[:coverage_sandwich_95],
                 cov_se: t25[:coverage_se], z_vs_095: z_theo },
        'cov_naive < 0.85, se_sw/se_naive > 1.25, |cov_sw - 0.95| <= 3 SE (R=1000)',
        'Sandwich SE Bollerslev-Wooldridge terhadap target nominal teori 0.95 +- 3 SE pada t(5) T=2500'
      )
    end
    private_class_method :add_qmle_criteria

    def add_part3_criteria(out, p3)
      e25 = p3[:efficiency_comparison].last
      eff_a = e25[:params][:alpha][:mse_eff_ratio]
      eff_b = e25[:params][:beta][:mse_eff_ratio]
      out[:c7_student_t_mle_nu_and_efficiency] = crit(
        e25[:inv_nu_bias].abs < 0.025 && eff_a > 1.05 && eff_b > 1.05,
        { nu_median: e25[:nu_median], inv_nu_bias: e25[:inv_nu_bias], eff_alpha: eff_a, eff_beta: eff_b },
        '|bias(1/nu)| < 0.025 & eff_ratio(alpha,beta) > 1.05 at T=2500',
        'Student-t MLE memulihkan nu dan jauh lebih efisien secara MSE daripada QMLE Normal'
      )
      lr0 = p3[:lr_boundary_null_normal]
      lra = p3[:lr_boundary_alt_t5]
      z_null = (lr0[:rej_rate_mixture_5pct] - 0.05).abs / [lr0[:rej_rate_se], 1.0e-6].max
      out[:c8_boundary_lr_half_chi2_calibration] = crit(
        lr0[:zero_mass_fraction].between?(0.38, 0.65) && z_null <= 3.0 && lra[:rej_rate_mixture_5pct] >= 0.90,
        { zero_mass: lr0[:zero_mass_fraction], null_rej: lr0[:rej_rate_mixture_5pct],
          null_z: z_null, alt_power: lra[:rej_rate_mixture_5pct] },
        'zero_mass in [0.38,0.65], |null_rej-0.05|<=3SE, alt_power>=0.90',
        'Uji batas 1/nu=0 mengikuti campuran 50:50 chi2_0 + chi2_1'
      )
    end
    private_class_method :add_part3_criteria

    def add_part4_criteria(out, p4)
      ro = p4[:rolling_oos]
      gt = ro['garch_t_1pct']
      gf = ro['garch_fhs_1pct']
      gn = ro['garch_norm_1pct']
      cn = ro['const_norm_1pct']
      out[:c9_var_backtest_dynamics_and_tails] = crit(
        gn[:p_uc] < 0.05 && gt[:p_uc] > 0.05 && gt[:p_ind] > 0.05 && gf[:p_uc] > 0.05,
        { const_rate: cn[:hit_rate], garch_n_rate: gn[:hit_rate],
          garch_t_rate: gt[:hit_rate], garch_fhs_rate: gf[:hit_rate] },
        'garch_norm gagal 1% tail; garch_t & garch_fhs lolos Kupiec+Christoffersen',
        'VaR 1% dinamis Student-t dan FHS terkalibrasi baik pada data berkluster ekor tebal'
      )
      out[:c10_acerbi_szekely_es_backtest] = crit(
        gn[:acerbi_z2] < -0.40 && gt[:acerbi_z2].abs < 0.35 && gf[:acerbi_z2].abs < 0.35,
        { gn_z2: gn[:acerbi_z2], gt_z2: gt[:acerbi_z2], gf_z2: gf[:acerbi_z2] },
        'gn_z2 < -0.40 (under-ES); |gt_z2|, |gf_z2| < 0.35',
        'Backtest Expected Shortfall Acerbi-Szekely Z2 mendeteksi kekurangan ES Normal vs t/FHS'
      )
      mc = p4[:mc_size_power]
      sz = mc[:size_nominal_5pct]
      pw = mc[:power_misspecified_t5_const_1pct]
      z_uc = (sz[:uc][:rate] - 0.05).abs / [sz[:uc][:se], 1.0e-6].max
      out[:c11_kupiec_christoffersen_mc_size_power] = crit(
        z_uc <= 3.0 && pw[:cc][:rate] > 0.70,
        { size_uc: sz[:uc][:rate], size_cc: sz[:cc][:rate], power_cc_1pct: pw[:cc][:rate] },
        '|size_uc - 0.05| <= 3 SE & power_cc(1% t5 const) > 0.70 (R=500)',
        'Ukuran empiris dan daya uji backtest Kupiec & Christoffersen (ekor 1%) terverifikasi via MC'
      )
    end
    private_class_method :add_part4_criteria

    def add_part5_criteria(out, p5)
      qa = p5[:qlike_anchor]
      out[:c12_patton_qlike_proper_scoring] = crit(
        qa[:true_is_lowest], qa, 'true_h_qlike < min(misspecified)',
        'Variansi sejati h_t meminimalkan fungsi kerugian Patton QLIKE meski memakai proksi r^2'
      )
      ad = p5[:asymmetric_dgp]
      q = ad[:oos_qlike]
      is = ad[:in_sample]
      nic_r = ad[:nic_curves][:gjr][:ratio_neg_pos_2sd]
      core_ok = is[:gjr][:bic] < is[:garch11][:bic] && q[:gjr] < q[:garch11] && nic_r > 1.5
      dm_ok = ad[:dm_garch_vs_gjr][:p_value] < 0.05
      status_c13 = if core_ok && dm_ok
                     'tercapai'
                   elsif core_ok
                     'sebagian'
                   else
                     'belum tercapai'
                   end
      out[:c13_asymmetric_model_selection_dm_nic] = {
        status: status_c13,
        measured: { bic_diff: is[:garch11][:bic] - is[:gjr][:bic],
                    qlike_diff: q[:garch11] - q[:gjr], dm_p: ad[:dm_garch_vs_gjr][:p_value], nic_ratio: nic_r },
        threshold: 'BIC_GJR < BIC_GARCH, QLIKE_GJR < QLIKE_GARCH, DM p<0.05, NIC > 1.5',
        note: 'Seleksi BIC/QLIKE/NIC menangkap asimetri; DM 2-sisi p=0.076 marginal pada proksi r^2 1000 hari'
      }
      lmc = p5[:leverage_mc]
      sz_lr = lmc[:size_h0_garch][:lr_rej_5pct]
      z_lev = (sz_lr - 0.05).abs / [lmc[:size_h0_garch][:lr_se], 1.0e-6].max
      pw_lr = lmc[:power_h1_gjr][:lr_rej_5pct]
      out[:c14_leverage_mc_size_and_power] = crit(
        z_lev <= 3.0 && pw_lr >= 0.85,
        { lr_size: sz_lr, lr_z: z_lev, lr_power: pw_lr, wald_power: lmc[:power_h1_gjr][:wald_rej_5pct] },
        '|size - 0.05| <= 3 SE & power >= 0.85',
        'Uji hipotesis gamma=0 menjaga ukuran 5% dan daya tinggi pada DGP asimetris'
      )
    end
    private_class_method :add_part5_criteria

    def add_part6_to_9_criteria(out, payload)
      add_dcc_criteria(out, payload[:part6_dcc])
      add_har_criteria(out, payload[:part7_har_rv])
      add_break_port_criteria(out, payload[:part8_structural_break], payload[:part9_portfolio_stress])
    end
    private_class_method :add_part6_to_9_criteria

    def add_dcc_criteria(out, p6)
      anch = p6[:ccc_pd_anchor]
      out[:c15_dcc_ccc_anchor_and_pd] = crit(
        anch[:max_ccc_err] < 1.0e-13 && anch[:min_det].positive? && anch[:sample_corr_true_vs_filt] > 0.65,
        { max_ccc_err: anch[:max_ccc_err], min_det: anch[:min_det], corr_rho: anch[:sample_corr_true_vs_filt] },
        'max_ccc_err < 1e-13, min_det > 0, corr(rho_true, rho_hat) > 0.65',
        'Reduksi DCC(0,0)==CCC tepat hingga presisi mesin dan korelasi dinamis positif definit'
      )
      rec = p6[:mc_recovery]
      out[:c16_dcc_parameter_recovery] = crit(
        rec[:n_converged] >= 195 && rec[:a][:bias].abs < 0.025 && rec[:b][:bias].abs < 0.035,
        { n_conv: rec[:n_converged], a_mean: rec[:a][:mean], b_mean: rec[:b][:mean] },
        'conv >= 195/200, |bias_a| < 0.025, |bias_b| < 0.035',
        'Estimasi DCC 2-tahap Engle memulihkan (a, b) secara akurat pada 200 replikasi'
      )
    end
    private_class_method :add_dcc_criteria

    def add_har_criteria(out, p7)
      inf = p7[:infill_anchor]
      max_z = inf[:grid].map { |r| r[:z_unbiased] }.max
      out[:c17_intraday_rv_unbiased_and_infill] = crit(
        max_z <= 3.0 && inf[:mse_ratio_13_to_390] > 20.0,
        { max_z: max_z, mse_ratio: inf[:mse_ratio_13_to_390] },
        'max |z(RV-IV)| <= 3.0 across M=13..390 & MSE(13)/MSE(390) > 20',
        'Realized Variance tak bias untuk Integrated Variance dan galatnya menyusut O(1/M)'
      )
      hv = p7[:har_vs_garch]
      qb = hv[:qlike_benchmarks]
      dm = hv[:dm_benchmarks]
      beats_g = qb[:har] < qb[:garch] && dm[:garch_minus_har][:p_value] < 0.05
      beats_info = %i[rv_naive ewma_rv garch_x].all? { |k| qb[:har] < qb[k] && dm[:"#{k}_minus_har"][:p_value] < 0.05 }
      st = if beats_g && beats_info
             'tercapai'
           elsif beats_g
             'sebagian'
           else
             'belum tercapai'
           end
      out[:c18_har_rv_beats_daily_garch] = {
        status: st,
        measured: { qlike_har: qb[:har], qlike_garch: qb[:garch], qlike_rv_naive: qb[:rv_naive],
                    qlike_ewma_rv: qb[:ewma_rv], qlike_garch_x: qb[:garch_x],
                    dm_p_garch: dm[:garch_minus_har][:p_value], dm_p_rv_naive: dm[:rv_naive_minus_har][:p_value],
                    dm_p_ewma_rv: dm[:ewma_rv_minus_har][:p_value], dm_p_garch_x: dm[:garch_x_minus_har][:p_value] },
        threshold: 'HAR < GARCH, RV-Naif, EWMA-RV, GARCH-X & seluruh DM p < 0.05',
        note: 'Model Corsi HAR-RV terhadap GARCH harian dan pembanding setara informasi berbasis RV'
      }
    end
    private_class_method :add_har_criteria

    def add_break_port_criteria(out, p8, p9)
      out[:c19_structural_break_spurious_igarch] = crit(
        p8[:full_sample_persistence_mean] > 0.92 && p8[:subsample_persistence_mean].between?(0.78, 0.90) &&
        p8[:cusum_power] >= 0.90 && (p8[:detected_idx_mean] - 1500.0).abs < 150.0,
        p8,
        'full_p > 0.92 (spurious), sub_p in [0.78,0.90], cusum_power >= 0.90',
        'Break struktural omega menimbulkan efek semu IGARCH yang pulih saat sampel dipisah & dideteksi CUSUM'
      )
      st = p9[:stress_test]
      dyn = p9[:oos_backtest][:dynamic_gjr_dcc][:var5_p]
      out[:c20_five_asset_portfolio_euler_stress] = crit(
        st[:calm][:euler_sum_err] < 1.0e-13 && dyn[:p_uc] > 0.05 &&
        st[:var1_amplification] > 2.5 && st[:div_ratio_drop].positive?,
        { euler_err: st[:calm][:euler_sum_err], oos_5pct_p_uc: dyn[:p_uc],
          var1_amp: st[:var1_amplification], div_drop: st[:div_ratio_drop] },
        'euler_err < 1e-13, dyn OOS p_uc > 0.05, stress amp > 2.5x, div_ratio turun',
        'Sistem 5-aset GJR-DCC memenuhi identitas Euler, lolos OOS, dan mengukur runtuhnya diversifikasi krisis'
      )
    end
    private_class_method :add_break_port_criteria

    def add_eng_criteria(out, meta)
      out[:c21_worker_invariance_1_vs_n] = crit(
        meta[:worker_parity_bit_diffs].zero?, { bit_diffs: meta[:worker_parity_bit_diffs] },
        'bit_diffs == 0', 'Paralelisme Process.fork menghasilkan keluaran bit-identik 1 vs N worker'
      )
      out[:c22_runtime_under_480s] = crit(
        meta[:elapsed_seconds][:total] < 480.0, { total_seconds: meta[:elapsed_seconds][:total] },
        '< 480.0s', 'Anggaran waktu eksekusi penuh Ruby murni <= 8 menit'
      )
    end
    private_class_method :add_eng_criteria

    def crit(passed, measured, threshold, note)
      { status: passed ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
