# frozen_string_literal: true

module BacktestOverfitting
  module Criteria
    module_function

    def evaluate(parts, bit_diffs, total_sec)
      c_first(parts).merge(c_second(parts, bit_diffs, total_sec))
    end

    def c_first(parts)
      {
        c1_engine_invariants_and_ar1_anchor: eval_c1(parts[:part1_anchors]),
        c2_sharpe_se_and_ci_coverage: eval_c2(parts[:part2_sharpe_stats]),
        c3_autocorr_lo_eta_and_stationary_bootstrap: eval_c3(parts[:part2_sharpe_stats][:autocorr_lo]),
        c4_selection_bias_expected_max: eval_c4(parts[:part3_selection_bias]),
        c5_dsr_calibration_vs_naive: eval_c5(parts[:part3_selection_bias][:dsr_calibration])
      }
    end
    private_class_method :c_first

    def c_second(parts, bit_diffs, total_sec)
      {
        c6_multiple_testing_fwer_and_bh_fdr: eval_c6(parts[:part4_multiple_testing]),
        c7_cscv_pbo_exact_and_planted: eval_c7(parts[:part5_cscv_pbo]),
        c8_transaction_cost_turnover_collapse: eval_c8(parts[:part5_cscv_pbo][:cost_study]),
        c9_purged_embargo_vs_naive_cv: eval_c9(parts[:part6_purged_cv]),
        c10_worker_invariance_1_vs_n: crit(bit_diffs.zero?, { bit_diffs: bit_diffs }, 'bit_diffs == 0',
                                           'Determinisme Process.fork bit-identik 1 vs N worker'),
        c11_runtime_under_360s: crit(total_sec < 360.0, { total_seconds: total_sec }, '< 360.0s',
                                     'Anggaran waktu eksekusi runner Ruby penuh')
      }
    end
    private_class_method :c_second

    def eval_c1(p1)
      ok = p1[:equity_identity_err] < 1.0e-13 && p1[:zero_cost_err] < 1.0e-14 &&
           p1[:no_lookahead_passed] && p1[:ar1_sign_z] <= 3.0 && p1[:gbm_bh_log_z] <= 3.0
      crit(ok, p1, 'eq_err < 1e-13, no_lookahead=true, ar1_z <= 3.0, gbm_z <= 3.0',
           'Identitas ekuitas/biaya/tanpa look-ahead dan anchor AR(1) & GBM B&H cocok teori')
    end
    private_class_method :eval_c1

    def eval_c2(p2)
      g_ok = p2[:gaussian_table].all? do |r|
        ((r[:emp_sd] - r[:theo_se]).abs / r[:theo_se]) < 0.03 && r[:coverage_95].between?(0.935, 0.965)
      end
      ng_ok = p2[:non_gaussian_table].all? do |r|
        ((r[:emp_sd] - r[:theo_se]).abs / r[:theo_se]) < 0.06 &&
          (r[:emp_sd] - r[:theo_se]).abs < (r[:emp_sd] - r[:naive_gauss_se]).abs &&
          r[:coverage_95_moment_adj].between?(0.92, 0.965)
      end
      crit(g_ok && ng_ok,
           { gauss_covs: p2[:gaussian_table].map { |r| r[:coverage_95] },
             nongauss_emp_vs_theo: p2[:non_gaussian_table].map { |r| [r[:dist], r[:emp_sd], r[:theo_se]] } },
           'rel_se_err < 3% Gauss & cov in [0.935,0.965]; non-Gauss adj lebih dekat ke theo_se drpd naif',
           'Metode delta SE Sharpe dan koreksi skewness/kurtosis cocok dengan MC')
    end
    private_class_method :eval_c2

    def eval_c3(ac)
      eta_ok = ac[:q_scaling].all? { |r| (r[:eta_mc] - r[:eta_theoretical]).abs < 0.04 }
      t20k = ac[:finite_sample_sweep].last[:by_q]
      t20k_ok = t20k.values.all? { |v| (v[:eta_mc_ratio] - v[:eta_asymp]).abs < 0.01 }
      cov = ac[:stationary_bootstrap_coverage_95]
      cov_ok = ac[:naive_coverage_95] < 0.88 && cov.between?(0.91, 0.97) && (cov - ac[:naive_coverage_95]) > 0.05
      status = if eta_ok && t20k_ok && cov_ok
                 ac[:margin_above_091_se] < 1.0 ? 'sebagian' : 'tercapai'
               else
                 'belum tercapai'
               end
      { status: status,
        measured: c3_measured(ac),
        threshold: '|eta_mc - eta_theor| < 0.04, |eta_T20k - asymp| < 0.01; naive < 0.88 & boot in [0.91,0.97]',
        note: 'eta(q) cocok setelah koreksi SR1_emp; cakupan bootstrap 0.914 +- 0.010 marginal (< 1 SE di atas 0.91)' }
    end
    private_class_method :eval_c3

    def c3_measured(ac)
      {
        eta_mc_vs_theo: ac[:q_scaling].map { |r| [r[:q], r[:eta_mc], r[:eta_theoretical], r[:eta_mc_pooled_true_sr1]] },
        eta_t20k_vs_asymp: ac[:finite_sample_sweep].last[:by_q].map { |q, v| [q, v[:eta_mc_ratio], v[:eta_asymp]] },
        naive_cov: ac[:naive_coverage_95], boot_cov: ac[:stationary_bootstrap_coverage_95],
        boot_se: ac[:boot_coverage_se], margin_above_091_se: ac[:margin_above_091_se]
      }
    end
    private_class_method :c3_measured

    def eval_c4(p3)
      t300 = p3[:expected_max_table].find { |r| r[:n_trials] == 300 }
      t1000 = p3[:expected_max_table].find { |r| r[:n_trials] == 1_000 }
      lib = p3[:library_noise_study]
      ok = t300[:abs_err] < 0.06 && t1000[:abs_err] < 0.04 &&
           lib[:avg_pairwise_corr] > 0.05 && lib[:n_eff_corr_formula] < 300 &&
           lib[:emp_mean_best_sharpe] < lib[:theo_max_iid_300]
      crit(ok,
           { exact_300: t300[:exact], approx_300: t300[:approx], avg_corr_lib: lib[:avg_pairwise_corr],
             n_eff: lib[:n_eff_corr_formula], lib_best_sr: lib[:emp_mean_best_sharpe],
             iid_bound_300: lib[:theo_max_iid_300] },
           '|exact-approx|<0.06 (N=300) & <0.04 (N=1k); korelasi pustaka -> n_eff < 300 & SR_best < iid',
           'Bias seleksi E[max] Gauss-Legendre vs aproksimasi & pengurangan N efektif pustaka')
    end
    private_class_method :eval_c4

    def eval_c5(dsr)
      se_best = (1.0 / Math.sqrt(dsr[:n_obs])) * 0.25 / Math.sqrt(dsr[:n_reps])
      iid_ok = (dsr[:emp_iid_mean_best_sharpe] - dsr[:theo_iid_mean_best_sharpe]).abs <= ((4.0 * se_best) + 0.004)
      cal_ok = dsr[:naive_pass_rate_95] > 0.95 && dsr[:dsr_pass_rate_95] <= 0.05
      crit(iid_ok && cal_ok,
           { emp_best: dsr[:emp_iid_mean_best_sharpe], theo_best: dsr[:theo_iid_mean_best_sharpe],
             naive_pass_rate_95: dsr[:naive_pass_rate_95], dsr_pass_rate_95: dsr[:dsr_pass_rate_95],
             dsr_median: dsr[:dsr_median], dsr_p95: dsr[:dsr_p95] },
           'naive_pass_rate > 0.95 & dsr_pass_rate <= 0.05; iid max SR cocok sigma_SR*E[max]',
           'Kalibrasi Deflated Sharpe Ratio vs ambang Sharpe naif pada 2000 noise nol')
    end
    private_class_method :eval_c5

    def eval_c6(p4)
      ind = p4[:independent]
      bh_fdr = ind[:methods][:bh][:fdr]
      bh_se_bound = 3.0 * Math.sqrt(bh_fdr * (1.0 - bh_fdr) / ind[:n_reps])
      ok_ind = ind[:methods][:bonferroni][:fwer] <= 0.05 &&
               ind[:methods][:holm][:fwer] <= 0.05 &&
               (bh_fdr - ind[:theoretical_bh_fdr]).abs <= [bh_se_bound, 0.008].max &&
               ind[:methods][:bh][:power] > ind[:methods][:bonferroni][:power]
      rw = p4[:romano_wolf_rho04]
      ok_rw = rw[:romano_wolf][:fwer] <= 0.075 && rw[:romano_wolf][:power] >= rw[:bonferroni][:power]
      crit(ok_ind && ok_rw,
           { bonf_fwer: ind[:methods][:bonferroni][:fwer], bh_fdr: bh_fdr,
             theo_bh_fdr: ind[:theoretical_bh_fdr], bh_power: ind[:methods][:bh][:power],
             rw_fwer_rho04: rw[:romano_wolf][:fwer], rw_power_rho04: rw[:romano_wolf][:power],
             bonf_power_rho04: rw[:bonferroni][:power] },
           'bonf_fwer<=0.05, |bh_fdr - pi0*alpha|<=3 SE, bh_power>bonf_power, RW_power>=bonf_power',
           'Pengujian ganda: Bonferroni/Holm/BH (indep & corr rho=0.4) dan Romano-Wolf')
    end
    private_class_method :eval_c6

    def eval_c7(p5)
      na = p5[:noise_anchor]
      z_noise = (na[:emp_mean_pbo] - na[:exact_theoretical_pbo]).abs / [na[:emp_se_pbo], 1.0e-12].max
      pl_pbo = p5[:planted_signal][:pbo]
      t_pbos = p5[:scaling][:vs_t_planted].map { |r| r[:pbo] }
      ok = z_noise <= 3.0 && pl_pbo < 0.10 && t_pbos.last <= t_pbos.first
      crit(ok,
           { noise_emp_pbo: na[:emp_mean_pbo], noise_exact_pbo: na[:exact_theoretical_pbo],
             noise_z: z_noise, planted_pbo: pl_pbo, t_scaling_pbos: t_pbos },
           'noise_z <= 3.0 vs 0.50; planted_pbo < 0.10; PBO turun saat T membesar',
           'CSCV C(16,8)=12870 pembagian: simetri noise 0.5 dan sinyal ditanam PBO << 0.5')
    end
    private_class_method :eval_c7

    def eval_c8(cs)
      r0 = cs[:rows][0]
      r20 = cs[:rows][2]
      hi_drop = r0[:hi_turnover_mean_sharpe] - r20[:hi_turnover_mean_sharpe]
      lo_drop = r0[:lo_turnover_mean_sharpe] - r20[:lo_turnover_mean_sharpe]
      ok = hi_drop > (lo_drop * 3.0) && r20[:hi_turnover_mean_sharpe].negative?
      crit(ok,
           { hi_to_0bps: r0[:hi_turnover_mean_sharpe], hi_to_20bps: r20[:hi_turnover_mean_sharpe],
             lo_to_0bps: r0[:lo_turnover_mean_sharpe], lo_to_20bps: r20[:lo_turnover_mean_sharpe],
             hi_drop: hi_drop, lo_drop: lo_drop },
           'hi_drop > 3 * lo_drop dan hi_turnover_sharpe(20 bps) < 0',
           'Strategi berfrekuensi/turnover tinggi runtuh pada biaya 5-20 bps')
    end
    private_class_method :eval_c8

    def eval_c9(p6)
      ok = p6[:naive_shuffled_corr_mean] > 0.30 && p6[:purged_zero_z] <= 3.0
      crit(ok,
           { naive_corr_mean: p6[:naive_shuffled_corr_mean],
             contig_unpurged_mean: p6[:contiguous_unpurged_corr_mean],
             purged_corr_mean: p6[:purged_embargo_corr_mean],
             purged_corr_se: p6[:purged_embargo_corr_se], purged_zero_z: p6[:purged_zero_z] },
           'naive_corr_mean > 0.30 (bocor) & |purged_corr_mean / se| <= 3.0 (anchor 0)',
           'Validasi K-fold naif membocorkan label tumpang-tindih h=10, purged+embargo bebas bias')
    end
    private_class_method :eval_c9

    def crit(ok, measured, threshold, note)
      { status: ok ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
