# frozen_string_literal: true

module CopulaTailDependence
  module Criteria
    module_function

    def evaluate(p1, p2, p3, p4, bit_diffs, total_sec)
      build_c1_c2(p1)
        .merge(build_c3_c4(p1, p2))
        .merge(c10_non_gaussian_euler_allocation: build_c10_data(p2))
        .merge(build_c5_c6_c11(p3))
        .merge(build_c7_c8_c9(p4, bit_diffs, total_sec))
    end

    def build_c1_c2(p1)
      max_tau_z = p1[:copulas].map { |c| c[:tau_z].abs }.max
      gauss_row = p1[:copulas].find { |c| c[:copula] == 'gaussian' }
      sp_err = (gauss_row[:spearman_emp] - gauss_row[:spearman_exact]).abs
      g_u_zs = gauss_row[:tail_by_u].map { |r| r[:z_finite_u].abs }
      g_grid = p1[:gaussian_chi_grid].map { |r| r[:chi] }
      g_mono = g_grid.each_cons(2).all? { |a, b| a > b } && g_grid.last < 0.005
      non_g_zs = p1[:copulas].reject { |c| c[:copula] == 'gaussian' }.map { |c| c[:tail_by_u].last[:z_asymp].abs }

      {
        c1_tau_and_spearman_vs_closed_form: crit(
          max_tau_z <= 3.0 && sp_err <= 0.015,
          { max_tau_z: max_tau_z, gaussian_spearman_err: sp_err },
          'max |tau_z| <= 3.0 dan |spearman_err| <= 0.015 (R=200k)',
          'Kendall tau (SE jackknife) dan Spearman cocok dengan rumus tertutup'
        ),
        c2_tail_dependence_convergence: crit(
          g_u_zs.max <= 3.0 && g_mono && non_g_zs.max <= 3.0,
          { gauss_max_finite_u_z: g_u_zs.max, gauss_chi_monotone_to_zero: g_mono,
            non_gauss_u999_max_z: non_g_zs.max },
          'Gaussian |z(u)| <= 3.0 di u in {0.9,0.99,0.999}, chi_Gauss turun monoton hingga 1-1e-10, ' \
          'dan t/Clayton/Gumbel |chi_hat(0.999)-lambda| <= 3 SE',
          'Konvergensi ketergantungan ekor: Gaussian cocok dengan chi(u) eksak berhingga ' \
          'dan kopula berekor cocok dengan lambda asimtotik'
        )
      }
    end
    private_class_method :build_c1_c2

    def build_c3_c4(p1, p2)
      ex_sum_err = (p2[:gaussian_exact][:es] - p2[:gaussian_exact][:es_contribs].sum(0.0)).abs
      {
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
    private_class_method :build_c3_c4

    def build_c10_data(p2)
      t4 = p2[:mc_t4]
      cl = p2[:mc_clayton]
      mg = p2[:mc_gaussian]
      t4_sum_err = (t4[:contribs].sum(0.0) - t4[:es]).abs
      cl_sum_err = (cl[:contribs].sum(0.0) - cl[:es]).abs
      t4_diff = t4[:contribs].each_with_index.map { |v, i| v - mg[:contribs][i] }
      cl_diff = cl[:contribs].each_with_index.map { |v, i| v - mg[:contribs][i] }
      t4_diff_se = t4[:contrib_ses].each_with_index.map { |s, i| Math.hypot(s, mg[:contrib_ses][i]) }
      cl_diff_se = cl[:contrib_ses].each_with_index.map { |s, i| Math.hypot(s, mg[:contrib_ses][i]) }

      crit(
        t4_sum_err < 1.0e-10 && cl_sum_err < 1.0e-10,
        {
          t4_sum_err: t4_sum_err, clayton_sum_err: cl_sum_err,
          t4_minus_gauss_contribs: t4_diff, t4_minus_gauss_ses: t4_diff_se,
          clayton_minus_gauss_contribs: cl_diff, clayton_minus_gauss_ses: cl_diff_se
        },
        '|sum(c_i) - ES| < 1e-10 (<= 3 SE) untuk t(4) dan Clayton serta selisih vs Gaussian dengan SE',
        'Alokasi Euler non-Gaussian tetap menjumlah ke ES total dan menggeser kontribusi aset'
      )
    end
    private_class_method :build_c10_data

    def build_c5_c6_c11(p3)
      n5000 = p3[:tau_vs_mle].last
      ms = p3[:misspecification]
      ss = ms[:single_sample]
      aic_ok = ss[:t_best][:aic] < ss[:gaussian][:aic] && ss[:t_best][:aic] < ss[:clayton][:aic]

      {
        c5_cml_vs_tau_inversion: crit(
          n5000[:rho_mle_bias].abs < 0.01 && n5000[:rho_tau_bias].abs < 0.01,
          { n5000_rho_mle_bias: n5000[:rho_mle_bias], n5000_rho_tau_bias: n5000[:rho_tau_bias] },
          '|rho_mle_bias| < 0.01 dan |rho_tau_bias| < 0.01 pada n=5000 (R=500)',
          'CML Gaussian dan inversi tau konsisten serta tak-bias'
        ),
        c6_misspecification_aic_lrt_es99: crit(
          ss[:gauss_rejected_at_5pct] && ms[:es99_fitted_gauss_mean] < ms[:es99_true_t4],
          { lrt_t_vs_gauss: ss[:lrt_t_vs_gauss], es99_true_t4: ms[:es99_true_t4],
            es99_fitted_gauss_mean: ms[:es99_fitted_gauss_mean], es99_rel_bias: ms[:es99_rel_bias] },
          'LRT > 3.841 (Gaussian ditolak) dan ES99_gauss < ES99_true_t4',
          'Kopula Gaussian ditolak oleh LRT dan meremehkan ES99 portofolio di bawah kebenaran t(4)'
        ),
        c11_aic_selects_t_copula: crit(
          aic_ok,
          { aic_t: ss[:t_best][:aic], aic_gaussian: ss[:gaussian][:aic], aic_clayton: ss[:clayton][:aic] },
          'AIC(t) < AIC(Gaussian) dan AIC(t) < AIC(Clayton)',
          'Kriteria informasi AIC memilih kopula Student-t di atas Gaussian dan Clayton'
        )
      }
    end
    private_class_method :build_c5_c6_c11

    def build_c7_c8_c9(p4, bit_diffs, total_sec)
      max_is_z = p4.map { |r| r[:is_z].abs }.max
      min_vr = p4.map { |r| r[:variance_reduction] }.min
      {
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
    private_class_method :build_c7_c8_c9

    def crit(ok, measured, threshold, note)
      { status: ok ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
