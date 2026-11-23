# frozen_string_literal: true

module CdoTranche
  module Criteria
    module_function

    def evaluate(p1, p2, p3, p4, p5, bit_diffs, total_sec)
      build_c1_c3(p1, p2).merge(build_c4_c6(p3, p4)).merge(build_c7_c9(p5, bit_diffs, total_sec))
    end

    def build_c1_c3(p1, p2)
      gran_errs = p1[:granularity].map { |r| r[:abs_err] }
      gran_mono = gran_errs.each_cons(2).all? { |a, b| a > b }
      {
        c1_asb_matches_exact_binomial: crit(
          p1[:max_asb_vs_binomial_err] < 1.0e-12,
          { max_asb_vs_binomial_err: p1[:max_asb_vs_binomial_err] },
          'max_asb_vs_binomial_err < 1e-12',
          'Rekursi ASB cocok dengan distribusi binomial eksak hingga presisi mesin'
        ),
        c2_lhp_closed_form_and_granularity: crit(
          p1[:max_lhp_bivariate_err] < 1.0e-9 && gran_mono,
          { max_lhp_bivariate_err: p1[:max_lhp_bivariate_err], granularity_monotone: gran_mono },
          'max_lhp_bivariate_err < 1e-9 dan galat LHP turun monoton terhadap N',
          'LHP kuadratur cocok dengan ortan bivariat dan konvergen O(1/N) ke ASB berhingga'
        ),
        c3_tranche_sum_conservation: crit(
          p2[:homog_pv_sum_err] < 1.0e-12 && p2[:hetero_pv_sum_err] < 1.0e-4,
          { homog_pv_sum_err: p2[:homog_pv_sum_err], hetero_pv_sum_err: p2[:hetero_pv_sum_err] },
          'homog_pv_sum_err < 1e-12 dan hetero_pv_sum_err < 1e-4',
          'Penjumlahan PV proteksi seluruh tranche sama dengan PV proteksi portofolio total'
        )
      }
    end
    private_class_method :build_c1_c3

    def build_c4_c6(p3, p4)
      bcs = p3[:base_correlations].map { |r| r[:base_corr] }
      bc_found = bcs.all?
      norm_els = p3[:normalized_tranche_els]
      norm_mono = norm_els.each_cons(2).all? { |a, b| a >= b }
      skew_dir = bcs.last - bcs.first
      max_is_z = p4.map { |r| r[:is_z].abs }.max
      vr_senior = p4.last[:variance_reduction]

      {
        c4_base_correlation_skew_and_arbitrage: crit(
          bc_found && norm_mono && p3[:quadrature_conv_diff] < 1.0e-8,
          { base_corrs: bcs, skew_30_minus_3: skew_dir, normalized_tranche_els: norm_els,
            norm_el_monotone_decreasing: norm_mono, quadrature_conv_diff: p3[:quadrature_conv_diff] },
          'base_corr di (0.01,0.95), EL ternormalisasi turun monoton, dan konvergensi kuadratur 2D < 1e-8',
          'Kalibrasi base correlation terhadap pasar t(4) deterministik & bebas arbitrase'
        ),
        c5_compound_correlation_non_uniqueness: crit(
          p3[:compound_scans].length == 5 && p3[:senior_compound_roots].empty?,
          { compound_scans: p3[:compound_scans] },
          'Pemindaian grid rho in (0.01,0.95) step 0.01 pada 5 tranche dan deteksi kegagalan/non-unik',
          'Menunjukkan non-monotonisitas dan kegagalan akar compound correlation pada tranche senior'
        ),
        c6_importance_sampling_vs_asb: crit(
          max_is_z <= 3.0 && vr_senior > 100.0,
          { max_is_z: max_is_z, super_senior_vr: vr_senior },
          'max |is_z| <= 3.0 dan super_senior_vr (30-100%) > 100.0 (R=100k)',
          'Conditional MC + Exponential Twisting tak-bias vs ASB dan hemat variansi besar di ekor'
        )
      }
    end
    private_class_method :build_c4_c6

    def build_c7_c9(p5, bit_diffs, total_sec)
      g_sr = p5[:gaussian][4][:fair_spread_bps]
      t_sr = p5[:t4][4][:fair_spread_bps]
      r_sr = p5[:random_recovery][4][:fair_spread_bps]
      {
        c7_model_risk_senior_spreads: crit(
          t_sr > g_sr && r_sr > g_sr,
          { gauss_15_30_bps: g_sr, t4_15_30_bps: t_sr, random_recovery_15_30_bps: r_sr },
          't4_15_30_bps > gauss_15_30_bps dan random_recovery_15_30_bps > gauss_15_30_bps',
          'Kopula berekor berat t(4) dan recovery prosiklis meningkatkan spread senior secara tajam'
        ),
        c8_worker_invariance_1_vs_n: crit(
          bit_diffs.zero?, { bit_diffs: bit_diffs }, 'bit_diffs == 0',
          'Determinisme paralel Process.fork bit-identik 1 vs N worker'
        ),
        c9_runtime_under_360s: crit(
          total_sec < 360.0, { total_seconds: total_sec }, '< 360.0s',
          'Anggaran waktu eksekusi runner Ruby penuh'
        )
      }
    end
    private_class_method :build_c7_c9

    def crit(ok, measured, threshold, note)
      { status: ok ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
