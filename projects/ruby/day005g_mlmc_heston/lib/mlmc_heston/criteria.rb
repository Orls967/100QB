# frozen_string_literal: true

module MlmcHeston
  module Criteria
    module_function

    def evaluate(p1, p2, p3, p4, p4_mm, bit_diffs, total_sec)
      build_c1_c4(p1, p2, p3).merge(build_c5_c8(p3, p4, p4_mm, bit_diffs, total_sec))
    end

    def build_c1_c4(p1, p2, p3)
      a = p1[:anchors]
      eu_b = p2[:gbm_euler][:beta_strong]
      mi_b = p2[:gbm_milstein][:beta_strong]
      mi_a = p2[:gbm_milstein][:alpha_weak]
      gbm_ok = (mi_b - 2.0).abs < 0.15 && (eu_b - 1.0).abs < 0.15 && (mi_a - 1.0).abs < 0.15
      fok_b = p2[:heston_feller_ok][:beta_strong]
      fvi_b = p2[:heston_feller_violated][:beta_strong]
      hest_ok = fok_b >= 1.0 && fvi_b.between?(0.5, 0.95) && (fok_b - fvi_b) >= 0.35
      g_row = p3[:gbm_milstein][:rows].last
      h_row = p3[:heston_european][:rows].last
      z_g = (g_row[:price] - a[:bs_exact]).abs / [g_row[:se], 1.0e-12].max
      z_h = (h_row[:price] - a[:heston_cf_exact]).abs / [h_row[:se], 1.0e-12].max

      {
        c1_anchors_bs_cir_and_heston_cf: crit(
          a[:xi_zero_bs_err] < 1.0e-10 && a[:parity_err] < 1.0e-12 && a[:cf_quad_err] < 1.0e-8,
          { xi_zero_bs_err: a[:xi_zero_bs_err], parity_err: a[:parity_err], cf_quad_err: a[:cf_quad_err] },
          'xi_zero_bs_err < 1e-10, parity_err < 1e-12, cf_quad_err < 1e-8',
          'Heston CF cocok dengan batas BS, paritas put-call, dan konvergensi kuadratur'
        ),
        c2_strong_weak_rates_gbm_milstein_vs_euler: crit(
          gbm_ok,
          { gbm_euler_beta: eu_b, gbm_milstein_beta: mi_b, gbm_milstein_alpha: mi_a },
          '|gbm_milstein_beta - 2.0| < 0.15, |gbm_euler_beta - 1.0| < 0.15, |gbm_milstein_alpha - 1.0| < 0.15',
          'Skema Milstein cocok dengan teori orde kuat beta=2 & alpha=1 vs Euler beta=1'
        ),
        c3_heston_full_truncation_rates_feller_effect: crit(
          hest_ok,
          { feller_ok_beta: fok_b, feller_violated_beta: fvi_b, beta_drop_feller_violation: fok_b - fvi_b },
          'feller_ok_beta >= 1.0, 0.5 <= feller_violated_beta <= 0.95, dan penurunan (ok - violated) >= 0.35',
          'Rezim Feller terpenuhi menjaga beta >= 1 sedangkan pelanggaran Feller menjatuhkan beta < 1'
        ),
        c4_mlmc_matches_exact_anchors: crit(
          z_g <= 3.5 && z_h <= 3.5,
          { gbm_mlmc_price: g_row[:price], bs_exact: a[:bs_exact], gbm_z: z_g,
            heston_mlmc_price: h_row[:price], heston_cf_exact: a[:heston_cf_exact], heston_z: z_h },
          'gbm_z <= 3.5 dan heston_z <= 3.5 pada eps = 0.01',
          'Estimator MLMC adaptif cocok dengan harga referensi tertutup GBM dan Heston'
        )
      }
    end
    private_class_method :build_c1_c4

    def build_c5_c8(p3, p4, p4_mm, bit_diffs, total_sec)
      {
        c5_complexity_scaling_eps2_vs_eps3: build_c5(p3[:gbm_milstein], p3[:heston_european]),
        c6_antithetic_level_variance_reduction: build_c6(p4, p4_mm),
        c7_worker_invariance_1_vs_n: crit(
          bit_diffs.zero?, { bit_diffs: bit_diffs }, 'bit_diffs == 0',
          'Determinisme paralel Process.fork bit-identik 1 vs N worker'
        ),
        c8_runtime_under_360s: crit(
          total_sec < 360.0, { total_seconds: total_sec }, '< 360.0s',
          'Anggaran waktu eksekusi runner Ruby penuh'
        )
      }
    end
    private_class_method :build_c5_c8

    def build_c5(gm, he)
      slope_ok = (gm[:mlmc_cost_slope] - gm[:giles_asymptotic_slope]).abs < 0.20 &&
                 (he[:mlmc_cost_slope] - he[:giles_preasymptotic_slope]).abs < 0.15 &&
                 (he[:mlmc_cost_slope] - he[:giles_asymptotic_slope]).abs < 0.35 &&
                 (he[:std_mc_cost_slope] - 3.0).abs < 0.10
      crit(
        slope_ok,
        c5_measured(gm, he),
        '|gbm_slope - 2| < 0.20, |heston_slope - preasymp| < 0.15, |heston_slope - asymp| < 0.35, |std - 3| < 0.10',
        'Kemiringan cocok dengan prediksi Giles; pada eps=0.01 MLMC Heston masih 5.4x lebih lambat dari MC standar'
      )
    end
    private_class_method :build_c5

    def c5_measured(gm, he)
      {
        gbm_milstein_slope: gm[:mlmc_cost_slope], gbm_giles_asymptotic_slope: gm[:giles_asymptotic_slope],
        gbm_eps001_speedup: gm[:rows].last[:speedup],
        heston_mlmc_slope: he[:mlmc_cost_slope], heston_giles_asymptotic_slope: he[:giles_asymptotic_slope],
        heston_giles_preasymptotic_slope: he[:giles_preasymptotic_slope],
        heston_beta_var: he[:beta_var_levels], heston_alpha_bias: he[:alpha_bias_part2],
        std_mc_slope: he[:std_mc_cost_slope], heston_eps001_speedup: he[:rows].last[:speedup],
        heston_v1_to_v0_ratio: he[:v1_to_v0_ratio],
        heston_crossover_eps_extrapolated: he[:crossover_eps_extrapolated]
      }
    end
    private_class_method :c5_measured

    def build_c6(p4, p4_mm)
      l3 = p4.find { |r| r[:level] == 3 }
      r_ceil = p4_mm.first[:heston_theoretical_ceiling]
      mm1 = p4_mm.first
      mm6 = p4_mm.last
      gbm_ok = (mm1[:gbm_euler_ratio] - 1.0).abs < 1.0e-9 && (mm6[:gbm_euler_ratio] - 1.0).abs < 1.0e-9
      cir_ok = mm1[:cir_euler_ratio] > mm6[:cir_euler_ratio] && mm6[:cir_euler_ratio].between?(1.0, 1.04)
      hest_ok = l3[:reduction_ratio].between?(1.15, r_ceil) &&
                mm6[:heston_euler_ratio].between?(1.10, mm1[:heston_euler_ratio])
      crit(
        gbm_ok && cir_ok && hest_ok,
        { level3_heston_ratio: l3[:reduction_ratio], heston_euler_ceiling: r_ceil,
          gbm_l1_l6: [mm1[:gbm_euler_ratio], mm6[:gbm_euler_ratio]],
          cir_l1_l6: [mm1[:cir_euler_ratio], mm6[:cir_euler_ratio]],
          heston_l1_l6: [mm1[:heston_euler_ratio], mm6[:heston_euler_ratio]] },
        'GBM_Euler == 1.0; CIR_Euler meluruh ke 1.0; Heston_Euler dalam [1.10, 2/(1+rho^2)=1.3423]',
        'Turunan teori Euler 1D (GBM=1, CIR->1) dan 2D Heston (batas atas 2/(1+rho^2)=1.3423) cocok di l=1..6'
      )
    end
    private_class_method :build_c6
    private_class_method :build_c5_c8

    def crit(ok, measured, threshold, note)
      { status: ok ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
