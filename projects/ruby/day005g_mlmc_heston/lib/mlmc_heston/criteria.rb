# frozen_string_literal: true

module MlmcHeston
  module Criteria
    module_function

    def evaluate(p1, p2, p3, p4, bit_diffs, total_sec)
      build_c1_c4(p1, p2, p3).merge(build_c5_c8(p3, p4, bit_diffs, total_sec))
    end

    def build_c1_c4(p1, p2, p3)
      a = p1[:anchors]
      gbm_ok = p2[:gbm_milstein][:beta_strong] > 1.5 &&
               p2[:gbm_milstein][:beta_strong] > p2[:gbm_euler][:beta_strong]
      hest_ok = p2[:heston_feller_ok][:beta_strong] > 0.7 && p2[:heston_feller_violated][:beta_strong] > 0.7
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
          { gbm_euler_beta: p2[:gbm_euler][:beta_strong], gbm_milstein_beta: p2[:gbm_milstein][:beta_strong],
            gbm_milstein_alpha: p2[:gbm_milstein][:alpha_weak] },
          'gbm_milstein_beta > 1.5 dan gbm_milstein_beta > gbm_euler_beta',
          'Skema Milstein mencapai peluruhan variansi orde ~2 vs orde ~1 Euler'
        ),
        c3_heston_full_truncation_rates_feller_effect: crit(
          hest_ok,
          { feller_ok_beta: p2[:heston_feller_ok][:beta_strong],
            feller_violated_beta: p2[:heston_feller_violated][:beta_strong] },
          'beta_strong > 0.7 pada rezim Feller terpenuhi dan Feller dilanggar',
          'Full-truncation Euler menjaga laju kuat pada rezim Feller terpenuhi maupun dilanggar'
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

    def build_c5_c8(p3, p4, bit_diffs, total_sec)
      gm_s = p3[:gbm_milstein][:mlmc_cost_slope]
      he_s = p3[:heston_european][:mlmc_cost_slope]
      std_s = p3[:heston_european][:std_mc_cost_slope]
      vr_l3 = p4.find { |r| r[:level] == 3 }[:reduction_ratio]

      {
        c5_complexity_scaling_eps2_vs_eps3: crit(
          gm_s < 2.4 && he_s < 2.6 && std_s > 2.7,
          { gbm_milstein_slope: gm_s, heston_mlmc_slope: he_s, std_mc_slope: std_s,
            heston_eps001_speedup: p3[:heston_european][:rows].last[:speedup] },
          'gbm_milstein_slope < 2.4, heston_mlmc_slope < 2.6, std_mc_slope > 2.7',
          'Hukum kompleksitas Giles: biaya MLMC ~ eps^-2..2.4 vs MC standar ~ eps^-3'
        ),
        c6_antithetic_level_variance_reduction: crit(
          vr_l3 > 1.0,
          { level3_var_standard: p4.find { |r| r[:level] == 3 }[:var_standard],
            level3_var_antithetic: p4.find { |r| r[:level] == 3 }[:var_antithetic],
            level3_reduction_ratio: vr_l3 },
          'level3_reduction_ratio > 1.0',
          'Koppel antitetik Giles-Szpruch menekan variansi level Heston'
        ),
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

    def crit(ok, measured, threshold, note)
      { status: ok ? 'tercapai' : 'belum tercapai', measured: measured, threshold: threshold, note: note }
    end
    private_class_method :crit
  end
end
