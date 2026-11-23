use std::path::Path;
use std::time::Instant;

use day005c_kelly_drawdown::estimation::{
    optimal_shrinkage_c, run_estimation_loss_study, run_shrinkage_study,
};
use day005c_kelly_drawdown::fractional::{
    build_fractional_curve, build_fractional_map, compute_fractional_roots,
};
use day005c_kelly_drawdown::gbm::GbmParams;
use day005c_kelly_drawdown::heavy_tail::{build_practical_summary, run_heavy_tail_study};
use day005c_kelly_drawdown::hitting::{infinite_horizon_table, run_hitting_grid};
use day005c_kelly_drawdown::jsonout::{obj, write_json_file, JsonVal};
use day005c_kelly_drawdown::kelly_binary::run_binary_study;
use day005c_kelly_drawdown::parallel::default_workers;
use day005c_kelly_drawdown::svg::{write_line_chart, Series};

const MAIN_SEED: u64 = 20261123;

fn criterion(status: &str, measured: JsonVal, threshold: &str, note: &str) -> JsonVal {
    obj([
        ("status", status.into()),
        ("measured", measured),
        ("threshold", threshold.into()),
        ("note", note.into()),
    ])
}

fn main() {
    let t_total = Instant::now();
    let n_workers = default_workers();
    println!("=== Day 005c: Kriteria Kelly & Risiko Drawdown (workers = {n_workers}) ===");

    // -------------------------------------------------------------------------
    // Bagian 1: Taruhan biner diskret
    // -------------------------------------------------------------------------
    let t1 = Instant::now();
    let fs_bin = [0.025, 0.05, 0.10, 0.15, 0.20, 0.30];
    let bin_res = run_binary_study(0.55, 1.0, 1000, 100_000, &fs_bin, MAIN_SEED, n_workers);
    let dt_p1 = t1.elapsed().as_secs_f64();

    let mut max_bin_z_log = 0.0_f64;
    let mut max_bin_med_rel = 0.0_f64;
    for pt in &bin_res.points {
        let z = (pt.e_ln_w_mc - pt.e_ln_w_exact).abs() / pt.e_ln_w_se;
        max_bin_z_log = max_bin_z_log.max(z);
        let rel = (pt.median_w_mc - pt.median_w_exact).abs() / pt.median_w_exact.max(1e-12);
        max_bin_med_rel = max_bin_med_rel.max(rel);
    }
    let f_star_bin_err = (bin_res.f_star_numeric - bin_res.f_star_exact).abs();
    println!(
        "[Bagian 1] selesai dalam {dt_p1:.2}s | f* = {:.6} (num diff = {f_star_bin_err:.2e}), f_be = {:.6}, max_z_log = {max_bin_z_log:.2}",
        bin_res.f_star_exact, bin_res.f_breakeven_numeric
    );

    // -------------------------------------------------------------------------
    // Bagian 2: Pasar kontinu GBM & peluang sentuh eksak + jembatan Brownian
    // -------------------------------------------------------------------------
    let t2 = Instant::now();
    let gbm = GbmParams::default();
    let dts = [1.0 / 12.0, 1.0 / 52.0, 1.0 / 252.0];
    let ts = [1.0, 5.0, 20.0];
    let xs = [0.5, 0.25];
    let fs_gbm = [0.5, 1.0, 1.25, 2.0];
    let hit_grid = run_hitting_grid(
        &gbm, &dts, &ts, &xs, &fs_gbm, 1_000_000, MAIN_SEED, n_workers,
    );
    let cs_inf = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0];
    let inf_table = infinite_horizon_table(&gbm, &cs_inf, 0.5);
    let dt_p2 = t2.elapsed().as_secs_f64();

    let n_within_3se = hit_grid.iter().filter(|r| r.bridge_z.abs() <= 3.0).count();
    let frac_within_3se = n_within_3se as f64 / hit_grid.len() as f64;
    let max_bridge_abs_err = hit_grid
        .iter()
        .map(|r| (r.bridge_mc - r.exact_prob).abs())
        .fold(0.0_f64, f64::max);
    let max_bridge_z = hit_grid
        .iter()
        .map(|r| r.bridge_z.abs())
        .fold(0.0_f64, f64::max);
    let min_naive_bias_monthly = hit_grid
        .iter()
        .filter(|r| (r.dt - 1.0 / 12.0).abs() < 1e-9)
        .map(|r| r.naive_bias)
        .fold(0.0_f64, f64::max);
    println!(
        "[Bagian 2] selesai dalam {dt_p2:.2}s | {n_within_3se}/{} bridge dalam 3 SE (max |z| = {max_bridge_z:.2}, max |err| = {max_bridge_abs_err:.5})",
        hit_grid.len()
    );

    // -------------------------------------------------------------------------
    // Bagian 3: Kelly fraksional
    // -------------------------------------------------------------------------
    let t3 = Instant::now();
    let frac_curve = build_fractional_curve(&gbm, &[0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]);
    let frac_map = build_fractional_map(&gbm, &cs_inf, &ts);
    let frac_roots = compute_fractional_roots(&gbm);
    let dt_p3 = t3.elapsed().as_secs_f64();

    let max_excess_ratio_err = frac_curve
        .iter()
        .map(|r| (r.excess_ratio_measured - r.excess_ratio_formula).abs())
        .fold(0.0_f64, f64::max);
    let root_max_err = (frac_roots.c_zero_excess_closed - frac_roots.c_zero_excess_numeric)
        .abs()
        .max((frac_roots.c_hit_half_eq_half_closed - frac_roots.c_hit_half_eq_half_numeric).abs())
        .max(
            (frac_roots.c_hit_half_eq_quarter_closed - frac_roots.c_hit_half_eq_quarter_numeric)
                .abs(),
        );
    println!(
        "[Bagian 3] selesai dalam {dt_p3:.2}s | max excess ratio diff = {max_excess_ratio_err:.2e}, root diff = {root_max_err:.2e}"
    );

    // -------------------------------------------------------------------------
    // Bagian 4: Galat estimasi, penyusutan, dan ekor berat Student-t
    // -------------------------------------------------------------------------
    let t4 = Instant::now();
    let alt_gbm = GbmParams {
        mu: 0.13,
        sigma: 0.35,
        r: 0.01,
    };
    let t_ests = [1.0, 5.0, 10.0, 20.0];
    let est_loss =
        run_estimation_loss_study(&gbm, &alt_gbm, &t_ests, 1_000_000, MAIN_SEED, n_workers);
    let shrink = run_shrinkage_study(&gbm, 5.0, 20.0, 200_000, MAIN_SEED, n_workers);
    let heavy_tail = run_heavy_tail_study(
        &gbm,
        4.0,
        &[0.5, 1.0, 2.0],
        20.0,
        100_000,
        MAIN_SEED,
        n_workers,
    );
    let practical_summary = build_practical_summary(
        &gbm,
        &[0.25, 0.5, 1.0, 1.5, 2.0],
        5.0,
        200_000,
        MAIN_SEED,
        n_workers,
    );
    let dt_p4 = t4.elapsed().as_secs_f64();

    let max_est_z = est_loss
        .iter()
        .map(|r| r.z_diff.abs())
        .fold(0.0_f64, f64::max);
    let max_est_alt_z = est_loss
        .iter()
        .map(|r| (r.loss_alt_params_mc - r.loss_exact).abs() / r.loss_alt_params_se)
        .fold(0.0_f64, f64::max);
    println!(
        "[Bagian 4] selesai dalam {dt_p4:.2}s | max est z = {max_est_z:.2} (alt = {max_est_alt_z:.2}), c_shrink = {:.4}",
        shrink.c_shrink_opt
    );

    // -------------------------------------------------------------------------
    // Worker invariance check (1 vs n_workers) on small deterministic slice
    // -------------------------------------------------------------------------
    let inv1 = run_binary_study(0.55, 1.0, 200, 1024, &[0.10], MAIN_SEED, 1);
    let inv_n = run_binary_study(0.55, 1.0, 200, 1024, &[0.10], MAIN_SEED, n_workers);
    let worker_bit_diff =
        (inv1.points[0].e_ln_w_mc.to_bits() != inv_n.points[0].e_ln_w_mc.to_bits()) as usize;

    // -------------------------------------------------------------------------
    // SVG Figures
    // -------------------------------------------------------------------------
    let fig_dir = Path::new("outputs/figures");
    // 1. Binary g(f)
    let f_grid: Vec<f64> = (1..=35).map(|i| i as f64 * 0.01).collect();
    let g_grid: Vec<f64> = f_grid
        .iter()
        .map(|&f| day005c_kelly_drawdown::kelly_binary::g_binary(0.55, 1.0, f))
        .collect();
    let zero_line: Vec<f64> = vec![0.0; f_grid.len()];
    write_line_chart(
        &fig_dir.join("day005c_binary_growth.svg"),
        "Laju Pertumbuhan Log Biner g(f) (p = 0.55, b = 1)",
        "Fraksi taruhan f",
        "g(f)",
        &[
            Series {
                label: "g(f)",
                x: &f_grid,
                y: &g_grid,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "g = 0",
                x: &f_grid,
                y: &zero_line,
                color: "#9ca3af",
                dashed: true,
            },
        ],
    )
    .expect("gagal menulis SVG binary_growth");

    // 2. Sample wealth paths (log10 W)
    let n_x: Vec<f64> = bin_res.sample_paths_n.iter().map(|&n| n as f64).collect();
    let log_half: Vec<f64> = bin_res
        .sample_path_half_kelly
        .iter()
        .map(|&w| w.log10())
        .collect();
    let log_full: Vec<f64> = bin_res
        .sample_path_full_kelly
        .iter()
        .map(|&w| w.log10())
        .collect();
    let log_over: Vec<f64> = bin_res
        .sample_path_over_kelly
        .iter()
        .map(|&w| w.log10())
        .collect();
    write_line_chart(
        &fig_dir.join("day005c_wealth_paths.svg"),
        "Lintasan Kekayaan Biner (log10 W_n, N = 1000)",
        "Langkah taruhan n",
        "log10(W_n / W_0)",
        &[
            Series {
                label: "0.5 f*",
                x: &n_x,
                y: &log_half,
                color: "#059669",
                dashed: false,
            },
            Series {
                label: "1.0 f*",
                x: &n_x,
                y: &log_full,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "2.0 f*",
                x: &n_x,
                y: &log_over,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .expect("gagal menulis SVG wealth_paths");

    // 3. Infinite-horizon hitting prob of 0.5 vs c
    let c_inf_x: Vec<f64> = inf_table.iter().map(|r| r.c).collect();
    let p_inf_y: Vec<f64> = inf_table.iter().map(|r| r.prob_hit_half).collect();
    write_line_chart(
        &fig_dir.join("day005c_hitting_vs_c.svg"),
        "Peluang Sentuh 0.5*W0 (Horizon Tak Hingga) vs c = f / f*",
        "Rasio fraksi Kelly c = f / f*",
        "P(inf W_t <= 0.5 W_0)",
        &[Series {
            label: "x^(2g/s^2)",
            x: &c_inf_x,
            y: &p_inf_y,
            color: "#dc2626",
            dashed: false,
        }],
    )
    .expect("gagal menulis SVG hitting_vs_c");

    // 4. Discretization bias vs dt for (T=5, x=0.5, f=1.25)
    let bias_pts: Vec<_> = hit_grid
        .iter()
        .filter(|r| {
            (r.t - 5.0).abs() < 1e-9 && (r.x - 0.5).abs() < 1e-9 && (r.f - 1.25).abs() < 1e-9
        })
        .collect();
    let dt_x: Vec<f64> = bias_pts.iter().map(|r| r.dt * 252.0).collect();
    let naive_y: Vec<f64> = bias_pts.iter().map(|r| r.naive_mc).collect();
    let bridge_y: Vec<f64> = bias_pts.iter().map(|r| r.bridge_mc).collect();
    let exact_y: Vec<f64> = bias_pts.iter().map(|r| r.exact_prob).collect();
    write_line_chart(
        &fig_dir.join("day005c_discretization_bias.svg"),
        "Bias Diskretisasi vs Jembatan Brownian (T=5, x=0.5, f=1.25)",
        "Langkah waktu (hari perdagangan)",
        "Peluang sentuh",
        &[
            Series {
                label: "Eksak tertutup",
                x: &dt_x,
                y: &exact_y,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "Jembatan Brownian",
                x: &dt_x,
                y: &bridge_y,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "Diskret naif",
                x: &dt_x,
                y: &naive_y,
                color: "#dc2626",
                dashed: false,
            },
        ],
    )
    .expect("gagal menulis SVG discretization_bias");

    // 5. Excess growth ratio c(2 - c)
    let fc_x: Vec<f64> = frac_curve.iter().map(|r| r.c).collect();
    let fc_y: Vec<f64> = frac_curve.iter().map(|r| r.excess_ratio_measured).collect();
    write_line_chart(
        &fig_dir.join("day005c_fractional_tradeoff.svg"),
        "Rasio Pertumbuhan Berlebih (g(c*f*) - r) / (g* - r) = c(2 - c)",
        "c = f / f*",
        "Rasio pertumbuhan berlebih",
        &[Series {
            label: "c(2 - c)",
            x: &fc_x,
            y: &fc_y,
            color: "#059669",
            dashed: false,
        }],
    )
    .expect("gagal menulis SVG fractional_tradeoff");

    // 6. Estimation loss vs T_est
    let te_x: Vec<f64> = est_loss.iter().map(|r| r.t_est).collect();
    let te_ex: Vec<f64> = est_loss.iter().map(|r| r.loss_exact).collect();
    let te_mc: Vec<f64> = est_loss.iter().map(|r| r.loss_mc).collect();
    write_line_chart(
        &fig_dir.join("day005c_estimation_loss.svg"),
        "Kerugian Pertumbuhan Galat Estimasi E[g(f*) - g(f_hat)] = 1/(2 T_est)",
        "Panjang sampel estimasi T_est (tahun)",
        "Kerugian laju pertumbuhan",
        &[
            Series {
                label: "Eksak 1/(2 T_est)",
                x: &te_x,
                y: &te_ex,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "MC (R = 1e6)",
                x: &te_x,
                y: &te_mc,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .expect("gagal menulis SVG estimation_loss");

    // 7. Heavy-tail Student-t(4) vs Gaussian 20y hitting prob
    let ht_x: Vec<f64> = heavy_tail.iter().map(|r| r.c).collect();
    let ht_g: Vec<f64> = heavy_tail
        .iter()
        .map(|r| r.prob_hit_half_gaussian_pred)
        .collect();
    let ht_t: Vec<f64> = heavy_tail.iter().map(|r| r.prob_hit_half_t4_mc).collect();
    write_line_chart(
        &fig_dir.join("day005c_heavy_tail.svg"),
        "Peluang Sentuh 0.5 (20 Tahun): Gaussian vs Student-t (nu = 4)",
        "c = f / f*",
        "P(min W_t <= 0.5 W_0)",
        &[
            Series {
                label: "Prediksi Gaussian",
                x: &ht_x,
                y: &ht_g,
                color: "#2563eb",
                dashed: true,
            },
            Series {
                label: "Student-t (nu=4)",
                x: &ht_x,
                y: &ht_t,
                color: "#dc2626",
                dashed: false,
            },
        ],
    )
    .expect("gagal menulis SVG heavy_tail");

    let total_elapsed = t_total.elapsed().as_secs_f64();
    println!("Total waktu eksekusi: {total_elapsed:.2}s");

    // -------------------------------------------------------------------------
    // Criteria evaluation (a-priori thresholds)
    // -------------------------------------------------------------------------
    let c1_ok = f_star_bin_err < 1e-10 && max_bin_z_log <= 3.0 && max_bin_med_rel < 0.05;
    let c2_ok = frac_within_3se >= 0.95 && min_naive_bias_monthly < 0.0;
    let c3_ok = max_excess_ratio_err < 1e-12 && root_max_err < 1e-10;
    let c4_ok = max_est_z <= 3.0 && max_est_alt_z <= 3.0;
    let c5_ok = shrink.diff_shrink_vs_full > 3.0 * shrink.se_diff_shrink_vs_full
        && shrink.prob_hit_half_20y_shrink < shrink.prob_hit_half_20y_full;
    let ht_full = &heavy_tail[1]; // c = 1.0
    let ht_over = &heavy_tail[2]; // c = 2.0
    let c6_ok = ht_full.prob_hit_half_t4_mc
        > ht_full.prob_hit_half_gaussian_pred + 3.0 * ht_full.prob_hit_half_t4_se
        && ht_over.g_t4_mc < ht_over.g_gaussian_pred;
    let c7_ok = worker_bit_diff == 0;
    let c8_ok = total_elapsed < 180.0;

    let criteria = obj([
        (
            "c1_binary_kelly_exact_and_mc",
            criterion(
                if c1_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("f_star_err", f_star_bin_err.into()),
                    ("max_z_log", max_bin_z_log.into()),
                    ("max_median_rel_err", max_bin_med_rel.into()),
                ]),
                "f_star_err < 1e-10, max_z_log <= 3.0, max_median_rel_err < 0.05",
                "Verifikasi fraksi Kelly biner, E[ln W_N], dan median W_N pada N=1000, R=100k",
            ),
        ),
        (
            "c2_gbm_hitting_bridge_unbiased_3se",
            criterion(
                if c2_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("frac_within_3se", frac_within_3se.into()),
                    ("max_bridge_z", max_bridge_z.into()),
                    ("max_bridge_abs_err", max_bridge_abs_err.into()),
                    ("max_naive_bias_monthly", min_naive_bias_monthly.into()),
                ]),
                "frac_within_3se >= 0.95 (72 kombinasi), max_naive_bias_monthly < 0",
                "Estimator jembatan Brownian tanpa bias vs rumus tertutup dan bias negatif naif pada R=1e6",
            ),
        ),
        (
            "c3_fractional_kelly_identity_and_roots",
            criterion(
                if c3_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("max_excess_ratio_err", max_excess_ratio_err.into()),
                    ("root_max_err", root_max_err.into()),
                ]),
                "max_excess_ratio_err < 1e-12, root_max_err < 1e-10",
                "Identitas c(2-c) dan kesesuaian akar numerik vs tertutup",
            ),
        ),
        (
            "c4_estimation_loss_1_over_2t_est",
            criterion(
                if c4_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("max_z_default", max_est_z.into()),
                    ("max_z_alt_params", max_est_alt_z.into()),
                ]),
                "max_z_default <= 3.0, max_z_alt_params <= 3.0",
                "Kerugian pertumbuhan 1/(2 T_est) pada R=1e6 dan invariansi terhadap (mu, sigma)",
            ),
        ),
        (
            "c5_shrinkage_improves_growth_and_drawdown",
            criterion(
                if c5_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("c_shrink_opt", shrink.c_shrink_opt.into()),
                    ("diff_shrink_vs_full", shrink.diff_shrink_vs_full.into()),
                    ("se_diff_shrink_vs_full", shrink.se_diff_shrink_vs_full.into()),
                    ("prob_hit_half_20y_full", shrink.prob_hit_half_20y_full.into()),
                    ("prob_hit_half_20y_shrink", shrink.prob_hit_half_20y_shrink.into()),
                ]),
                "diff_shrink_vs_full > 3 SE dan prob_hit_half_20y_shrink < prob_hit_half_20y_full",
                "Penyusutan optimal mengungguli plug-in penuh pada pertumbuhan dan drawdown 20 tahun",
            ),
        ),
        (
            "c6_heavy_tail_t4_exceeds_gaussian_drawdown",
            criterion(
                if c6_ok { "tercapai" } else { "belum tercapai" },
                obj([
                    ("prob_hit_half_t4_fstar", ht_full.prob_hit_half_t4_mc.into()),
                    ("prob_hit_half_gauss_fstar", ht_full.prob_hit_half_gaussian_pred.into()),
                    ("g_t4_2fstar", ht_over.g_t4_mc.into()),
                    ("g_gauss_2fstar", ht_over.g_gaussian_pred.into()),
                ]),
                "P_t4(hit 0.5 | f*) > P_gauss + 3 SE dan g_t4(2f*) < g_gauss(2f*)",
                "Ekor berat Student-t(4) mematahkan prediksi risiko dan pertumbuhan Gaussian",
            ),
        ),
        (
            "c7_worker_invariance_1_vs_n",
            criterion(
                if c7_ok { "tercapai" } else { "belum tercapai" },
                obj([("bit_diffs", worker_bit_diff.into())]),
                "bit_diffs == 0",
                "Reduksi blok tetap menghasilkan keluaran bit-identik 1 vs N worker",
            ),
        ),
        (
            "c8_runtime_under_180s",
            criterion(
                if c8_ok { "tercapai" } else { "belum tercapai" },
                obj([("total_seconds", total_elapsed.into())]),
                "< 180.0s",
                "Anggaran waktu eksekusi runner penuh",
            ),
        ),
    ]);

    // -------------------------------------------------------------------------
    // Write results.json
    // -------------------------------------------------------------------------
    let _ = optimal_shrinkage_c(&gbm, 5.0);
    let root_json = obj([
        (
            "metadata",
            obj([
                ("day", 5_usize.into()),
                ("session", "c".into()),
                ("language", "rust".into()),
                ("main_seed", (MAIN_SEED as usize).into()),
                ("workers_used", n_workers.into()),
                ("run_completed_without_errors", true.into()),
                (
                    "elapsed_seconds",
                    obj([
                        ("part1_binary", dt_p1.into()),
                        ("part2_gbm_hitting", dt_p2.into()),
                        ("part3_fractional", dt_p3.into()),
                        ("part4_estimation_heavy_tail", dt_p4.into()),
                        ("total", total_elapsed.into()),
                    ]),
                ),
            ]),
        ),
        ("criteria", criteria),
        (
            "part1_binary_kelly",
            obj([
                ("p", bin_res.p.into()),
                ("b", bin_res.b.into()),
                ("f_star_exact", bin_res.f_star_exact.into()),
                ("f_star_numeric", bin_res.f_star_numeric.into()),
                ("g_star_exact", bin_res.g_star_exact.into()),
                ("f_breakeven_numeric", bin_res.f_breakeven_numeric.into()),
                (
                    "points",
                    JsonVal::Arr(
                        bin_res
                            .points
                            .iter()
                            .map(|pt| {
                                obj([
                                    ("f", pt.f.into()),
                                    ("g_exact", pt.g_exact.into()),
                                    ("e_ln_w_mc", pt.e_ln_w_mc.into()),
                                    ("e_ln_w_se", pt.e_ln_w_se.into()),
                                    ("e_ln_w_exact", pt.e_ln_w_exact.into()),
                                    ("median_w_mc", pt.median_w_mc.into()),
                                    ("median_w_exact", pt.median_w_exact.into()),
                                    ("arith_mean_mc", pt.arith_mean_mc.into()),
                                    ("arith_mean_se", pt.arith_mean_se.into()),
                                    ("arith_mean_exact", pt.arith_mean_exact.into()),
                                    ("prob_loss_mc", pt.prob_loss_mc.into()),
                                    ("prob_loss_se", pt.prob_loss_se.into()),
                                    ("prob_loss_cont_approx", pt.prob_loss_cont_approx.into()),
                                    ("prob_touch_half_mc", pt.prob_touch_half_mc.into()),
                                    ("prob_touch_half_se", pt.prob_touch_half_se.into()),
                                    (
                                        "prob_touch_half_cont_approx",
                                        pt.prob_touch_half_cont_approx.into(),
                                    ),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "part2_gbm_hitting",
            obj([
                ("mu", gbm.mu.into()),
                ("sigma", gbm.sigma.into()),
                ("r", gbm.r.into()),
                ("f_star", gbm.kelly_fraction().into()),
                ("g_star", gbm.kelly_growth().into()),
                (
                    "infinite_horizon_half_table",
                    JsonVal::Arr(
                        inf_table
                            .iter()
                            .map(|r| {
                                obj([
                                    ("c", r.c.into()),
                                    ("f", r.f.into()),
                                    ("g", r.g.into()),
                                    ("s", r.s.into()),
                                    ("two_g_over_s2", r.two_g_over_s2.into()),
                                    ("prob_hit_half", r.prob_hit_half.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "grid",
                    JsonVal::Arr(
                        hit_grid
                            .iter()
                            .map(|r| {
                                obj([
                                    ("dt", r.dt.into()),
                                    ("t", r.t.into()),
                                    ("x", r.x.into()),
                                    ("f", r.f.into()),
                                    ("exact_prob", r.exact_prob.into()),
                                    ("naive_mc", r.naive_mc.into()),
                                    ("naive_se", r.naive_se.into()),
                                    ("naive_bias", r.naive_bias.into()),
                                    ("bridge_mc", r.bridge_mc.into()),
                                    ("bridge_se", r.bridge_se.into()),
                                    ("bridge_z", r.bridge_z.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "part3_fractional_kelly",
            obj([
                (
                    "roots",
                    obj([
                        (
                            "c_zero_excess_closed",
                            frac_roots.c_zero_excess_closed.into(),
                        ),
                        (
                            "c_zero_excess_numeric",
                            frac_roots.c_zero_excess_numeric.into(),
                        ),
                        (
                            "c_hit_half_eq_half_closed",
                            frac_roots.c_hit_half_eq_half_closed.into(),
                        ),
                        (
                            "c_hit_half_eq_half_numeric",
                            frac_roots.c_hit_half_eq_half_numeric.into(),
                        ),
                        (
                            "c_hit_half_eq_quarter_closed",
                            frac_roots.c_hit_half_eq_quarter_closed.into(),
                        ),
                        (
                            "c_hit_half_eq_quarter_numeric",
                            frac_roots.c_hit_half_eq_quarter_numeric.into(),
                        ),
                    ]),
                ),
                (
                    "curve",
                    JsonVal::Arr(
                        frac_curve
                            .iter()
                            .map(|r| {
                                obj([
                                    ("c", r.c.into()),
                                    ("f", r.f.into()),
                                    ("g", r.g.into()),
                                    ("s", r.s.into()),
                                    ("excess_ratio_measured", r.excess_ratio_measured.into()),
                                    ("excess_ratio_formula", r.excess_ratio_formula.into()),
                                    ("prob_hit_half_inf", r.prob_hit_half_inf.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "map_c_t",
                    JsonVal::Arr(
                        frac_map
                            .iter()
                            .map(|r| {
                                obj([
                                    ("c", r.c.into()),
                                    ("t", r.t.into()),
                                    ("median_wealth", r.median_wealth.into()),
                                    ("p05_wealth", r.p05_wealth.into()),
                                    ("sd_log_wealth", r.sd_log_wealth.into()),
                                    ("prob_loss", r.prob_loss.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
        (
            "part4_estimation_and_heavy_tail",
            obj([
                (
                    "estimation_loss",
                    JsonVal::Arr(
                        est_loss
                            .iter()
                            .map(|r| {
                                obj([
                                    ("t_est", r.t_est.into()),
                                    ("loss_exact", r.loss_exact.into()),
                                    ("loss_mc", r.loss_mc.into()),
                                    ("loss_se", r.loss_se.into()),
                                    ("z_diff", r.z_diff.into()),
                                    ("loss_alt_params_mc", r.loss_alt_params_mc.into()),
                                    ("loss_alt_params_se", r.loss_alt_params_se.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "shrinkage",
                    obj([
                        ("t_est", shrink.t_est.into()),
                        ("c_shrink_opt", shrink.c_shrink_opt.into()),
                        ("horizon_years", shrink.horizon_years.into()),
                        ("g_full_plugin", shrink.g_full_plugin.into()),
                        ("se_g_full_plugin", shrink.se_g_full_plugin.into()),
                        ("g_half_kelly", shrink.g_half_kelly.into()),
                        ("se_g_half_kelly", shrink.se_g_half_kelly.into()),
                        ("g_shrink_opt", shrink.g_shrink_opt.into()),
                        ("se_g_shrink_opt", shrink.se_g_shrink_opt.into()),
                        ("diff_shrink_vs_full", shrink.diff_shrink_vs_full.into()),
                        (
                            "se_diff_shrink_vs_full",
                            shrink.se_diff_shrink_vs_full.into(),
                        ),
                        ("diff_half_vs_full", shrink.diff_half_vs_full.into()),
                        ("se_diff_half_vs_full", shrink.se_diff_half_vs_full.into()),
                        (
                            "prob_hit_half_20y_full",
                            shrink.prob_hit_half_20y_full.into(),
                        ),
                        (
                            "prob_hit_half_20y_half",
                            shrink.prob_hit_half_20y_half.into(),
                        ),
                        (
                            "prob_hit_half_20y_shrink",
                            shrink.prob_hit_half_20y_shrink.into(),
                        ),
                    ]),
                ),
                (
                    "heavy_tail_t4",
                    JsonVal::Arr(
                        heavy_tail
                            .iter()
                            .map(|r| {
                                obj([
                                    ("c", r.c.into()),
                                    ("f", r.f.into()),
                                    ("g_gaussian_pred", r.g_gaussian_pred.into()),
                                    ("g_t4_mc", r.g_t4_mc.into()),
                                    ("g_t4_se", r.g_t4_se.into()),
                                    (
                                        "prob_hit_half_gaussian_pred",
                                        r.prob_hit_half_gaussian_pred.into(),
                                    ),
                                    ("prob_hit_half_t4_mc", r.prob_hit_half_t4_mc.into()),
                                    ("prob_hit_half_t4_se", r.prob_hit_half_t4_se.into()),
                                    ("ruin_prob_t4_mc", r.ruin_prob_t4_mc.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
                (
                    "practical_summary",
                    JsonVal::Arr(
                        practical_summary
                            .iter()
                            .map(|r| {
                                obj([
                                    ("c", r.c.into()),
                                    ("f", r.f.into()),
                                    ("g_no_err", r.g_no_err.into()),
                                    ("g_with_err", r.g_with_err.into()),
                                    (
                                        "prob_hit_half_inf_no_err",
                                        r.prob_hit_half_inf_no_err.into(),
                                    ),
                                    (
                                        "prob_hit_half_20y_no_err",
                                        r.prob_hit_half_20y_no_err.into(),
                                    ),
                                    (
                                        "prob_hit_half_20y_with_err",
                                        r.prob_hit_half_20y_with_err.into(),
                                    ),
                                    ("prob_loss_5y_no_err", r.prob_loss_5y_no_err.into()),
                                    ("prob_loss_5y_with_err", r.prob_loss_5y_with_err.into()),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]),
        ),
    ]);

    write_json_file(Path::new("outputs/day005c/results.json"), &root_json)
        .expect("gagal menulis outputs/day005c/results.json");
}
