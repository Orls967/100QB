use std::path::Path;
use std::time::Instant;

use day005a_sabr::black::{black_call, black_implied_vol, black_vega};
use day005a_sabr::calib::calibrate_sabr_smile;
use day005a_sabr::cev::{
    cev_absorb_prob, cev_call_exact, cev_call_pde_cn, cev_mc_euler, cev_mc_exact,
};
use day005a_sabr::density::{hagan_digital_and_butterfly, scan_hagan_negative_density};
use day005a_sabr::hagan::hagan_implied_vol;
use day005a_sabr::hedge::{run_sabr_hedging_crn, sabr_bartlett_delta, sabr_naive_delta};
use day005a_sabr::jsonout::{obj, write_json_file, JsonVal};
use day005a_sabr::sabr_mc::{default_sabr_workers, sabr_cond_rho0_mc, sabr_euler_mc, EulerScheme};
use day005a_sabr::svg::{write_bar_chart, write_line_chart, Series};

const MAIN_SEED: u64 = 20261123;

fn repo_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

fn main() {
    let t_total = Instant::now();
    let n_workers = default_sabr_workers();
    let root = repo_root();
    let out_json = root.join("outputs/day005a/results.json");
    let fig_dir = root.join("outputs/figures");

    let f0 = 0.04_f64;
    let alpha0 = 0.0325_f64;
    let beta0 = 0.5_f64;
    let nu0 = 0.4_f64;
    let rho0 = -0.3_f64;
    let ts = [1.0_f64, 5.0, 10.0];
    let ks = [0.005_f64, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.10];

    // =========================================================================
    // PART 1: CEV (nu = 0) — Exact vs PDE vs MC + Absorption & Martingale
    // =========================================================================
    let t_p1 = Instant::now();
    let mut cev_3way_rows = Vec::new();
    let mut max_pde_rel_err = 0.0_f64;
    let mut max_mc_exact_z = 0.0_f64;
    let mut max_mc_euler_z = 0.0_f64;

    for &b in &[0.3_f64, 0.5, 0.7] {
        for &k in &[0.02_f64, 0.04, 0.06] {
            let t_h = 5.0_f64;
            let c_ex = cev_call_exact(f0, k, alpha0, b, t_h);
            let p1 = cev_call_pde_cn(f0, k, alpha0, b, t_h, 1000, 600);
            let p2 = cev_call_pde_cn(f0, k, alpha0, b, t_h, 2000, 1200);
            let c_pde = (4.0 * p2 - p1) / 3.0;
            let rel_pde = (c_ex - c_pde).abs() / c_ex.max(1e-12);
            max_pde_rel_err = max_pde_rel_err.max(rel_pde);

            let (c_mce, se_mce, _, _, _) =
                cev_mc_exact(f0, k, alpha0, b, t_h, 200_000, n_workers, MAIN_SEED);
            let z_mce = (c_mce - c_ex).abs() / se_mce.max(1e-15);
            max_mc_exact_z = max_mc_exact_z.max(z_mce);

            let (c_eul, se_eul, _, _) =
                cev_mc_euler(f0, k, alpha0, b, t_h, 1000, 100_000, n_workers, MAIN_SEED);
            let z_eul = (c_eul - c_ex).abs() / se_eul.max(1e-15);
            max_mc_euler_z = max_mc_euler_z.max(z_eul);

            cev_3way_rows.push((
                b, k, c_ex, c_pde, rel_pde, c_mce, se_mce, z_mce, c_eul, z_eul,
            ));
        }
    }

    let p_abs_t10_exact = cev_absorb_prob(f0, alpha0, beta0, 10.0);
    let p_abs_t10_feller = (-2.0 * f0 / (alpha0 * alpha0 * 10.0)).exp();
    let (_, _, p_abs_t10_mc, mean_ft_t10_mc, _) =
        cev_mc_exact(f0, 0.04, alpha0, beta0, 10.0, 200_000, n_workers, MAIN_SEED);
    let p1_secs = t_p1.elapsed().as_secs_f64();

    // =========================================================================
    // PART 2: Hagan (2002) & nu -> 0 anchor vs CEV exact
    // =========================================================================
    let t_p2 = Instant::now();
    let mut hagan_nu0_rows = Vec::new();
    for &t_h in &ts {
        for &k in &ks {
            let c_cev = cev_call_exact(f0, k, alpha0, beta0, t_h);
            let iv_cev = black_implied_vol(c_cev, f0, k, t_h);
            let iv_hag = hagan_implied_vol(f0, k, t_h, alpha0, beta0, rho0, 1e-4);
            let err_bps = iv_cev.map(|v| (iv_hag - v) * 1e4);
            hagan_nu0_rows.push((t_h, k, iv_cev, iv_hag, err_bps));
        }
    }
    let p2_secs = t_p2.elapsed().as_secs_f64();

    // =========================================================================
    // PART 3: Reference Pricing for nu > 0 (Conditional rho=0 + Euler MC)
    // =========================================================================
    let t_p3 = Instant::now();

    // 3.1 Conditional CEV MC for rho = 0 (R = 10^6)
    let mut cond_rho0_rows = Vec::new();
    for &t_h in &ts {
        let n_simp = (200.0 * t_h) as usize;
        let res = sabr_cond_rho0_mc(
            f0,
            &ks,
            t_h,
            alpha0,
            beta0,
            nu0,
            n_simp,
            1_000_000,
            n_workers,
            MAIN_SEED ^ (t_h.to_bits()),
        );
        cond_rho0_rows.push((t_h, res));
    }

    // 3.2 Full Euler convergence & scheme comparison at rho = 0, T = 1 (n in {100, 200, 400, 800}, R = 2e6)
    let mut euler_conv_rows = Vec::new();
    let cond_t1_atm = cond_rho0_rows[0].1.prices[4];
    let cond_t1_atm_se = cond_rho0_rows[0].1.ses[4];
    for &n_yr in &[100_usize, 200, 400, 800] {
        let e_abs = sabr_euler_mc(
            f0,
            &[0.04],
            1.0,
            alpha0,
            beta0,
            0.0,
            nu0,
            n_yr,
            2_000_000,
            EulerScheme::Absorbing,
            n_workers,
            MAIN_SEED,
        );
        let e_ft = sabr_euler_mc(
            f0,
            &[0.04],
            1.0,
            alpha0,
            beta0,
            0.0,
            nu0,
            n_yr,
            2_000_000,
            EulerScheme::FullTruncation,
            n_workers,
            MAIN_SEED,
        );
        let z_abs = (e_abs.prices[0] - cond_t1_atm).abs()
            / (e_abs.ses[0].powi(2) + cond_t1_atm_se.powi(2)).sqrt();
        let z_ft = (e_ft.prices[0] - cond_t1_atm).abs()
            / (e_ft.ses[0].powi(2) + cond_t1_atm_se.powi(2)).sqrt();
        euler_conv_rows.push((
            n_yr,
            e_abs.prices[0],
            e_abs.ses[0],
            z_abs,
            e_ft.prices[0],
            e_ft.ses[0],
            z_ft,
        ));
    }

    // 3.3 Reference smile at rho = -0.3 (R = 2e6, n = 400 steps/yr)
    let mut ref_smiles = Vec::new();
    for &t_h in &ts {
        let n_st = (400.0 * t_h).round() as usize;
        let mc_res = sabr_euler_mc(
            f0,
            &ks,
            t_h,
            alpha0,
            beta0,
            rho0,
            nu0,
            n_st,
            2_000_000,
            EulerScheme::Absorbing,
            n_workers,
            MAIN_SEED ^ t_h.to_bits(),
        );
        let mut pts = Vec::new();
        for (j, &k) in ks.iter().enumerate() {
            let p_mc = mc_res.prices[j];
            let se_p = mc_res.ses[j];
            let iv_mc = black_implied_vol(p_mc, f0, k, t_h);
            let se_iv_bps = iv_mc.map(|v| {
                let vega = black_vega(f0, k, v, t_h).max(1e-12);
                (se_p / vega) * 1e4
            });
            let iv_hag = hagan_implied_vol(f0, k, t_h, alpha0, beta0, rho0, nu0);
            let p_hag = black_call(f0, k, iv_hag, t_h);
            let err_iv_bps = iv_mc.map(|v| (iv_hag - v) * 1e4);
            pts.push((
                k,
                p_mc * 1e4,
                se_p * 1e4,
                p_hag * 1e4,
                iv_mc,
                se_iv_bps,
                iv_hag,
                err_iv_bps,
            ));
        }
        ref_smiles.push((t_h, mc_res.p_absorb, mc_res.mean_ft, mc_res.se_ft, pts));
    }

    // 3.4 Parameter variation map at T = 5 (nu in {0.2,0.4,0.8}, rho in {-0.6,-0.3,0}, beta in {0.3,0.5,0.8})
    let mut error_map_rows = Vec::new();
    let ks_map = [0.01_f64, 0.02, 0.04, 0.06, 0.08];
    for &nu_v in &[0.2_f64, 0.4, 0.8] {
        for &rho_v in &[-0.6_f64, -0.3, 0.0] {
            for &beta_v in &[0.3_f64, 0.5, 0.8] {
                // Match initial ATM local vol: alpha * F0^beta = alpha0 * F0^beta0
                let a_v = alpha0 * f0.powf(beta0 - beta_v);
                let mc_m = sabr_euler_mc(
                    f0,
                    &ks_map,
                    5.0,
                    a_v,
                    beta_v,
                    rho_v,
                    nu_v,
                    1000,
                    200_000,
                    EulerScheme::Absorbing,
                    n_workers,
                    MAIN_SEED ^ nu_v.to_bits() ^ rho_v.to_bits() ^ beta_v.to_bits(),
                );
                let mut sq_err = 0.0_f64;
                let mut cnt = 0_usize;
                let mut max_abs_bps = 0.0_f64;
                for (j, &k) in ks_map.iter().enumerate() {
                    if let Some(iv_m) = black_implied_vol(mc_m.prices[j], f0, k, 5.0) {
                        let iv_h = hagan_implied_vol(f0, k, 5.0, a_v, beta_v, rho_v, nu_v);
                        let diff_bps = (iv_h - iv_m) * 1e4;
                        sq_err += diff_bps * diff_bps;
                        max_abs_bps = max_abs_bps.max(diff_bps.abs());
                        cnt += 1;
                    }
                }
                let rmse_bps = if cnt > 0 {
                    (sq_err / cnt as f64).sqrt()
                } else {
                    0.0
                };
                error_map_rows.push((nu_v, rho_v, beta_v, rmse_bps, max_abs_bps));
            }
        }
    }
    let p3_secs = t_p3.elapsed().as_secs_f64();

    // =========================================================================
    // PART 4: Negative Implied Density & Butterfly Arbitrage
    // =========================================================================
    let t_p4 = Instant::now();
    let mut density_rows = Vec::new();
    let mut pricing_impact_rows = Vec::new();
    for (ti, &t_h) in ts.iter().enumerate() {
        let scan =
            scan_hagan_negative_density(f0, t_h, alpha0, beta0, rho0, nu0, 0.0005, 0.08, 400);
        density_rows.push((t_h, ref_smiles[ti].1, scan));

        // Compare digital & butterfly at K = 0.005 (floor 0.5%) and K = 0.01
        let dk = 0.001;
        for &k_tgt in &[0.005_f64, 0.01] {
            let (dig_h, but_h) =
                hagan_digital_and_butterfly(f0, k_tgt, dk, t_h, alpha0, beta0, rho0, nu0);
            let mc_3 = sabr_euler_mc(
                f0,
                &[k_tgt - dk, k_tgt, k_tgt + dk],
                t_h,
                alpha0,
                beta0,
                rho0,
                nu0,
                (400.0 * t_h).round() as usize,
                500_000,
                EulerScheme::Absorbing,
                n_workers,
                MAIN_SEED ^ (t_h.to_bits() ^ k_tgt.to_bits()),
            );
            let dig_mc = (mc_3.prices[0] - mc_3.prices[2]) / (2.0 * dk);
            let but_mc = mc_3.prices[2] - 2.0 * mc_3.prices[1] + mc_3.prices[0];
            pricing_impact_rows.push((t_h, k_tgt, dig_h, dig_mc, but_h * 1e4, but_mc * 1e4));
        }
    }
    let p4_secs = t_p4.elapsed().as_secs_f64();

    // =========================================================================
    // PART 5: Calibration & Bartlett Hedging
    // =========================================================================
    let t_p5 = Instant::now();
    // Calibrate on strikes where MC IV is well-resolved (se_iv_bps < 50 bps)
    let mut calib_fixed_beta = Vec::new();
    let init_guess = [0.045, -0.1, 0.25];
    for (ti, &t_h) in ts.iter().enumerate() {
        let mut ks_fit = Vec::new();
        let mut ivs_fit = Vec::new();
        for pt in &ref_smiles[ti].4 {
            if let (Some(iv), Some(se_bps)) = (pt.4, pt.5) {
                if se_bps < 50.0 {
                    ks_fit.push(pt.0);
                    ivs_fit.push(iv);
                }
            }
        }
        let fit = calibrate_sabr_smile(f0, t_h, beta0, &ks_fit, &ivs_fit, init_guess);
        calib_fixed_beta.push((t_h, fit));
    }

    // Beta non-identifiability at T = 5 across beta in {0.0, 0.5, 1.0}
    let mut calib_beta_sweep_t5 = Vec::new();
    let mut ks_t5 = Vec::new();
    let mut ivs_t5 = Vec::new();
    for pt in &ref_smiles[1].4 {
        if let (Some(iv), Some(se_bps)) = (pt.4, pt.5) {
            if se_bps < 50.0 {
                ks_t5.push(pt.0);
                ivs_t5.push(iv);
            }
        }
    }
    for &b_val in &[0.0_f64, 0.5, 1.0] {
        let a_init = alpha0 * f0.powf(beta0 - b_val) * 1.2;
        let fit = calibrate_sabr_smile(f0, 5.0, b_val, &ks_t5, &ivs_t5, [a_init, -0.1, 0.25]);
        calib_beta_sweep_t5.push(fit);
    }

    // Bartlett hedging at T = 1 (K = 0.04, 252 daily steps, R = 20,000 CRN)
    let fit_t1 = &calib_fixed_beta[0].1;
    let hedge_main = run_sabr_hedging_crn(
        f0,
        f0,
        1.0,
        alpha0,
        beta0,
        rho0,
        nu0,
        fit_t1.alpha,
        beta0,
        fit_t1.rho,
        fit_t1.nu,
        252,
        20_000,
        n_workers,
        MAIN_SEED,
    );

    // Misspecified beta hedging: calibrate beta in {0.3, 0.5, 0.8} on T=1 smile, hedge true beta=0.5 world
    let mut ks_t1 = Vec::new();
    let mut ivs_t1 = Vec::new();
    for pt in &ref_smiles[0].4 {
        if let (Some(iv), Some(se_bps)) = (pt.4, pt.5) {
            if se_bps < 50.0 {
                ks_t1.push(pt.0);
                ivs_t1.push(iv);
            }
        }
    }
    let mut hedge_misspec_beta = Vec::new();
    for &b_val in &[0.3_f64, 0.5, 0.8] {
        let a_init = alpha0 * f0.powf(beta0 - b_val) * 1.2;
        let fit_b = calibrate_sabr_smile(f0, 1.0, b_val, &ks_t1, &ivs_t1, [a_init, -0.1, 0.25]);
        let h_b = run_sabr_hedging_crn(
            f0,
            f0,
            1.0,
            alpha0,
            beta0,
            rho0,
            nu0,
            fit_b.alpha,
            b_val,
            fit_b.rho,
            fit_b.nu,
            252,
            20_000,
            n_workers,
            MAIN_SEED,
        );
        let d_bart = sabr_bartlett_delta(f0, f0, 1.0, fit_b.alpha, b_val, fit_b.rho, fit_b.nu);
        let d_naiv = sabr_naive_delta(f0, f0, 1.0, fit_b.alpha, b_val, fit_b.rho, fit_b.nu);
        hedge_misspec_beta.push((
            b_val,
            fit_b.rmse_iv_bps,
            d_bart,
            d_naiv,
            h_b.sd_pnl_bps[0],
            h_b.sd_pnl_bps[1],
        ));
    }
    let p5_secs = t_p5.elapsed().as_secs_f64();

    // =========================================================================
    // SVG Figures (7 charts)
    // =========================================================================
    let k_cev_x: Vec<f64> = cev_3way_rows
        .iter()
        .filter(|r| (r.0 - 0.5).abs() < 1e-9)
        .map(|r| r.1 * 100.0)
        .collect();
    let y_cev_ex: Vec<f64> = cev_3way_rows
        .iter()
        .filter(|r| (r.0 - 0.5).abs() < 1e-9)
        .map(|r| r.2 * 1e4)
        .collect();
    let y_cev_pde: Vec<f64> = cev_3way_rows
        .iter()
        .filter(|r| (r.0 - 0.5).abs() < 1e-9)
        .map(|r| r.3 * 1e4)
        .collect();
    let y_cev_mc: Vec<f64> = cev_3way_rows
        .iter()
        .filter(|r| (r.0 - 0.5).abs() < 1e-9)
        .map(|r| r.5 * 1e4)
        .collect();
    write_line_chart(
        &fig_dir.join("day005a_cev-three-ways.svg"),
        "CEV (beta=0.5, T=5): Exact vs Crank-Nicolson PDE vs Exact MC (bps)",
        "Strike K (%)",
        "Call Price (bps)",
        &[
            Series {
                label: "Exact BESQ",
                x: &k_cev_x,
                y: &y_cev_ex,
                color: "#111827",
                dashed: false,
            },
            Series {
                label: "PDE Crank-Nicolson",
                x: &k_cev_x,
                y: &y_cev_pde,
                color: "#2563eb",
                dashed: true,
            },
            Series {
                label: "Exact MC",
                x: &k_cev_x,
                y: &y_cev_mc,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let k_pct: Vec<f64> = ks.iter().map(|&k| k * 100.0).collect();
    let nu0_err_t1: Vec<f64> = hagan_nu0_rows
        .iter()
        .filter(|r| (r.0 - 1.0).abs() < 1e-9)
        .map(|r| r.4.unwrap_or(0.0))
        .collect();
    let nu0_err_t5: Vec<f64> = hagan_nu0_rows
        .iter()
        .filter(|r| (r.0 - 5.0).abs() < 1e-9)
        .map(|r| r.4.unwrap_or(0.0))
        .collect();
    let nu0_err_t10: Vec<f64> = hagan_nu0_rows
        .iter()
        .filter(|r| (r.0 - 10.0).abs() < 1e-9)
        .map(|r| r.4.unwrap_or(0.0))
        .collect();
    write_line_chart(
        &fig_dir.join("day005a_hagan-nu0-error.svg"),
        "Hagan (nu->0) Implied Vol Error vs Exact CEV (bps)",
        "Strike K (%)",
        "Vol Error (Hagan - CEV, bps)",
        &[
            Series {
                label: "T = 1y",
                x: &k_pct,
                y: &nu0_err_t1,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "T = 5y",
                x: &k_pct,
                y: &nu0_err_t5,
                color: "#f59e0b",
                dashed: false,
            },
            Series {
                label: "T = 10y",
                x: &k_pct,
                y: &nu0_err_t10,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let smile_t1_mc: Vec<f64> = ref_smiles[0]
        .4
        .iter()
        .map(|p| p.4.unwrap_or(0.0) * 100.0)
        .collect();
    let smile_t1_hag: Vec<f64> = ref_smiles[0].4.iter().map(|p| p.6 * 100.0).collect();
    let smile_t10_mc: Vec<f64> = ref_smiles[2]
        .4
        .iter()
        .map(|p| p.4.unwrap_or(0.0) * 100.0)
        .collect();
    let smile_t10_hag: Vec<f64> = ref_smiles[2].4.iter().map(|p| p.6 * 100.0).collect();
    write_line_chart(
        &fig_dir.join("day005a_smile-ref-vs-hagan.svg"),
        "SABR Implied Vol Smile: MC Reference (R=2M) vs Hagan (T=1y & T=10y)",
        "Strike K (%)",
        "Black Implied Vol (%)",
        &[
            Series {
                label: "MC Ref T=1y",
                x: &k_pct,
                y: &smile_t1_mc,
                color: "#111827",
                dashed: false,
            },
            Series {
                label: "Hagan T=1y",
                x: &k_pct,
                y: &smile_t1_hag,
                color: "#2563eb",
                dashed: true,
            },
            Series {
                label: "MC Ref T=10y",
                x: &k_pct,
                y: &smile_t10_mc,
                color: "#dc2626",
                dashed: false,
            },
            Series {
                label: "Hagan T=10y",
                x: &k_pct,
                y: &smile_t10_hag,
                color: "#f59e0b",
                dashed: true,
            },
        ],
    )
    .unwrap();

    // Hagan error map bar chart across nu at rho=-0.3, beta=0.5, T=5
    let nu_map_vals: Vec<f64> = error_map_rows
        .iter()
        .filter(|r| (r.1 - (-0.3)).abs() < 1e-9 && (r.2 - 0.5).abs() < 1e-9)
        .map(|r| r.3)
        .collect();
    write_bar_chart(
        &fig_dir.join("day005a_hagan-error-map.svg"),
        "Hagan IV RMSE vs nu (T=5y, rho=-0.3, beta=0.5, bps)",
        "IV RMSE (bps)",
        &["nu=0.2", "nu=0.4", "nu=0.8"],
        &nu_map_vals,
        &["#16a34a", "#f59e0b", "#dc2626"],
    )
    .unwrap();

    let dens_x: Vec<f64> = density_rows[0].2.ks.iter().map(|&k| k * 100.0).collect();
    let dens_t1: Vec<f64> = density_rows[0].2.qs.clone();
    let dens_t5: Vec<f64> = density_rows[1].2.qs.clone();
    let dens_t10: Vec<f64> = density_rows[2].2.qs.clone();
    write_line_chart(
        &fig_dir.join("day005a_implied-density.svg"),
        "Hagan Breeden-Litzenberger Implied Density q(K) across T",
        "Strike K (%)",
        "Implied Density q(K)",
        &[
            Series {
                label: "T = 1y",
                x: &dens_x,
                y: &dens_t1,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "T = 5y",
                x: &dens_x,
                y: &dens_t5,
                color: "#f59e0b",
                dashed: false,
            },
            Series {
                label: "T = 10y",
                x: &dens_x,
                y: &dens_t10,
                color: "#dc2626",
                dashed: false,
            },
        ],
    )
    .unwrap();

    let beta_rmses: Vec<f64> = calib_beta_sweep_t5.iter().map(|f| f.rmse_iv_bps).collect();
    write_bar_chart(
        &fig_dir.join("day005a_calib-beta-equiv.svg"),
        "Non-Identifiability of beta from Smile Alone (T=5y Fit RMSE, bps)",
        "Fit IV RMSE (bps)",
        &["beta=0.0", "beta=0.5", "beta=1.0"],
        &beta_rmses,
        &["#2563eb", "#16a34a", "#7c3aed"],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day005a_bartlett-hedge-pnl.svg"),
        "Daily Delta-Hedging PnL SD (T=1y ATM Call, R=20k CRN, bps)",
        "PnL Std Dev (bps)",
        &["Bartlett", "Naive", "Black ATM", "CEV"],
        &hedge_main.sd_pnl_bps,
        &["#16a34a", "#dc2626", "#f59e0b", "#6b7280"],
    )
    .unwrap();

    let total_secs = t_total.elapsed().as_secs_f64();

    // =========================================================================
    // Criteria evaluation (established a priori)
    // =========================================================================
    let c1_pass = max_pde_rel_err < 1e-6 && max_mc_exact_z <= 3.0 && max_mc_euler_z <= 3.0;
    let c2_pass = (p_abs_t10_exact - p_abs_t10_feller).abs() < 1e-12
        && (p_abs_t10_mc - p_abs_t10_exact).abs() < 0.005
        && (mean_ft_t10_mc - f0).abs() < 0.002;
    let c3_pass = euler_conv_rows.iter().any(|r| r.3 <= 3.0 && r.6 <= 3.0);
    let max_abs_err_t1_liquid_bps = ref_smiles[0]
        .4
        .iter()
        .filter(|p| p.5.is_some_and(|se| se < 50.0))
        .filter_map(|p| p.7)
        .map(f64::abs)
        .fold(0.0_f64, f64::max);
    let c4_pass = max_abs_err_t1_liquid_bps < 30.0;
    let c5_pass = density_rows[2].2.total_neg_mass > 0.0 && density_rows[2].2.min_q < 0.0;
    let b0_rmse = calib_beta_sweep_t5[0].rmse_iv_bps;
    let b1_rmse = calib_beta_sweep_t5[2].rmse_iv_bps;
    let c6_pass = (b0_rmse - b1_rmse).abs() < 20.0 && calib_fixed_beta[0].1.rmse_iv_bps < 15.0;
    let c7_pass = hedge_main.sd_pnl_bps[0] < hedge_main.sd_pnl_bps[1]
        && hedge_main.paired_var_diff_bps2 < 0.0
        && hedge_main.paired_var_diff_bps2.abs() > 3.0 * hedge_main.paired_var_diff_se_bps2;
    let c8_pass = total_secs < 480.0;

    let status_str = |b: bool| if b { "tercapai" } else { "belum tercapai" };

    let doc = obj([
        (
            "metadata",
            obj([
                ("session", "day005a".into()),
                ("date_wita", "2026-11-23".into()),
                ("language", "rust".into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as i64).into()),
                (
                    "elapsed_seconds",
                    obj([
                        ("part1_cev", p1_secs.into()),
                        ("part2_hagan_nu0", p2_secs.into()),
                        ("part3_sabr_mc", p3_secs.into()),
                        ("part4_density", p4_secs.into()),
                        ("part5_calib_hedge", p5_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part1_cev",
            obj([
                ("max_pde_rel_err", max_pde_rel_err.into()),
                ("max_mc_exact_z", max_mc_exact_z.into()),
                ("max_mc_euler_z", max_mc_euler_z.into()),
                ("p_absorb_t10_exact", p_abs_t10_exact.into()),
                ("p_absorb_t10_feller", p_abs_t10_feller.into()),
                ("p_absorb_t10_mc", p_abs_t10_mc.into()),
                ("mean_ft_t10_mc", mean_ft_t10_mc.into()),
                (
                    "three_way_grid",
                    cev_3way_rows
                        .iter()
                        .map(|r| {
                            obj([
                                ("beta", r.0.into()),
                                ("k", r.1.into()),
                                ("exact_bps", (r.2 * 1e4).into()),
                                ("pde_bps", (r.3 * 1e4).into()),
                                ("pde_rel_err", r.4.into()),
                                ("mc_exact_bps", (r.5 * 1e4).into()),
                                ("mc_se_bps", (r.6 * 1e4).into()),
                                ("z_mc_exact", r.7.into()),
                                ("mc_euler_bps", (r.8 * 1e4).into()),
                                ("z_mc_euler", r.9.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part2_hagan_nu0",
            obj([(
                "nu0_vs_cev",
                hagan_nu0_rows
                    .iter()
                    .map(|r| {
                        obj([
                            ("t", r.0.into()),
                            ("k", r.1.into()),
                            ("iv_cev_exact", r.2.map_or(JsonVal::Null, Into::into)),
                            ("iv_hagan_nu1e4", r.3.into()),
                            ("err_bps", r.4.map_or(JsonVal::Null, Into::into)),
                        ])
                    })
                    .collect::<Vec<_>>()
                    .into(),
            )]),
        ),
        (
            "part3_sabr_reference",
            obj([
                (
                    "cond_rho0_mc",
                    cond_rho0_rows
                        .iter()
                        .map(|(t_h, res)| {
                            obj([
                                ("t", (*t_h).into()),
                                (
                                    "prices_bps",
                                    res.prices
                                        .iter()
                                        .map(|&p| JsonVal::from(p * 1e4))
                                        .collect::<Vec<_>>()
                                        .into(),
                                ),
                                (
                                    "ses_bps",
                                    res.ses
                                        .iter()
                                        .map(|&s| JsonVal::from(s * 1e4))
                                        .collect::<Vec<_>>()
                                        .into(),
                                ),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "euler_rho0_convergence_t1",
                    euler_conv_rows
                        .iter()
                        .map(|r| {
                            obj([
                                ("n_steps_yr", r.0.into()),
                                ("cond_exact_bps", (cond_t1_atm * 1e4).into()),
                                ("cond_se_bps", (cond_t1_atm_se * 1e4).into()),
                                ("absorbing_bps", (r.1 * 1e4).into()),
                                ("absorbing_se_bps", (r.2 * 1e4).into()),
                                ("z_absorbing", r.3.into()),
                                ("full_trunc_bps", (r.4 * 1e4).into()),
                                ("full_trunc_se_bps", (r.5 * 1e4).into()),
                                ("z_full_trunc", r.6.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "ref_smiles_rho_neg03",
                    ref_smiles
                        .iter()
                        .map(|(t_h, p_abs, m_ft, se_ft, pts)| {
                            obj([
                                ("t", (*t_h).into()),
                                ("p_absorb", (*p_abs).into()),
                                ("mean_ft", (*m_ft).into()),
                                ("se_ft", (*se_ft).into()),
                                (
                                    "points",
                                    pts.iter()
                                        .map(|p| {
                                            obj([
                                                ("k", p.0.into()),
                                                ("mc_price_bps", p.1.into()),
                                                ("mc_se_bps", p.2.into()),
                                                ("hagan_price_bps", p.3.into()),
                                                ("mc_iv", p.4.map_or(JsonVal::Null, Into::into)),
                                                (
                                                    "mc_iv_se_bps",
                                                    p.5.map_or(JsonVal::Null, Into::into),
                                                ),
                                                ("hagan_iv", p.6.into()),
                                                (
                                                    "iv_err_bps",
                                                    p.7.map_or(JsonVal::Null, Into::into),
                                                ),
                                            ])
                                        })
                                        .collect::<Vec<_>>()
                                        .into(),
                                ),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "hagan_error_map_t5",
                    error_map_rows
                        .iter()
                        .map(|r| {
                            obj([
                                ("nu", r.0.into()),
                                ("rho", r.1.into()),
                                ("beta", r.2.into()),
                                ("rmse_iv_bps", r.3.into()),
                                ("max_abs_err_bps", r.4.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part4_density_and_arbitrage",
            obj([
                (
                    "density_scans",
                    density_rows
                        .iter()
                        .map(|(t_h, p_abs, sc)| {
                            obj([
                                ("t", (*t_h).into()),
                                ("ref_p_absorb", (*p_abs).into()),
                                ("k_neg_min", sc.k_neg_min.map_or(JsonVal::Null, Into::into)),
                                ("k_neg_max", sc.k_neg_max.map_or(JsonVal::Null, Into::into)),
                                ("min_q", sc.min_q.into()),
                                ("total_neg_mass", sc.total_neg_mass.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "pricing_impact",
                    pricing_impact_rows
                        .iter()
                        .map(|r| {
                            obj([
                                ("t", r.0.into()),
                                ("k", r.1.into()),
                                ("digital_hagan", r.2.into()),
                                ("digital_mc", r.3.into()),
                                ("butterfly_hagan_bps", r.4.into()),
                                ("butterfly_mc_bps", r.5.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part5_calibration_and_hedging",
            obj([
                (
                    "calib_fixed_beta05",
                    calib_fixed_beta
                        .iter()
                        .map(|(t_h, fit)| {
                            obj([
                                ("t", (*t_h).into()),
                                ("alpha_hat", fit.alpha.into()),
                                ("rho_hat", fit.rho.into()),
                                ("nu_hat", fit.nu.into()),
                                ("rmse_iv_bps", fit.rmse_iv_bps.into()),
                                ("corr_alpha_rho", fit.corr_matrix[0][1].into()),
                                ("corr_alpha_nu", fit.corr_matrix[0][2].into()),
                                ("corr_rho_nu", fit.corr_matrix[1][2].into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "calib_beta_sweep_t5",
                    calib_beta_sweep_t5
                        .iter()
                        .map(|fit| {
                            obj([
                                ("beta", fit.beta.into()),
                                ("alpha_hat", fit.alpha.into()),
                                ("rho_hat", fit.rho.into()),
                                ("nu_hat", fit.nu.into()),
                                ("rmse_iv_bps", fit.rmse_iv_bps.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "bartlett_hedging_t1",
                    obj([
                        ("initial_price_bps", hedge_main.initial_price_bps.into()),
                        ("bartlett_mean_bps", hedge_main.mean_pnl_bps[0].into()),
                        ("bartlett_sd_bps", hedge_main.sd_pnl_bps[0].into()),
                        ("naive_mean_bps", hedge_main.mean_pnl_bps[1].into()),
                        ("naive_sd_bps", hedge_main.sd_pnl_bps[1].into()),
                        ("black_atm_mean_bps", hedge_main.mean_pnl_bps[2].into()),
                        ("black_atm_sd_bps", hedge_main.sd_pnl_bps[2].into()),
                        ("cev_mean_bps", hedge_main.mean_pnl_bps[3].into()),
                        ("cev_sd_bps", hedge_main.sd_pnl_bps[3].into()),
                        (
                            "paired_var_diff_bps2",
                            hedge_main.paired_var_diff_bps2.into(),
                        ),
                        (
                            "paired_var_diff_se_bps2",
                            hedge_main.paired_var_diff_se_bps2.into(),
                        ),
                        (
                            "paired_se_pnl_diff_bps",
                            hedge_main.paired_se_diff_bps.into(),
                        ),
                    ]),
                ),
                (
                    "misspec_beta_hedging_t1",
                    hedge_misspec_beta
                        .iter()
                        .map(|r| {
                            obj([
                                ("beta_calib", r.0.into()),
                                ("smile_rmse_bps", r.1.into()),
                                ("delta_bartlett", r.2.into()),
                                ("delta_naive", r.3.into()),
                                ("sd_pnl_bartlett_bps", r.4.into()),
                                ("sd_pnl_naive_bps", r.5.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            obj([
                (
                    "c1_cev_three_ways_match",
                    obj([
                        ("status", status_str(c1_pass).into()),
                        ("max_pde_rel_err", max_pde_rel_err.into()),
                        ("max_mc_exact_z", max_mc_exact_z.into()),
                        ("max_mc_euler_z", max_mc_euler_z.into()),
                    ]),
                ),
                (
                    "c2_cev_absorption_and_martingale",
                    obj([
                        ("status", status_str(c2_pass).into()),
                        ("p_abs_exact", p_abs_t10_exact.into()),
                        ("p_abs_mc", p_abs_t10_mc.into()),
                        ("mean_ft_mc", mean_ft_t10_mc.into()),
                    ]),
                ),
                (
                    "c3_sabr_rho0_cond_matches_euler_3se",
                    obj([
                        ("status", status_str(c3_pass).into()),
                        (
                            "z_abs_n800",
                            euler_conv_rows.last().map_or(0.0, |r| r.3).into(),
                        ),
                        (
                            "z_ft_n800",
                            euler_conv_rows.last().map_or(0.0, |r| r.6).into(),
                        ),
                    ]),
                ),
                (
                    "c4_hagan_t1_accurate_under_30bps",
                    obj([
                        ("status", status_str(c4_pass).into()),
                        (
                            "max_abs_err_t1_liquid_bps",
                            max_abs_err_t1_liquid_bps.into(),
                        ),
                    ]),
                ),
                (
                    "c5_hagan_negative_density_detected",
                    obj([
                        ("status", status_str(c5_pass).into()),
                        ("t10_neg_mass", density_rows[2].2.total_neg_mass.into()),
                        ("t10_min_q", density_rows[2].2.min_q.into()),
                    ]),
                ),
                (
                    "c6_beta_non_identifiable_from_smile",
                    obj([
                        ("status", status_str(c6_pass).into()),
                        ("t5_beta0_rmse_bps", b0_rmse.into()),
                        ("t5_beta1_rmse_bps", b1_rmse.into()),
                    ]),
                ),
                (
                    "c7_bartlett_delta_reduces_pnl_var_3se",
                    obj([
                        ("status", status_str(c7_pass).into()),
                        ("bartlett_sd_bps", hedge_main.sd_pnl_bps[0].into()),
                        ("naive_sd_bps", hedge_main.sd_pnl_bps[1].into()),
                        (
                            "paired_var_diff_bps2",
                            hedge_main.paired_var_diff_bps2.into(),
                        ),
                        (
                            "paired_var_diff_se_bps2",
                            hedge_main.paired_var_diff_se_bps2.into(),
                        ),
                    ]),
                ),
                (
                    "c8_runtime_under_480s",
                    obj([
                        ("status", status_str(c8_pass).into()),
                        ("total_seconds", total_secs.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &doc).expect("write results.json");
    println!(
        "day005a selesai dalam {:.2}s | workers={} | c1..c8={}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        n_workers,
        c1_pass,
        c2_pass,
        c3_pass,
        c4_pass,
        c5_pass,
        c6_pass,
        c7_pass,
        c8_pass
    );
}
