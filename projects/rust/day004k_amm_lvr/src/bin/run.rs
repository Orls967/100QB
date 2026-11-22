use std::path::Path;
use std::time::Instant;

use day004k_amm_lvr::blocks::{
    block_arb_outcome, exact_block_expectations_gl, mc_block_annualized,
};
use day004k_amm_lvr::conc::{expected_exit_time_gbm, mc_conc_metrics, ConcRange};
use day004k_amm_lvr::cpmm::{
    expected_pool_value, impermanent_loss, lvr_rate_analytic, mc_pool_moments, pathwise_lvr_ito,
    pool_value, variance_pool_value,
};
use day004k_amm_lvr::jsonout::{obj, write_json_file};
use day004k_amm_lvr::misspec::{misspec_drift, misspec_merton_jumps, misspec_stoch_vol};
use day004k_amm_lvr::parallel::default_workers;
use day004k_amm_lvr::reflect::{
    mc_noise_trader_breakeven, mc_reflected_block_sim, reflected_limit_rates,
};
use day004k_amm_lvr::stats;
use day004k_amm_lvr::strategies::compare_lp_strategies;
use day004k_amm_lvr::svg::{write_bar_chart, write_line_chart, Series};

const MAIN_SEED: u64 = 20261122;

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
    let n_workers = default_workers().min(8);
    let root = repo_root();
    let out_json = root.join("outputs/day004k/results.json");
    let fig_dir = root.join("outputs/figures");

    let sigma = 0.8_f64;
    let v0 = 1.0_f64;
    let l = 0.5_f64;
    let sec_yr = 365.0 * 86_400.0;

    // -------------------------------------------------------------------------
    // PART 1: Closed-form CPMM, GBM moments (R=10^6), and pathwise Ito LVR
    // -------------------------------------------------------------------------
    let t_p1 = Instant::now();
    let horizons = [1.0 / 365.0, 7.0 / 365.0, 30.0 / 365.0, 1.0];
    let mut moment_rows = Vec::new();
    let mut all_moments_3se = true;
    for (idx, &t_h) in horizons.iter().enumerate() {
        let (m_mc, se_mc, var_mc) = mc_pool_moments(
            v0,
            sigma,
            t_h,
            1_000_000,
            n_workers,
            MAIN_SEED ^ (idx as u64),
        );
        let m_ex = expected_pool_value(v0, sigma, t_h);
        let var_ex = variance_pool_value(v0, sigma, t_h);
        let z = (m_mc - m_ex).abs() / se_mc;
        if z > 3.0 {
            all_moments_3se = false;
        }
        moment_rows.push((t_h, m_ex, m_mc, se_mc, z, var_ex, var_mc));
    }

    let (lvr_rep, lvr_int, lvr_rmse_path) =
        pathwise_lvr_ito(v0, sigma, 30.0 / 365.0, 1000, 50_000, n_workers, MAIN_SEED);

    let mut sigma_rows = Vec::new();
    for &sig in &[0.3_f64, 0.5, 0.8, 1.2] {
        let th_rate = lvr_rate_analytic(v0, sig);
        let (rep_s, int_s, _) = pathwise_lvr_ito(
            v0,
            sig,
            1.0 / 365.0,
            500,
            20_000,
            n_workers,
            MAIN_SEED ^ sig.to_bits(),
        );
        let rate_mc = rep_s * 365.0;
        let rate_int = int_s * 365.0;
        sigma_rows.push((sig, th_rate, rate_mc, rate_int, (rate_mc - th_rate).abs()));
    }
    let p1_secs = t_p1.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 2: Discrete blocks, fees, reflected Brownian limit & breakeven
    // -------------------------------------------------------------------------
    let t_p2 = Instant::now();
    let deltas_sec = [12.0_f64, 60.0, 600.0];
    let phis = [0.0_f64, 0.0005, 0.003, 0.01];
    let mut block_rows = Vec::new();
    let mut all_block_gl_3se = true;
    for &d_sec in &deltas_sec {
        let dt = d_sec / sec_yr;
        for &phi in &phis {
            let (gp, gf, gv) = exact_block_expectations_gl(sigma, dt, phi, l);
            let gp_a = gp / dt;
            let gf_a = gf / dt;
            let gv_a = gv / dt;
            let (mp_a, se_a, mf_a, mv_a, _, _) = mc_block_annualized(
                sigma,
                dt,
                phi,
                500_000,
                n_workers,
                MAIN_SEED ^ (d_sec as u64) ^ ((phi * 1e5) as u64),
            );
            let z = if se_a > 1e-12 {
                (gp_a - mp_a).abs() / se_a
            } else {
                0.0
            };
            if z > 3.0 {
                all_block_gl_3se = false;
            }
            block_rows.push((d_sec, phi, gp_a, mp_a, se_a, z, gf_a, mf_a, gv_a, mv_a));
        }
    }

    // Continuous reflected limit (Delta -> 0) vs Delta = 1s, 12s, 60s at phi = 0.003
    let (delta_band_003, vol_lim_003, fee_lim_003) = reflected_limit_rates(sigma, 0.003);
    let mut reflect_rows = Vec::new();
    for &(d_s, t_days, np) in &[
        (1.0_f64, 1.0_f64, 20_000_usize),
        (12.0, 3.0, 30_000),
        (60.0, 7.0, 40_000),
    ] {
        let (pr, se, fe, vo) =
            mc_reflected_block_sim(sigma, 0.003, d_s, t_days, np, n_workers, MAIN_SEED);
        reflect_rows.push((d_s, pr, se, fe, vo));
    }

    // Noise-trader breakeven volume/TVL at phi = 0.003
    let mut breakeven_rows = Vec::new();
    let req_daily_vol_tvl = (sigma * sigma / 8.0) / (0.003 * 365.0);
    for &v_day in &[0.0_f64, 0.05, 0.1, 0.3, 1.0] {
        let (net_lp, se_net, fee_yr) =
            mc_noise_trader_breakeven(sigma, 0.003, v_day, 20_000, n_workers, MAIN_SEED);
        breakeven_rows.push((v_day, net_lp, se_net, fee_yr));
    }
    let p2_secs = t_p2.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 3: Concentrated liquidity (Uniswap v3 style) & LP strategies
    // -------------------------------------------------------------------------
    let t_p3 = Instant::now();
    let ranges = [
        ("pm_5pct", ConcRange::normalized(0.95, 1.05)),
        ("pm_10pct", ConcRange::normalized(0.90, 1.10)),
        ("pm_30pct", ConcRange::normalized(0.70, 1.30)),
        ("full_range", ConcRange::full_range()),
    ];
    let mut conc_rows = Vec::new();
    let mut all_conc_3se = true;
    for (nm, rg) in &ranges {
        let ce = rg.capital_efficiency();
        let lvr0_over_v0 = rg.lvr_rate(1.0, sigma);
        let (m_rep, se_rep, m_a, m_in) =
            mc_conc_metrics(*rg, sigma, 14.0, 40_000, n_workers, MAIN_SEED);
        let z = (m_rep - m_a).abs() / se_rep;
        if z > 3.0 {
            all_conc_3se = false;
        }
        let exp_exit_yr = expected_exit_time_gbm(rg.pa, rg.pb, sigma);
        let req_vol_tvl_day = lvr0_over_v0 / (0.003 * 365.0);
        conc_rows.push((
            (*nm).to_string(),
            ce,
            lvr0_over_v0,
            m_rep,
            se_rep,
            m_a,
            z,
            m_in,
            exp_exit_yr,
            req_vol_tvl_day,
        ));
    }

    let strat_res = compare_lp_strategies(
        sigma, 0.003, 0.15, 0.0005, 3, 30, 50_000, n_workers, MAIN_SEED,
    );
    let hodl_pnl = &strat_res[0].pnl;
    let pas_pnl = &strat_res[2].pnl;
    let hedged_pnl = &strat_res[4].pnl;
    let diff_hedged_vs_pas: Vec<f64> = hedged_pnl
        .iter()
        .zip(pas_pnl.iter())
        .map(|(a, b)| a - b)
        .collect();
    let diff_pas_vs_hodl: Vec<f64> = pas_pnl
        .iter()
        .zip(hodl_pnl.iter())
        .map(|(a, b)| a - b)
        .collect();
    let p3_secs = t_p3.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // PART 4: Misspecification (Jumps, Stoch Vol, Drift, Block Size)
    // -------------------------------------------------------------------------
    let t_p4 = Instant::now();
    let (jmp_real, jmp_diff, jmp_se) = misspec_merton_jumps(
        sigma,
        20.0,
        -0.02,
        0.05,
        30.0 / 365.0,
        1500,
        40_000,
        n_workers,
        MAIN_SEED,
    );
    let (sv_real, sv_pred, sv_se, sv_corr) = misspec_stoch_vol(
        sigma * sigma,
        3.0,
        1.0,
        30.0 / 365.0,
        1500,
        40_000,
        n_workers,
        MAIN_SEED,
    );
    let dr_neg = misspec_drift(
        sigma,
        -0.5,
        30.0 / 365.0,
        1000,
        40_000,
        n_workers,
        MAIN_SEED,
    );
    let dr_pos = misspec_drift(sigma, 0.5, 30.0 / 365.0, 1000, 40_000, n_workers, MAIN_SEED);

    let blk_600_phi0 = block_rows
        .iter()
        .find(|x| (x.0 - 600.0).abs() < 1e-6 && x.1 == 0.0)
        .unwrap()
        .2;
    let blk_600_phi003 = block_rows
        .iter()
        .find(|x| (x.0 - 600.0).abs() < 1e-6 && (x.1 - 0.003).abs() < 1e-9)
        .unwrap()
        .2;

    // -------------------------------------------------------------------------
    // SVG Figures (8 charts)
    // -------------------------------------------------------------------------
    let rs_grid: Vec<f64> = (0..81).map(|i| 0.25 + 3.75 * i as f64 / 80.0).collect();
    let vp_ys: Vec<f64> = rs_grid.iter().map(|&r| pool_value(l, r)).collect();
    let il_ys: Vec<f64> = rs_grid.iter().map(|&r| impermanent_loss(r)).collect();
    write_line_chart(
        &fig_dir.join("day004k_vp-and-il.svg"),
        "CPMM Value V(P) and Impermanent Loss IL(r)",
        "Price Ratio r = P / P0",
        "Value / IL",
        &[
            Series {
                label: "Pool Value V(r)",
                x: &rs_grid,
                y: &vp_ys,
                color: "#2563eb",
                dashed: false,
            },
            Series {
                label: "Impermanent Loss IL(r)",
                x: &rs_grid,
                y: &il_ys,
                color: "#dc2626",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let t_days_x: Vec<f64> = moment_rows.iter().map(|x| x.0 * 365.0).collect();
    let ev_ex_y: Vec<f64> = moment_rows.iter().map(|x| x.1).collect();
    let ev_mc_y: Vec<f64> = moment_rows.iter().map(|x| x.2).collect();
    write_line_chart(
        &fig_dir.join("day004k_ev-vs-mc.svg"),
        "Martingale Anchor: E[V(P_T)] = V0 * exp(-sigma^2 T / 8) vs MC (R=10^6)",
        "Horizon T (days)",
        "Expected Pool Value",
        &[
            Series {
                label: "Analytic",
                x: &t_days_x,
                y: &ev_ex_y,
                color: "#111827",
                dashed: false,
            },
            Series {
                label: "MC (R=10^6)",
                x: &t_days_x,
                y: &ev_mc_y,
                color: "#2563eb",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let d_x: Vec<f64> = deltas_sec.to_vec();
    let lvr_phi0: Vec<f64> = block_rows
        .iter()
        .filter(|x| x.1 == 0.0)
        .map(|x| x.2)
        .collect();
    let lvr_phi003: Vec<f64> = block_rows
        .iter()
        .filter(|x| (x.1 - 0.003).abs() < 1e-9)
        .map(|x| x.2)
        .collect();
    write_line_chart(
        &fig_dir.join("day004k_lvr-vs-delta.svg"),
        "Annualized Net LVR vs Block Time Delta (sigma = 0.8)",
        "Block Interval Delta (seconds)",
        "Annualized LVR / V0",
        &[
            Series {
                label: "phi = 0 (no fee)",
                x: &d_x,
                y: &lvr_phi0,
                color: "#111827",
                dashed: false,
            },
            Series {
                label: "phi = 0.003 (30 bps)",
                x: &d_x,
                y: &lvr_phi003,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .unwrap();

    let zs_grid: Vec<f64> = (0..81).map(|i| -0.02 + 0.04 * i as f64 / 80.0).collect();
    let blk_pr_y: Vec<f64> = zs_grid
        .iter()
        .map(|&z| block_arb_outcome(z, 0.003, l, 1.0).0)
        .collect();
    let blk_fe_y: Vec<f64> = zs_grid
        .iter()
        .map(|&z| block_arb_outcome(z, 0.003, l, 1.0).1)
        .collect();
    write_line_chart(
        &fig_dir.join("day004k_block-profit-vs-z.svg"),
        "Per-Block Arbitrage Profit & LP Fee vs Log-Price Shock z (phi = 0.003)",
        "Log-Price Shock z = ln(P'/p)",
        "Fraction of Pool Value",
        &[
            Series {
                label: "Arb Profit",
                x: &zs_grid,
                y: &blk_pr_y,
                color: "#dc2626",
                dashed: false,
            },
            Series {
                label: "LP Fee",
                x: &zs_grid,
                y: &blk_fe_y,
                color: "#16a34a",
                dashed: true,
            },
        ],
    )
    .unwrap();

    let ref_dx: Vec<f64> = reflect_rows.iter().map(|x| x.0).collect();
    let ref_fe: Vec<f64> = reflect_rows.iter().map(|x| x.3).collect();
    let ref_lim: Vec<f64> = vec![fee_lim_003; ref_dx.len()];
    write_line_chart(
        &fig_dir.join("day004k_reflect-vs-block.svg"),
        "Reflected Brownian Limit vs Discrete-Block LP Fee Rate (phi = 0.003)",
        "Block Interval Delta (seconds)",
        "Annualized LP Fee Rate / V0",
        &[
            Series {
                label: "Reflected Limit (Delta->0)",
                x: &ref_dx,
                y: &ref_lim,
                color: "#111827",
                dashed: true,
            },
            Series {
                label: "Discrete Block MC",
                x: &ref_dx,
                y: &ref_fe,
                color: "#2563eb",
                dashed: false,
            },
        ],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day004k_conc-range-lvr.svg"),
        "Concentrated Liquidity: Instantaneous LVR/V0 by Range (sigma = 0.8)",
        "Annualized LVR / V0",
        &["+-5%", "+-10%", "+-30%", "Full Range"],
        &[
            conc_rows[0].2,
            conc_rows[1].2,
            conc_rows[2].2,
            conc_rows[3].2,
        ],
        &["#dc2626", "#f59e0b", "#2563eb", "#16a34a"],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day004k_lp-strategies.svg"),
        "30-Day LP & Hedging Strategies (CRN, R=50k): Mean PnL",
        "Mean 30-Day PnL",
        &[
            "Hodl",
            "Full LP",
            "Passive +-10%",
            "Rebal +-10%",
            "Hedged +-10%",
        ],
        &[
            strat_res[0].mean_pnl,
            strat_res[1].mean_pnl,
            strat_res[2].mean_pnl,
            strat_res[3].mean_pnl,
            strat_res[4].mean_pnl,
        ],
        &["#6b7280", "#16a34a", "#f59e0b", "#7c3aed", "#2563eb"],
    )
    .unwrap();

    write_bar_chart(
        &fig_dir.join("day004k_misspec-lvr.svg"),
        "Misspecification: Realized Annualized LVR Across Regimes",
        "Annualized LVR / V0",
        &[
            "Diffusion (0.08)",
            "Merton Jumps",
            "Stoch Vol (CIR)",
            "Drift mu=-0.5",
            "Drift mu=+0.5",
        ],
        &[0.08, jmp_real, sv_real, dr_neg.2, dr_pos.2],
        &["#111827", "#dc2626", "#2563eb", "#f59e0b", "#16a34a"],
    )
    .unwrap();

    let p4_secs = t_p4.elapsed().as_secs_f64();
    let total_secs = t_total.elapsed().as_secs_f64();

    // -------------------------------------------------------------------------
    // Criteria evaluation (a priori)
    // -------------------------------------------------------------------------
    let c1_pass = all_moments_3se;
    let c2_pass = (lvr_rep - lvr_int).abs() < 2e-4;
    let c3_pass = all_block_gl_3se;
    let c4_pass = (reflect_rows[0].3 - fee_lim_003).abs() / fee_lim_003 < 0.25;
    let c5_pass =
        all_conc_3se && conc_rows[0].2 > conc_rows[1].2 && conc_rows[1].2 > conc_rows[3].2;
    let c6_pass = jmp_real > jmp_diff;
    let c7_pass = (sv_real - sv_pred).abs() <= 3.0 * sv_se;
    let c8_pass = (dr_neg.3 - 0.08).abs() < 0.005 && (dr_pos.3 - 0.08).abs() < 0.005;
    let c9_pass = total_secs < 480.0;
    let status_str = |b: bool| if b { "tercapai" } else { "belum tercapai" };

    let root_json = obj([
        (
            "metadata",
            obj([
                ("session", "day004k".into()),
                ("date_wita", "2026-11-22".into()),
                ("language", "rust".into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as i64).into()),
                (
                    "elapsed_seconds",
                    obj([
                        ("part1_closed_form", p1_secs.into()),
                        ("part2_blocks_reflect", p2_secs.into()),
                        ("part3_conc_strategies", p3_secs.into()),
                        ("part4_misspec_svg", p4_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part1_closed_form",
            obj([
                (
                    "gbm_moments",
                    moment_rows
                        .iter()
                        .map(|&(t, me, mm, se, z, ve, vm)| {
                            obj([
                                ("t_years", t.into()),
                                ("e_v_exact", me.into()),
                                ("e_v_mc", mm.into()),
                                ("se_mc", se.into()),
                                ("z_diff", z.into()),
                                ("var_exact", ve.into()),
                                ("var_mc", vm.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("pathwise_lvr_rebal_30d", lvr_rep.into()),
                ("pathwise_lvr_ito_30d", lvr_int.into()),
                ("pathwise_lvr_rmse", lvr_rmse_path.into()),
                (
                    "sigma_sweep",
                    sigma_rows
                        .iter()
                        .map(|&(s, th, mc, it, err)| {
                            obj([
                                ("sigma", s.into()),
                                ("analytic_lvr_rate", th.into()),
                                ("mc_rebal_rate", mc.into()),
                                ("ito_int_rate", it.into()),
                                ("abs_err", err.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part2_blocks_and_reflect",
            obj([
                (
                    "block_grid",
                    block_rows
                        .iter()
                        .map(|&(ds, ph, gp, mp, se, z, gf, mf, gv, mv)| {
                            obj([
                                ("delta_sec", ds.into()),
                                ("phi", ph.into()),
                                ("gl_lvr_net_ann", gp.into()),
                                ("mc_lvr_net_ann", mp.into()),
                                ("se_mc", se.into()),
                                ("z_diff", z.into()),
                                ("gl_fee_ann", gf.into()),
                                ("mc_fee_ann", mf.into()),
                                ("gl_vol_ann", gv.into()),
                                ("mc_vol_ann", mv.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("reflect_band_delta_phi003", delta_band_003.into()),
                ("reflect_vol_rate_limit_phi003", vol_lim_003.into()),
                ("reflect_fee_rate_limit_phi003", fee_lim_003.into()),
                (
                    "reflect_delta_convergence",
                    reflect_rows
                        .iter()
                        .map(|&(ds, pr, se, fe, vo)| {
                            obj([
                                ("delta_sec", ds.into()),
                                ("net_lvr_ann", pr.into()),
                                ("se", se.into()),
                                ("fee_rate_ann", fe.into()),
                                ("arb_vol_ann", vo.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("breakeven_daily_vol_tvl_phi003", req_daily_vol_tvl.into()),
                (
                    "noise_trader_breakeven",
                    breakeven_rows
                        .iter()
                        .map(|&(vd, net, se, fe)| {
                            obj([
                                ("daily_vol_tvl", vd.into()),
                                ("net_lp_ann", net.into()),
                                ("se", se.into()),
                                ("total_fee_ann", fe.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part3_concentrated_and_strategies",
            obj([
                (
                    "ranges",
                    conc_rows
                        .iter()
                        .map(|r| {
                            obj([
                                ("name", r.0.clone().into()),
                                ("capital_efficiency", r.1.into()),
                                ("analytic_lvr_v0", r.2.into()),
                                ("mc_rebal_lvr_v0", r.3.into()),
                                ("se_mc", r.4.into()),
                                ("ito_int_lvr_v0", r.5.into()),
                                ("z_diff", r.6.into()),
                                ("frac_time_in_range", r.7.into()),
                                ("expected_exit_years", r.8.into()),
                                ("breakeven_daily_vol_tvl", r.9.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "strategies_30d",
                    strat_res
                        .iter()
                        .map(|s| {
                            obj([
                                ("name", s.name.clone().into()),
                                ("mean_pnl", s.mean_pnl.into()),
                                ("sd_pnl", s.sd_pnl.into()),
                                ("se_pnl", s.se_pnl.into()),
                                ("sharpe", s.sharpe.into()),
                                ("frac_in_range", s.frac_in_range.into()),
                                ("mean_gas", s.mean_gas.into()),
                                ("mean_rebalances", s.mean_rebalances.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "paired_hedged_minus_passive_mean",
                    stats::mean(&diff_hedged_vs_pas).into(),
                ),
                (
                    "paired_hedged_minus_passive_se",
                    stats::se(&diff_hedged_vs_pas).into(),
                ),
                (
                    "paired_passive_minus_hodl_mean",
                    stats::mean(&diff_pas_vs_hodl).into(),
                ),
                (
                    "paired_passive_minus_hodl_se",
                    stats::se(&diff_pas_vs_hodl).into(),
                ),
            ]),
        ),
        (
            "part4_misspecification",
            obj([
                ("merton_lvr_realized_ann", jmp_real.into()),
                ("merton_lvr_diffusion_ann", jmp_diff.into()),
                ("merton_jump_contribution", (jmp_real - jmp_diff).into()),
                ("merton_se", jmp_se.into()),
                ("stoch_vol_lvr_realized_ann", sv_real.into()),
                ("stoch_vol_lvr_expected_v_over_8", sv_pred.into()),
                ("stoch_vol_se", sv_se.into()),
                ("stoch_vol_corr_lvr_v", sv_corr.into()),
                ("drift_neg05_lp_ret", dr_neg.0.into()),
                ("drift_neg05_il", dr_neg.1.into()),
                ("drift_neg05_lvr_per_v_yr", dr_neg.3.into()),
                ("drift_pos05_lp_ret", dr_pos.0.into()),
                ("drift_pos05_il", dr_pos.1.into()),
                ("drift_pos05_lvr_per_v_yr", dr_pos.3.into()),
                ("block_600s_phi0_lvr_ann", blk_600_phi0.into()),
                ("block_600s_phi003_lvr_ann", blk_600_phi003.into()),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            obj([
                (
                    "c1_gbm_pool_moments_3se",
                    obj([
                        ("status", status_str(c1_pass).into()),
                        (
                            "max_z",
                            moment_rows
                                .iter()
                                .map(|x| x.4)
                                .fold(0.0_f64, f64::max)
                                .into(),
                        ),
                    ]),
                ),
                (
                    "c2_pathwise_lvr_ito_matches",
                    obj([
                        ("status", status_str(c2_pass).into()),
                        ("abs_diff", (lvr_rep - lvr_int).abs().into()),
                    ]),
                ),
                (
                    "c3_block_gl_matches_mc_3se",
                    obj([
                        ("status", status_str(c3_pass).into()),
                        (
                            "max_z",
                            block_rows
                                .iter()
                                .map(|x| x.5)
                                .fold(0.0_f64, f64::max)
                                .into(),
                        ),
                    ]),
                ),
                (
                    "c4_reflected_limit_matches_1s_block",
                    obj([
                        ("status", status_str(c4_pass).into()),
                        ("fee_limit", fee_lim_003.into()),
                        ("fee_1s_mc", reflect_rows[0].3.into()),
                    ]),
                ),
                (
                    "c5_concentrated_lvr_scales_with_ce",
                    obj([
                        ("status", status_str(c5_pass).into()),
                        ("lvr_5pct", conc_rows[0].2.into()),
                        ("lvr_10pct", conc_rows[1].2.into()),
                        ("lvr_full", conc_rows[3].2.into()),
                    ]),
                ),
                (
                    "c6_merton_jumps_increase_lvr",
                    obj([
                        ("status", status_str(c6_pass).into()),
                        ("lvr_jumps", jmp_real.into()),
                        ("lvr_diffusion", jmp_diff.into()),
                    ]),
                ),
                (
                    "c7_stoch_vol_lvr_matches_ev_over_8",
                    obj([
                        ("status", status_str(c7_pass).into()),
                        ("z_diff", ((sv_real - sv_pred).abs() / sv_se).into()),
                        ("corr_lvr_v", sv_corr.into()),
                    ]),
                ),
                (
                    "c8_drift_invariance_lvr_per_v",
                    obj([
                        ("status", status_str(c8_pass).into()),
                        ("lvr_per_v_neg05", dr_neg.3.into()),
                        ("lvr_per_v_pos05", dr_pos.3.into()),
                    ]),
                ),
                (
                    "c9_runtime_under_480s",
                    obj([
                        ("status", status_str(c9_pass).into()),
                        ("total_seconds", total_secs.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &root_json).expect("write results.json");
    println!(
        "day004k selesai dalam {:.2}s | workers={} | c1..c9={}/{}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        n_workers,
        c1_pass,
        c2_pass,
        c3_pass,
        c4_pass,
        c5_pass,
        c6_pass,
        c7_pass,
        c8_pass,
        c9_pass
    );
}
