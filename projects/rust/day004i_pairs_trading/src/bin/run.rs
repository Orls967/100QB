//! Unified CLI runner for day004i:
//! Part 1: OU exact vs Euler transition, closed-form vs numerical MLE, R=20,000 replication bias (Kendall).
//! Part 2: Bertram FPT closed-form vs BVP vs MC (R=200k), optimal threshold a*, 200-yr MC check, sensitivity.
//! Part 3: Dickey-Fuller (T=1000, R=10^6) + Brownian limit, Engle-Granger (R=20k), Spurious, Johansen (R=20k).
//! Part 4: Rolling backtest (W in {1,2,5}, R=2000 10-yr OOS), EG filter, structural break + stop-loss, cost curve.

use std::path::PathBuf;
use std::time::Instant;

use day004i_pairs_trading::backtest::run_backtest_study;
use day004i_pairs_trading::fpt::{expected_fpt_exact, simulate_fpt_mc, solve_fpt_bvp};
use day004i_pairs_trading::johansen::run_johansen_study;
use day004i_pairs_trading::jsonout::{write_json_file, JsonVal};
use day004i_pairs_trading::ou::{
    euler_loglik, exact_loglik, mle_closed_form, mle_numerical, run_ou_replication_study,
    simulate_ou_path, OuParams,
};
use day004i_pairs_trading::parallel;
use day004i_pairs_trading::rng::{derive_seed, Rng, MAIN_SEED};
use day004i_pairs_trading::stats;
use day004i_pairs_trading::strategy::{
    expected_return_rate, optimal_threshold_direct, optimal_threshold_foc, theoretical_sharpe_at_a,
    verify_strategy_mc_200yr,
};
use day004i_pairs_trading::svg::{self, LineSeries, PALETTE};
use day004i_pairs_trading::unitroot::{
    run_eg_power_grid, simulate_dickey_fuller_dist, simulate_engle_granger_dist,
    simulate_spurious_regression,
};

fn find_repo_root() -> PathBuf {
    let mut dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    for _ in 0..8 {
        if dir.join("outputs").is_dir() && dir.join("projects").join("rust").is_dir() {
            return dir;
        }
        if !dir.pop() {
            break;
        }
    }
    PathBuf::from(".")
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let quick = args.iter().any(|a| a == "--quick" || a == "-quick");
    let n_workers = parallel::default_workers();
    let root = find_repo_root();
    let out_json = if quick {
        std::env::temp_dir().join("day004i_quick_results.json")
    } else {
        root.join("outputs").join("day004i").join("results.json")
    };
    let fig_dir = if quick {
        std::env::temp_dir().join("day004i_quick_figures")
    } else {
        root.join("outputs").join("figures")
    };

    let t_total = Instant::now();
    let p = OuParams::default_params();
    let c_cost = 0.002_f64;

    // ==================== BAGIAN 1: OU & MLE ====================
    let t1 = Instant::now();
    let mut rg0 = Rng::new(MAIN_SEED);
    let sample_path = simulate_ou_path(p, 252 * 5, true, &mut rg0);
    let cf_mle = mle_closed_form(&sample_path, p.dt);
    let num_mle = mle_numerical(
        &sample_path,
        p.dt,
        [cf_mle.theta * 0.9, cf_mle.mu + 0.01, cf_mle.sigma * 1.1],
    );
    let mle_max_diff = (cf_mle.theta - num_mle[0])
        .abs()
        .max((cf_mle.mu - num_mle[1]).abs())
        .max((cf_mle.sigma - num_mle[2]).abs());
    let ll_exact = exact_loglik(&sample_path, p.dt, p.theta, p.mu, p.sigma);
    let ll_euler = euler_loglik(&sample_path, p.dt, p.theta, p.mu, p.sigma);

    // Stationary var & lag-1 ACF check
    let acf1_exact = (-p.theta * p.dt).exp();
    let var_stat_exact = p.stationary_sd().powi(2);
    let acf_samples = parallel::run_chunks(if quick { 2000 } else { 20_000 }, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(MAIN_SEED, 0x1050, r as u64));
        let path = simulate_ou_path(p, 252 * 10, true, &mut rg);
        let m = stats::mean(&path);
        let mut s0 = 0.0_f64;
        let mut s1 = 0.0_f64;
        for i in 0..(path.len() - 1) {
            s0 += (path[i] - m) * (path[i] - m);
            s1 += (path[i] - m) * (path[i + 1] - m);
        }
        (s0 / path.len() as f64, s1 / s0)
    });
    let emp_var_stat = stats::mean(&acf_samples.iter().map(|x| x.0).collect::<Vec<_>>());
    let emp_acf1 = stats::mean(&acf_samples.iter().map(|x| x.1).collect::<Vec<_>>());

    let r_rep_ou = if quick { 1000 } else { 20_000 };
    let ou_rep =
        run_ou_replication_study(p, &[1.0, 2.0, 5.0, 10.0], r_rep_ou, n_workers, MAIN_SEED);
    let part1_secs = t1.elapsed().as_secs_f64();

    // ==================== BAGIAN 2: FPT & Bertram Threshold ====================
    let t2 = Instant::now();
    let sd_stat = p.stationary_sd();
    let fpt_pairs = [
        (-sd_stat, sd_stat),
        (-1.5 * sd_stat, 1.5 * sd_stat),
        (-2.0 * sd_stat, 2.0 * sd_stat),
        (0.0, sd_stat),
    ];
    let r_fpt = if quick { 5_000 } else { 200_000 };
    let mut fpt_rows = Vec::new();
    for (idx, &(a, m)) in fpt_pairs.iter().enumerate() {
        let ex = expected_fpt_exact(a, m, p.theta, p.sigma);
        let (bvp, var_bvp) = solve_fpt_bvp(a, m, p.theta, p.sigma, 2000);
        let mc = simulate_fpt_mc(
            a,
            m,
            p.theta,
            p.sigma,
            1.0 / 2520.0,
            r_fpt,
            n_workers,
            MAIN_SEED,
            0x2100 + idx as u64,
        );
        let z_mc = (mc.mean_tau - ex).abs() / mc.se_tau.max(1e-12);
        fpt_rows.push((a, m, ex, bvp, var_bvp.sqrt(), mc.mean_tau, mc.se_tau, z_mc));
    }

    let (a_dir, r_dir) = optimal_threshold_direct(p.theta, p.sigma, c_cost);
    let a_foc = optimal_threshold_foc(p.theta, p.sigma, c_cost);
    let sharpe_a_star = theoretical_sharpe_at_a(a_foc, p.theta, p.sigma, c_cost);
    let mc_200yr = verify_strategy_mc_200yr(p.theta, p.mu, p.sigma, a_foc, c_cost, MAIN_SEED);

    // Sensitivity map a*(theta, c)
    let mut sens_map = Vec::new();
    for &th in &[8.0_f64, 16.0, 32.0] {
        for &cv in &[0.0_f64, 0.001, 0.002, 0.005] {
            let (a_opt, r_opt) = if cv == 0.0 {
                let a0 = -0.5 * (p.sigma / (2.0 * th).sqrt());
                (a0, expected_return_rate(a0, th, p.sigma, 0.0))
            } else {
                optimal_threshold_direct(th, p.sigma, cv)
            };
            sens_map.push((th, cv, a_opt, r_opt));
        }
    }
    let part2_secs = t2.elapsed().as_secs_f64();

    // ==================== BAGIAN 3: Unit Root & Cointegration ====================
    let t3 = Instant::now();
    let r_df = if quick { 20_000 } else { 1_000_000 };
    let r_eg = if quick { 2_000 } else { 20_000 };
    let df_dist = simulate_dickey_fuller_dist(1000, r_df, n_workers, MAIN_SEED);
    let eg_dist = simulate_engle_granger_dist(1000, r_eg, df_dist.q05, n_workers, MAIN_SEED);
    let spurious = simulate_spurious_regression(1000, r_eg, n_workers, MAIN_SEED);
    let eg_power = run_eg_power_grid(
        &[2.0, 8.0, 16.0, 32.0],
        &[2.0, 5.0, 10.0],
        eg_dist.q05,
        r_eg / 4,
        n_workers,
        MAIN_SEED,
    );
    let johansen = run_johansen_study(500, r_eg, n_workers, MAIN_SEED);
    let part3_secs = t3.elapsed().as_secs_f64();

    // ==================== BAGIAN 4: Trading Backtest ====================
    let t4 = Instant::now();
    let r_bt = if quick { 200 } else { 2_000 };
    let bt = run_backtest_study(p, c_cost, eg_dist.q05, r_bt, n_workers, MAIN_SEED);
    let part4_secs = t4.elapsed().as_secs_f64();

    // ==================== BAGIAN 5: 7+ SVG Figures ====================
    let t5 = Instant::now();
    // 1. Sample spread & optimal bands +/- a*
    let n_plot = 504_usize;
    let t_days: Vec<f64> = (0..=n_plot).map(|i| i as f64).collect();
    let sp_plot: Vec<f64> = sample_path[..=n_plot].to_vec();
    let upper_b = vec![-a_foc; n_plot + 1];
    let lower_b = vec![a_foc; n_plot + 1];
    svg::write_line_chart(
        &fig_dir.join("day004i_spread-and-thresholds.svg"),
        "Spread OU Contoh (theta=16, sigma=0.20) dan Ambang Optimal +/- a*",
        "Hari Perdagangan",
        "Spread X_t",
        &[
            LineSeries {
                label: "Spread X_t",
                x: &t_days,
                y: &sp_plot,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "+a* (Jual)",
                x: &t_days,
                y: &upper_b,
                color: PALETTE[1],
                dashed: true,
            },
            LineSeries {
                label: "-a* (Beli)",
                x: &t_days,
                y: &lower_b,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 2. Bias theta_hat vs T vs Kendall 4/T
    let ts_f: Vec<f64> = ou_rep.iter().map(|r| r.t_years).collect();
    let bias_obs: Vec<f64> = ou_rep.iter().map(|r| r.bias_theta).collect();
    let bias_knd: Vec<f64> = ou_rep.iter().map(|r| r.kendall_approx_bias).collect();
    svg::write_line_chart(
        &fig_dir.join("day004i_ou-theta-bias-vs-t.svg"),
        "Bias MLE theta_hat vs T (R=20.000) vs Aproksimasi Kendall",
        "Horizon T (Tahun)",
        "Bias E[theta_hat] - theta",
        &[
            LineSeries {
                label: "Bias Empiris (MLE)",
                x: &ts_f,
                y: &bias_obs,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Kendall (1+3*phi)/(T*phi)",
                x: &ts_f,
                y: &bias_knd,
                color: PALETTE[1],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 3. R(a) curve around a* (+/- 50%)
    let a_grid: Vec<f64> = (50..=150)
        .step_by(5)
        .map(|pct| a_foc * (pct as f64 / 100.0))
        .collect();
    let r_curve: Vec<f64> = a_grid
        .iter()
        .map(|&a| expected_return_rate(a, p.theta, p.sigma, c_cost))
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004i_bertram-return-rate.svg"),
        "Laju Laba Teoretis Bertram R(a) di Sekitar Ambang Optimal a*",
        "Ambang Bawah a",
        "Laba per Tahun R(a)",
        &[LineSeries {
            label: "R(a) (c=0.002)",
            x: &a_grid,
            y: &r_curve,
            color: PALETTE[0],
            dashed: false,
        }],
    )
    .unwrap();

    // 4. EG Power vs theta for T in {2, 5, 10}
    let ths_plot = [2.0_f64, 8.0, 16.0, 32.0];
    let pwr_2: Vec<f64> = eg_power
        .iter()
        .filter(|r| (r.t_years - 2.0).abs() < 1e-9)
        .map(|r| r.power)
        .collect();
    let pwr_5: Vec<f64> = eg_power
        .iter()
        .filter(|r| (r.t_years - 5.0).abs() < 1e-9)
        .map(|r| r.power)
        .collect();
    let pwr_10: Vec<f64> = eg_power
        .iter()
        .filter(|r| (r.t_years - 10.0).abs() < 1e-9)
        .map(|r| r.power)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004i_eg-power-vs-theta.svg"),
        "Power Uji Engle-Granger (alpha=5%) terhadap Mean-Reversion theta",
        "Kecepatan Mean-Reversion theta",
        "Power Uji EG",
        &[
            LineSeries {
                label: "T = 2 yr",
                x: &ths_plot,
                y: &pwr_2,
                color: PALETTE[1],
                dashed: false,
            },
            LineSeries {
                label: "T = 5 yr",
                x: &ths_plot,
                y: &pwr_5,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "T = 10 yr",
                x: &ths_plot,
                y: &pwr_10,
                color: PALETTE[2],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 5. Rolling profit & penalty vs window W
    let ws_f: Vec<f64> = bt
        .window_rows
        .iter()
        .map(|r| r.window_years as f64)
        .collect();
    let prof_un: Vec<f64> = bt.window_rows.iter().map(|r| r.profit_per_yr).collect();
    let prof_fl: Vec<f64> = bt
        .window_rows
        .iter()
        .map(|r| r.filtered_profit_per_yr)
        .collect();
    let ora_line = vec![bt.oracle_profit_per_yr; ws_f.len()];
    svg::write_line_chart(
        &fig_dir.join("day004i_rolling-profit-vs-window.svg"),
        "Laba OOS per Tahun vs Jendela Estimasi W (R=2.000 x 10 yr)",
        "Jendela Estimasi W (Tahun)",
        "Laba per Tahun",
        &[
            LineSeries {
                label: "Oracle (Parameter Benar)",
                x: &ws_f,
                y: &ora_line,
                color: PALETTE[2],
                dashed: true,
            },
            LineSeries {
                label: "Rolling Tanpa Filter",
                x: &ws_f,
                y: &prof_un,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Rolling + Filter EG (p<0.05)",
                x: &ws_f,
                y: &prof_fl,
                color: PALETTE[1],
                dashed: false,
            },
        ],
    )
    .unwrap();

    // 6. Cost curve profit vs c
    let cs_f: Vec<f64> = bt.cost_curve.iter().map(|r| r.c).collect();
    let c_ora: Vec<f64> = bt
        .cost_curve
        .iter()
        .map(|r| r.profit_per_yr_oracle)
        .collect();
    let c_est: Vec<f64> = bt
        .cost_curve
        .iter()
        .map(|r| r.profit_per_yr_est_w2)
        .collect();
    svg::write_line_chart(
        &fig_dir.join("day004i_cost-sensitivity.svg"),
        "Sensitivitas Laba per Tahun terhadap Biaya Transaksi c",
        "Biaya per Siklus c",
        "Laba per Tahun",
        &[
            LineSeries {
                label: "Oracle",
                x: &cs_f,
                y: &c_ora,
                color: PALETTE[0],
                dashed: false,
            },
            LineSeries {
                label: "Estimasi W=2 yr",
                x: &cs_f,
                y: &c_est,
                color: PALETTE[1],
                dashed: true,
            },
        ],
    )
    .unwrap();

    // 7. FPT Closed-form vs BVP vs MC bar chart
    svg::write_bar_chart(
        &fig_dir.join("day004i_fpt-comparison.svg"),
        "Perbandingan E[tau(-SD -> +SD)]: Tertutup vs BVP vs MC (R=200k)",
        "E[tau] (Tahun)",
        &["Tertutup (erfi)", "BVP Beda Hingga", "MC Fine-Step"],
        &[fpt_rows[0].2, fpt_rows[0].3, fpt_rows[0].5],
        PALETTE[0],
    )
    .unwrap();

    let part5_secs = t5.elapsed().as_secs_f64();
    let total_secs = t_total.elapsed().as_secs_f64();

    // ==================== Criteria Evaluation (pre-set thresholds) ====================
    let c1_mle_exact = mle_max_diff < 1e-4 && (emp_acf1 - acf1_exact).abs() < 0.01;
    let c2_fpt_matches = fpt_rows
        .iter()
        .all(|r| (r.2 - r.3).abs() / r.2 < 1e-3 && r.7 <= 3.5);
    let c3_opt_threshold_and_mc =
        (a_dir - a_foc).abs() < 1e-6 && mc_200yr.z_trades <= 3.0 && mc_200yr.z_profit <= 3.0;
    let c4_df_and_eg_size =
        (eg_dist.size_at_own_q05 - 0.05).abs() < 0.015 && spurious.frac_abs_t_gt_196 > 0.50;
    let c5_johansen_power = (johansen.spurious_size_at_q95 - 0.05).abs() < 0.02
        && johansen.coint_power_theta16_2yr > 0.90;
    let c6_rolling_est_penalty =
        bt.window_rows[2].est_penalty_vs_oracle < bt.window_rows[0].est_penalty_vs_oracle;
    let c7_spurious_filter_nonpos =
        bt.spurious_w2.filtered_profit_per_yr <= 3.0 * bt.spurious_w2.se_filtered_profit;
    let c8_runtime = total_secs < 300.0 || quick;

    let crit_str = |ok: bool| if ok { "tercapai" } else { "belum tercapai" };

    let doc = JsonVal::Obj(vec![
        (
            "metadata",
            JsonVal::Obj(vec![
                ("session", "day004i".into()),
                ("date_wita", "2026-11-22".into()),
                ("language", "rust".into()),
                ("quick_mode", quick.into()),
                ("workers_used", n_workers.into()),
                ("main_seed", (MAIN_SEED as usize).into()),
                (
                    "elapsed_seconds",
                    JsonVal::Obj(vec![
                        ("part1_ou_mle", part1_secs.into()),
                        ("part2_fpt_bertram", part2_secs.into()),
                        ("part3_unitroot_coint", part3_secs.into()),
                        ("part4_backtest", part4_secs.into()),
                        ("part5_figures", part5_secs.into()),
                        ("total", total_secs.into()),
                    ]),
                ),
            ]),
        ),
        (
            "part1_ou_estimation",
            JsonVal::Obj(vec![
                ("mle_closed_vs_num_max_diff", mle_max_diff.into()),
                ("ll_exact", ll_exact.into()),
                ("ll_euler", ll_euler.into()),
                ("stationary_var_exact", var_stat_exact.into()),
                ("stationary_var_emp", emp_var_stat.into()),
                ("acf1_exact", acf1_exact.into()),
                ("acf1_emp", emp_acf1.into()),
                (
                    "replication",
                    ou_rep
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("t_years", r.t_years.into()),
                                ("n_rep", r.n_rep.into()),
                                ("bias_theta", r.bias_theta.into()),
                                ("sd_theta", r.sd_theta.into()),
                                ("rmse_theta", r.rmse_theta.into()),
                                ("cov95_theta", r.cov95_theta.into()),
                                ("kendall_approx_bias", r.kendall_approx_bias.into()),
                                ("bias_ratio_to_kendall", r.bias_ratio_to_kendall.into()),
                                ("bias_mu", r.bias_mu.into()),
                                ("sd_mu", r.sd_mu.into()),
                                ("cov95_mu", r.cov95_mu.into()),
                                ("bias_sigma", r.bias_sigma.into()),
                                ("sd_sigma", r.sd_sigma.into()),
                                ("cov95_sigma", r.cov95_sigma.into()),
                                ("half_life_mean", r.half_life_mean.into()),
                                ("half_life_skewness", r.half_life_skewness.into()),
                                ("prob_theta_nonpos", r.prob_theta_nonpos.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part2_fpt_and_strategy",
            JsonVal::Obj(vec![
                (
                    "fpt_benchmarks",
                    fpt_rows
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("a", r.0.into()),
                                ("m", r.1.into()),
                                ("exact", r.2.into()),
                                ("bvp", r.3.into()),
                                ("sd_bvp", r.4.into()),
                                ("mc", r.5.into()),
                                ("mc_se", r.6.into()),
                                ("z_mc", r.7.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("a_star_direct", a_dir.into()),
                ("a_star_foc", a_foc.into()),
                ("a_star_in_stat_sd", (-a_foc / sd_stat).into()),
                ("expected_return_rate_a_star", r_dir.into()),
                ("theoretical_sharpe_a_star", sharpe_a_star.into()),
                ("mc_200yr_trades_per_yr", mc_200yr.trades_per_year_mc.into()),
                ("mc_200yr_se_trades", mc_200yr.se_trades_per_year.into()),
                (
                    "mc_200yr_exact_trades",
                    mc_200yr.trades_per_year_exact.into(),
                ),
                ("mc_200yr_z_trades", mc_200yr.z_trades.into()),
                ("mc_200yr_profit_per_yr", mc_200yr.profit_per_year_mc.into()),
                ("mc_200yr_se_profit", mc_200yr.se_profit_per_year.into()),
                (
                    "mc_200yr_exact_profit",
                    mc_200yr.profit_per_year_exact.into(),
                ),
                ("mc_200yr_z_profit", mc_200yr.z_profit.into()),
                (
                    "sensitivity_theta_c",
                    sens_map
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("theta", r.0.into()),
                                ("c", r.1.into()),
                                ("a_star", r.2.into()),
                                ("return_rate", r.3.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
            ]),
        ),
        (
            "part3_unitroot_cointegration",
            JsonVal::Obj(vec![
                (
                    "dickey_fuller",
                    JsonVal::Obj(vec![
                        ("n_paths", df_dist.n_paths.into()),
                        ("q01", df_dist.q01.into()),
                        ("se_q01", df_dist.se_q01.into()),
                        ("q05", df_dist.q05.into()),
                        ("se_q05", df_dist.se_q05.into()),
                        ("q10", df_dist.q10.into()),
                        ("se_q10", df_dist.se_q10.into()),
                        ("brownian_q05", df_dist.brownian_q05.into()),
                    ]),
                ),
                (
                    "engle_granger",
                    JsonVal::Obj(vec![
                        ("n_paths", eg_dist.n_paths.into()),
                        ("q01", eg_dist.q01.into()),
                        ("se_q01", eg_dist.se_q01.into()),
                        ("q05", eg_dist.q05.into()),
                        ("se_q05", eg_dist.se_q05.into()),
                        ("q10", eg_dist.q10.into()),
                        ("se_q10", eg_dist.se_q10.into()),
                        ("size_at_own_q05", eg_dist.size_at_own_q05.into()),
                        ("size_if_using_df_q05", eg_dist.size_if_using_df_q05.into()),
                    ]),
                ),
                (
                    "spurious_regression",
                    JsonVal::Obj(vec![
                        ("frac_abs_t_gt_196", spurious.frac_abs_t_gt_196.into()),
                        ("mean_r2", spurious.mean_r2.into()),
                        ("mean_dw", spurious.mean_dw.into()),
                    ]),
                ),
                (
                    "eg_power_grid",
                    eg_power
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("theta", r.theta.into()),
                                ("t_years", r.t_years.into()),
                                ("power", r.power.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "johansen_trace",
                    JsonVal::Obj(vec![
                        ("n_paths", johansen.n_paths.into()),
                        ("t_len", johansen.t_len.into()),
                        ("q90", johansen.q90.into()),
                        ("q95", johansen.q95.into()),
                        ("se_q95", johansen.se_q95.into()),
                        ("q99", johansen.q99.into()),
                        ("spurious_size_at_q95", johansen.spurious_size_at_q95.into()),
                        (
                            "coint_power_theta16_2yr",
                            johansen.coint_power_theta16_2yr.into(),
                        ),
                    ]),
                ),
            ]),
        ),
        (
            "part4_backtest",
            JsonVal::Obj(vec![
                ("oracle_profit_per_yr", bt.oracle_profit_per_yr.into()),
                ("oracle_se", bt.oracle_se.into()),
                ("oracle_sharpe", bt.oracle_sharpe.into()),
                ("oracle_trades_per_yr", bt.oracle_trades_per_yr.into()),
                ("oracle_max_dd", bt.oracle_max_dd.into()),
                (
                    "window_rows",
                    bt.window_rows
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("window_years", r.window_years.into()),
                                ("profit_per_yr", r.profit_per_yr.into()),
                                ("se_profit", r.se_profit.into()),
                                ("sharpe", r.sharpe.into()),
                                ("trades_per_yr", r.trades_per_yr.into()),
                                ("mean_max_dd", r.mean_max_dd.into()),
                                ("est_penalty_vs_oracle", r.est_penalty_vs_oracle.into()),
                                ("se_penalty", r.se_penalty.into()),
                                ("filtered_profit_per_yr", r.filtered_profit_per_yr.into()),
                                ("filtered_sharpe", r.filtered_sharpe.into()),
                                ("filtered_trades_per_yr", r.filtered_trades_per_yr.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                (
                    "spurious_filter_w2",
                    JsonVal::Obj(vec![
                        (
                            "unfiltered_profit_per_yr",
                            bt.spurious_w2.unfiltered_profit_per_yr.into(),
                        ),
                        (
                            "filtered_profit_per_yr",
                            bt.spurious_w2.filtered_profit_per_yr.into(),
                        ),
                        (
                            "se_filtered_profit",
                            bt.spurious_w2.se_filtered_profit.into(),
                        ),
                        (
                            "filtered_trade_frac",
                            bt.spurious_w2.filtered_trade_frac.into(),
                        ),
                    ]),
                ),
                (
                    "structural_break",
                    JsonVal::Obj(vec![
                        ("profit_no_stop", bt.structural_break.profit_no_stop.into()),
                        (
                            "profit_with_stop",
                            bt.structural_break.profit_with_stop.into(),
                        ),
                        ("max_dd_no_stop", bt.structural_break.max_dd_no_stop.into()),
                        (
                            "max_dd_with_stop",
                            bt.structural_break.max_dd_with_stop.into(),
                        ),
                        (
                            "mean_days_to_eg_detect",
                            bt.structural_break.mean_days_to_eg_detect.into(),
                        ),
                        ("frac_detected", bt.structural_break.frac_detected.into()),
                    ]),
                ),
                (
                    "cost_curve",
                    bt.cost_curve
                        .iter()
                        .map(|r| {
                            JsonVal::Obj(vec![
                                ("c", r.c.into()),
                                ("a_star", r.a_star.into()),
                                ("profit_per_yr_oracle", r.profit_per_yr_oracle.into()),
                                ("profit_per_yr_est_w2", r.profit_per_yr_est_w2.into()),
                            ])
                        })
                        .collect::<Vec<_>>()
                        .into(),
                ),
                ("breakeven_cost_approx", bt.breakeven_cost_approx.into()),
            ]),
        ),
        ("run_completed_without_errors", true.into()),
        (
            "criteria",
            JsonVal::Obj(vec![
                (
                    "c1_ou_mle_and_stationary_moments",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c1_mle_exact).into()),
                        ("mle_max_diff", mle_max_diff.into()),
                        ("acf1_diff", (emp_acf1 - acf1_exact).abs().into()),
                    ]),
                ),
                (
                    "c2_fpt_closed_bvp_mc_agree",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c2_fpt_matches).into()),
                        ("sym_exact", fpt_rows[0].2.into()),
                        ("sym_bvp", fpt_rows[0].3.into()),
                        ("sym_mc", fpt_rows[0].5.into()),
                        ("sym_z_mc", fpt_rows[0].7.into()),
                    ]),
                ),
                (
                    "c3_optimal_threshold_and_200yr_mc",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c3_opt_threshold_and_mc).into()),
                        ("a_star_diff", (a_dir - a_foc).abs().into()),
                        ("z_trades", mc_200yr.z_trades.into()),
                        ("z_profit", mc_200yr.z_profit.into()),
                    ]),
                ),
                (
                    "c4_df_eg_size_and_spurious_rejection",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c4_df_and_eg_size).into()),
                        ("eg_size_at_own_q05", eg_dist.size_at_own_q05.into()),
                        ("spurious_abs_t_gt_196", spurious.frac_abs_t_gt_196.into()),
                    ]),
                ),
                (
                    "c5_johansen_size_and_power",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c5_johansen_power).into()),
                        ("johansen_q95", johansen.q95.into()),
                        ("spurious_size", johansen.spurious_size_at_q95.into()),
                        ("coint_power_2yr", johansen.coint_power_theta16_2yr.into()),
                    ]),
                ),
                (
                    "c6_rolling_penalty_shrinks_with_w",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c6_rolling_est_penalty).into()),
                        ("penalty_w1", bt.window_rows[0].est_penalty_vs_oracle.into()),
                        ("penalty_w5", bt.window_rows[2].est_penalty_vs_oracle.into()),
                    ]),
                ),
                (
                    "c7_spurious_filtered_profit_nonpos",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c7_spurious_filter_nonpos).into()),
                        (
                            "spurious_filtered_profit",
                            bt.spurious_w2.filtered_profit_per_yr.into(),
                        ),
                        ("spurious_se", bt.spurious_w2.se_filtered_profit.into()),
                    ]),
                ),
                (
                    "c8_runtime_under_300s",
                    JsonVal::Obj(vec![
                        ("status", crit_str(c8_runtime).into()),
                        ("total_seconds", total_secs.into()),
                        ("budget_seconds", 300.0_f64.into()),
                    ]),
                ),
            ]),
        ),
    ]);

    write_json_file(&out_json, &doc).expect("failed to write results.json");
    println!(
        "day004i selesai dalam {:.2}s | workers={n_workers} | a*={a_foc:.5} | R(a*)={r_dir:.4} | DF_q05={:.3} | EG_q05={:.3} | Joh_q95={:.3} | c1..c8={}/{}/{}/{}/{}/{}/{}/{}",
        total_secs,
        df_dist.q05,
        eg_dist.q05,
        johansen.q95,
        c1_mle_exact,
        c2_fpt_matches,
        c3_opt_threshold_and_mc,
        c4_df_and_eg_size,
        c5_johansen_power,
        c6_rolling_est_penalty,
        c7_spurious_filter_nonpos,
        c8_runtime
    );
}
