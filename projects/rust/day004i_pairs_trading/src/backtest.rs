//! Bagian 4: Rolling-window trading backtest (W in {1, 2, 5} yr, R=2000 10-yr OOS paths),
//! estimation-error penalty vs oracle, EG test filter (cointegrated vs spurious),
//! structural break detection + stop-loss, and cost sensitivity c in {0, 0.002, 0.005}.

use crate::ou::{mle_closed_form, simulate_ou_path, OuParams};
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;
use crate::strategy::optimal_threshold_foc;
use crate::unitroot::dickey_fuller_t_stat;

#[derive(Clone, Debug)]
pub struct BacktestWindowRow {
    pub window_years: usize,
    pub profit_per_yr: f64,
    pub se_profit: f64,
    pub sharpe: f64,
    pub trades_per_yr: f64,
    pub mean_max_dd: f64,
    pub est_penalty_vs_oracle: f64,
    pub se_penalty: f64,
    pub filtered_profit_per_yr: f64,
    pub filtered_sharpe: f64,
    pub filtered_trades_per_yr: f64,
}

#[derive(Clone, Debug)]
pub struct SpuriousFilterSummary {
    pub unfiltered_profit_per_yr: f64,
    pub filtered_profit_per_yr: f64,
    pub se_filtered_profit: f64,
    pub filtered_trade_frac: f64,
}

#[derive(Clone, Debug)]
pub struct StructuralBreakSummary {
    pub profit_no_stop: f64,
    pub profit_with_stop: f64,
    pub max_dd_no_stop: f64,
    pub max_dd_with_stop: f64,
    pub mean_days_to_eg_detect: f64,
    pub frac_detected: f64,
}

#[derive(Clone, Debug)]
pub struct CostRow {
    pub c: f64,
    pub a_star: f64,
    pub profit_per_yr_oracle: f64,
    pub profit_per_yr_est_w2: f64,
}

#[derive(Clone, Debug)]
pub struct BacktestStudySummary {
    pub oracle_profit_per_yr: f64,
    pub oracle_se: f64,
    pub oracle_sharpe: f64,
    pub oracle_trades_per_yr: f64,
    pub oracle_max_dd: f64,
    pub window_rows: Vec<BacktestWindowRow>,
    pub spurious_w2: SpuriousFilterSummary,
    pub structural_break: StructuralBreakSummary,
    pub cost_curve: Vec<CostRow>,
    pub breakeven_cost_approx: f64,
}

fn trade_one_year(
    spread_yr: &[f64],
    mu_hat: f64,
    a_star: f64,
    c: f64,
    stop_dist: Option<f64>,
    state_pos: &mut i8, // 0 = flat, +1 = long spread, -1 = short spread
) -> (f64, usize, f64) {
    let mut pnl = 0.0_f64;
    let mut trades = 0_usize;
    let mut peak = 0.0_f64;
    let mut max_dd = 0.0_f64;
    let mut prev_x = spread_yr[0];

    for &x in spread_yr.iter().skip(1) {
        // Daily mark-to-market PnL from previous position
        if *state_pos != 0 {
            pnl += (*state_pos as f64) * (x - prev_x);
        }
        prev_x = x;

        if let Some(sd_lim) = stop_dist {
            if *state_pos != 0 && (x - mu_hat).abs() > sd_lim {
                pnl -= 0.5 * c;
                trades += 1;
                *state_pos = 0;
            }
        }
        if *state_pos == 0 {
            if x <= mu_hat + a_star {
                *state_pos = 1;
                pnl -= 0.5 * c;
            } else if x >= mu_hat - a_star {
                *state_pos = -1;
                pnl -= 0.5 * c;
            }
        } else if *state_pos == 1 && x >= mu_hat - a_star {
            pnl -= c; // close long + open short
            trades += 1;
            *state_pos = -1;
        } else if *state_pos == -1 && x <= mu_hat + a_star {
            pnl -= c; // close short + open long
            trades += 1;
            *state_pos = 1;
        }
        if pnl > peak {
            peak = pnl;
        }
        let dd = peak - pnl;
        if dd > max_dd {
            max_dd = dd;
        }
    }
    (pnl, trades, max_dd)
}

pub fn run_backtest_study(
    p: OuParams,
    c_cost: f64,
    eg_q05: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> BacktestStudySummary {
    let burn_yrs = 5_usize;
    let oos_yrs = 10_usize;
    let days_per_yr = 252_usize;
    let total_steps = (burn_yrs + oos_yrs) * days_per_yr;
    let a_oracle = optimal_threshold_foc(p.theta, p.sigma, c_cost);

    // Pre-run R=n_paths 15-year paths and compute oracle + W in {1, 2, 5}
    struct PathMetrics {
        ora_pnl: f64,
        ora_tr: f64,
        ora_dd: f64,
        ora_sh: f64,
        w_pnl: [f64; 3],
        w_tr: [f64; 3],
        w_dd: [f64; 3],
        w_sh: [f64; 3],
        wf_pnl: [f64; 3],
        wf_tr: [f64; 3],
        wf_sh: [f64; 3],
    }

    let ws = [1_usize, 2, 5];
    let rows = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x4100, r as u64));
        let spread = simulate_ou_path(p, total_steps, true, &mut rg);
        let oos_start = burn_yrs * days_per_yr;

        // Oracle
        let mut st_ora = 0_i8;
        let mut ora_yr_pnls = Vec::with_capacity(oos_yrs);
        let mut ora_tr_tot = 0.0_f64;
        let mut ora_dd_max = 0.0_f64;
        for y in 0..oos_yrs {
            let s0 = oos_start + y * days_per_yr;
            let (pnl, tr, dd) = trade_one_year(
                &spread[s0..s0 + days_per_yr],
                p.mu,
                a_oracle,
                c_cost,
                None,
                &mut st_ora,
            );
            ora_yr_pnls.push(pnl);
            ora_tr_tot += tr as f64;
            ora_dd_max = ora_dd_max.max(dd);
        }
        let ora_m = stats::mean(&ora_yr_pnls);
        let ora_s = stats::sd(&ora_yr_pnls);

        let mut w_pnl = [0.0_f64; 3];
        let mut w_tr = [0.0_f64; 3];
        let mut w_dd = [0.0_f64; 3];
        let mut w_sh = [0.0_f64; 3];
        let mut wf_pnl = [0.0_f64; 3];
        let mut wf_tr = [0.0_f64; 3];
        let mut wf_sh = [0.0_f64; 3];

        for (wi, &w_yr) in ws.iter().enumerate() {
            let mut st_est = 0_i8;
            let mut st_filt = 0_i8;
            let mut pnls = Vec::with_capacity(oos_yrs);
            let mut f_pnls = Vec::with_capacity(oos_yrs);
            let mut tr_tot = 0.0_f64;
            let mut f_tr_tot = 0.0_f64;
            let mut dd_max = 0.0_f64;
            for y in 0..oos_yrs {
                let s0 = oos_start + y * days_per_yr;
                let win = &spread[s0 - w_yr * days_per_yr..s0];
                let mle = mle_closed_form(win, p.dt);
                let (pnl, tr, dd) = if mle.theta.is_finite() && mle.theta > 0.1 && mle.sigma > 1e-4
                {
                    let a_hat = optimal_threshold_foc(mle.theta, mle.sigma, c_cost);
                    trade_one_year(
                        &spread[s0..s0 + days_per_yr],
                        mle.mu,
                        a_hat,
                        c_cost,
                        None,
                        &mut st_est,
                    )
                } else {
                    st_est = 0;
                    (0.0, 0, 0.0)
                };
                pnls.push(pnl);
                tr_tot += tr as f64;
                dd_max = dd_max.max(dd);

                // EG filter on trailing window
                let eg_t = dickey_fuller_t_stat(win);
                if eg_t < eg_q05 && mle.theta.is_finite() && mle.theta > 0.1 && mle.sigma > 1e-4 {
                    let a_hat = optimal_threshold_foc(mle.theta, mle.sigma, c_cost);
                    let (fp, ft, _) = trade_one_year(
                        &spread[s0..s0 + days_per_yr],
                        mle.mu,
                        a_hat,
                        c_cost,
                        None,
                        &mut st_filt,
                    );
                    f_pnls.push(fp);
                    f_tr_tot += ft as f64;
                } else {
                    st_filt = 0;
                    f_pnls.push(0.0);
                }
            }
            w_pnl[wi] = stats::mean(&pnls);
            let sd_p = stats::sd(&pnls);
            w_sh[wi] = if sd_p > 0.0 { w_pnl[wi] / sd_p } else { 0.0 };
            w_tr[wi] = tr_tot / oos_yrs as f64;
            w_dd[wi] = dd_max;

            wf_pnl[wi] = stats::mean(&f_pnls);
            let sd_fp = stats::sd(&f_pnls);
            wf_sh[wi] = if sd_fp > 0.0 { wf_pnl[wi] / sd_fp } else { 0.0 };
            wf_tr[wi] = f_tr_tot / oos_yrs as f64;
        }

        PathMetrics {
            ora_pnl: ora_m,
            ora_tr: ora_tr_tot / oos_yrs as f64,
            ora_dd: ora_dd_max,
            ora_sh: if ora_s > 0.0 { ora_m / ora_s } else { 0.0 },
            w_pnl,
            w_tr,
            w_dd,
            w_sh,
            wf_pnl,
            wf_tr,
            wf_sh,
        }
    });

    let ora_pnls: Vec<f64> = rows.iter().map(|r| r.ora_pnl).collect();
    let ora_trs: Vec<f64> = rows.iter().map(|r| r.ora_tr).collect();
    let ora_dds: Vec<f64> = rows.iter().map(|r| r.ora_dd).collect();
    let ora_shs: Vec<f64> = rows.iter().map(|r| r.ora_sh).collect();

    let mut window_rows = Vec::with_capacity(3);
    for (wi, &w_yr) in ws.iter().enumerate() {
        let pnls: Vec<f64> = rows.iter().map(|r| r.w_pnl[wi]).collect();
        let pens: Vec<f64> = rows.iter().map(|r| r.ora_pnl - r.w_pnl[wi]).collect();
        let shs: Vec<f64> = rows.iter().map(|r| r.w_sh[wi]).collect();
        let trs: Vec<f64> = rows.iter().map(|r| r.w_tr[wi]).collect();
        let dds: Vec<f64> = rows.iter().map(|r| r.w_dd[wi]).collect();
        let f_pnls: Vec<f64> = rows.iter().map(|r| r.wf_pnl[wi]).collect();
        let f_shs: Vec<f64> = rows.iter().map(|r| r.wf_sh[wi]).collect();
        let f_trs: Vec<f64> = rows.iter().map(|r| r.wf_tr[wi]).collect();
        window_rows.push(BacktestWindowRow {
            window_years: w_yr,
            profit_per_yr: stats::mean(&pnls),
            se_profit: stats::se(&pnls),
            sharpe: stats::mean(&shs),
            trades_per_yr: stats::mean(&trs),
            mean_max_dd: stats::mean(&dds),
            est_penalty_vs_oracle: stats::mean(&pens),
            se_penalty: stats::se(&pens),
            filtered_profit_per_yr: stats::mean(&f_pnls),
            filtered_sharpe: stats::mean(&f_shs),
            filtered_trades_per_yr: stats::mean(&f_trs),
        });
    }

    // Spurious pairs (two independent random walks, residual spread, W=2 yr)
    let sp_rows = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x4200, r as u64));
        let dt = p.dt;
        let sqrt_dt = dt.sqrt();
        let mut x = Vec::with_capacity(total_steps + 1);
        let mut y = Vec::with_capacity(total_steps + 1);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..=total_steps {
            cx += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            cy += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            x.push(cx);
            y.push(cy);
        }
        let oos_start = burn_yrs * days_per_yr;
        let mut st_un = 0_i8;
        let mut st_fl = 0_i8;
        let mut un_pnls = Vec::with_capacity(oos_yrs);
        let mut fl_pnls = Vec::with_capacity(oos_yrs);
        let mut fl_active = 0_usize;
        for yr in 0..oos_yrs {
            let s0 = oos_start + yr * days_per_yr;
            let w0 = s0 - 2 * days_per_yr;
            let ols = stats::ols_simple(&x[w0..s0], &y[w0..s0]);
            let mle = mle_closed_form(&ols.residuals, dt);
            let oos_spread: Vec<f64> = (s0..s0 + days_per_yr)
                .map(|t| y[t] - (ols.alpha + ols.beta * x[t]))
                .collect();
            let a_hat = if mle.theta.is_finite() && mle.theta > 0.1 && mle.sigma > 1e-4 {
                optimal_threshold_foc(mle.theta, mle.sigma, c_cost)
            } else {
                -0.05
            };
            let (up, _, _) = trade_one_year(&oos_spread, mle.mu, a_hat, c_cost, None, &mut st_un);
            un_pnls.push(up);

            let eg_t = dickey_fuller_t_stat(&ols.residuals);
            if eg_t < eg_q05 && mle.theta.is_finite() && mle.theta > 0.1 {
                let (fp, _, _) =
                    trade_one_year(&oos_spread, mle.mu, a_hat, c_cost, None, &mut st_fl);
                fl_pnls.push(fp);
                fl_active += 1;
            } else {
                st_fl = 0;
                fl_pnls.push(0.0);
            }
        }
        (
            stats::mean(&un_pnls),
            stats::mean(&fl_pnls),
            fl_active as f64 / oos_yrs as f64,
        )
    });
    let sp_un: Vec<f64> = sp_rows.iter().map(|r| r.0).collect();
    let sp_fl: Vec<f64> = sp_rows.iter().map(|r| r.1).collect();
    let sp_fr: Vec<f64> = sp_rows.iter().map(|r| r.2).collect();
    let spurious_w2 = SpuriousFilterSummary {
        unfiltered_profit_per_yr: stats::mean(&sp_un),
        filtered_profit_per_yr: stats::mean(&sp_fl),
        se_filtered_profit: stats::se(&sp_fl),
        filtered_trade_frac: stats::mean(&sp_fr),
    };

    // Structural break: theta = 16 for years 0..5, then drops to theta = 2 for years 5..10
    let p_pre = p;
    let p_post = OuParams {
        theta: 2.0,
        mu: p.mu,
        sigma: p.sigma,
        dt: p.dt,
    };
    let sb_rows = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x4300, r as u64));
        let pre = simulate_ou_path(p_pre, 5 * days_per_yr, true, &mut rg);
        let mut post = Vec::with_capacity(5 * days_per_yr);
        let phi2 = (-p_post.theta * p_post.dt).exp();
        let sd2 = (p_post.sigma * p_post.sigma * (1.0 - phi2 * phi2) / (2.0 * p_post.theta)).sqrt();
        let mut cur = *pre.last().unwrap();
        for _ in 0..(5 * days_per_yr) {
            cur = cur * phi2 + sd2 * rg.normal();
            post.push(cur);
        }
        let mut full = pre;
        full.extend_from_slice(&post);
        let brk = 5 * days_per_yr;

        let mut st_ns = 0_i8;
        let mut st_sl = 0_i8;
        let (pnl_ns, _, dd_ns) =
            trade_one_year(&full[brk..], p.mu, a_oracle, c_cost, None, &mut st_ns);
        let stop_lim = 4.0 * p.stationary_sd();
        let (pnl_sl, _, dd_sl) = trade_one_year(
            &full[brk..],
            p.mu,
            a_oracle,
            c_cost,
            Some(stop_lim),
            &mut st_sl,
        );

        // Days after break until rolling 2-year EG test fails to reject (eg_t >= eg_q05)
        let mut detect_day = f64::NAN;
        for d in (0..(5 * days_per_yr)).step_by(5) {
            let end_idx = brk + d;
            let win = &full[end_idx - 2 * days_per_yr..end_idx];
            if dickey_fuller_t_stat(win) >= eg_q05 {
                detect_day = d as f64;
                break;
            }
        }
        (pnl_ns / 5.0, pnl_sl / 5.0, dd_ns, dd_sl, detect_day)
    });
    let det_days: Vec<f64> = sb_rows
        .iter()
        .filter(|r| r.4.is_finite())
        .map(|r| r.4)
        .collect();
    let structural_break = StructuralBreakSummary {
        profit_no_stop: stats::mean(&sb_rows.iter().map(|r| r.0).collect::<Vec<_>>()),
        profit_with_stop: stats::mean(&sb_rows.iter().map(|r| r.1).collect::<Vec<_>>()),
        max_dd_no_stop: stats::mean(&sb_rows.iter().map(|r| r.2).collect::<Vec<_>>()),
        max_dd_with_stop: stats::mean(&sb_rows.iter().map(|r| r.3).collect::<Vec<_>>()),
        mean_days_to_eg_detect: if det_days.is_empty() {
            f64::NAN
        } else {
            stats::mean(&det_days)
        },
        frac_detected: det_days.len() as f64 / n_paths as f64,
    };

    // Cost sensitivity c in {0, 0.002, 0.005}
    let mut cost_curve = Vec::new();
    for &c_val in &[0.0, 0.002, 0.005] {
        let a_c = if c_val == 0.0 {
            -0.5 * p.stationary_sd()
        } else {
            optimal_threshold_foc(p.theta, p.sigma, c_val)
        };
        let sub_rows = parallel::run_chunks(n_paths / 2, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x4400 ^ c_val.to_bits(), r as u64));
            let spread = simulate_ou_path(p, (2 + 5) * days_per_yr, true, &mut rg);
            let s0 = 2 * days_per_yr;
            let mut st_o = 0_i8;
            let mut st_e = 0_i8;
            let (po, _, _) = trade_one_year(&spread[s0..], p.mu, a_c, c_val, None, &mut st_o);
            let mle = mle_closed_form(&spread[..s0], p.dt);
            let ae = if mle.theta.is_finite() && mle.theta > 0.1 && mle.sigma > 1e-4 {
                if c_val == 0.0 {
                    -0.5 * (mle.sigma / (2.0 * mle.theta).sqrt())
                } else {
                    optimal_threshold_foc(mle.theta, mle.sigma, c_val)
                }
            } else {
                a_c
            };
            let (pe, _, _) = trade_one_year(&spread[s0..], mle.mu, ae, c_val, None, &mut st_e);
            (po / 5.0, pe / 5.0)
        });
        cost_curve.push(CostRow {
            c: c_val,
            a_star: a_c,
            profit_per_yr_oracle: stats::mean(&sub_rows.iter().map(|r| r.0).collect::<Vec<_>>()),
            profit_per_yr_est_w2: stats::mean(&sub_rows.iter().map(|r| r.1).collect::<Vec<_>>()),
        });
    }
    let breakeven_cost_approx = 6.0 * p.stationary_sd();

    BacktestStudySummary {
        oracle_profit_per_yr: stats::mean(&ora_pnls),
        oracle_se: stats::se(&ora_pnls),
        oracle_sharpe: stats::mean(&ora_shs),
        oracle_trades_per_yr: stats::mean(&ora_trs),
        oracle_max_dd: stats::mean(&ora_dds),
        window_rows,
        spurious_w2,
        structural_break,
        cost_curve,
        breakeven_cost_approx,
    }
}
