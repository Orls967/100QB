//! Bagian 4: Synthetic market generation, Nelder-Mead calibration of (H, eta, rho),
//! H-profile likelihood/objective, 3x3 numerical Hessian & H-eta correlation,
//! 30-replication identifiability study, and T >= 0.5 vs T <= 0.1 misspecification.

#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

use crate::bs;
use crate::linalg;
use crate::parallel;
use crate::rbergomi::RBergomiParams;
use crate::rng::{derive_seed, Rng};
use crate::smile::{compute_smile, compute_smile_with_scheme, SmilePoint};
use crate::stats;
use crate::volterra::{exact_var, HybridScheme, HybridWeights};

/// Pre-generated CRN normal drivers `(u, v)` on `n_steps` for fast repeated calibration evaluations.
#[derive(Clone)]
pub struct CalibCrnBank {
    pub n_paths: usize,
    pub n_steps: usize,
    pub u_flat: Vec<f64>,
    pub v_flat: Vec<f64>,
}

impl CalibCrnBank {
    pub fn build(n_paths: usize, n_steps: usize, n_workers: usize, seed: u64, salt: u64) -> Self {
        let rows = parallel::run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, salt, r as u64));
            let mut u = vec![0.0_f64; n_steps];
            let mut v = vec![0.0_f64; n_steps];
            for i in 0..n_steps {
                u[i] = rg.normal();
                v[i] = rg.normal();
            }
            (u, v)
        });
        let mut u_flat = Vec::with_capacity(n_paths * n_steps);
        let mut v_flat = Vec::with_capacity(n_paths * n_steps);
        for (u, v) in rows {
            u_flat.extend_from_slice(&u);
            v_flat.extend_from_slice(&v);
        }
        Self {
            n_paths,
            n_steps,
            u_flat,
            v_flat,
        }
    }
}

/// Evaluates implied vol RMSE (in vol points) against `market_smile` for `(h, eta, rho)` on `T_max = 1.0`.
/// Maturities `ts` must be multiples of `dt = 1.0 / n_steps` (e.g. `{0.1, 0.25, 0.5, 1.0}` for `n_steps=100`).
pub fn eval_smile_sse(
    bank: &CalibCrnBank,
    xi0: f64,
    s0: f64,
    h: f64,
    eta: f64,
    rho: f64,
    market: &[SmilePoint],
    ts: &[f64],
    log_ks: &[f64],
    n_workers: usize,
) -> f64 {
    eval_smile_sse_with_scheme(
        bank,
        xi0,
        s0,
        h,
        eta,
        rho,
        market,
        ts,
        log_ks,
        n_workers,
        HybridScheme::BlpOptimal,
    )
}

pub fn eval_smile_sse_with_scheme(
    bank: &CalibCrnBank,
    xi0: f64,
    s0: f64,
    h: f64,
    eta: f64,
    rho: f64,
    market: &[SmilePoint],
    ts: &[f64],
    log_ks: &[f64],
    n_workers: usize,
    scheme: HybridScheme,
) -> f64 {
    if !(0.01..=0.495).contains(&h)
        || !(0.05..=5.0).contains(&eta)
        || !(-0.99..=0.99).contains(&rho)
    {
        return 1e6;
    }
    let n = bank.n_steps;
    let hw = HybridWeights::with_scheme(n, 1.0, h, scheme);
    let dt = hw.dt;
    let half_eta2 = 0.5 * eta * eta;
    let one_m_rho2 = (1.0 - rho * rho).max(0.0);
    let step_targets: Vec<usize> = ts.iter().map(|&t| (t / dt).round() as usize).collect();
    let n_t = ts.len();
    let n_k = log_ks.len();

    let prices_sum = parallel::run_chunk_reduce(
        bank.n_paths,
        n_workers,
        |start, end| {
            let mut local_sums = vec![0.0_f64; n_t * n_k];
            let mut wt = vec![0.0_f64; n];
            let mut db = vec![0.0_f64; n];
            for r in start..end {
                let off = r * n;
                let u = &bank.u_flat[off..off + n];
                let v = &bank.v_flat[off..off + n];
                hw.generate_path(u, v, &mut wt, &mut db);

                let mut int_v = 0.0_f64;
                let mut stoch_int = 0.0_f64;
                let mut t_ptr = 0_usize;
                for i in 0..n {
                    let v_i = if i == 0 {
                        xi0
                    } else {
                        let t_i = i as f64 * dt;
                        xi0 * (eta * wt[i - 1] - half_eta2 * exact_var(t_i, h)).exp()
                    };
                    int_v += v_i * dt;
                    stoch_int += v_i.sqrt() * db[i];

                    if t_ptr < n_t && (i + 1) == step_targets[t_ptr] {
                        let s_eff = s0 * (rho * stoch_int - 0.5 * rho * rho * int_v).exp();
                        let vc = one_m_rho2 * int_v;
                        let row_off = t_ptr * n_k;
                        for (ki, &lk) in log_ks.iter().enumerate() {
                            let strike = s0 * lk.exp();
                            local_sums[row_off + ki] += bs::bs_call_total_var(s_eff, strike, vc);
                        }
                        t_ptr += 1;
                    }
                }
            }
            local_sums
        },
        |parts| {
            let mut tot = vec![0.0_f64; n_t * n_k];
            for p in parts {
                for (a, b) in tot.iter_mut().zip(p.iter()) {
                    *a += *b;
                }
            }
            tot
        },
    );

    let nf = bank.n_paths as f64;
    let mut sse = 0.0_f64;
    for (ti, &t_h) in ts.iter().enumerate() {
        for (ki, &lk) in log_ks.iter().enumerate() {
            let call_p = prices_sum[ti * n_k + ki] / nf;
            let strike = s0 * lk.exp();
            let iv_mod = bs::implied_vol(s0, strike, t_h, call_p);
            let iv_mkt = market
                .iter()
                .find(|pt| (pt.t_horizon - t_h).abs() < 1e-12 && (pt.log_k - lk).abs() < 1e-9)
                .map(|pt| pt.iv)
                .unwrap_or(iv_mod);
            if iv_mod.is_finite() && iv_mkt.is_finite() {
                let d = iv_mod - iv_mkt;
                sse += d * d;
            } else {
                sse += 1.0;
            }
        }
    }
    sse
}

/// Nelder-Mead simplex minimizer for `f: &[f64] -> f64`.
pub fn nelder_mead<F>(
    init: &[f64],
    steps: &[f64],
    max_evals: usize,
    tol: f64,
    mut f: F,
) -> (Vec<f64>, f64, usize)
where
    F: FnMut(&[f64]) -> f64,
{
    let n = init.len();
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    simplex.push((init.to_vec(), f(init)));
    let mut evals = 1;
    for i in 0..n {
        let mut pt = init.to_vec();
        pt[i] += steps[i];
        let val = f(&pt);
        evals += 1;
        simplex.push((pt, val));
    }

    while evals < max_evals {
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let f_best = simplex[0].1;
        let f_worst = simplex[n].1;
        let mut max_dx = 0.0_f64;
        for i in 1..=n {
            for j in 0..n {
                max_dx = max_dx.max((simplex[i].0[j] - simplex[0].0[j]).abs());
            }
        }
        if (f_worst - f_best).abs() < tol && max_dx < 1e-4 {
            break;
        }

        let mut centroid = vec![0.0_f64; n];
        for item in simplex.iter().take(n) {
            for j in 0..n {
                centroid[j] += item.0[j] / n as f64;
            }
        }

        let mut xr = vec![0.0_f64; n];
        for j in 0..n {
            xr[j] = centroid[j] + (centroid[j] - simplex[n].0[j]);
        }
        let fr = f(&xr);
        evals += 1;

        if fr < simplex[0].1 {
            let mut xe = vec![0.0_f64; n];
            for j in 0..n {
                xe[j] = centroid[j] + 2.0 * (xr[j] - centroid[j]);
            }
            let fe = f(&xe);
            evals += 1;
            if fe < fr {
                simplex[n] = (xe, fe);
            } else {
                simplex[n] = (xr, fr);
            }
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (xr, fr);
        } else {
            let mut xc = vec![0.0_f64; n];
            if fr < simplex[n].1 {
                for j in 0..n {
                    xc[j] = centroid[j] + 0.5 * (xr[j] - centroid[j]);
                }
            } else {
                for j in 0..n {
                    xc[j] = centroid[j] - 0.5 * (centroid[j] - simplex[n].0[j]);
                }
            }
            let fc = f(&xc);
            evals += 1;
            if fc < simplex[n].1.min(fr) {
                simplex[n] = (xc, fc);
            } else {
                let x0 = simplex[0].0.clone();
                for item in simplex.iter_mut().skip(1) {
                    for j in 0..n {
                        item.0[j] = x0[j] + 0.5 * (item.0[j] - x0[j]);
                    }
                    item.1 = f(&item.0);
                    evals += 1;
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    (simplex[0].0.clone(), simplex[0].1, evals)
}

#[derive(Clone, Debug)]
pub struct CalibResult {
    pub h_hat: f64,
    pub eta_hat: f64,
    pub rho_hat: f64,
    pub rmse_opt: f64,
    pub rmse_at_true: f64,
    pub evals: usize,
}

#[derive(Clone, Debug)]
pub struct HProfilePoint {
    pub h: f64,
    pub eta_opt: f64,
    pub rho_opt: f64,
    pub rmse: f64,
}

#[derive(Clone, Debug)]
pub struct ReplicationSummary {
    pub n_rep: usize,
    pub bias_h: f64,
    pub sd_h: f64,
    pub bias_eta: f64,
    pub sd_eta: f64,
    pub bias_rho: f64,
    pub sd_rho: f64,
}

#[derive(Clone, Debug)]
pub struct CalibCoreSummary {
    pub main_fit: CalibResult,
    pub h_profile: Vec<HProfilePoint>,
    pub corr_h_eta: f64,
    pub replication: ReplicationSummary,
}

#[derive(Clone, Debug)]
pub struct RmseGridStats {
    pub n_steps: usize,
    pub n_paths: usize,
    pub n_seeds: usize,
    pub median_rmse: f64,
    pub iqr_rmse: f64,
    pub min_rmse: f64,
    pub max_rmse: f64,
}

#[derive(Clone, Debug)]
pub struct HillEstimate {
    pub k: usize,
    pub alpha_st: f64,
    pub alpha_payoff2: f64,
}

#[derive(Clone, Debug)]
pub struct CalibA6ExtrasSummary {
    pub self_consistency_fit: CalibResult,
    pub fit_r80k: CalibResult,
    pub fit_vega_price: CalibResult,
    pub true_rmse_grid_stats: Vec<RmseGridStats>,
    pub hill_estimates: Vec<HillEstimate>,
    pub rep_diff_crn: ReplicationSummary,
    pub rep_old_diff_n: ReplicationSummary,
    pub rep_blp_diff_n: ReplicationSummary,
    pub rep_cond_bs_fixed_crn: ReplicationSummary,
    pub rep_cond_bs_diff_crn: ReplicationSummary,
    pub fit_long_only: CalibResult,
    pub fit_short_only: CalibResult,
}

pub fn run_calibration_core(
    true_p: RBergomiParams,
    n_workers: usize,
    seed: u64,
) -> CalibCoreSummary {
    let ts_main = [0.1, 0.25, 0.5, 1.0];
    let ks_calib = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3];
    let n_pts_main = ts_main.len() * ks_calib.len();

    // 1. Independent synthetic market (R=200,000, n=500)
    let market_seed = seed ^ 0xCAFE_BABE_1234_5678;
    let market_smile = compute_smile(
        true_p,
        &ts_main,
        &ks_calib,
        500,
        200_000,
        n_workers,
        market_seed,
        0x4100,
    );

    // Fixed CRN bank for calibration (R=20,000, n=100)
    let bank = CalibCrnBank::build(20_000, 100, n_workers, seed, 0x4200);
    let sse_true = eval_smile_sse(
        &bank,
        true_p.xi0,
        true_p.s0,
        true_p.h,
        true_p.eta,
        true_p.rho,
        &market_smile,
        &ts_main,
        &ks_calib,
        n_workers,
    );

    let init = [0.22_f64, 1.45_f64, -0.55_f64];
    let steps = [-0.06_f64, 0.25_f64, -0.1_f64];
    let (opt_x, opt_sse, evals) = nelder_mead(&init, &steps, 110, 1e-10, |x| {
        eval_smile_sse(
            &bank,
            true_p.xi0,
            true_p.s0,
            x[0],
            x[1],
            x[2],
            &market_smile,
            &ts_main,
            &ks_calib,
            n_workers,
        )
    });
    let main_fit = CalibResult {
        h_hat: opt_x[0],
        eta_hat: opt_x[1],
        rho_hat: opt_x[2],
        rmse_opt: (opt_sse / n_pts_main as f64).sqrt(),
        rmse_at_true: (sse_true / n_pts_main as f64).sqrt(),
        evals,
    };

    // 2. H-profile (fixing H, optimizing eta & rho)
    let mut h_profile = Vec::new();
    for &h_grid in &[0.04, 0.06, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30] {
        let (x2, sse2, _) = nelder_mead(&[opt_x[1], opt_x[2]], &[0.15, -0.08], 55, 1e-9, |y| {
            eval_smile_sse(
                &bank,
                true_p.xi0,
                true_p.s0,
                h_grid,
                y[0],
                y[1],
                &market_smile,
                &ts_main,
                &ks_calib,
                n_workers,
            )
        });
        h_profile.push(HProfilePoint {
            h: h_grid,
            eta_opt: x2[0],
            rho_opt: x2[1],
            rmse: (sse2 / n_pts_main as f64).sqrt(),
        });
    }

    // 3. Numerical 3x3 Hessian at optimum -> H-eta correlation
    let eps_h = [0.004_f64, 0.03_f64, 0.015_f64];
    let mut hess = vec![vec![0.0_f64; 3]; 3];
    let f0 = opt_sse;
    for i in 0..3 {
        let mut xp = opt_x.clone();
        let mut xm = opt_x.clone();
        xp[i] += eps_h[i];
        xm[i] -= eps_h[i];
        let fp = eval_smile_sse(
            &bank,
            true_p.xi0,
            true_p.s0,
            xp[0],
            xp[1],
            xp[2],
            &market_smile,
            &ts_main,
            &ks_calib,
            n_workers,
        );
        let fm = eval_smile_sse(
            &bank,
            true_p.xi0,
            true_p.s0,
            xm[0],
            xm[1],
            xm[2],
            &market_smile,
            &ts_main,
            &ks_calib,
            n_workers,
        );
        hess[i][i] = (fp - 2.0 * f0 + fm) / (eps_h[i] * eps_h[i]);
        for j in (i + 1)..3 {
            let mut xpp = opt_x.clone();
            let mut xpm = opt_x.clone();
            let mut xmp = opt_x.clone();
            let mut xmm = opt_x.clone();
            xpp[i] += eps_h[i];
            xpp[j] += eps_h[j];
            xpm[i] += eps_h[i];
            xpm[j] -= eps_h[j];
            xmp[i] -= eps_h[i];
            xmp[j] += eps_h[j];
            xmm[i] -= eps_h[i];
            xmm[j] -= eps_h[j];
            let fpp = eval_smile_sse(
                &bank,
                true_p.xi0,
                true_p.s0,
                xpp[0],
                xpp[1],
                xpp[2],
                &market_smile,
                &ts_main,
                &ks_calib,
                n_workers,
            );
            let fpm = eval_smile_sse(
                &bank,
                true_p.xi0,
                true_p.s0,
                xpm[0],
                xpm[1],
                xpm[2],
                &market_smile,
                &ts_main,
                &ks_calib,
                n_workers,
            );
            let fmp = eval_smile_sse(
                &bank,
                true_p.xi0,
                true_p.s0,
                xmp[0],
                xmp[1],
                xmp[2],
                &market_smile,
                &ts_main,
                &ks_calib,
                n_workers,
            );
            let fmm = eval_smile_sse(
                &bank,
                true_p.xi0,
                true_p.s0,
                xmm[0],
                xmm[1],
                xmm[2],
                &market_smile,
                &ts_main,
                &ks_calib,
                n_workers,
            );
            let hij = (fpp - fpm - fmp + fmm) / (4.0 * eps_h[i] * eps_h[j]);
            hess[i][j] = hij;
            hess[j][i] = hij;
        }
    }
    let corr_h_eta = match linalg::invert_dense(&hess) {
        Ok(cov_est) => {
            let denom = (cov_est[0][0].abs() * cov_est[1][1].abs()).sqrt();
            if denom > 1e-18 {
                cov_est[0][1] / denom
            } else {
                0.0
            }
        }
        Err(_) => 0.0,
    };

    // 4. Three-condition replication study to separate scheme bias from pure identifiability:
    //    (a) Old scheme (RightEndpoint), diff n (n_mkt=250, n_calib=100, R_rep=15)
    //    (b) New scheme (BlpOptimal), diff n (n_mkt=250, n_calib=100, R_rep=15)
    //    (c) New scheme (BlpOptimal), SAME n (n_mkt=100, n_calib=100, R_rep=30) -> pure identifiability
    let run_rep =
        |n_rep: usize, n_mkt: usize, scheme: HybridScheme, salt_base: u64| -> ReplicationSummary {
            let mut hs_rep = Vec::with_capacity(n_rep);
            let mut etas_rep = Vec::with_capacity(n_rep);
            let mut rhos_rep = Vec::with_capacity(n_rep);
            for rep in 0..n_rep {
                let rep_seed = seed ^ ((rep as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                let rep_mkt = compute_smile_with_scheme(
                    true_p, &ts_main, &ks_calib, n_mkt, 20_000, n_workers, rep_seed, salt_base,
                    scheme,
                );
                let (xr, _, _) =
                    nelder_mead(&[0.18, 1.6, -0.6], &[-0.04, 0.2, -0.08], 105, 1e-8, |x| {
                        eval_smile_sse_with_scheme(
                            &bank, true_p.xi0, true_p.s0, x[0], x[1], x[2], &rep_mkt, &ts_main,
                            &ks_calib, n_workers, scheme,
                        )
                    });
                hs_rep.push(xr[0]);
                etas_rep.push(xr[1]);
                rhos_rep.push(xr[2]);
            }
            ReplicationSummary {
                n_rep,
                bias_h: stats::mean(&hs_rep) - true_p.h,
                sd_h: stats::sd(&hs_rep),
                bias_eta: stats::mean(&etas_rep) - true_p.eta,
                sd_eta: stats::sd(&etas_rep),
                bias_rho: stats::mean(&rhos_rep) - true_p.rho,
                sd_rho: stats::sd(&rhos_rep),
            }
        };
    let replication = run_rep(30, 100, HybridScheme::BlpOptimal, 0x4300);

    CalibCoreSummary {
        main_fit,
        h_profile,
        corr_h_eta,
        replication,
    }
}

/// Conditional Black-Scholes smile builder (high VRF) for A6 calibration check.
pub fn compute_smile_cond_bs(
    bank: &CalibCrnBank,
    p: RBergomiParams,
    ts: &[f64],
    log_ks: &[f64],
    n_workers: usize,
) -> Vec<SmilePoint> {
    let n = bank.n_steps;
    let hw = HybridWeights::with_scheme(n, 1.0, p.h, HybridScheme::BlpOptimal);
    let dt = hw.dt;
    let half_eta2 = 0.5 * p.eta * p.eta;
    let one_m_rho2 = (1.0 - p.rho * p.rho).max(0.0);
    let step_targets: Vec<usize> = ts.iter().map(|&t| (t / dt).round() as usize).collect();
    let n_t = ts.len();
    let n_k = log_ks.len();
    let sums = parallel::run_chunk_reduce(
        bank.n_paths,
        n_workers,
        |start, end| {
            let mut local = vec![0.0_f64; n_t * n_k];
            let mut wt = vec![0.0_f64; n];
            let mut db = vec![0.0_f64; n];
            for r in start..end {
                let off = r * n;
                hw.generate_path(
                    &bank.u_flat[off..off + n],
                    &bank.v_flat[off..off + n],
                    &mut wt,
                    &mut db,
                );
                let mut int_v = 0.0_f64;
                let mut stoch_int = 0.0_f64;
                let mut tp = 0_usize;
                for i in 0..n {
                    let v_i = if i == 0 {
                        p.xi0
                    } else {
                        p.xi0
                            * (p.eta * wt[i - 1] - half_eta2 * exact_var(i as f64 * dt, p.h)).exp()
                    };
                    int_v += v_i * dt;
                    stoch_int += v_i.sqrt() * db[i];
                    if tp < n_t && (i + 1) == step_targets[tp] {
                        let s_eff = p.s0 * (p.rho * stoch_int - 0.5 * p.rho * p.rho * int_v).exp();
                        let vc = one_m_rho2 * int_v;
                        for (ki, &lk) in log_ks.iter().enumerate() {
                            local[tp * n_k + ki] +=
                                bs::bs_call_total_var(s_eff, p.s0 * lk.exp(), vc);
                        }
                        tp += 1;
                    }
                }
            }
            local
        },
        |parts| {
            let mut tot = vec![0.0_f64; n_t * n_k];
            for part in parts {
                for (a, b) in tot.iter_mut().zip(part.iter()) {
                    *a += *b;
                }
            }
            tot
        },
    );
    let mut pts = Vec::with_capacity(n_t * n_k);
    for (ti, &t_h) in ts.iter().enumerate() {
        for (ki, &lk) in log_ks.iter().enumerate() {
            let cp = sums[ti * n_k + ki] / bank.n_paths as f64;
            let iv = bs::implied_vol(p.s0, p.s0 * lk.exp(), t_h, cp);
            pts.push(SmilePoint {
                t_horizon: t_h,
                log_k: lk,
                call_price: cp,
                call_se: 0.0,
                iv,
                iv_se: 0.0,
            });
        }
    }
    pts
}

fn hill_alpha(sorted_asc: &[f64], k: usize) -> f64 {
    let n = sorted_asc.len();
    if k == 0 || k >= n {
        return 0.0;
    }
    let thresh = sorted_asc[n - 1 - k].max(1e-30);
    let mut sum_log = 0.0_f64;
    for &x in &sorted_asc[(n - k)..] {
        sum_log += (x.max(thresh) / thresh).ln();
    }
    if sum_log > 1e-15 {
        k as f64 / sum_log
    } else {
        0.0
    }
}

pub fn run_calibration_a6_extras(
    true_p: RBergomiParams,
    n_workers: usize,
    seed: u64,
) -> CalibA6ExtrasSummary {
    let ts_main = [0.1, 0.25, 0.5, 1.0];
    let ks_calib = [-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3];
    let n_pts_main = ts_main.len() * ks_calib.len();

    let market_seed = seed ^ 0xCAFE_BABE_1234_5678;
    let market_smile = compute_smile(
        true_p,
        &ts_main,
        &ks_calib,
        500,
        200_000,
        n_workers,
        market_seed,
        0x4100,
    );
    let bank = CalibCrnBank::build(20_000, 100, n_workers, seed, 0x4200);
    let init = [0.18, 1.6, -0.6];
    let steps = [-0.04, 0.2, -0.08];

    let run_rep =
        |n_rep: usize, n_mkt: usize, scheme: HybridScheme, salt_base: u64| -> ReplicationSummary {
            let mut hs_rep = Vec::with_capacity(n_rep);
            let mut etas_rep = Vec::with_capacity(n_rep);
            let mut rhos_rep = Vec::with_capacity(n_rep);
            for rep in 0..n_rep {
                let rep_seed = seed ^ ((rep as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                let rep_mkt = compute_smile_with_scheme(
                    true_p, &ts_main, &ks_calib, n_mkt, 20_000, n_workers, rep_seed, salt_base,
                    scheme,
                );
                let (xr, _, _) =
                    nelder_mead(&[0.18, 1.6, -0.6], &[-0.04, 0.2, -0.08], 105, 1e-8, |x| {
                        eval_smile_sse_with_scheme(
                            &bank, true_p.xi0, true_p.s0, x[0], x[1], x[2], &rep_mkt, &ts_main,
                            &ks_calib, n_workers, scheme,
                        )
                    });
                hs_rep.push(xr[0]);
                etas_rep.push(xr[1]);
                rhos_rep.push(xr[2]);
            }
            ReplicationSummary {
                n_rep,
                bias_h: stats::mean(&hs_rep) - true_p.h,
                sd_h: stats::sd(&hs_rep),
                bias_eta: stats::mean(&etas_rep) - true_p.eta,
                sd_eta: stats::sd(&etas_rep),
                bias_rho: stats::mean(&rhos_rep) - true_p.rho,
                sd_rho: stats::sd(&rhos_rep),
            }
        };
    let rep_old_diff_n = run_rep(10, 250, HybridScheme::RightEndpoint, 0x4310);
    let rep_blp_diff_n = run_rep(10, 250, HybridScheme::BlpOptimal, 0x4320);

    // 10-seed median & IQR for true_rmse_grid (n in {100,250,500} x R in {20k,80k,200k})
    let mut true_rmse_grid_stats = Vec::new();
    for &n_g in &[100_usize, 250, 500] {
        for &r_g in &[20_000_usize, 80_000, 200_000] {
            let mut rmses = Vec::with_capacity(10);
            for s_idx in 0..10_u64 {
                let bank_g = CalibCrnBank::build(
                    r_g,
                    n_g,
                    n_workers,
                    seed ^ 0x7777 ^ (s_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15)),
                    (n_g * 10 + r_g / 1000) as u64,
                );
                let sse_g = eval_smile_sse(
                    &bank_g,
                    true_p.xi0,
                    true_p.s0,
                    true_p.h,
                    true_p.eta,
                    true_p.rho,
                    &market_smile,
                    &ts_main,
                    &ks_calib,
                    n_workers,
                );
                rmses.push((sse_g / n_pts_main as f64).sqrt());
            }
            rmses.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let median_rmse = 0.5 * (rmses[4] + rmses[5]);
            let q25 = rmses[2];
            let q75 = rmses[7];
            true_rmse_grid_stats.push(RmseGridStats {
                n_steps: n_g,
                n_paths: r_g,
                n_seeds: 10,
                median_rmse,
                iqr_rmse: q75 - q25,
                min_rmse: rmses[0],
                max_rmse: rmses[9],
            });
        }
    }

    // Hill tail index estimator on S_T and payoff^2 (N=200k, n=500, T=1.0)
    let hill_estimates = {
        let n_h = 500_usize;
        let r_h = 200_000_usize;
        let hw = HybridWeights::with_scheme(n_h, 1.0, true_p.h, HybridScheme::BlpOptimal);
        let dt = hw.dt;
        let half_eta2 = 0.5 * true_p.eta * true_p.eta;
        let one_m_rho2 = (1.0 - true_p.rho * true_p.rho).max(0.0).sqrt();
        let pairs = parallel::run_chunks(r_h, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x8801, r as u64));
            let mut u = vec![0.0_f64; n_h];
            let mut v = vec![0.0_f64; n_h];
            for i in 0..n_h {
                u[i] = rg.normal();
                v[i] = rg.normal();
            }
            let mut wt = vec![0.0_f64; n_h];
            let mut db = vec![0.0_f64; n_h];
            hw.generate_path(&u, &v, &mut wt, &mut db);
            let mut log_s = true_p.s0.ln();
            for i in 0..n_h {
                let v_i = if i == 0 {
                    true_p.xi0
                } else {
                    true_p.xi0
                        * (true_p.eta * wt[i - 1] - half_eta2 * exact_var(i as f64 * dt, true_p.h))
                            .exp()
                };
                let z2 = rg.normal();
                let dw = true_p.rho * db[i] + one_m_rho2 * dt.sqrt() * z2;
                log_s += -0.5 * v_i * dt + v_i.sqrt() * dw;
            }
            let st = log_s.exp();
            let pay2 = (st - true_p.s0).max(0.0).powi(2);
            (st, pay2)
        });
        let mut sts: Vec<f64> = pairs.iter().map(|x| x.0).collect();
        let mut pay2s: Vec<f64> = pairs.iter().map(|x| x.1).filter(|&p| p > 0.0).collect();
        sts.sort_by(|a, b| a.partial_cmp(b).unwrap());
        pay2s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        [500_usize, 1000, 2000, 5000]
            .iter()
            .map(|&k| HillEstimate {
                k,
                alpha_st: hill_alpha(&sts, k),
                alpha_payoff2: hill_alpha(&pay2s, k),
            })
            .collect()
    };

    // Conditional-BS market + Conditional-BS calibration (2 replication sets, R_rep=30)
    let (rep_cond_bs_fixed_crn, rep_cond_bs_diff_crn) = {
        let n_rep = 30_usize;
        let mut hs_f = Vec::with_capacity(n_rep);
        let mut etas_f = Vec::with_capacity(n_rep);
        let mut rhos_f = Vec::with_capacity(n_rep);
        let mut hs_d = Vec::with_capacity(n_rep);
        let mut etas_d = Vec::with_capacity(n_rep);
        let mut rhos_d = Vec::with_capacity(n_rep);
        for rep in 0..n_rep {
            let rep_seed = seed ^ ((rep as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let mkt_bank = CalibCrnBank::build(20_000, 100, n_workers, rep_seed, 0x4900);
            let rep_mkt_cbs =
                compute_smile_cond_bs(&mkt_bank, true_p, &ts_main, &ks_calib, n_workers);
            let (xf, _, _) =
                nelder_mead(&[0.18, 1.6, -0.6], &[-0.04, 0.2, -0.08], 105, 1e-8, |x| {
                    eval_smile_sse(
                        &bank,
                        true_p.xi0,
                        true_p.s0,
                        x[0],
                        x[1],
                        x[2],
                        &rep_mkt_cbs,
                        &ts_main,
                        &ks_calib,
                        n_workers,
                    )
                });
            hs_f.push(xf[0]);
            etas_f.push(xf[1]);
            rhos_f.push(xf[2]);

            let bank_rep = CalibCrnBank::build(20_000, 100, n_workers, rep_seed ^ 0x5555, 0x4910);
            let (xd, _, _) =
                nelder_mead(&[0.18, 1.6, -0.6], &[-0.04, 0.2, -0.08], 105, 1e-8, |x| {
                    eval_smile_sse(
                        &bank_rep,
                        true_p.xi0,
                        true_p.s0,
                        x[0],
                        x[1],
                        x[2],
                        &rep_mkt_cbs,
                        &ts_main,
                        &ks_calib,
                        n_workers,
                    )
                });
            hs_d.push(xd[0]);
            etas_d.push(xd[1]);
            rhos_d.push(xd[2]);
        }
        (
            ReplicationSummary {
                n_rep,
                bias_h: stats::mean(&hs_f) - true_p.h,
                sd_h: stats::sd(&hs_f),
                bias_eta: stats::mean(&etas_f) - true_p.eta,
                sd_eta: stats::sd(&etas_f),
                bias_rho: stats::mean(&rhos_f) - true_p.rho,
                sd_rho: stats::sd(&rhos_f),
            },
            ReplicationSummary {
                n_rep,
                bias_h: stats::mean(&hs_d) - true_p.h,
                sd_h: stats::sd(&hs_d),
                bias_eta: stats::mean(&etas_d) - true_p.eta,
                sd_eta: stats::sd(&etas_d),
                bias_rho: stats::mean(&rhos_d) - true_p.rho,
                sd_rho: stats::sd(&rhos_d),
            },
        )
    };

    // Self-consistency test & remaining moved fits
    let self_mkt = {
        let n = bank.n_steps;
        let hw = HybridWeights::with_scheme(n, 1.0, true_p.h, HybridScheme::BlpOptimal);
        let dt = hw.dt;
        let half_eta2 = 0.5 * true_p.eta * true_p.eta;
        let one_m_rho2 = (1.0 - true_p.rho * true_p.rho).max(0.0);
        let step_targets: Vec<usize> = ts_main.iter().map(|&t| (t / dt).round() as usize).collect();
        let n_t = ts_main.len();
        let n_k = ks_calib.len();
        let sums = parallel::run_chunk_reduce(
            bank.n_paths,
            n_workers,
            |start, end| {
                let mut local = vec![0.0_f64; n_t * n_k];
                let mut wt = vec![0.0_f64; n];
                let mut db = vec![0.0_f64; n];
                for r in start..end {
                    let off = r * n;
                    hw.generate_path(
                        &bank.u_flat[off..off + n],
                        &bank.v_flat[off..off + n],
                        &mut wt,
                        &mut db,
                    );
                    let mut int_v = 0.0_f64;
                    let mut stoch_int = 0.0_f64;
                    let mut tp = 0_usize;
                    for i in 0..n {
                        let v_i = if i == 0 {
                            true_p.xi0
                        } else {
                            true_p.xi0
                                * (true_p.eta * wt[i - 1]
                                    - half_eta2 * exact_var(i as f64 * dt, true_p.h))
                                .exp()
                        };
                        int_v += v_i * dt;
                        stoch_int += v_i.sqrt() * db[i];
                        if tp < n_t && (i + 1) == step_targets[tp] {
                            let s_eff = true_p.s0
                                * (true_p.rho * stoch_int - 0.5 * true_p.rho * true_p.rho * int_v)
                                    .exp();
                            let vc = one_m_rho2 * int_v;
                            for (ki, &lk) in ks_calib.iter().enumerate() {
                                local[tp * n_k + ki] +=
                                    bs::bs_call_total_var(s_eff, true_p.s0 * lk.exp(), vc);
                            }
                            tp += 1;
                        }
                    }
                }
                local
            },
            |parts| {
                let mut tot = vec![0.0_f64; n_t * n_k];
                for part in parts {
                    for (a, b) in tot.iter_mut().zip(part.iter()) {
                        *a += *b;
                    }
                }
                tot
            },
        );
        let mut pts = Vec::with_capacity(n_t * n_k);
        for (ti, &t_h) in ts_main.iter().enumerate() {
            for (ki, &lk) in ks_calib.iter().enumerate() {
                let cp = sums[ti * n_k + ki] / bank.n_paths as f64;
                let iv = bs::implied_vol(true_p.s0, true_p.s0 * lk.exp(), t_h, cp);
                pts.push(SmilePoint {
                    t_horizon: t_h,
                    log_k: lk,
                    call_price: cp,
                    call_se: 0.0,
                    iv,
                    iv_se: 0.0,
                });
            }
        }
        pts
    };
    let (x_sc, sse_sc, ev_sc) = nelder_mead(&init, &steps, 140, 1e-12, |x| {
        eval_smile_sse(
            &bank, true_p.xi0, true_p.s0, x[0], x[1], x[2], &self_mkt, &ts_main, &ks_calib,
            n_workers,
        )
    });
    let self_consistency_fit = CalibResult {
        h_hat: x_sc[0],
        eta_hat: x_sc[1],
        rho_hat: x_sc[2],
        rmse_opt: (sse_sc / n_pts_main as f64).sqrt(),
        rmse_at_true: 0.0,
        evals: ev_sc,
    };

    // (ii) Per-replication DIFFERENT CRN bank (same-n=100, R_rep=15)
    let rep_diff_crn = {
        let n_rep = 15;
        let mut hs_rep = Vec::with_capacity(n_rep);
        let mut etas_rep = Vec::with_capacity(n_rep);
        let mut rhos_rep = Vec::with_capacity(n_rep);
        for rep in 0..n_rep {
            let rep_seed = seed ^ ((rep as u64 + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15));
            let rep_mkt = compute_smile_with_scheme(
                true_p,
                &ts_main,
                &ks_calib,
                100,
                20_000,
                n_workers,
                rep_seed,
                0x4350,
                HybridScheme::BlpOptimal,
            );
            let bank_rep = CalibCrnBank::build(20_000, 100, n_workers, rep_seed ^ 0x5555, 0x4360);
            let (xr, _, _) =
                nelder_mead(&[0.18, 1.6, -0.6], &[-0.04, 0.2, -0.08], 105, 1e-8, |x| {
                    eval_smile_sse(
                        &bank_rep, true_p.xi0, true_p.s0, x[0], x[1], x[2], &rep_mkt, &ts_main,
                        &ks_calib, n_workers,
                    )
                });
            hs_rep.push(xr[0]);
            etas_rep.push(xr[1]);
            rhos_rep.push(xr[2]);
        }
        ReplicationSummary {
            n_rep,
            bias_h: stats::mean(&hs_rep) - true_p.h,
            sd_h: stats::sd(&hs_rep),
            bias_eta: stats::mean(&etas_rep) - true_p.eta,
            sd_eta: stats::sd(&etas_rep),
            bias_rho: stats::mean(&rhos_rep) - true_p.rho,
            sd_rho: stats::sd(&rhos_rep),
        }
    };

    // (iii) R_calib = 80,000 fit
    let bank_80k = CalibCrnBank::build(80_000, 100, n_workers, seed, 0x4280);
    let (x_80k, sse_80k, ev_80k) = nelder_mead(&init, &steps, 105, 1e-10, |x| {
        eval_smile_sse(
            &bank_80k,
            true_p.xi0,
            true_p.s0,
            x[0],
            x[1],
            x[2],
            &market_smile,
            &ts_main,
            &ks_calib,
            n_workers,
        )
    });
    let fit_r80k = CalibResult {
        h_hat: x_80k[0],
        eta_hat: x_80k[1],
        rho_hat: x_80k[2],
        rmse_opt: (sse_80k / n_pts_main as f64).sqrt(),
        rmse_at_true: 0.0,
        evals: ev_80k,
    };

    // (iv) Moved to 10-seed true_rmse_grid_stats above.

    // (v) Vega-weighted conditional-BS price-space calibration (avoids IV inversion wing noise)
    let eval_vega_price_sse = |h: f64, eta: f64, rho: f64| -> f64 {
        if !(0.01..=0.495).contains(&h)
            || !(0.05..=5.0).contains(&eta)
            || !(-0.99..=0.99).contains(&rho)
        {
            return 1e6;
        }
        let n = bank.n_steps;
        let hw = HybridWeights::with_scheme(n, 1.0, h, HybridScheme::BlpOptimal);
        let dt = hw.dt;
        let half_eta2 = 0.5 * eta * eta;
        let one_m_rho2 = (1.0 - rho * rho).max(0.0);
        let step_targets: Vec<usize> = ts_main.iter().map(|&t| (t / dt).round() as usize).collect();
        let n_t = ts_main.len();
        let n_k = ks_calib.len();
        let sums = parallel::run_chunk_reduce(
            bank.n_paths,
            n_workers,
            |start, end| {
                let mut local = vec![0.0_f64; n_t * n_k];
                let mut wt = vec![0.0_f64; n];
                let mut db = vec![0.0_f64; n];
                for r in start..end {
                    let off = r * n;
                    hw.generate_path(
                        &bank.u_flat[off..off + n],
                        &bank.v_flat[off..off + n],
                        &mut wt,
                        &mut db,
                    );
                    let mut int_v = 0.0_f64;
                    let mut stoch_int = 0.0_f64;
                    let mut tp = 0_usize;
                    for i in 0..n {
                        let v_i = if i == 0 {
                            true_p.xi0
                        } else {
                            true_p.xi0
                                * (eta * wt[i - 1] - half_eta2 * exact_var(i as f64 * dt, h)).exp()
                        };
                        int_v += v_i * dt;
                        stoch_int += v_i.sqrt() * db[i];
                        if tp < n_t && (i + 1) == step_targets[tp] {
                            let s_eff =
                                true_p.s0 * (rho * stoch_int - 0.5 * rho * rho * int_v).exp();
                            let vc = one_m_rho2 * int_v;
                            for (ki, &lk) in ks_calib.iter().enumerate() {
                                local[tp * n_k + ki] +=
                                    bs::bs_call_total_var(s_eff, true_p.s0 * lk.exp(), vc);
                            }
                            tp += 1;
                        }
                    }
                }
                local
            },
            |parts| {
                let mut tot = vec![0.0_f64; n_t * n_k];
                for p in parts {
                    for (a, b) in tot.iter_mut().zip(p.iter()) {
                        *a += *b;
                    }
                }
                tot
            },
        );
        let nf = bank.n_paths as f64;
        let mut sse = 0.0_f64;
        for (ti, &t_h) in ts_main.iter().enumerate() {
            for (ki, &lk) in ks_calib.iter().enumerate() {
                let c_mod = sums[ti * n_k + ki] / nf;
                if let Some(pt) = market_smile
                    .iter()
                    .find(|pt| (pt.t_horizon - t_h).abs() < 1e-12 && (pt.log_k - lk).abs() < 1e-9)
                {
                    let strike = true_p.s0 * lk.exp();
                    let vega = bs::bs_vega(true_p.s0, strike, t_h, pt.iv).max(1e-4);
                    let d = (c_mod - pt.call_price) / vega;
                    sse += d * d;
                }
            }
        }
        sse
    };
    let (x_vp, sse_vp, ev_vp) = nelder_mead(&init, &steps, 105, 1e-10, |x| {
        eval_vega_price_sse(x[0], x[1], x[2])
    });
    let fit_vega_price = CalibResult {
        h_hat: x_vp[0],
        eta_hat: x_vp[1],
        rho_hat: x_vp[2],
        rmse_opt: (sse_vp / n_pts_main as f64).sqrt(),
        rmse_at_true: 0.0,
        evals: ev_vp,
    };

    // 5. Misspecification: T >= 0.5 only vs T <= 0.1 only
    let ts_long = [0.5, 1.0];
    let (xl, ssel, evl) = nelder_mead(&init, &steps, 90, 1e-9, |x| {
        eval_smile_sse(
            &bank,
            true_p.xi0,
            true_p.s0,
            x[0],
            x[1],
            x[2],
            &market_smile,
            &ts_long,
            &ks_calib,
            n_workers,
        )
    });
    let fit_long_only = CalibResult {
        h_hat: xl[0],
        eta_hat: xl[1],
        rho_hat: xl[2],
        rmse_opt: (ssel / (ts_long.len() * ks_calib.len()) as f64).sqrt(),
        rmse_at_true: 0.0,
        evals: evl,
    };

    let ts_short = [0.05, 0.1];
    let mkt_short = compute_smile(
        true_p,
        &ts_short,
        &ks_calib,
        500,
        200_000,
        n_workers,
        market_seed,
        0x4400,
    );
    let bank_short = CalibCrnBank::build(20_000, 100, n_workers, seed, 0x4500);
    // Use n=100 on T=1.0 -> step 5 is T=0.05, step 10 is T=0.10!
    let (xs, sses, evs) = nelder_mead(&init, &steps, 90, 1e-9, |x| {
        eval_smile_sse(
            &bank_short,
            true_p.xi0,
            true_p.s0,
            x[0],
            x[1],
            x[2],
            &mkt_short,
            &ts_short,
            &ks_calib,
            n_workers,
        )
    });
    let fit_short_only = CalibResult {
        h_hat: xs[0],
        eta_hat: xs[1],
        rho_hat: xs[2],
        rmse_opt: (sses / (ts_short.len() * ks_calib.len()) as f64).sqrt(),
        rmse_at_true: 0.0,
        evals: evs,
    };

    CalibA6ExtrasSummary {
        self_consistency_fit,
        fit_r80k,
        fit_vega_price,
        true_rmse_grid_stats,
        hill_estimates,
        rep_diff_crn,
        rep_old_diff_n,
        rep_blp_diff_n,
        rep_cond_bs_fixed_crn,
        rep_cond_bs_diff_crn,
        fit_long_only,
        fit_short_only,
    }
}
