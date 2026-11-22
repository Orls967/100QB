//! Bagian 3: Implied volatility smile, delta-method SE, ATM skew scaling psi(T) ~ T^{H - 1/2},
//! and Markov (H=0.5) comparison at T = 0.1.

use crate::bs;
use crate::parallel;
use crate::rbergomi::RBergomiParams;
use crate::rng::{derive_seed, Rng};
use crate::stats;
use crate::volterra::{exact_var, HybridScheme, HybridWeights};

#[derive(Clone, Debug)]
pub struct SmilePoint {
    pub t_horizon: f64,
    pub log_k: f64,
    pub call_price: f64,
    pub call_se: f64,
    pub iv: f64,
    pub iv_se: f64,
}

#[derive(Clone, Debug)]
pub struct AtmSkewRow {
    pub h: f64,
    pub t_horizons: Vec<f64>,
    pub atm_skews: Vec<f64>,
    pub slope: f64,
    pub se_slope: f64,
    pub theoretical_slope: f64,
    pub r2: f64,
}

#[derive(Clone, Debug)]
pub struct MarkovCompareSummary {
    pub t_horizon: f64,
    pub eta_markov: f64,
    pub rbergomi_atm_iv: f64,
    pub markov_atm_iv: f64,
    pub rbergomi_skew: f64,
    pub markov_skew: f64,
    pub rbergomi_smile: Vec<SmilePoint>,
    pub markov_smile: Vec<SmilePoint>,
}

/// Computes implied volatility smile on `(ts, log_ks)` using conditional BS estimator (`R` hybrid paths, `n_steps`).
#[allow(clippy::too_many_arguments)]
pub fn compute_smile(
    p: RBergomiParams,
    ts: &[f64],
    log_ks: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    salt: u64,
) -> Vec<SmilePoint> {
    compute_smile_with_scheme(
        p,
        ts,
        log_ks,
        n_steps,
        n_paths,
        n_workers,
        seed,
        salt,
        HybridScheme::BlpOptimal,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn compute_smile_with_scheme(
    p: RBergomiParams,
    ts: &[f64],
    log_ks: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
    salt: u64,
    scheme: HybridScheme,
) -> Vec<SmilePoint> {
    let half_eta2 = 0.5 * p.eta * p.eta;
    let one_m_rho2 = (1.0 - p.rho * p.rho).max(0.0);
    let mut out = Vec::with_capacity(ts.len() * log_ks.len());

    for &t_h in ts {
        let hw = HybridWeights::with_scheme(n_steps, t_h, p.h, scheme);
        let dt = hw.dt;
        let cond_paths = parallel::run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, salt ^ t_h.to_bits(), r as u64));
            let mut u = vec![0.0_f64; n_steps];
            let mut v = vec![0.0_f64; n_steps];
            let mut wt = vec![0.0_f64; n_steps];
            let mut db = vec![0.0_f64; n_steps];
            for i in 0..n_steps {
                u[i] = rg.normal();
                v[i] = rg.normal();
            }
            hw.generate_path(&u, &v, &mut wt, &mut db);

            let mut int_v = 0.0_f64;
            let mut stoch_int = 0.0_f64;
            for i in 0..n_steps {
                let v_i = if i == 0 {
                    p.xi0
                } else {
                    let t_i = i as f64 * dt;
                    p.xi0 * (p.eta * wt[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
                };
                int_v += v_i * dt;
                stoch_int += v_i.sqrt() * db[i];
            }
            let s_eff = p.s0 * (p.rho * stoch_int - 0.5 * p.rho * p.rho * int_v).exp();
            (s_eff, one_m_rho2 * int_v)
        });

        for &lk in log_ks {
            let strike = p.s0 * lk.exp();
            let prices: Vec<f64> = cond_paths
                .iter()
                .map(|&(seff, vc)| bs::bs_call_total_var(seff, strike, vc))
                .collect();
            let m_p = stats::mean(&prices);
            let se_p = stats::se(&prices);
            let iv = bs::implied_vol(p.s0, strike, t_h, m_p);
            let vega = bs::bs_vega(p.s0, strike, t_h, iv);
            let iv_se = if vega > 1e-12 && iv.is_finite() {
                se_p / vega
            } else {
                f64::NAN
            };
            out.push(SmilePoint {
                t_horizon: t_h,
                log_k: lk,
                call_price: m_p,
                call_se: se_p,
                iv,
                iv_se,
            });
        }
    }
    out
}

/// Computes ATM skew `|d(sigma_imp)/dk|_{k=0}` via symmetric finite difference at `k = +/- dk`.
pub fn atm_skew_from_smile(smile: &[SmilePoint], t_h: f64, dk: f64) -> f64 {
    let mut iv_minus = f64::NAN;
    let mut iv_plus = f64::NAN;
    for pt in smile {
        if (pt.t_horizon - t_h).abs() < 1e-12 {
            if (pt.log_k + dk).abs() < 1e-9 {
                iv_minus = pt.iv;
            }
            if (pt.log_k - dk).abs() < 1e-9 {
                iv_plus = pt.iv;
            }
        }
    }
    ((iv_plus - iv_minus) / (2.0 * dk)).abs()
}

pub fn run_atm_skew_scaling_study(
    base_p: RBergomiParams,
    hs: &[f64],
    ts: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<AtmSkewRow> {
    let dk = 0.05;
    let ks = [-dk, 0.0, dk];
    let log_ts: Vec<f64> = ts.iter().map(|t| t.ln()).collect();
    let mut rows = Vec::with_capacity(hs.len());

    for &h in hs {
        let mut ph = base_p;
        ph.h = h;
        let sm = compute_smile(ph, ts, &ks, n_steps, n_paths, n_workers, seed, 0x3100);
        let skews: Vec<f64> = ts
            .iter()
            .map(|&t| atm_skew_from_smile(&sm, t, dk))
            .collect();
        let log_skews: Vec<f64> = skews.iter().map(|s| s.max(1e-12).ln()).collect();
        let (_, slope, se_slope, r2) = stats::ols_line(&log_ts, &log_skews);
        rows.push(AtmSkewRow {
            h,
            t_horizons: ts.to_vec(),
            atm_skews: skews,
            slope,
            se_slope,
            theoretical_slope: h - 0.5,
            r2,
        });
    }
    rows
}

/// Compares rough Bergomi (`H=0.1`) vs 1-factor Markov Bergomi (`H=0.5`) at `T=0.1`
/// with `eta_markov` tuned so ATM implied vol matches (cumulative expected variance is already `xi0*T`).
pub fn run_markov_comparison(
    base_p: RBergomiParams,
    rb_smile_t01: &[SmilePoint],
    log_ks: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> MarkovCompareSummary {
    let t_h = 0.1;
    let target_atm_iv = rb_smile_t01
        .iter()
        .find(|pt| (pt.t_horizon - t_h).abs() < 1e-12 && pt.log_k.abs() < 1e-9)
        .map(|pt| pt.iv)
        .unwrap_or(0.19);

    let mut p_mk = base_p;
    p_mk.h = 0.5;
    // Bisection for eta_markov in [0.1, 6.0] to match ATM IV at T=0.1
    let mut lo = 0.1_f64;
    let mut hi = 6.0_f64;
    for _ in 0..22 {
        let mid = 0.5 * (lo + hi);
        p_mk.eta = mid;
        let sm = compute_smile(
            p_mk,
            &[t_h],
            &[0.0],
            n_steps,
            n_paths / 4,
            n_workers,
            seed,
            0x3200,
        );
        let iv = sm[0].iv;
        // Larger eta_markov -> stronger vol-of-vol convexity -> lower ATM IV
        if iv > target_atm_iv {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    p_mk.eta = 0.5 * (lo + hi);
    let mk_smile = compute_smile(
        p_mk,
        &[t_h],
        log_ks,
        n_steps,
        n_paths,
        n_workers,
        seed,
        0x3200,
    );
    let mk_atm_iv = mk_smile
        .iter()
        .find(|pt| pt.log_k.abs() < 1e-9)
        .map(|pt| pt.iv)
        .unwrap_or(f64::NAN);
    let rb_skew = atm_skew_from_smile(rb_smile_t01, t_h, 0.05);
    let mk_skew = atm_skew_from_smile(&mk_smile, t_h, 0.05);

    MarkovCompareSummary {
        t_horizon: t_h,
        eta_markov: p_mk.eta,
        rbergomi_atm_iv: target_atm_iv,
        markov_atm_iv: mk_atm_iv,
        rbergomi_skew: rb_skew,
        markov_skew: mk_skew,
        rbergomi_smile: rb_smile_t01.to_vec(),
        markov_smile: mk_smile,
    }
}
