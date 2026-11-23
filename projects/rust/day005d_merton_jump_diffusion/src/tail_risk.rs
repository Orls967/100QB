use crate::bs::{bs_put, implied_vol_call};
use crate::merton::{merton_call_series, merton_put_series, MertonParams};
use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::stats::{mean, quantile};

#[derive(Clone, Debug)]
pub struct OtmPutComparisonRow {
    pub k: f64,
    pub t: f64,
    pub merton_put: f64,
    pub bs_atm_put: f64,
    pub ratio_merton_to_bs: f64,
}

#[derive(Clone, Debug)]
pub struct ShortPutTailRisk {
    pub t: f64,
    pub k: f64,
    pub r_paths: usize,
    pub bs_atm_vol: f64,
    pub var99_merton: f64,
    pub es99_merton: f64,
    pub var99_bs: f64,
    pub es99_bs: f64,
}

pub fn compare_otm_puts(p: &MertonParams, ks: &[f64], t: f64) -> Vec<OtmPutComparisonRow> {
    let atm_call = merton_call_series(p, p.s0, t);
    let atm_iv = implied_vol_call(atm_call, p.s0, p.s0, p.r, p.q, t);
    ks.iter()
        .map(|&k| {
            let merton_put = merton_put_series(p, k, t);
            let bs_atm_put = bs_put(p.s0, k, p.r, p.q, atm_iv, t);
            OtmPutComparisonRow {
                k,
                t,
                merton_put,
                bs_atm_put,
                ratio_merton_to_bs: merton_put / bs_atm_put.max(1e-12),
            }
        })
        .collect()
}

/// Computes 99% VaR and Expected Shortfall of an unhedged short OTM put (`K = 90, T = 1/12`)
/// under Merton jump-diffusion vs Black-Scholes matched to ATM implied volatility.
pub fn simulate_short_put_tail_risk(
    p: &MertonParams,
    k: f64,
    t: f64,
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> ShortPutTailRisk {
    let atm_call = merton_call_series(p, p.s0, t);
    let bs_atm_vol = implied_vol_call(atm_call, p.s0, p.s0, p.r, p.q, t);

    let p0_m = merton_put_series(p, k, t);
    let p0_b = bs_put(p.s0, k, p.r, p.q, bs_atm_vol, t);

    let kappa = p.kappa();
    let lam_t = p.lambda * t;
    let drift_m = (p.r - p.q - 0.5 * p.sigma * p.sigma - p.lambda * kappa) * t;
    let vol_m = p.sigma * t.sqrt();
    let drift_b = (p.r - p.q - 0.5 * bs_atm_vol * bs_atm_vol) * t;
    let vol_b = bs_atm_vol * t.sqrt();

    let losses: Vec<(f64, f64)> = run_chunks(r_paths, n_workers, |idx| {
        let mut rng = Rng::new(derive_seed(main_seed, 400, idx as u64));
        let z = rng.normal();
        let n_jumps = rng.poisson(lam_t);
        let jump_sum = if n_jumps == 0 {
            0.0
        } else {
            let nf = n_jumps as f64;
            nf * p.mu_j + (nf.sqrt() * p.delta_j) * rng.normal()
        };
        let s_m = p.s0 * (drift_m + vol_m * z + jump_sum).exp();
        let s_b = p.s0 * (drift_b + vol_b * z).exp();
        let loss_m = (k - s_m).max(0.0) - p0_m * (p.r * t).exp();
        let loss_b = (k - s_b).max(0.0) - p0_b * (p.r * t).exp();
        (loss_m, loss_b)
    });

    let mut lm: Vec<f64> = losses.iter().map(|x| x.0).collect();
    let mut lb: Vec<f64> = losses.iter().map(|x| x.1).collect();
    let var99_merton = quantile(&lm, 0.99);
    let var99_bs = quantile(&lb, 0.99);

    lm.sort_by(|a, b| a.total_cmp(b));
    lb.sort_by(|a, b| a.total_cmp(b));
    let tail_start = ((0.99 * r_paths as f64).floor() as usize).min(r_paths - 1);
    let es99_merton = mean(&lm[tail_start..]);
    let es99_bs = mean(&lb[tail_start..]);

    ShortPutTailRisk {
        t,
        k,
        r_paths,
        bs_atm_vol,
        var99_merton,
        es99_merton,
        var99_bs,
        es99_bs,
    }
}
