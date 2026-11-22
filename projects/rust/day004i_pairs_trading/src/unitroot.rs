//! Unit-root and cointegration tests:
//! (1) Dickey-Fuller with intercept (T=1000, R=10^6) & Brownian functional limit,
//! (2) Engle-Granger cointegration test (residuals of y = alpha + beta*x) critical quantiles & size,
//! (3) Spurious regression diagnostics (|t|>1.96 fraction, R^2, Durbin-Watson) & EG power vs (theta, T).

use crate::ou::{simulate_ou_path, OuParams};
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

/// Computes Dickey-Fuller t-statistic (with intercept, no lag) on `y[0..=T]`:
/// `Delta y_t = alpha + gamma * y_{t-1} + eps_t`, returning `t_gamma = gamma_hat / SE(gamma_hat)`.
pub fn dickey_fuller_t_stat(y: &[f64]) -> f64 {
    let n = y.len().saturating_sub(1);
    if n < 3 {
        return 0.0;
    }
    let mut sum_x = 0.0_f64;
    let mut sum_dy = 0.0_f64;
    for t in 0..n {
        sum_x += y[t];
        sum_dy += y[t + 1] - y[t];
    }
    let nf = n as f64;
    let mx = sum_x / nf;
    let mdy = sum_dy / nf;

    let mut sxx = 0.0_f64;
    let mut sxy = 0.0_f64;
    let mut syy = 0.0_f64;
    for t in 0..n {
        let dx = y[t] - mx;
        let ddy = (y[t + 1] - y[t]) - mdy;
        sxx += dx * dx;
        sxy += dx * ddy;
        syy += ddy * ddy;
    }
    if sxx <= 1e-18 {
        return 0.0;
    }
    let gamma = sxy / sxx;
    let rss = (syy - gamma * sxy).max(0.0);
    let s2 = rss / (nf - 2.0);
    let se_gamma = (s2 / sxx).sqrt();
    if se_gamma > 0.0 {
        gamma / se_gamma
    } else {
        0.0
    }
}

/// Computes `dickey_fuller_t_stat` directly on a single random walk realization of length `n_steps`
/// without allocating a `Vec<f64>`.
#[inline]
fn df_random_walk_stat_inplace(n_steps: usize, rg: &mut Rng) -> f64 {
    let mut y = 0.0_f64;
    let mut sum_x = 0.0_f64;
    let mut sum_x2 = 0.0_f64;
    let mut sum_dy = 0.0_f64;
    let mut sum_dy2 = 0.0_f64;
    let mut sum_x_dy = 0.0_f64;
    for _ in 0..n_steps {
        let dy = rg.normal();
        sum_x += y;
        sum_x2 += y * y;
        sum_dy += dy;
        sum_dy2 += dy * dy;
        sum_x_dy += y * dy;
        y += dy;
    }
    let nf = n_steps as f64;
    let sxx = sum_x2 - sum_x * sum_x / nf;
    let syy = sum_dy2 - sum_dy * sum_dy / nf;
    let sxy = sum_x_dy - sum_x * sum_dy / nf;
    if sxx <= 1e-18 {
        return 0.0;
    }
    let gamma = sxy / sxx;
    let rss = (syy - gamma * sxy).max(0.0);
    let s2 = rss / (nf - 2.0);
    gamma / (s2 / sxx).sqrt()
}

/// Computes the demeaned Brownian functional limit of Dickey-Fuller with intercept:
/// `T_limit = (0.5*(W_tilde(1)^2 - W_tilde(0)^2 - 1) ...)` = `(int_0^1 W_tilde dW) / sqrt(int_0^1 W_tilde^2 dt)`
/// where `W_tilde(t) = W(t) - int_0^1 W(s) ds`.
#[inline]
fn df_brownian_functional_inplace(n_steps: usize, rg: &mut Rng) -> f64 {
    let dt = 1.0 / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let mut w = 0.0_f64;
    let mut sum_w = 0.0_f64;
    let mut sum_w2 = 0.0_f64;
    let mut sum_w_dw = 0.0_f64;
    for _ in 0..n_steps {
        let dw = sqrt_dt * rg.normal();
        sum_w += w;
        sum_w2 += w * w;
        sum_w_dw += w * dw;
        w += dw;
    }
    let mean_w = sum_w * dt;
    let int_w_tilde2 = (sum_w2 * dt - mean_w * mean_w).max(1e-18);
    // int W_tilde dW = int W dW - mean_w * W(1)
    let num = sum_w_dw - mean_w * w;
    num / int_w_tilde2.sqrt()
}

#[derive(Clone, Debug)]
pub struct DfDistSummary {
    pub n_paths: usize,
    pub t_len: usize,
    pub q01: f64,
    pub se_q01: f64,
    pub q05: f64,
    pub se_q05: f64,
    pub q10: f64,
    pub se_q10: f64,
    pub brownian_q05: f64,
}

pub fn simulate_dickey_fuller_dist(
    t_len: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> DfDistSummary {
    let mut stats_vec = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3100, r as u64));
        df_random_walk_stat_inplace(t_len, &mut rg)
    });
    stats_vec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (q01, se01) = stats::quantile_and_se(&stats_vec, 0.01);
    let (q05, se05) = stats::quantile_and_se(&stats_vec, 0.05);
    let (q10, se10) = stats::quantile_and_se(&stats_vec, 0.10);

    // Brownian functional cross-check on fine grid (n_steps=2000, R=50_000)
    let n_bm = (n_paths / 20).clamp(2000, 50_000);
    let mut bm_vec = parallel::run_chunks(n_bm, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3150, r as u64));
        df_brownian_functional_inplace(2000, &mut rg)
    });
    bm_vec.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (bm_q05, _) = stats::quantile_and_se(&bm_vec, 0.05);

    DfDistSummary {
        n_paths,
        t_len,
        q01,
        se_q01: se01,
        q05,
        se_q05: se05,
        q10,
        se_q10: se10,
        brownian_q05: bm_q05,
    }
}

#[derive(Clone, Debug)]
pub struct EgDistSummary {
    pub n_paths: usize,
    pub t_len: usize,
    pub q01: f64,
    pub se_q01: f64,
    pub q05: f64,
    pub se_q05: f64,
    pub q10: f64,
    pub se_q10: f64,
    pub size_at_own_q05: f64,
    pub size_if_using_df_q05: f64,
    pub size_at_lit_eg_q05: f64,
}

/// Simulates Engle-Granger residual ADF t-statistic distribution under H0 (independent random walks),
/// plus empirical size on independent pairs.
pub fn simulate_engle_granger_dist(
    t_len: usize,
    n_paths: usize,
    df_q05: f64,
    n_workers: usize,
    seed: u64,
) -> EgDistSummary {
    let mut eg_stats = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3200, r as u64));
        let mut x = Vec::with_capacity(t_len);
        let mut y = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += rg.normal();
            cy += rg.normal();
            x.push(cx);
            y.push(cy);
        }
        let ols = stats::ols_simple(&x, &y);
        dickey_fuller_t_stat(&ols.residuals)
    });
    eg_stats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let (q01, se01) = stats::quantile_and_se(&eg_stats, 0.01);
    let (q05, se05) = stats::quantile_and_se(&eg_stats, 0.05);
    let (q10, se10) = stats::quantile_and_se(&eg_stats, 0.10);

    // Independent out-of-sample size check on fresh seed R=n_paths independent pairs
    let oos = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed ^ 0x7777_0002, 0x3250, r as u64));
        let mut x = Vec::with_capacity(t_len);
        let mut y = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += rg.normal();
            cy += rg.normal();
            x.push(cx);
            y.push(cy);
        }
        let ols = stats::ols_simple(&x, &y);
        let tstat = dickey_fuller_t_stat(&ols.residuals);
        (tstat < q05, tstat < df_q05, tstat < -3.34)
    });
    let rej_own = oos.iter().filter(|x| x.0).count() as f64 / n_paths as f64;
    let rej_df = oos.iter().filter(|x| x.1).count() as f64 / n_paths as f64;
    let rej_lit = oos.iter().filter(|x| x.2).count() as f64 / n_paths as f64;

    EgDistSummary {
        n_paths,
        t_len,
        q01,
        se_q01: se01,
        q05,
        se_q05: se05,
        q10,
        se_q10: se10,
        size_at_own_q05: rej_own,
        size_if_using_df_q05: rej_df,
        size_at_lit_eg_q05: rej_lit,
    }
}

#[derive(Clone, Debug)]
pub struct SpuriousSummary {
    pub n_paths: usize,
    pub t_len: usize,
    pub frac_abs_t_gt_196: f64,
    pub mean_r2: f64,
    pub mean_dw: f64,
}

pub fn simulate_spurious_regression(
    t_len: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> SpuriousSummary {
    let dt = 1.0_f64 / 252.0;
    let sqrt_dt = dt.sqrt();
    let rows = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x3300, r as u64));
        let mut x = Vec::with_capacity(t_len);
        let mut y = Vec::with_capacity(t_len);
        let mut cx = 0.0_f64;
        let mut cy = 0.0_f64;
        for _ in 0..t_len {
            cx += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            cy += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
            x.push(cx);
            y.push(cy);
        }
        let ols = stats::ols_simple(&x, &y);
        (ols.t_beta.abs() > 1.96, ols.r2, ols.durbin_watson)
    });
    let frac = rows.iter().filter(|r| r.0).count() as f64 / n_paths as f64;
    let r2s: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let dws: Vec<f64> = rows.iter().map(|r| r.2).collect();
    SpuriousSummary {
        n_paths,
        t_len,
        frac_abs_t_gt_196: frac,
        mean_r2: stats::mean(&r2s),
        mean_dw: stats::mean(&dws),
    }
}

#[derive(Clone, Debug)]
pub struct EgPowerRow {
    pub theta: f64,
    pub t_years: f64,
    pub power: f64,
}

pub fn run_eg_power_grid(
    thetas: &[f64],
    t_years_list: &[f64],
    eg_q05: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<EgPowerRow> {
    let dt = 1.0_f64 / 252.0;
    let sqrt_dt = dt.sqrt();
    let mut out = Vec::new();
    for &theta in thetas {
        for &t_yr in t_years_list {
            let n_steps = (t_yr * 252.0).round() as usize;
            let ou_p = OuParams {
                theta,
                mu: 0.0,
                sigma: 0.20,
                dt,
            };
            let hits = parallel::run_chunks(n_paths, n_workers, |r| {
                let mut rg = Rng::new(derive_seed(
                    seed,
                    0x3400 ^ theta.to_bits() ^ t_yr.to_bits(),
                    r as u64,
                ));
                let spread = simulate_ou_path(ou_p, n_steps - 1, true, &mut rg);
                let mut x = Vec::with_capacity(n_steps);
                let mut y = Vec::with_capacity(n_steps);
                let mut cx = 0.0_f64;
                for &sp in &spread {
                    cx += 0.05 * dt + 0.20 * sqrt_dt * rg.normal();
                    x.push(cx);
                    y.push(0.1 + 1.2 * cx + sp);
                }
                let ols = stats::ols_simple(&x, &y);
                dickey_fuller_t_stat(&ols.residuals) < eg_q05
            });
            let pwr = hits.iter().filter(|&&b| b).count() as f64 / n_paths as f64;
            out.push(EgPowerRow {
                theta,
                t_years: t_yr,
                power: pwr,
            });
        }
    }
    out
}
