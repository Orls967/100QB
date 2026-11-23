use crate::parallel::run_chunks;
use crate::rng::{derive_seed, Rng};
use crate::special::norm_cdf;
use crate::stats::{mean, quantile, se};

#[derive(Clone, Debug)]
pub struct BinaryKellyPoint {
    pub f: f64,
    pub g_exact: f64,
    pub e_ln_w_mc: f64,
    pub e_ln_w_se: f64,
    pub e_ln_w_exact: f64,
    pub median_w_mc: f64,
    pub median_w_exact: f64,
    pub arith_mean_mc: f64,
    pub arith_mean_se: f64,
    pub arith_mean_exact: f64,
    pub prob_loss_mc: f64,
    pub prob_loss_se: f64,
    pub prob_loss_cont_approx: f64,
    pub prob_touch_half_mc: f64,
    pub prob_touch_half_se: f64,
    pub prob_touch_half_cont_approx: f64,
}

#[derive(Clone, Debug)]
pub struct BinaryKellySummary {
    pub p: f64,
    pub b: f64,
    pub f_star_exact: f64,
    pub f_star_numeric: f64,
    pub g_star_exact: f64,
    pub f_breakeven_numeric: f64,
    pub points: Vec<BinaryKellyPoint>,
    pub sample_paths_n: Vec<usize>,
    pub sample_path_half_kelly: Vec<f64>,
    pub sample_path_full_kelly: Vec<f64>,
    pub sample_path_over_kelly: Vec<f64>,
}

#[inline]
pub fn g_binary(p: f64, b: f64, f: f64) -> f64 {
    p * (1.0 + b * f).ln() + (1.0 - p) * (1.0 - f).ln()
}

#[inline]
pub fn f_star_binary(p: f64, b: f64) -> f64 {
    p - (1.0 - p) / b
}

/// Numerical maximization of `g_binary(p, b, f)` on `(0, 1)` via golden-section + Newton polish.
pub fn maximize_g_binary(p: f64, b: f64) -> f64 {
    let phi = 0.5 * (5.0_f64.sqrt() - 1.0);
    let mut a = 1e-9_f64;
    let mut c = 1.0 - 1e-9_f64;
    let mut x1 = c - phi * (c - a);
    let mut x2 = a + phi * (c - a);
    let mut f1 = g_binary(p, b, x1);
    let mut f2 = g_binary(p, b, x2);
    for _ in 0..50 {
        if f1 < f2 {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (c - a);
            f2 = g_binary(p, b, x2);
        } else {
            c = x2;
            x2 = x1;
            f2 = f1;
            x1 = c - phi * (c - a);
            f1 = g_binary(p, b, x1);
        }
    }
    let mut x = 0.5 * (a + c);
    for _ in 0..6 {
        let gp = p * b / (1.0 + b * x) - (1.0 - p) / (1.0 - x);
        let gpp = -p * b * b / (1.0 + b * x).powi(2) - (1.0 - p) / (1.0 - x).powi(2);
        x = (x - gp / gpp).clamp(1e-12, 1.0 - 1e-12);
    }
    x
}

/// Numerical root of `g_binary(p, b, f) = 0` for `f > f*`.
pub fn breakeven_binary(p: f64, b: f64) -> f64 {
    let mut lo = f_star_binary(p, b) + 1e-9;
    let mut hi = 1.0 - 1e-9;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if g_binary(p, b, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

pub fn run_binary_study(
    p: f64,
    b: f64,
    n_bets: usize,
    r_paths: usize,
    fs: &[f64],
    main_seed: u64,
    n_workers: usize,
) -> BinaryKellySummary {
    let f_star_exact = f_star_binary(p, b);
    let f_star_numeric = maximize_g_binary(p, b);
    let g_star_exact = g_binary(p, b, f_star_exact);
    let f_breakeven_numeric = breakeven_binary(p, b);

    let mut points = Vec::with_capacity(fs.len());
    for (fi, &f) in fs.iter().enumerate() {
        let g = g_binary(p, b, f);
        let up = (1.0 + b * f).ln();
        let dn = (1.0 - f).ln();

        let res: Vec<(f64, f64, f64, f64)> = run_chunks(r_paths, n_workers, |idx| {
            let mut rng = Rng::new(derive_seed(main_seed, 100 + fi as u64, idx as u64));
            let mut log_w = 0.0_f64;
            let mut touched_half = 0.0_f64;
            let half_ln = 0.5_f64.ln();
            for _ in 0..n_bets {
                if rng.bernoulli(p) {
                    log_w += up;
                } else {
                    log_w += dn;
                    if log_w <= half_ln {
                        touched_half = 1.0;
                    }
                }
            }
            let w = log_w.exp();
            let loss = if log_w < 0.0 { 1.0 } else { 0.0 };
            (log_w, w, loss, touched_half)
        });

        let mut log_ws = Vec::with_capacity(r_paths);
        let mut ws = Vec::with_capacity(r_paths);
        let mut losses = Vec::with_capacity(r_paths);
        let mut touches = Vec::with_capacity(r_paths);
        for (lw, w, l, t) in res {
            log_ws.push(lw);
            ws.push(w);
            losses.push(l);
            touches.push(t);
        }

        let e_ln_w_mc = mean(&log_ws);
        let e_ln_w_se = se(&log_ws);
        let e_ln_w_exact = (n_bets as f64) * g;
        let median_ln_w = quantile(&log_ws, 0.5);
        let median_w_mc = median_ln_w.exp();
        let median_w_exact = e_ln_w_exact.exp();
        let arith_mean_mc = mean(&ws);
        let arith_mean_se = se(&ws);
        let arith_mean_exact = (1.0 + f * (b * p - (1.0 - p))).powi(n_bets as i32);
        let prob_loss_mc = mean(&losses);
        let prob_loss_se = se(&losses);
        let prob_touch_half_mc = mean(&touches);
        let prob_touch_half_se = se(&touches);

        // Continuous GBM approximation matching 1-step log-return mean `g` and variance `v`.
        let s_step = (p * (1.0 - p)).sqrt() * (up - dn);
        let t_tot = n_bets as f64;
        let prob_loss_cont_approx = norm_cdf(-g * t_tot.sqrt() / s_step);
        let a = 0.5_f64.ln();
        let sqrt_t = t_tot.sqrt();
        let term1 = norm_cdf((a - g * t_tot) / (s_step * sqrt_t));
        let exp_arg = (2.0 * g * a) / (s_step * s_step);
        let term2 = if exp_arg < -700.0 {
            0.0
        } else {
            exp_arg.exp() * norm_cdf((a + g * t_tot) / (s_step * sqrt_t))
        };
        let prob_touch_half_cont_approx = (term1 + term2).clamp(0.0, 1.0);

        points.push(BinaryKellyPoint {
            f,
            g_exact: g,
            e_ln_w_mc,
            e_ln_w_se,
            e_ln_w_exact,
            median_w_mc,
            median_w_exact,
            arith_mean_mc,
            arith_mean_se,
            arith_mean_exact,
            prob_loss_mc,
            prob_loss_se,
            prob_loss_cont_approx,
            prob_touch_half_mc,
            prob_touch_half_se,
            prob_touch_half_cont_approx,
        });
    }

    // Generate one deterministic sample trajectory (CRN across 0.5*f*, f*, 2*f*) for plotting.
    let mut rng_plot = Rng::new(derive_seed(main_seed, 199, 0));
    let mut sample_paths_n = Vec::with_capacity(n_bets + 1);
    let mut sample_path_half_kelly = Vec::with_capacity(n_bets + 1);
    let mut sample_path_full_kelly = Vec::with_capacity(n_bets + 1);
    let mut sample_path_over_kelly = Vec::with_capacity(n_bets + 1);
    let mut w_half = 1.0_f64;
    let mut w_full = 1.0_f64;
    let mut w_over = 1.0_f64;
    sample_paths_n.push(0);
    sample_path_half_kelly.push(w_half);
    sample_path_full_kelly.push(w_full);
    sample_path_over_kelly.push(w_over);
    let f_half = 0.5 * f_star_exact;
    let f_over = 2.0 * f_star_exact;
    for step in 1..=n_bets {
        if rng_plot.bernoulli(p) {
            w_half *= 1.0 + b * f_half;
            w_full *= 1.0 + b * f_star_exact;
            w_over *= 1.0 + b * f_over;
        } else {
            w_half *= 1.0 - f_half;
            w_full *= 1.0 - f_star_exact;
            w_over *= 1.0 - f_over;
        }
        sample_paths_n.push(step);
        sample_path_half_kelly.push(w_half);
        sample_path_full_kelly.push(w_full);
        sample_path_over_kelly.push(w_over);
    }

    BinaryKellySummary {
        p,
        b,
        f_star_exact,
        f_star_numeric,
        g_star_exact,
        f_breakeven_numeric,
        points,
        sample_paths_n,
        sample_path_half_kelly,
        sample_path_full_kelly,
        sample_path_over_kelly,
    }
}
