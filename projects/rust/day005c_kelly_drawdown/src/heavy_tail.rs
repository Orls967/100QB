use crate::estimation::optimal_shrinkage_c;
use crate::gbm::GbmParams;
use crate::hitting::{prob_hit_finite, prob_hit_infinite};
use crate::parallel::reduce_vec_fixed_blocks;
use crate::rng::{derive_seed, Rng};
use crate::special::norm_cdf;

#[derive(Clone, Debug)]
pub struct HeavyTailRow {
    pub c: f64,
    pub f: f64,
    pub g_gaussian_pred: f64,
    pub g_t4_mc: f64,
    pub g_t4_se: f64,
    pub prob_hit_half_gaussian_pred: f64,
    pub prob_hit_half_t4_mc: f64,
    pub prob_hit_half_t4_se: f64,
    pub ruin_prob_t4_mc: f64,
}

#[derive(Clone, Debug)]
pub struct PracticalSummaryRow {
    pub c: f64,
    pub f: f64,
    pub g_no_err: f64,
    pub g_with_err: f64,
    pub prob_hit_half_inf_no_err: f64,
    pub prob_hit_half_20y_no_err: f64,
    pub prob_hit_half_20y_with_err: f64,
    pub prob_loss_5y_no_err: f64,
    pub prob_loss_5y_with_err: f64,
}

/// Daily Student-t (nu = 4) returns scaled to match GBM daily variance (`1 / sqrt(nu / (nu - 2))`),
/// compared against Gaussian GBM predictions over `horizon_years` (`252` steps/year).
pub fn run_heavy_tail_study(
    gbm: &GbmParams,
    nu: f64,
    cs: &[f64],
    horizon_years: f64,
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> Vec<HeavyTailRow> {
    let dt = 1.0_f64 / 252.0;
    let sqrt_dt = dt.sqrt();
    let n_steps = (horizon_years * 252.0).round() as usize;
    let t_scale = ((nu - 2.0) / nu).sqrt(); // Var(t_nu) = nu/(nu-2) => unit variance
    let f_star = gbm.kelly_fraction();
    let n_c = cs.len();
    let fs: Vec<f64> = cs.iter().map(|&c| c * f_star).collect();
    let half_ln = 0.5_f64.ln();
    let ruin_floor = 1e-12_f64;

    // 5 accumulators per c: [sum_g, sum_g_sq, sum_hit, sum_hit_sq, sum_ruin]
    let sums = reduce_vec_fixed_blocks(r_paths, n_workers, n_c * 5, |idx, acc| {
        let mut rng = Rng::new(derive_seed(main_seed, 500, idx as u64));
        let mut log_w = vec![0.0_f64; n_c];
        let mut hit_half = vec![false; n_c];
        let mut ruined = vec![false; n_c];

        for _ in 0..n_steps {
            let z_t = rng.student_t(nu) * t_scale;
            let r_asset_excess = (gbm.mu - gbm.r) * dt + gbm.sigma * sqrt_dt * z_t;
            for ci in 0..n_c {
                if ruined[ci] {
                    continue;
                }
                let ret = gbm.r * dt + fs[ci] * r_asset_excess;
                let gross = 1.0 + ret;
                if gross <= ruin_floor {
                    ruined[ci] = true;
                    hit_half[ci] = true;
                    log_w[ci] = ruin_floor.ln();
                } else {
                    log_w[ci] += gross.ln();
                    if log_w[ci] <= half_ln {
                        hit_half[ci] = true;
                    }
                }
            }
        }

        for ci in 0..n_c {
            let g_real = log_w[ci] / horizon_years;
            let h = if hit_half[ci] { 1.0 } else { 0.0 };
            let r = if ruined[ci] { 1.0 } else { 0.0 };
            let off = ci * 5;
            acc[off] += g_real;
            acc[off + 1] += g_real * g_real;
            acc[off + 2] += h;
            acc[off + 3] += h * h;
            acc[off + 4] += r;
        }
    });

    let rf = r_paths as f64;
    cs.iter()
        .enumerate()
        .map(|(ci, &c)| {
            let f = fs[ci];
            let off = ci * 5;
            let g_t4_mc = sums[off] / rf;
            let g_t4_se = ((sums[off + 1] / rf - g_t4_mc * g_t4_mc).max(0.0) / rf).sqrt();
            let prob_hit_half_t4_mc = sums[off + 2] / rf;
            let prob_hit_half_t4_se =
                ((sums[off + 3] / rf - prob_hit_half_t4_mc * prob_hit_half_t4_mc).max(0.0) / rf)
                    .sqrt();
            let ruin_prob_t4_mc = sums[off + 4] / rf;

            let g_gaussian_pred = gbm.drift_log(f);
            let s_pred = gbm.vol_log(f);
            let prob_hit_half_gaussian_pred =
                prob_hit_finite(g_gaussian_pred, s_pred, horizon_years, 0.5);

            HeavyTailRow {
                c,
                f,
                g_gaussian_pred,
                g_t4_mc,
                g_t4_se,
                prob_hit_half_gaussian_pred,
                prob_hit_half_t4_mc,
                prob_hit_half_t4_se,
                ruin_prob_t4_mc,
            }
        })
        .collect()
}

/// Builds the practical summary table for `f / f* in {0.25, 0.5, 1.0, 1.5, 2.0}`
/// with and without estimation error (`T_est = 5` years, `r_paths = 200_000` MC for plug-in).
pub fn build_practical_summary(
    gbm: &GbmParams,
    cs: &[f64],
    t_est: f64,
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> Vec<PracticalSummaryRow> {
    let f_star = gbm.kelly_fraction();
    let sd_mu = gbm.sigma / t_est.sqrt();
    let s2 = gbm.sigma * gbm.sigma;
    let n_c = cs.len();

    // 3 accumulators per c: [g_with_err, hit_half_20y_with_err, loss_5y_with_err]
    let sums = reduce_vec_fixed_blocks(r_paths, n_workers, n_c * 3, |idx, acc| {
        let mut rng = Rng::new(derive_seed(main_seed, 550, idx as u64));
        let mu_hat = gbm.mu + sd_mu * rng.normal();
        let f_hat = (mu_hat - gbm.r) / s2;

        for (ci, &c) in cs.iter().enumerate() {
            let f = c * f_hat;
            let g_true = gbm.drift_log(f);
            let s_true = gbm.vol_log(f);
            let p_hit_20 = prob_hit_finite(g_true, s_true, 20.0, 0.5);
            let p_loss_5 = if s_true > 1e-12 {
                norm_cdf(-g_true * 5.0_f64.sqrt() / s_true)
            } else if g_true < 0.0 {
                1.0
            } else {
                0.0
            };
            let off = ci * 3;
            acc[off] += g_true;
            acc[off + 1] += p_hit_20;
            acc[off + 2] += p_loss_5;
        }
    });

    let _ = optimal_shrinkage_c(gbm, t_est);
    let rf = r_paths as f64;
    cs.iter()
        .enumerate()
        .map(|(ci, &c)| {
            let f = c * f_star;
            let g_no_err = gbm.drift_log(f);
            let s_no_err = gbm.vol_log(f);
            let off = ci * 3;
            PracticalSummaryRow {
                c,
                f,
                g_no_err,
                g_with_err: sums[off] / rf,
                prob_hit_half_inf_no_err: prob_hit_infinite(g_no_err, s_no_err, 0.5),
                prob_hit_half_20y_no_err: prob_hit_finite(g_no_err, s_no_err, 20.0, 0.5),
                prob_hit_half_20y_with_err: sums[off + 1] / rf,
                prob_loss_5y_no_err: norm_cdf(-g_no_err * 5.0_f64.sqrt() / s_no_err),
                prob_loss_5y_with_err: sums[off + 2] / rf,
            }
        })
        .collect()
}
