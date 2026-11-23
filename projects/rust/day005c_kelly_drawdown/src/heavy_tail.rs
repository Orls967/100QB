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
    pub g_gauss_daily_mc: f64,
    pub g_gauss_daily_se: f64,
    pub g_t4_mc: f64,
    pub g_t4_se: f64,
    pub diff_g_t4_vs_gauss_daily: f64,
    pub se_diff_g_t4_vs_gauss_daily: f64,
    pub prob_hit_half_gaussian_pred: f64,
    pub prob_hit_half_gauss_daily_mc: f64,
    pub prob_hit_half_gauss_daily_se: f64,
    pub prob_hit_half_t4_mc: f64,
    pub prob_hit_half_t4_se: f64,
    pub diff_hit_t4_vs_gauss_daily: f64,
    pub se_diff_hit_t4_vs_gauss_daily: f64,
    pub effect_discrete_hit: f64,
    pub effect_heavy_tail_hit: f64,
    pub effect_total_hit: f64,
    pub ruin_prob_gauss_daily_mc: f64,
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

    // 14 accumulators per c:
    // 0,1: g_t4 sum, sq | 2,3: hit_t4 sum, sq | 4: ruin_t4 sum
    // 5,6: g_gauss_daily sum, sq | 7,8: hit_gauss_daily sum, sq | 9: ruin_gauss_daily sum
    // 10,11: (g_t4 - g_gauss_daily) sum, sq | 12,13: (hit_t4 - hit_gauss_daily) sum, sq
    let sums = reduce_vec_fixed_blocks(r_paths, n_workers, n_c * 14, |idx, acc| {
        let mut rng_t = Rng::new(derive_seed(main_seed, 500, idx as u64));
        let mut rng_g = Rng::new(derive_seed(main_seed, 501, idx as u64));
        let mut log_w_t = vec![0.0_f64; n_c];
        let mut hit_t = vec![false; n_c];
        let mut ruin_t = vec![false; n_c];
        let mut log_w_g = vec![0.0_f64; n_c];
        let mut hit_g = vec![false; n_c];
        let mut ruin_g = vec![false; n_c];

        for _ in 0..n_steps {
            let z_t = rng_t.student_t(nu) * t_scale;
            let z_g = rng_g.normal();
            let r_ex_t = (gbm.mu - gbm.r) * dt + gbm.sigma * sqrt_dt * z_t;
            let r_ex_g = (gbm.mu - gbm.r) * dt + gbm.sigma * sqrt_dt * z_g;
            for ci in 0..n_c {
                if !ruin_t[ci] {
                    let gross = 1.0 + gbm.r * dt + fs[ci] * r_ex_t;
                    if gross <= ruin_floor {
                        ruin_t[ci] = true;
                        hit_t[ci] = true;
                        log_w_t[ci] = ruin_floor.ln();
                    } else {
                        log_w_t[ci] += gross.ln();
                        if log_w_t[ci] <= half_ln {
                            hit_t[ci] = true;
                        }
                    }
                }
                if !ruin_g[ci] {
                    let gross = 1.0 + gbm.r * dt + fs[ci] * r_ex_g;
                    if gross <= ruin_floor {
                        ruin_g[ci] = true;
                        hit_g[ci] = true;
                        log_w_g[ci] = ruin_floor.ln();
                    } else {
                        log_w_g[ci] += gross.ln();
                        if log_w_g[ci] <= half_ln {
                            hit_g[ci] = true;
                        }
                    }
                }
            }
        }

        for ci in 0..n_c {
            let gt = log_w_t[ci] / horizon_years;
            let gg = log_w_g[ci] / horizon_years;
            let ht = if hit_t[ci] { 1.0 } else { 0.0 };
            let hg = if hit_g[ci] { 1.0 } else { 0.0 };
            let rt = if ruin_t[ci] { 1.0 } else { 0.0 };
            let rg = if ruin_g[ci] { 1.0 } else { 0.0 };
            let dg = gt - gg;
            let dh = ht - hg;
            let off = ci * 14;
            acc[off] += gt;
            acc[off + 1] += gt * gt;
            acc[off + 2] += ht;
            acc[off + 3] += ht * ht;
            acc[off + 4] += rt;
            acc[off + 5] += gg;
            acc[off + 6] += gg * gg;
            acc[off + 7] += hg;
            acc[off + 8] += hg * hg;
            acc[off + 9] += rg;
            acc[off + 10] += dg;
            acc[off + 11] += dg * dg;
            acc[off + 12] += dh;
            acc[off + 13] += dh * dh;
        }
    });

    let rf = r_paths as f64;
    let stat = |sum: f64, sum_sq: f64| -> (f64, f64) {
        let m = sum / rf;
        let se = ((sum_sq / rf - m * m).max(0.0) / rf).sqrt();
        (m, se)
    };
    cs.iter()
        .enumerate()
        .map(|(ci, &c)| {
            let f = fs[ci];
            let off = ci * 14;
            let (g_t4_mc, g_t4_se) = stat(sums[off], sums[off + 1]);
            let (prob_hit_half_t4_mc, prob_hit_half_t4_se) = stat(sums[off + 2], sums[off + 3]);
            let ruin_prob_t4_mc = sums[off + 4] / rf;
            let (g_gauss_daily_mc, g_gauss_daily_se) = stat(sums[off + 5], sums[off + 6]);
            let (prob_hit_half_gauss_daily_mc, prob_hit_half_gauss_daily_se) =
                stat(sums[off + 7], sums[off + 8]);
            let ruin_prob_gauss_daily_mc = sums[off + 9] / rf;
            let (diff_g_t4_vs_gauss_daily, se_diff_g_t4_vs_gauss_daily) =
                stat(sums[off + 10], sums[off + 11]);
            let (diff_hit_t4_vs_gauss_daily, se_diff_hit_t4_vs_gauss_daily) =
                stat(sums[off + 12], sums[off + 13]);

            let g_gaussian_pred = gbm.drift_log(f);
            let s_pred = gbm.vol_log(f);
            let prob_hit_half_gaussian_pred =
                prob_hit_finite(g_gaussian_pred, s_pred, horizon_years, 0.5);

            HeavyTailRow {
                c,
                f,
                g_gaussian_pred,
                g_gauss_daily_mc,
                g_gauss_daily_se,
                g_t4_mc,
                g_t4_se,
                diff_g_t4_vs_gauss_daily,
                se_diff_g_t4_vs_gauss_daily,
                prob_hit_half_gaussian_pred,
                prob_hit_half_gauss_daily_mc,
                prob_hit_half_gauss_daily_se,
                prob_hit_half_t4_mc,
                prob_hit_half_t4_se,
                diff_hit_t4_vs_gauss_daily,
                se_diff_hit_t4_vs_gauss_daily,
                effect_discrete_hit: prob_hit_half_gauss_daily_mc - prob_hit_half_gaussian_pred,
                effect_heavy_tail_hit: diff_hit_t4_vs_gauss_daily,
                effect_total_hit: prob_hit_half_t4_mc - prob_hit_half_gaussian_pred,
                ruin_prob_gauss_daily_mc,
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
