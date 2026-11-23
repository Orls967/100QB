use crate::gbm::GbmParams;
use crate::hitting::prob_hit_finite;
use crate::parallel::reduce_vec_fixed_blocks;
use crate::rng::{derive_seed, Rng};

#[derive(Clone, Debug)]
pub struct EstimationLossRow {
    pub t_est: f64,
    pub loss_exact: f64,
    pub loss_mc: f64,
    pub loss_se: f64,
    pub z_diff: f64,
    pub loss_alt_params_mc: f64,
    pub loss_alt_params_se: f64,
}

#[derive(Clone, Debug)]
pub struct ShrinkageSummary {
    pub t_est: f64,
    pub c_shrink_opt: f64,
    pub horizon_years: f64,
    pub g_full_plugin: f64,
    pub se_g_full_plugin: f64,
    pub g_half_kelly: f64,
    pub se_g_half_kelly: f64,
    pub g_shrink_opt: f64,
    pub se_g_shrink_opt: f64,
    pub diff_shrink_vs_full: f64,
    pub se_diff_shrink_vs_full: f64,
    pub diff_half_vs_full: f64,
    pub se_diff_half_vs_full: f64,
    pub prob_hit_half_20y_full: f64,
    pub prob_hit_half_20y_half: f64,
    pub prob_hit_half_20y_shrink: f64,
}

/// Tests `E[g(f*) - g(f_hat)] = 1 / (2 * T_est)` via MC (`R = 1_000_000`)
/// and demonstrates invariance to `(mu, sigma)` using an alternative pair.
pub fn run_estimation_loss_study(
    gbm: &GbmParams,
    alt_gbm: &GbmParams,
    t_ests: &[f64],
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> Vec<EstimationLossRow> {
    let g_star = gbm.kelly_growth();
    let g_star_alt = alt_gbm.kelly_growth();

    t_ests
        .iter()
        .enumerate()
        .map(|(ti, &t_est)| {
            let sd_mu = gbm.sigma / t_est.sqrt();
            let sd_mu_alt = alt_gbm.sigma / t_est.sqrt();
            let sums = reduce_vec_fixed_blocks(r_paths, n_workers, 4, |idx, acc| {
                let mut rng = Rng::new(derive_seed(main_seed, 400 + ti as u64, idx as u64));
                let z = rng.normal();
                let mu_hat = gbm.mu + sd_mu * z;
                let f_hat = (mu_hat - gbm.r) / (gbm.sigma * gbm.sigma);
                let loss = g_star - gbm.drift_log(f_hat);

                let mu_hat_alt = alt_gbm.mu + sd_mu_alt * z;
                let f_hat_alt = (mu_hat_alt - alt_gbm.r) / (alt_gbm.sigma * alt_gbm.sigma);
                let loss_alt = g_star_alt - alt_gbm.drift_log(f_hat_alt);

                acc[0] += loss;
                acc[1] += loss * loss;
                acc[2] += loss_alt;
                acc[3] += loss_alt * loss_alt;
            });

            let rf = r_paths as f64;
            let loss_mc = sums[0] / rf;
            let loss_se = ((sums[1] / rf - loss_mc * loss_mc).max(0.0) / rf).sqrt();
            let loss_alt_params_mc = sums[2] / rf;
            let loss_alt_params_se =
                ((sums[3] / rf - loss_alt_params_mc * loss_alt_params_mc).max(0.0) / rf).sqrt();
            let loss_exact = 0.5 / t_est;
            let z_diff = if loss_se > 1e-15 {
                (loss_mc - loss_exact) / loss_se
            } else {
                0.0
            };

            EstimationLossRow {
                t_est,
                loss_exact,
                loss_mc,
                loss_se,
                z_diff,
                loss_alt_params_mc,
                loss_alt_params_se,
            }
        })
        .collect()
}

/// Optimal scalar shrinkage `c* = f*^2 / (f*^2 + Var(f_hat))` maximizing expected log growth
/// under estimation error `mu_hat ~ N(mu, sigma^2 / T_est)`.
pub fn optimal_shrinkage_c(gbm: &GbmParams, t_est: f64) -> f64 {
    let f_star = gbm.kelly_fraction();
    let var_f_hat = 1.0 / (gbm.sigma * gbm.sigma * t_est);
    (f_star * f_star) / (f_star * f_star + var_f_hat)
}

/// CRN comparison of full plug-in `f_hat`, half-Kelly `0.5 * f_hat`, and optimal shrinkage `c* * f_hat`
/// over a 20-year horizon (realized growth + 20y hitting probability of 0.5*W0).
pub fn run_shrinkage_study(
    gbm: &GbmParams,
    t_est: f64,
    horizon_years: f64,
    r_paths: usize,
    main_seed: u64,
    n_workers: usize,
) -> ShrinkageSummary {
    let c_shrink_opt = optimal_shrinkage_c(gbm, t_est);
    let sd_mu = gbm.sigma / t_est.sqrt();
    let s2 = gbm.sigma * gbm.sigma;

    // 10 accumulators:
    // 0,1: g_full sum, sq
    // 2,3: g_half sum, sq
    // 4,5: g_shrink sum, sq
    // 6,7: (g_shrink - g_full) sum, sq
    // 8,9: (g_half - g_full) sum, sq
    // 10,11,12: prob_hit_half_20y sums (full, half, shrink)
    let sums = reduce_vec_fixed_blocks(r_paths, n_workers, 13, |idx, acc| {
        let mut rng = Rng::new(derive_seed(main_seed, 450, idx as u64));
        let z_est = rng.normal();
        let mu_hat = gbm.mu + sd_mu * z_est;
        let f_hat = (mu_hat - gbm.r) / s2;

        let f_full = f_hat;
        let f_half = 0.5 * f_hat;
        let f_shr = c_shrink_opt * f_hat;

        // Realized 20-year log return under true GBM parameters given estimated f
        let z_mkt = rng.normal();
        let sqrt_h = horizon_years.sqrt();

        let g_full_true = gbm.drift_log(f_full);
        let s_full = gbm.vol_log(f_full);
        let r_full = g_full_true + s_full * z_mkt / sqrt_h;

        let g_half_true = gbm.drift_log(f_half);
        let s_half = gbm.vol_log(f_half);
        let r_half = g_half_true + s_half * z_mkt / sqrt_h;

        let g_shr_true = gbm.drift_log(f_shr);
        let s_shr = gbm.vol_log(f_shr);
        let r_shr = g_shr_true + s_shr * z_mkt / sqrt_h;

        let d_sf = r_shr - r_full;
        let d_hf = r_half - r_full;

        acc[0] += r_full;
        acc[1] += r_full * r_full;
        acc[2] += r_half;
        acc[3] += r_half * r_half;
        acc[4] += r_shr;
        acc[5] += r_shr * r_shr;
        acc[6] += d_sf;
        acc[7] += d_sf * d_sf;
        acc[8] += d_hf;
        acc[9] += d_hf * d_hf;

        // Conditional on f, exact 20y hitting prob of 0.5*W0 (Rao-Blackwellized)
        acc[10] += prob_hit_finite(g_full_true, s_full, horizon_years, 0.5);
        acc[11] += prob_hit_finite(g_half_true, s_half, horizon_years, 0.5);
        acc[12] += prob_hit_finite(g_shr_true, s_shr, horizon_years, 0.5);
    });

    let rf = r_paths as f64;
    let stat = |sum: f64, sum_sq: f64| -> (f64, f64) {
        let m = sum / rf;
        let se = ((sum_sq / rf - m * m).max(0.0) / rf).sqrt();
        (m, se)
    };

    let (g_full_plugin, se_g_full_plugin) = stat(sums[0], sums[1]);
    let (g_half_kelly, se_g_half_kelly) = stat(sums[2], sums[3]);
    let (g_shrink_opt, se_g_shrink_opt) = stat(sums[4], sums[5]);
    let (diff_shrink_vs_full, se_diff_shrink_vs_full) = stat(sums[6], sums[7]);
    let (diff_half_vs_full, se_diff_half_vs_full) = stat(sums[8], sums[9]);

    ShrinkageSummary {
        t_est,
        c_shrink_opt,
        horizon_years,
        g_full_plugin,
        se_g_full_plugin,
        g_half_kelly,
        se_g_half_kelly,
        g_shrink_opt,
        se_g_shrink_opt,
        diff_shrink_vs_full,
        se_diff_shrink_vs_full,
        diff_half_vs_full,
        se_diff_half_vs_full,
        prob_hit_half_20y_full: sums[10] / rf,
        prob_hit_half_20y_half: sums[11] / rf,
        prob_hit_half_20y_shrink: sums[12] / rf,
    }
}
