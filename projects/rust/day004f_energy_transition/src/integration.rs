//! Panel F: 10-year battery storage project discounted with the AFNS zero-coupon curve,
//! AFNS parameter-uncertainty sensitivity (4b/4c/4d finite fits), duration/convexity/KRD,
//! parallel PV01, and two-bond KRD hedge.

use crate::afns::AfnsParams;
use crate::kalman::MleResult;
use crate::linalg;
use crate::ns::NsParams;
use crate::risk::{self, KrdReport, STANDARD_KRD_NODES};
use crate::rng::{derive_seed, Rng, SEED_DAY004F_INTEGRATION};
use crate::stats;

#[derive(Clone, Debug)]
pub struct AfnsCurveFitDiag {
    pub source: &'static str,
    pub total_fits: usize,
    pub finite_fits: usize,
    pub non_converged_count: usize,
    pub bound_hit_count: usize,
    pub npv_mean: f64,
    pub npv_sd: f64,
    pub npv_q05: f64,
    pub npv_median: f64,
    pub npv_q95: f64,
}

#[derive(Clone, Debug)]
pub struct CapexScenario {
    pub capex_eur: f64,
    pub annual_om_eur: f64,
    pub post_capex_npv_eur: f64,
}

#[derive(Clone, Debug)]
pub struct ProjectTwoBondHedge {
    pub node1: f64,
    pub node2: f64,
    pub weight_bond_5y_4pct: f64,
    pub weight_bond_10y_2pct: f64,
    pub residual_node1: f64,
    pub residual_node2: f64,
    pub target_total_dollar_dur: f64,
    pub hedged_total_dollar_dur: f64,
}

#[derive(Clone, Debug)]
pub struct IntegrationResult {
    pub elapsed_seconds: f64,
    pub project_years: usize,
    pub days_per_year: usize,
    pub num_resamples: usize,
    pub mean_annual_cf_eur: f64,
    pub sd_annual_cf_eur: f64,
    pub base_afns_expected_annual_cfs: Vec<f64>,
    pub base_npv_pre_capex_eur: f64,
    pub break_even_capex_eur: f64,
    pub resampled_npv_mean_eur: f64,
    pub resampled_npv_sd_eur: f64,
    pub resampled_npv_mcse_eur: f64,
    pub resampled_npv_q05_eur: f64,
    pub resampled_npv_median_eur: f64,
    pub resampled_npv_q95_eur: f64,
    pub resampled_npvs: Vec<f64>,
    pub afns_param_sensitivities: Vec<AfnsCurveFitDiag>,
    pub eff_duration_closed: f64,
    pub eff_duration_fd: f64,
    pub duration_rel_err: f64,
    pub eff_convexity_closed: f64,
    pub eff_convexity_fd: f64,
    pub convexity_rel_err: f64,
    pub parallel_pv01_eur: f64,
    pub krd: KrdReport,
    pub hedge: ProjectTwoBondHedge,
    pub capex_scenarios: Vec<CapexScenario>,
}

pub fn run_integration(
    dp33_daily_profits: &[f64],
    fits_4b: &[MleResult],
    fits_4c: &[MleResult],
    fits_4d: &[MleResult],
    quick: bool,
) -> IntegrationResult {
    let t0 = std::time::Instant::now();
    let project_years = 10usize;
    let days_per_year = 365usize;
    let num_resamples = if quick { 200 } else { 1000 };
    let n_days_src = dp33_daily_profits.len().max(1);

    let base_afns = AfnsParams::default_params();
    let base_yield = |tau: f64| base_afns.yield_at(base_afns.theta_p, tau);
    let times: Vec<f64> = (1..=project_years).map(|y| y as f64).collect();
    let dfs_base: Vec<f64> = times.iter().map(|&t| (-t * base_yield(t)).exp()).collect();

    let mean_daily = stats::mean(dp33_daily_profits);
    let expected_annual_cf = mean_daily * days_per_year as f64;
    let expected_cfs = vec![expected_annual_cf; project_years];
    let base_npv_pre_capex: f64 = expected_cfs
        .iter()
        .zip(dfs_base.iter())
        .map(|(c, d)| c * d)
        .sum();

    let mut rg = Rng::new(derive_seed(SEED_DAY004F_INTEGRATION, 0xF001, 0));
    let mut resampled_npvs = Vec::with_capacity(num_resamples);
    let mut all_annual_cfs = Vec::with_capacity(num_resamples * project_years);
    for _ in 0..num_resamples {
        let mut npv = 0.0;
        for y in 0..project_years {
            let mut ann_cf = 0.0;
            for _ in 0..days_per_year {
                let idx = (rg.next_u64() as usize) % n_days_src;
                ann_cf += dp33_daily_profits[idx];
            }
            all_annual_cfs.push(ann_cf);
            npv += ann_cf * dfs_base[y];
        }
        resampled_npvs.push(npv);
    }

    let eval_fit_group = |source: &'static str, fits: &[MleResult]| -> AfnsCurveFitDiag {
        let mut npvs = Vec::new();
        let mut non_conv = 0usize;
        let mut bnd_hit = 0usize;
        for f in fits {
            if !f.converged {
                non_conv += 1;
            }
            if f.hit_bound {
                bnd_hit += 1;
            }
            if !f.finite {
                continue;
            }
            let p = AfnsParams {
                lambda: f.est_params[0],
                kp: [
                    [f.est_params[1], 0.0, 0.0],
                    [0.0, f.est_params[2], 0.0],
                    [0.0, 0.0, f.est_params[3]],
                ],
                theta_p: [f.est_params[4], f.est_params[5], f.est_params[6]],
                sigma: [f.est_params[7], f.est_params[8], f.est_params[9]],
            };
            let npv: f64 = times
                .iter()
                .map(|&t| expected_annual_cf * (-t * p.yield_at(p.theta_p, t)).exp())
                .sum();
            if npv.is_finite() {
                npvs.push(npv);
            }
        }
        AfnsCurveFitDiag {
            source,
            total_fits: fits.len(),
            finite_fits: npvs.len(),
            non_converged_count: non_conv,
            bound_hit_count: bnd_hit,
            npv_mean: stats::mean(&npvs),
            npv_sd: stats::sd(&npvs),
            npv_q05: stats::quantile(&npvs, 0.05),
            npv_median: stats::quantile(&npvs, 0.50),
            npv_q95: stats::quantile(&npvs, 0.95),
        }
    };

    let afns_param_sensitivities = vec![
        eval_fit_group("day004b_baseline_1000w", fits_4b),
        eval_fit_group("day004c_outer_20y", fits_4c),
        eval_fit_group("day004d_outer_20y", fits_4d),
    ];

    // Duration, convexity, PV01, KRD, and two-bond hedge on base AFNS curve
    let bump_d = 1e-6;
    let bump_c = 1e-4;
    let mut sum_t = 0.0;
    let mut sum_t2 = 0.0;
    for (&t, &cf) in times.iter().zip(expected_cfs.iter()) {
        let pv_cf = cf * (-t * base_yield(t)).exp();
        sum_t += t * pv_cf;
        sum_t2 += t * t * pv_cf;
    }
    let eff_dur_closed = sum_t / base_npv_pre_capex;
    let eff_conv_closed = sum_t2 / base_npv_pre_capex;

    let price_shift = |shift: f64| -> f64 {
        times
            .iter()
            .zip(expected_cfs.iter())
            .map(|(&t, &cf)| cf * (-t * (base_yield(t) + shift)).exp())
            .sum()
    };
    let eff_dur_fd =
        -(price_shift(bump_d) - price_shift(-bump_d)) / (2.0 * bump_d * base_npv_pre_capex);
    let eff_conv_fd = (price_shift(bump_c) - 2.0 * base_npv_pre_capex + price_shift(-bump_c))
        / (bump_c * bump_c * base_npv_pre_capex);

    let krd = risk::compute_krd_generic(
        &base_yield,
        &times,
        &expected_cfs,
        &STANDARD_KRD_NODES,
        bump_d,
    );
    let parallel_pv01 = sum_t * 1e-4;

    // Hedge KRD at 5y (idx 3) and 10y (idx 5) using 5y 4% and 10y 2% bonds priced on the same AFNS curve
    let (t_h1, cf_h1) = risk::bond_cashflows(5.0, 0.04, 2.0);
    let (t_h2, cf_h2) = risk::bond_cashflows(10.0, 0.02, 2.0);
    let k_h1 = risk::compute_krd_generic(&base_yield, &t_h1, &cf_h1, &STANDARD_KRD_NODES, bump_d);
    let k_h2 = risk::compute_krd_generic(&base_yield, &t_h2, &cf_h2, &STANDARD_KRD_NODES, bump_d);
    let idx1 = 3usize;
    let idx2 = 5usize;
    let mat = vec![
        vec![k_h1.dollar_krd_closed[idx1], k_h2.dollar_krd_closed[idx1]],
        vec![k_h1.dollar_krd_closed[idx2], k_h2.dollar_krd_closed[idx2]],
    ];
    let rhs = vec![-krd.dollar_krd_closed[idx1], -krd.dollar_krd_closed[idx2]];
    let w = linalg::solve_lu(&mat, &rhs).unwrap_or_else(|_| vec![0.0, 0.0]);
    let res1 = krd.dollar_krd_closed[idx1]
        + w[0] * k_h1.dollar_krd_closed[idx1]
        + w[1] * k_h2.dollar_krd_closed[idx1];
    let res2 = krd.dollar_krd_closed[idx2]
        + w[0] * k_h1.dollar_krd_closed[idx2]
        + w[1] * k_h2.dollar_krd_closed[idx2];
    let tgt_tot: f64 = krd.dollar_krd_closed.iter().sum();
    let hedged_tot: f64 = (0..STANDARD_KRD_NODES.len())
        .map(|i| {
            krd.dollar_krd_closed[i]
                + w[0] * k_h1.dollar_krd_closed[i]
                + w[1] * k_h2.dollar_krd_closed[i]
        })
        .sum();

    let pv_annuity: f64 = dfs_base.iter().sum();
    let capex_scenarios = vec![
        CapexScenario {
            capex_eur: 400_000.0,
            annual_om_eur: 15_000.0,
            post_capex_npv_eur: base_npv_pre_capex - 400_000.0 - 15_000.0 * pv_annuity,
        },
        CapexScenario {
            capex_eur: 600_000.0,
            annual_om_eur: 15_000.0,
            post_capex_npv_eur: base_npv_pre_capex - 600_000.0 - 15_000.0 * pv_annuity,
        },
        CapexScenario {
            capex_eur: 800_000.0,
            annual_om_eur: 15_000.0,
            post_capex_npv_eur: base_npv_pre_capex - 800_000.0 - 15_000.0 * pv_annuity,
        },
    ];

    let _ = NsParams::default_ns();
    IntegrationResult {
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        project_years,
        days_per_year,
        num_resamples,
        mean_annual_cf_eur: stats::mean(&all_annual_cfs),
        sd_annual_cf_eur: stats::sd(&all_annual_cfs),
        base_afns_expected_annual_cfs: expected_cfs,
        base_npv_pre_capex_eur: base_npv_pre_capex,
        break_even_capex_eur: base_npv_pre_capex,
        resampled_npv_mean_eur: stats::mean(&resampled_npvs),
        resampled_npv_sd_eur: stats::sd(&resampled_npvs),
        resampled_npv_mcse_eur: stats::se(&resampled_npvs),
        resampled_npv_q05_eur: stats::quantile(&resampled_npvs, 0.05),
        resampled_npv_median_eur: stats::quantile(&resampled_npvs, 0.50),
        resampled_npv_q95_eur: stats::quantile(&resampled_npvs, 0.95),
        resampled_npvs,
        afns_param_sensitivities,
        eff_duration_closed: eff_dur_closed,
        eff_duration_fd: eff_dur_fd,
        duration_rel_err: (eff_dur_fd - eff_dur_closed).abs() / eff_dur_closed.abs().max(1e-12),
        eff_convexity_closed: eff_conv_closed,
        eff_convexity_fd: eff_conv_fd,
        convexity_rel_err: (eff_conv_fd - eff_conv_closed).abs() / eff_conv_closed.abs().max(1e-12),
        parallel_pv01_eur: parallel_pv01,
        krd,
        hedge: ProjectTwoBondHedge {
            node1: STANDARD_KRD_NODES[idx1],
            node2: STANDARD_KRD_NODES[idx2],
            weight_bond_5y_4pct: w[0],
            weight_bond_10y_2pct: w[1],
            residual_node1: res1,
            residual_node2: res2,
            target_total_dollar_dur: tgt_tot,
            hedged_total_dollar_dur: hedged_tot,
        },
        capex_scenarios,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_integration_krd_and_hedge_anchors() {
        let daily = vec![235.0; 100];
        let res = run_integration(&daily, &[], &[], &[], true);
        assert!(res.duration_rel_err < 1e-6);
        assert!(res.convexity_rel_err < 1e-6);
        assert!(res.krd.sum_vs_par_closed_err < 1e-10);
        assert!(res.hedge.residual_node1.abs() < 1e-9);
        assert!(res.hedge.residual_node2.abs() < 1e-9);
    }
}
