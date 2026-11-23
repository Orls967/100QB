pub mod estimation;
pub mod fractional;
pub mod gbm;
pub mod heavy_tail;
pub mod hitting;
pub mod jsonout;
pub mod kelly_binary;
pub mod parallel;
pub mod rng;
pub mod special;
pub mod stats;
pub mod svg;

#[cfg(test)]
mod tests {
    use super::estimation::{optimal_shrinkage_c, run_estimation_loss_study};
    use super::fractional::compute_fractional_roots;
    use super::gbm::GbmParams;
    use super::hitting::{bridge_cross_prob, prob_hit_finite, prob_hit_infinite, run_hitting_grid};
    use super::kelly_binary::{f_star_binary, g_binary, maximize_g_binary};
    use super::special::norm_cdf;

    #[test]
    fn test_binary_kelly_exact_vs_numeric() {
        let p = 0.55;
        let b = 1.0;
        let f_exact = f_star_binary(p, b);
        let f_num = maximize_g_binary(p, b);
        assert!((f_exact - 0.10).abs() < 1e-15);
        assert!((f_num - f_exact).abs() < 1e-10);
        assert!(g_binary(p, b, f_exact) > g_binary(p, b, 0.05));
        assert!(g_binary(p, b, f_exact) > g_binary(p, b, 0.15));
    }

    #[test]
    fn test_gbm_kelly_and_hitting_formulas() {
        let gbm = GbmParams::default();
        assert!((gbm.kelly_fraction() - 1.25).abs() < 1e-12);
        assert!((gbm.kelly_growth() - 0.06125).abs() < 1e-12);

        // r = 0 => P(inf hit x) = x at f = f*
        let gbm_r0 = GbmParams {
            mu: 0.08,
            sigma: 0.20,
            r: 0.0,
        };
        let f_star_r0 = gbm_r0.kelly_fraction();
        let p_inf = prob_hit_infinite(gbm_r0.drift_log(f_star_r0), gbm_r0.vol_log(f_star_r0), 0.5);
        assert!((p_inf - 0.5).abs() < 1e-12);

        // Finite-horizon converges to infinite-horizon as T -> large
        let f_star = gbm.kelly_fraction();
        let p_fin_large =
            prob_hit_finite(gbm.drift_log(f_star), gbm.vol_log(f_star), 10_000.0, 0.5);
        let p_inf_default = prob_hit_infinite(gbm.drift_log(f_star), gbm.vol_log(f_star), 0.5);
        assert!((p_fin_large - p_inf_default).abs() < 1e-8);
    }

    #[test]
    fn test_norm_cdf_and_fractional_roots() {
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_cdf(1.959963984540054) - 0.975).abs() < 1e-12);

        let gbm = GbmParams::default();
        let roots = compute_fractional_roots(&gbm);
        assert!((roots.c_zero_excess_closed - roots.c_zero_excess_numeric).abs() < 1e-10);
        assert!((roots.c_hit_half_eq_half_closed - roots.c_hit_half_eq_half_numeric).abs() < 1e-10);
        assert!(
            (roots.c_hit_half_eq_quarter_closed - roots.c_hit_half_eq_quarter_numeric).abs()
                < 1e-10
        );
    }

    #[test]
    fn test_bridge_unbiased_and_estimation_loss_small() {
        let gbm = GbmParams::default();
        let rows = run_hitting_grid(&gbm, &[1.0 / 12.0], &[5.0], &[0.5], &[1.25], 20_000, 42, 2);
        assert_eq!(rows.len(), 1);
        assert!(rows[0].bridge_z.abs() < 3.0);

        let alt = GbmParams {
            mu: 0.12,
            sigma: 0.35,
            r: 0.01,
        };
        let est = run_estimation_loss_study(&gbm, &alt, &[5.0], 20_000, 42, 2);
        assert!(est[0].z_diff.abs() < 3.0);
        assert!(optimal_shrinkage_c(&gbm, 5.0) < 1.0);
    }

    #[test]
    fn test_binary_kelly_fine_grid_peak() {
        let p = 0.55;
        let b = 1.0;
        let f_star = f_star_binary(p, b);
        let g_star = g_binary(p, b, f_star);
        for k in 1..200 {
            let f = k as f64 * 0.002;
            if (f - f_star).abs() > 1e-12 {
                assert!(g_binary(p, b, f) < g_star);
            }
        }
    }

    #[test]
    fn test_hitting_r_zero_identity_multiple_barriers() {
        let gbm_r0 = GbmParams {
            mu: 0.09,
            sigma: 0.25,
            r: 0.0,
        };
        let f_star = gbm_r0.kelly_fraction();
        let g = gbm_r0.drift_log(f_star);
        let s = gbm_r0.vol_log(f_star);
        for x in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let p_inf = prob_hit_infinite(g, s, x);
            assert!((p_inf - x).abs() < 1e-12, "x={x}, p_inf={p_inf}");
        }
    }

    #[test]
    fn test_bridge_cross_prob_at_or_below_barrier() {
        let a = 0.5_f64.ln();
        assert_eq!(bridge_cross_prob(a, a + 0.1, a, 0.01), 1.0);
        assert_eq!(bridge_cross_prob(a + 0.1, a, a, 0.01), 1.0);
        assert_eq!(bridge_cross_prob(a - 0.05, a + 0.2, a, 0.01), 1.0);
        let p_above = bridge_cross_prob(a + 0.1, a + 0.1, a, 0.01);
        assert!(p_above > 0.0 && p_above < 1.0);
    }

    #[test]
    fn test_shrinkage_monotonic_and_bounds() {
        let gbm = GbmParams::default();
        let c1 = optimal_shrinkage_c(&gbm, 1.0);
        let c5 = optimal_shrinkage_c(&gbm, 5.0);
        let c20 = optimal_shrinkage_c(&gbm, 20.0);
        let c_large = optimal_shrinkage_c(&gbm, 100_000.0);
        assert!(0.0 < c1 && c1 < c5 && c5 < c20 && c20 < c_large && c_large < 1.0);
        assert!(c_large > 0.999);
    }
}
