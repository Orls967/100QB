pub mod bs;
pub mod calib;
pub mod jsonout;
pub mod linalg;
pub mod merton;
pub mod parallel;
pub mod pide;
pub mod rng;
pub mod special;
pub mod stats;
pub mod svg;
pub mod tail_risk;

#[cfg(test)]
mod tests {
    use super::bs::{bs_call, bs_put, implied_vol_call};
    use super::calib::calibrate_merton_nelder_mead;
    use super::linalg::solve_tridiagonal;
    use super::merton::{merton_call_series, merton_put_series, simulate_merton_mc, MertonParams};
    use super::pide::price_call_pide_cn;

    #[test]
    fn test_merton_reduces_to_bs_when_lambda_zero() {
        let p = MertonParams {
            lambda: 0.0,
            ..MertonParams::default()
        };
        for k in [80.0, 100.0, 120.0] {
            let c_m = merton_call_series(&p, k, 1.0);
            let c_b = bs_call(p.s0, k, p.r, p.q, p.sigma, 1.0);
            assert!((c_m - c_b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_merton_put_call_parity() {
        let p = MertonParams::default();
        for t in [0.25, 1.0, 2.0] {
            for k in [70.0, 100.0, 130.0] {
                let c = merton_call_series(&p, k, t);
                let pt = merton_put_series(&p, k, t);
                let rhs = p.s0 * (-p.q * t).exp() - k * (-p.r * t).exp();
                assert!((c - pt - rhs).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_bs_implied_vol_roundtrip() {
        let s = 100.0;
        let r = 0.05;
        let q = 0.01;
        let t = 1.0;
        for sigma in [0.10, 0.25, 0.60] {
            for k in [80.0, 100.0, 125.0] {
                let px = bs_call(s, k, r, q, sigma, t);
                let iv = implied_vol_call(px, s, k, r, q, t);
                assert!((iv - sigma).abs() < 1e-9);
            }
        }
        assert!(
            (bs_put(s, 100.0, r, q, 0.2, t)
                - (bs_call(s, 100.0, r, q, 0.2, t) - s * (-q * t).exp() + 100.0 * (-r * t).exp()))
            .abs()
                < 1e-12
        );
    }

    #[test]
    fn test_thomas_tridiagonal_solver() {
        let a = vec![0.0, -1.0, -1.0, -1.0];
        let b = vec![4.0, 4.0, 4.0, 4.0];
        let c = vec![-1.0, -1.0, -1.0, 0.0];
        let x_true = [1.0, 2.0, -1.0, 3.0];
        let mut d = vec![0.0; 4];
        for i in 0..4 {
            d[i] = b[i] * x_true[i];
            if i > 0 {
                d[i] += a[i] * x_true[i - 1];
            }
            if i + 1 < 4 {
                d[i] += c[i] * x_true[i + 1];
            }
        }
        let x_sol = solve_tridiagonal(&a, &b, &c, &d);
        for i in 0..4 {
            assert!((x_sol[i] - x_true[i]).abs() < 1e-13);
        }
    }

    #[test]
    fn test_pide_cn_matches_bs_at_zero_lambda() {
        let p = MertonParams {
            lambda: 0.0,
            ..MertonParams::default()
        };
        let exact = bs_call(p.s0, 100.0, p.r, p.q, p.sigma, 1.0);
        let pide = price_call_pide_cn(&p, 100.0, 1.0, 401, 200);
        assert!((pide - exact).abs() < 5e-3, "pide={pide}, exact={exact}");
    }

    #[test]
    fn test_pide_cn_matches_merton_series() {
        let p = MertonParams::default();
        let exact = merton_call_series(&p, 100.0, 1.0);
        let pide = price_call_pide_cn(&p, 100.0, 1.0, 401, 200);
        assert!((pide - exact).abs() < 2e-2, "pide={pide}, exact={exact}");
    }

    #[test]
    fn test_merton_mc_martingale_and_price_small() {
        let p = MertonParams::default();
        let mc = simulate_merton_mc(&p, &[100.0], 1.0, 30_000, 42, 2);
        assert!(
            mc.martingale_z.abs() < 3.0,
            "martingale_z={}",
            mc.martingale_z
        );
        assert!(
            mc.rows[0].call_z.abs() < 3.0,
            "call_z={}",
            mc.rows[0].call_z
        );
    }

    #[test]
    fn test_merton_calibration_recovers_prices() {
        let p_true = MertonParams::default();
        let p_init = MertonParams {
            sigma: 0.18,
            lambda: 0.50,
            mu_j: -0.05,
            delta_j: 0.20,
            ..p_true
        };
        let cal = calibrate_merton_nelder_mead(
            &p_true,
            &p_init,
            &[0.25, 1.0],
            &[85.0, 100.0, 115.0],
            400,
        );
        assert!(cal.price_rmse < 1e-2, "price_rmse={}", cal.price_rmse);
    }
}
