use crate::black::black_call;

/// Hagan et al. (2002) asymptotic Black-76 implied volatility for SABR `(alpha, beta, rho, nu)`.
/// Handles `K -> F` and `nu -> 0` (`|z| < 1e-5`) smoothly via Taylor expansion of `z / x(z)`.
pub fn hagan_implied_vol(f: f64, k: f64, t: f64, alpha: f64, beta: f64, rho: f64, nu: f64) -> f64 {
    if f <= 0.0 || k <= 0.0 || t <= 0.0 || alpha <= 0.0 {
        return 0.0;
    }
    let one_m_b = 1.0 - beta;
    let fk = f * k;
    let fk_pow_half = fk.powf(0.5 * one_m_b);
    let fk_pow = fk_pow_half * fk_pow_half;

    let log_fk = (f / k).ln();
    let log2 = log_fk * log_fk;
    let denom = fk_pow_half
        * (1.0 + (one_m_b * one_m_b / 24.0) * log2 + (one_m_b.powi(4) / 1920.0) * log2 * log2);

    let z = (nu / alpha) * fk_pow_half * log_fk;
    let z_over_x = if z.abs() < 1e-5 {
        1.0 - 0.5 * rho * z + ((2.0 - 3.0 * rho * rho) / 12.0) * z * z
    } else {
        let disc = (1.0 - 2.0 * rho * z + z * z).max(0.0).sqrt();
        let arg = (disc + z - rho) / (1.0 - rho);
        if arg <= 0.0 {
            1.0
        } else {
            z / arg.ln()
        }
    };

    let t_corr = 1.0
        + ((one_m_b * one_m_b / 24.0) * (alpha * alpha / fk_pow)
            + (rho * beta * nu * alpha) / (4.0 * fk_pow_half)
            + ((2.0 - 3.0 * rho * rho) / 24.0) * nu * nu)
            * t;

    (alpha / denom) * z_over_x * t_corr
}

/// Hagan SABR forward call price via Black-76 formula with `hagan_implied_vol`.
#[inline]
pub fn hagan_call(f: f64, k: f64, t: f64, alpha: f64, beta: f64, rho: f64, nu: f64) -> f64 {
    let iv = hagan_implied_vol(f, k, t, alpha, beta, rho, nu);
    black_call(f, k, iv, t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hagan_limits() {
        let f = 0.04;
        let alpha = 0.0325;
        let t = 5.0;
        // beta = 1, nu -> 0 => sigma_B = alpha
        let v_ln = hagan_implied_vol(f, 0.06, t, alpha, 1.0, -0.3, 0.0);
        assert!((v_ln - alpha).abs() < 1e-14, "v_ln={v_ln}");

        // ATM limit (K = F): z/x(z) = 1
        let v_atm = hagan_implied_vol(f, f, t, alpha, 0.5, -0.3, 0.4);
        let v_near = hagan_implied_vol(f, f * (1.0 + 1e-8), t, alpha, 0.5, -0.3, 0.4);
        assert!((v_atm - v_near).abs() < 1e-7, "atm={v_atm} near={v_near}");

        // nu -> 0 limit matches nu = 0 CEV expansion smoothly
        let v_nu0 = hagan_implied_vol(f, 0.02, t, alpha, 0.5, -0.3, 0.0);
        let v_small = hagan_implied_vol(f, 0.02, t, alpha, 0.5, -0.3, 1e-6);
        assert!(
            (v_nu0 - v_small).abs() < 1e-5,
            "nu0={v_nu0} small={v_small}"
        );
    }
}
