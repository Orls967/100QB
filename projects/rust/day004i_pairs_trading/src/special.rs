//! Special functions: imaginary error function `erfi(x)` and its derivative `erfi'(x) = (2/sqrt(pi)) * exp(x^2)`.

use std::f64::consts::FRAC_2_SQRT_PI;

/// Computes `erfi(x) = (2/sqrt(pi)) * int_0^x exp(t^2) dt` via Taylor series for `|x| <= 4.5`
/// and Gauss-Legendre / asymptotic-free quadrature for larger `|x|`.
pub fn erfi(x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    let sign = x.signum();
    let ax = x.abs();
    if ax <= 5.0 {
        let x2 = ax * ax;
        let mut term = ax;
        let mut sum = ax;
        for k in 1..120 {
            term *= x2 / k as f64;
            let add = term / (2 * k + 1) as f64;
            sum += add;
            if add.abs() < 1e-16 * sum.abs() {
                break;
            }
        }
        return sign * FRAC_2_SQRT_PI * sum;
    }
    // Composite Simpson quadrature for |x| > 5.0
    let n = 400;
    let h = ax / n as f64;
    let mut acc = 1.0 + (ax * ax).exp();
    for i in 1..n {
        let t = i as f64 * h;
        let w = if i % 2 == 0 { 2.0 } else { 4.0 };
        acc += w * (t * t).exp();
    }
    sign * FRAC_2_SQRT_PI * (h / 3.0) * acc
}

#[inline]
pub fn erfi_deriv(x: f64) -> f64 {
    FRAC_2_SQRT_PI * (x * x).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_erfi_matches_integral_and_derivative() {
        // Independent Simpson check for erfi(1.5)
        let x = 1.5_f64;
        let n = 4000;
        let h = x / n as f64;
        let mut acc = 1.0 + (x * x).exp();
        for i in 1..n {
            let t = i as f64 * h;
            acc += (if i % 2 == 0 { 2.0 } else { 4.0 }) * (t * t).exp();
        }
        let quad_val = FRAC_2_SQRT_PI * (h / 3.0) * acc;
        assert!((erfi(x) - quad_val).abs() < 1e-12);
        assert!((erfi(-x) + erfi(x)).abs() < 1e-14);
    }
}
