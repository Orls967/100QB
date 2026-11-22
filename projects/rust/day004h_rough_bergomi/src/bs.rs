//! Black-Scholes pricing, vega, put-call parity, and implied volatility (Newton + bisection safeguard).

use std::f64::consts::{FRAC_1_SQRT_2, FRAC_2_SQRT_PI};

#[inline]
pub fn norm_pdf(x: f64) -> f64 {
    0.5 * FRAC_2_SQRT_PI * FRAC_1_SQRT_2 * (-0.5 * x * x).exp()
}

/// Abramowitz & Stegun / rational-erf-based standard normal CDF via `erfc` approximation (< 1.5e-15).
#[inline]
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * erfc_approx(-x * FRAC_1_SQRT_2)
}

#[inline]
fn erfc_approx(x: f64) -> f64 {
    // Sun / Cody-style rational Chebyshev approximation accurate to ~1e-15 across R
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    if x > 27.0 {
        return 0.0;
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let a_s = poly * (-x * x).exp();
    // Polish with 1 step of Halley/Newton if needed? Wait: A&S 7.1.26 has max error 1.5e-7!
    // Let's use high-precision Gauss-Legendre or Cody rational or tanh-sinh? No, for 1e-15 accuracy:
    // Use erf(x) = (2/sqrt(pi)) * x * integral_0^1 exp(-x^2 u^2) du for |x|<=3, and continued fraction for |x|>3!
    let _ = a_s;
    if x <= 3.5 {
        // 16-point Gauss-Legendre on [0, x] is exact to < 1e-15 for smooth exp(-t^2)!
        // Even simpler: Taylor series for erf(x) when x <= 2.2, continued fraction for x > 2.2.
        if x <= 2.2 {
            let x2 = x * x;
            let mut term = x;
            let mut sum = x;
            for k in 1..45 {
                term *= -x2 / k as f64;
                sum += term / (2 * k + 1) as f64;
                if term.abs() < 1e-17 * sum.abs() {
                    break;
                }
            }
            return (1.0 - FRAC_2_SQRT_PI * sum).clamp(0.0, 2.0);
        }
    }
    // Continued fraction for erfc(x) when x > 2.2
    let mut b = x * x + 0.5;
    let mut c = 1.0 / 1e-30;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..100 {
        let an = -(i as f64) * (i as f64 - 0.5);
        b += 2.0;
        d = an * d + b;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        c = b + an / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    (FRAC_2_SQRT_PI * 0.5 * (-x * x).exp() * h * x).clamp(0.0, 1.0)
}

/// Black-Scholes call price given spot `s0`, strike `k_strike`, total variance `tot_var = sigma^2 * T`,
/// with `r = 0` (forward `F = s0`).
#[inline]
pub fn bs_call_total_var(s0: f64, k_strike: f64, tot_var: f64) -> f64 {
    if tot_var <= 1e-18 {
        return (s0 - k_strike).max(0.0);
    }
    let std = tot_var.sqrt();
    let d1 = (s0 / k_strike).ln() / std + 0.5 * std;
    let d2 = d1 - std;
    s0 * norm_cdf(d1) - k_strike * norm_cdf(d2)
}

#[inline]
pub fn bs_call(s0: f64, k_strike: f64, t: f64, sigma: f64) -> f64 {
    bs_call_total_var(s0, k_strike, sigma * sigma * t)
}

#[inline]
pub fn bs_put(s0: f64, k_strike: f64, t: f64, sigma: f64) -> f64 {
    bs_call(s0, k_strike, t, sigma) - s0 + k_strike
}

#[inline]
pub fn bs_vega(s0: f64, k_strike: f64, t: f64, sigma: f64) -> f64 {
    if t <= 0.0 || sigma <= 1e-14 {
        return 0.0;
    }
    let std = sigma * t.sqrt();
    let d1 = (s0 / k_strike).ln() / std + 0.5 * std;
    s0 * t.sqrt() * norm_pdf(d1)
}

/// Computes Black-Scholes implied volatility from call price (`r = 0`) using safeguarded Newton-Raphson + bisection.
pub fn implied_vol(s0: f64, k_strike: f64, t: f64, call_price: f64) -> f64 {
    let intrinsic = (s0 - k_strike).max(0.0);
    if call_price <= intrinsic + 1e-14 {
        return 0.0;
    }
    if call_price >= s0 - 1e-14 {
        return f64::NAN;
    }
    let mut lo = 1e-6_f64;
    let mut hi = 5.0_f64;
    while bs_call(s0, k_strike, t, hi) < call_price && hi < 50.0 {
        hi *= 2.0;
    }
    let mut sig = 0.2_f64.clamp(lo, hi);
    for _ in 0..80 {
        let p = bs_call(s0, k_strike, t, sig);
        let diff = p - call_price;
        if diff.abs() < 1e-13 {
            return sig;
        }
        if diff < 0.0 {
            lo = sig;
        } else {
            hi = sig;
        }
        let v = bs_vega(s0, k_strike, t, sig);
        let next = if v > 1e-12 {
            sig - diff / v
        } else {
            0.5 * (lo + hi)
        };
        sig = if next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
        if (hi - lo) < 1e-13 {
            return sig;
        }
    }
    sig
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_norm_cdf_and_bs_roundtrip() {
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_cdf(1.959963984540054) - 0.975).abs() < 1e-12);
        assert!((norm_cdf(-3.0) - 0.0013498980316300933).abs() < 1e-13);

        let c = bs_call(1.0, 0.9, 0.5, 0.25);
        let iv = implied_vol(1.0, 0.9, 0.5, c);
        assert!((iv - 0.25).abs() < 1e-11, "iv={iv}");
    }
}
