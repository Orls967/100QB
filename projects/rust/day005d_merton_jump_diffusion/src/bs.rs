use crate::special::{norm_cdf, norm_pdf};

#[inline]
pub fn bs_call(s: f64, k: f64, r: f64, q: f64, sigma: f64, t: f64) -> f64 {
    if t <= 0.0 {
        return (s - k).max(0.0);
    }
    let disc_s = s * (-q * t).exp();
    let disc_k = k * (-r * t).exp();
    if sigma <= 0.0 {
        return (disc_s - disc_k).max(0.0);
    }
    let vs = sigma * t.sqrt();
    let d1 = ((s / k).ln() + (r - q + 0.5 * sigma * sigma) * t) / vs;
    let d2 = d1 - vs;
    disc_s * norm_cdf(d1) - disc_k * norm_cdf(d2)
}

#[inline]
pub fn bs_put(s: f64, k: f64, r: f64, q: f64, sigma: f64, t: f64) -> f64 {
    let call = bs_call(s, k, r, q, sigma, t);
    call - s * (-q * t).exp() + k * (-r * t).exp()
}

#[inline]
pub fn bs_vega(s: f64, k: f64, r: f64, q: f64, sigma: f64, t: f64) -> f64 {
    if t <= 0.0 || sigma <= 0.0 {
        return 0.0;
    }
    let vs = sigma * t.sqrt();
    let d1 = ((s / k).ln() + (r - q + 0.5 * sigma * sigma) * t) / vs;
    s * (-q * t).exp() * norm_pdf(d1) * t.sqrt()
}

/// Implied volatility from European call price via Newton-Raphson + bisection fallback.
pub fn implied_vol_call(target_price: f64, s: f64, k: f64, r: f64, q: f64, t: f64) -> f64 {
    let intrinsic = (s * (-q * t).exp() - k * (-r * t).exp()).max(0.0);
    let upper = s * (-q * t).exp();
    if target_price <= intrinsic + 1e-12 {
        return 0.0;
    }
    if target_price >= upper - 1e-12 {
        return 5.0;
    }
    let mut lo = 1e-6_f64;
    let mut hi = 5.0_f64;
    let mut vol = 0.25_f64;
    for _ in 0..40 {
        let p = bs_call(s, k, r, q, vol, t);
        let diff = p - target_price;
        if diff.abs() < 1e-11 {
            return vol;
        }
        if diff < 0.0 {
            lo = lo.max(vol);
        } else {
            hi = hi.min(vol);
        }
        let vega = bs_vega(s, k, r, q, vol, t);
        if vega > 1e-12 {
            let next = vol - diff / vega;
            if next > lo && next < hi {
                vol = next;
                continue;
            }
        }
        vol = 0.5 * (lo + hi);
    }
    vol
}
