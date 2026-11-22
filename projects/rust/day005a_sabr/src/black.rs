use crate::special::{norm_cdf, norm_pdf};

/// Black-76 forward call option price (`r = 0`).
#[inline]
pub fn black_call(f: f64, k: f64, vol: f64, t: f64) -> f64 {
    if t <= 0.0 || vol <= 0.0 {
        return (f - k).max(0.0);
    }
    if k <= 0.0 {
        return f;
    }
    if f <= 0.0 {
        return 0.0;
    }
    let std_dev = vol * t.sqrt();
    let d1 = ((f / k).ln() + 0.5 * std_dev * std_dev) / std_dev;
    let d2 = d1 - std_dev;
    (f * norm_cdf(d1) - k * norm_cdf(d2)).max(0.0)
}

/// Black-76 forward put option price (`r = 0`).
#[inline]
pub fn black_put(f: f64, k: f64, vol: f64, t: f64) -> f64 {
    black_call(f, k, vol, t) - (f - k)
}

/// Black-76 forward call delta `dC/dF`.
#[inline]
pub fn black_delta(f: f64, k: f64, vol: f64, t: f64) -> f64 {
    if t <= 0.0 || vol <= 0.0 {
        return if f > k { 1.0 } else { 0.0 };
    }
    if k <= 0.0 {
        return 1.0;
    }
    if f <= 0.0 {
        return 0.0;
    }
    let std_dev = vol * t.sqrt();
    let d1 = ((f / k).ln() + 0.5 * std_dev * std_dev) / std_dev;
    norm_cdf(d1)
}

/// Black-76 forward vega `dC/d(vol)`.
#[inline]
pub fn black_vega(f: f64, k: f64, vol: f64, t: f64) -> f64 {
    if t <= 0.0 || vol <= 0.0 || f <= 0.0 || k <= 0.0 {
        return 0.0;
    }
    let sqrt_t = t.sqrt();
    let std_dev = vol * sqrt_t;
    let d1 = ((f / k).ln() + 0.5 * std_dev * std_dev) / std_dev;
    f * sqrt_t * norm_pdf(d1)
}

/// Implied Black-76 volatility from a forward call price via safeguarded Newton + Brent.
pub fn black_implied_vol(price: f64, f: f64, k: f64, t: f64) -> Option<f64> {
    if t <= 0.0 || f <= 0.0 || k <= 0.0 {
        return None;
    }
    let intrinsic = (f - k).max(0.0);
    if price <= intrinsic + 1e-15 || price >= f - 1e-15 {
        return None;
    }
    let mut lo = 1e-6_f64;
    let mut hi = 5.0_f64;
    while black_call(f, k, hi, t) < price && hi < 50.0 {
        hi *= 2.0;
    }
    if black_call(f, k, hi, t) < price {
        return None;
    }
    let mut vol = 0.2_f64.clamp(lo, hi);
    for _ in 0..100 {
        let c = black_call(f, k, vol, t);
        let diff = c - price;
        if diff.abs() < 1e-15 * f.max(1e-8) {
            return Some(vol);
        }
        if diff < 0.0 {
            lo = vol;
        } else {
            hi = vol;
        }
        let vega = black_vega(f, k, vol, t);
        let next = if vega > 1e-14 * f {
            let step = diff / vega;
            if step.abs() < 1e-12 * (1.0 + vol) && diff.abs() < 1e-12 * f {
                return Some(vol - step);
            }
            vol - step
        } else {
            0.5 * (lo + hi)
        };
        vol = if next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
        if (hi - lo) < 1e-12 * (1.0 + vol) {
            return Some(vol);
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_black_iv_roundtrip() {
        for &f in &[0.04, 0.02, 0.08] {
            for &k in &[0.01, 0.04, 0.08] {
                for &vol in &[0.15, 0.35, 0.8] {
                    let c = black_call(f, k, vol, 5.0);
                    if let Some(iv) = black_implied_vol(c, f, k, 5.0) {
                        assert!((iv - vol).abs() < 1e-6, "vol={vol} iv={iv}");
                    }
                }
            }
        }
    }
}
