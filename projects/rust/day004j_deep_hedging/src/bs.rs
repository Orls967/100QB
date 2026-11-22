use std::f64::consts::{FRAC_2_SQRT_PI, PI, SQRT_2};

#[inline]
pub fn norm_pdf(x: f64) -> f64 {
    (-0.5 * x * x).exp() / (2.0 * PI).sqrt()
}

/// Standard normal CDF via Taylor + continued-fraction `erfc` (~1e-15 accuracy).
#[inline]
pub fn norm_cdf(x: f64) -> f64 {
    0.5 * erfc_approx(-x / SQRT_2)
}

#[inline]
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    if x > 27.0 {
        return 0.0;
    }
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

pub fn bs_call(s: f64, k: f64, r: f64, sigma: f64, tau: f64) -> f64 {
    if tau <= 1e-14 || sigma <= 1e-14 {
        return (s - k * (-r * tau).exp()).max(0.0);
    }
    let sqrt_t = tau.sqrt();
    let d1 = ((s / k).ln() + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqrt_t);
    let d2 = d1 - sigma * sqrt_t;
    s * norm_cdf(d1) - k * (-r * tau).exp() * norm_cdf(d2)
}

pub fn bs_delta(s: f64, k: f64, r: f64, sigma: f64, tau: f64) -> f64 {
    if tau <= 1e-14 || sigma <= 1e-14 {
        return if s > k { 1.0 } else { 0.0 };
    }
    let sqrt_t = tau.sqrt();
    let d1 = ((s / k).ln() + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqrt_t);
    norm_cdf(d1)
}

pub fn bs_gamma(s: f64, k: f64, r: f64, sigma: f64, tau: f64) -> f64 {
    if tau <= 1e-14 || sigma <= 1e-14 || s <= 0.0 {
        return 0.0;
    }
    let sqrt_t = tau.sqrt();
    let d1 = ((s / k).ln() + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqrt_t);
    norm_pdf(d1) / (s * sigma * sqrt_t)
}

pub fn bs_vega(s: f64, k: f64, r: f64, sigma: f64, tau: f64) -> f64 {
    if tau <= 1e-14 || sigma <= 1e-14 || s <= 0.0 {
        return 0.0;
    }
    let sqrt_t = tau.sqrt();
    let d1 = ((s / k).ln() + (r + 0.5 * sigma * sigma) * tau) / (sigma * sqrt_t);
    s * norm_pdf(d1) * sqrt_t
}

pub fn implied_vol_call(price: f64, s: f64, k: f64, r: f64, tau: f64) -> f64 {
    let intrinsic = (s - k * (-r * tau).exp()).max(0.0);
    if price <= intrinsic + 1e-12 {
        return 1e-6;
    }
    let mut lo = 1e-5_f64;
    let mut hi = 3.0_f64;
    let mut sig = 0.2_f64;
    for _ in 0..60 {
        let p = bs_call(s, k, r, sig, tau);
        let diff = p - price;
        if diff.abs() < 1e-12 {
            return sig;
        }
        if diff < 0.0 {
            lo = lo.max(sig);
        } else {
            hi = hi.min(sig);
        }
        let v = bs_vega(s, k, r, sig, tau);
        if v > 1e-12 {
            let next = sig - diff / v;
            if next > lo && next < hi {
                sig = next;
                continue;
            }
        }
        sig = 0.5 * (lo + hi);
    }
    sig
}

/// Leland (1985) adjusted volatility for one-way proportional transaction cost `c`
/// (`cost = c * |dDelta| * S`) and rebalancing step `dt`.
pub fn leland_sigma(sigma: f64, c: f64, dt: f64) -> f64 {
    let le = (8.0 / PI).sqrt() * c / (sigma * dt.sqrt());
    sigma * (1.0 + le).sqrt()
}

/// Whalley-Wilmott (1997) asymptotic no-trade half-width around BS delta:
/// `H = (3 * exp(-r * tau) * c * S * Gamma^2 / (2 * lambda))^(1/3)`.
pub fn whalley_wilmott_half_width(
    s: f64,
    k: f64,
    r: f64,
    sigma: f64,
    tau: f64,
    c: f64,
    lambda: f64,
) -> f64 {
    let gamma = bs_gamma(s, k, r, sigma, tau);
    let num = 3.0 * (-r * tau).exp() * c * s * gamma * gamma;
    (num / (2.0 * lambda)).max(0.0).cbrt()
}

/// Bertsimas-Kogan-Lo (2000) / Boyle-Emanuel (1980) asymptotic discrete-hedging SD:
/// `SD(PnL) approx sqrt( (T / (2*n)) * int_0^T E[(S_t^2 * Gamma_t * sigma^2)^2] dt )`.
pub fn discrete_hedging_sd_approx(
    s0: f64,
    k: f64,
    r: f64,
    sigma: f64,
    t_mat: f64,
    n_steps: usize,
) -> f64 {
    let nt = 40_usize;
    let nz = 61_usize;
    let dt_q = t_mat / nt as f64;
    let z_max = 5.5_f64;
    let dz = (2.0 * z_max) / (nz - 1) as f64;

    let mut integral_t = 0.0_f64;
    for it in 0..nt {
        let t = (it as f64 + 0.5) * dt_q;
        let tau = (t_mat - t).max(1e-6);
        let sqrt_t = t.sqrt();
        let mut exp_sq = 0.0_f64;
        for iz in 0..nz {
            let z = -z_max + iz as f64 * dz;
            let w = if iz == 0 || iz == nz - 1 { 0.5 } else { 1.0 } * dz * norm_pdf(z);
            let s = s0 * ((r - 0.5 * sigma * sigma) * t + sigma * sqrt_t * z).exp();
            let g = bs_gamma(s, k, r, sigma, tau);
            let term = s * s * g * sigma * sigma;
            exp_sq += w * term * term;
        }
        integral_t += exp_sq * dt_q;
    }
    ((t_mat / (2.0 * n_steps as f64)) * integral_t).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bs_and_iv_roundtrip() {
        let p = bs_call(1.0, 1.0, 0.0, 0.2, 0.25);
        assert!((p - 0.03987761167674492).abs() < 1e-12);
        let iv = implied_vol_call(p, 1.0, 1.0, 0.0, 0.25);
        assert!((iv - 0.2).abs() < 1e-10);
    }
}
