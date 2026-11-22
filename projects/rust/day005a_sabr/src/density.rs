#![allow(clippy::too_many_arguments)]

use crate::hagan::hagan_call;

/// Breeden-Litzenberger implied density `q(K) = d^2 C / dK^2` from Hagan call prices
/// via central finite difference `(C(K+h) - 2 C(K) + C(K-h)) / h^2`.
#[inline]
pub fn hagan_density_fd(
    f: f64,
    k: f64,
    t: f64,
    alpha: f64,
    beta: f64,
    rho: f64,
    nu: f64,
    h: f64,
) -> f64 {
    if k <= h {
        return 0.0;
    }
    let cp = hagan_call(f, k + h, t, alpha, beta, rho, nu);
    let c0 = hagan_call(f, k, t, alpha, beta, rho, nu);
    let cm = hagan_call(f, k - h, t, alpha, beta, rho, nu);
    (cp - 2.0 * c0 + cm) / (h * h)
}

#[derive(Clone, Debug)]
pub struct DensityScanResult {
    pub ks: Vec<f64>,
    pub qs: Vec<f64>,
    pub k_neg_min: Option<f64>,
    pub k_neg_max: Option<f64>,
    pub min_q: f64,
    pub total_neg_mass: f64,
}

/// Scans `q_Hagan(K)` on `[k_min, k_max]` with `n_pts` points to measure negative-density
/// butterfly-arbitrage region and total negative probability mass `int q^-(K) dK`.
pub fn scan_hagan_negative_density(
    f: f64,
    t: f64,
    alpha: f64,
    beta: f64,
    rho: f64,
    nu: f64,
    k_min: f64,
    k_max: f64,
    n_pts: usize,
) -> DensityScanResult {
    let dk = (k_max - k_min) / (n_pts - 1) as f64;
    let h = (0.25 * dk).min(1e-5);
    let mut ks = Vec::with_capacity(n_pts);
    let mut qs = Vec::with_capacity(n_pts);
    let mut k_neg_min: Option<f64> = None;
    let mut k_neg_max: Option<f64> = None;
    let mut min_q = f64::INFINITY;
    let mut neg_mass = 0.0_f64;

    for i in 0..n_pts {
        let k = k_min + dk * i as f64;
        let q = hagan_density_fd(f, k, t, alpha, beta, rho, nu, h);
        ks.push(k);
        qs.push(q);
        if q < min_q {
            min_q = q;
        }
        if q < 0.0 {
            if k_neg_min.is_none() {
                k_neg_min = Some(k);
            }
            k_neg_max = Some(k);
            let w = if i == 0 || i + 1 == n_pts { 0.5 } else { 1.0 };
            neg_mass += w * (-q) * dk;
        }
    }

    DensityScanResult {
        ks,
        qs,
        k_neg_min,
        k_neg_max,
        min_q,
        total_neg_mass: neg_mass,
    }
}

/// Digital call `P(F_T > K) = -dC/dK` (approximated via central difference) and butterfly spread
/// `C(K - dK) - 2 C(K) + C(K + dK)` for Hagan.
pub fn hagan_digital_and_butterfly(
    f: f64,
    k: f64,
    dk: f64,
    t: f64,
    alpha: f64,
    beta: f64,
    rho: f64,
    nu: f64,
) -> (f64, f64) {
    let cp = hagan_call(f, k + dk, t, alpha, beta, rho, nu);
    let c0 = hagan_call(f, k, t, alpha, beta, rho, nu);
    let cm = hagan_call(f, k - dk, t, alpha, beta, rho, nu);
    let digital = (cm - cp) / (2.0 * dk);
    let butterfly = cp - 2.0 * c0 + cm;
    (digital, butterfly)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hagan_negative_density_at_long_t() {
        let scan =
            scan_hagan_negative_density(0.04, 10.0, 0.0325, 0.5, -0.3, 0.4, 0.001, 0.08, 400);
        assert!(
            scan.total_neg_mass > 0.0,
            "expected negative density at T=10"
        );
        assert!(scan.min_q < 0.0);
    }
}
