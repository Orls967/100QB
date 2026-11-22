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

/// Computes `n`-point Gauss-Legendre nodes and weights on `[-1, 1]`.
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut x = vec![0.0_f64; n];
    let mut w = vec![0.0_f64; n];
    let m = n.div_ceil(2);
    for i in 0..m {
        let mut z = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut pp = 0.0_f64;
        for _ in 0..40 {
            let mut p1 = 1.0_f64;
            let mut p2 = 0.0_f64;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp;
            if (z - z1).abs() < 1e-15 {
                break;
            }
        }
        x[i] = -z;
        x[n - 1 - i] = z;
        let wi = 2.0 / ((1.0 - z * z) * pp * pp);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
    (x, w)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gauss_legendre_and_norm_cdf() {
        let (x, w) = gauss_legendre(32);
        let int_x4: f64 = x
            .iter()
            .zip(w.iter())
            .map(|(&xi, &wi)| wi * xi.powi(4))
            .sum();
        assert!((int_x4 - 0.4).abs() < 1e-14);
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_cdf(1.959963984540054) - 0.975).abs() < 1e-13);
    }
}
