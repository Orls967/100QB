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

/// Log-Gamma function via Lanczos approximation (`g = 7`, `n = 9`).
pub fn lgamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        PI.ln() - (PI * z).sin().abs().ln() - lgamma(1.0 - z)
    } else {
        let z1 = z - 1.0;
        let mut x = COEFFS[0];
        for (i, &c) in COEFFS.iter().enumerate().skip(1) {
            x += c / (z1 + i as f64);
        }
        let t = z1 + 7.5;
        0.5 * (2.0 * PI).ln() + (z1 + 0.5) * t.ln() - t + x.ln()
    }
}

/// Regularized lower incomplete gamma `P(a, x)` (`a > 0, x >= 0`).
pub fn gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_p_series(a, x)
    } else {
        1.0 - gamma_q_cf(a, x)
    }
}

/// Regularized upper incomplete gamma `Q(a, x) = 1 - P(a, x)` (`a > 0, x >= 0`).
pub fn gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        1.0 - gamma_p_series(a, x)
    } else {
        gamma_q_cf(a, x)
    }
}

fn gamma_p_series(a: f64, x: f64) -> f64 {
    let gln = lgamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..400 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    (sum * (-x + a * x.ln() - gln).exp()).clamp(0.0, 1.0)
}

fn gamma_q_cf(a: f64, x: f64) -> f64 {
    let gln = lgamma(a);
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / 1e-30;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..400 {
        let an = -(i as f64) * (i as f64 - a);
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
    ((-x + a * x.ln() - gln).exp() * h).clamp(0.0, 1.0)
}

/// Non-central chi-squared CDF `P(chi'^2(df, ncp) <= x)` via mode-centered Poisson-Gamma series.
/// Supports `df >= 0.0` (at `df == 0.0`, has point mass `exp(-ncp/2)` at `x = 0`).
pub fn ncx2_cdf(x: f64, df: f64, ncp: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    }
    let lam = 0.5 * ncp.max(0.0);
    if x == 0.0 {
        return if df <= 0.0 { (-lam).exp() } else { 0.0 };
    }
    if lam == 0.0 {
        return if df <= 0.0 {
            1.0
        } else {
            gamma_p(0.5 * df, 0.5 * x)
        };
    }
    let k0 = lam.floor() as usize;
    let log_w0 = -lam + (k0 as f64) * lam.ln() - lgamma(k0 as f64 + 1.0);
    let w0 = log_w0.exp();
    let half_x = 0.5 * x;
    let half_df = 0.5 * df.max(0.0);

    // Recurrence: P(a+1, z) = P(a, z) - t_a, where t_a = z^a e^{-z} / Gamma(a+1)
    // and t_{a+1} = t_a * z / (a + 1).
    let a0 = half_df + k0 as f64;
    let p0 = if a0 <= 0.0 { 1.0 } else { gamma_p(a0, half_x) };
    let log_hx = half_x.ln();
    let mut t_up = (-half_x + a0 * log_hx - lgamma(a0 + 1.0)).exp();
    let mut p_up = p0;
    let mut sum = w0 * p0;
    let mut w_up = w0;
    for j in 1..1000 {
        let k = k0 + j;
        w_up *= lam / k as f64;
        p_up = (p_up - t_up).clamp(0.0, 1.0);
        let a_cur = half_df + k as f64;
        t_up *= half_x / a_cur;
        if w_up < 1e-17 * (sum.max(1e-30)) {
            break;
        }
        sum += w_up * p_up;
    }
    let mut w_dn = w0;
    let mut p_dn = p0;
    let mut t_dn = t_up; // not needed, recompute from a0
    let _ = &mut t_dn;
    let mut t_cur = (-half_x + a0 * log_hx - lgamma(a0 + 1.0)).exp();
    for k in (1..=k0).rev() {
        let a_prev = half_df + (k - 1) as f64;
        t_cur *= (a_prev + 1.0) / half_x;
        p_dn = if a_prev <= 0.0 {
            1.0
        } else {
            (p_dn + t_cur).clamp(0.0, 1.0)
        };
        w_dn *= (k as f64) / lam;
        if w_dn < 1e-17 * (sum.max(1e-30)) {
            break;
        }
        sum += w_dn * p_dn;
    }
    sum.clamp(0.0, 1.0)
}

/// Non-central chi-squared survival function `P(chi'^2(df, ncp) > x)`.
pub fn ncx2_sf(x: f64, df: f64, ncp: f64) -> f64 {
    if x < 0.0 {
        return 1.0;
    }
    let lam = 0.5 * ncp.max(0.0);
    if x == 0.0 {
        return if df <= 0.0 {
            (1.0 - (-lam).exp()).clamp(0.0, 1.0)
        } else {
            1.0
        };
    }
    if lam == 0.0 {
        return if df <= 0.0 {
            0.0
        } else {
            gamma_q(0.5 * df, 0.5 * x)
        };
    }
    let k0 = lam.floor() as usize;
    let log_w0 = -lam + (k0 as f64) * lam.ln() - lgamma(k0 as f64 + 1.0);
    let w0 = log_w0.exp();
    let half_x = 0.5 * x;
    let half_df = 0.5 * df.max(0.0);

    // Recurrence: Q(a+1, z) = Q(a, z) + t_a, where t_a = z^a e^{-z} / Gamma(a+1)
    let a0 = half_df + k0 as f64;
    let q0 = if a0 <= 0.0 { 0.0 } else { gamma_q(a0, half_x) };
    let log_hx = half_x.ln();
    let mut t_up = (-half_x + a0 * log_hx - lgamma(a0 + 1.0)).exp();
    let mut q_up = q0;
    let mut sum = w0 * q0;
    let mut w_up = w0;
    for j in 1..1000 {
        let k = k0 + j;
        w_up *= lam / k as f64;
        q_up = (q_up + t_up).clamp(0.0, 1.0);
        let a_cur = half_df + k as f64;
        t_up *= half_x / a_cur;
        if w_up < 1e-17 * (sum.max(1e-30)) {
            break;
        }
        sum += w_up * q_up;
    }
    let mut w_dn = w0;
    let mut q_dn = q0;
    let mut t_cur = (-half_x + a0 * log_hx - lgamma(a0 + 1.0)).exp();
    for k in (1..=k0).rev() {
        let a_prev = half_df + (k - 1) as f64;
        t_cur *= (a_prev + 1.0) / half_x;
        q_dn = if a_prev <= 0.0 {
            0.0
        } else {
            (q_dn - t_cur).clamp(0.0, 1.0)
        };
        w_dn *= (k as f64) / lam;
        if w_dn < 1e-17 * (sum.max(1e-30)) {
            break;
        }
        sum += w_dn * q_dn;
    }
    sum.clamp(0.0, 1.0)
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
    fn test_gauss_legendre_norm_and_ncx2() {
        let (x, w) = gauss_legendre(32);
        let int_x4: f64 = x
            .iter()
            .zip(w.iter())
            .map(|(&xi, &wi)| wi * xi.powi(4))
            .sum();
        assert!((int_x4 - 0.4).abs() < 1e-14);
        assert!((norm_cdf(0.0) - 0.5).abs() < 1e-15);
        assert!((norm_cdf(1.959963984540054) - 0.975).abs() < 1e-13);

        // ncx2_cdf + ncx2_sf = 1
        for &(xv, df, ncp) in &[(0.5, 0.0, 4.0), (5.0, 2.0, 3.5), (120.0, 1.5, 100.0)] {
            let c = ncx2_cdf(xv, df, ncp);
            let s = ncx2_sf(xv, df, ncp);
            assert!((c + s - 1.0).abs() < 1e-13, "c={c} s={s}");
        }
        // Central chi2(df=2) = Exp(1/2), CDF(x) = 1 - exp(-x/2)
        assert!((ncx2_cdf(3.0, 2.0, 0.0) - (1.0 - (-1.5_f64).exp())).abs() < 1e-14);
        // Zero df at x=0 has atom exp(-ncp/2)
        assert!((ncx2_cdf(0.0, 0.0, 4.0) - (-2.0_f64).exp()).abs() < 1e-15);
    }
}
