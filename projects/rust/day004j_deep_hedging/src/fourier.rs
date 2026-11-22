use std::f64::consts::PI;

#[derive(Clone, Copy, Debug)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    #[inline]
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    #[inline]
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self::new(self.re - o.re, self.im - o.im)
    }
    #[inline]
    fn mul(self, o: Self) -> Self {
        Self::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
    #[inline]
    fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Self::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
    #[inline]
    fn exp(self) -> Self {
        let r = self.re.exp();
        Self::new(r * self.im.cos(), r * self.im.sin())
    }
    #[inline]
    fn ln(self) -> Self {
        let r = (self.re * self.re + self.im * self.im).sqrt().ln();
        let th = self.im.atan2(self.re);
        Self::new(r, th)
    }
    #[inline]
    fn sqrt(self) -> Self {
        let r = (self.re * self.re + self.im * self.im).sqrt();
        let th = self.im.atan2(self.re) * 0.5;
        let sr = r.sqrt();
        Self::new(sr * th.cos(), sr * th.sin())
    }
}

/// Normalized forward Heston characteristic function `E[exp(i * z * ln(S_T / (S0 * exp(r*T))))]`
/// using the Albrecher et al. (2007) "Little Heston Trap" formulation (branch-cut free).
#[allow(clippy::too_many_arguments)]
fn heston_cf_norm(
    z: Complex,
    v0: f64,
    kappa: f64,
    theta: f64,
    xi: f64,
    rho: f64,
    tau: f64,
) -> Complex {
    let i = Complex::new(0.0, 1.0);
    let xi2 = xi * xi;
    let b = Complex::new(kappa, 0.0).sub(i.mul(z).mul(Complex::new(rho * xi, 0.0)));
    let u_term = z.mul(z).add(i.mul(z)).mul(Complex::new(xi2, 0.0));
    let d = b.mul(b).add(u_term).sqrt();

    // Little Heston Trap: g = (b - d) / (b + d) so Re(d) >= 0 => |exp(-d*tau)| <= 1
    let g = b.sub(d).div(b.add(d));
    let exp_dt = d.mul(Complex::new(-tau, 0.0)).exp();
    let one = Complex::new(1.0, 0.0);
    let one_m_g_exp = one.sub(g.mul(exp_dt));
    let one_m_g = one.sub(g);

    let b_m_d = b.sub(d);
    let c_term = b_m_d
        .mul(Complex::new(tau, 0.0))
        .sub(Complex::new(2.0, 0.0).mul(one_m_g_exp.div(one_m_g).ln()))
        .mul(Complex::new(kappa * theta / xi2, 0.0));
    let d_term = b_m_d
        .div(Complex::new(xi2, 0.0))
        .mul(one.sub(exp_dt).div(one_m_g_exp))
        .mul(Complex::new(v0, 0.0));

    c_term.add(d_term).exp()
}

/// Heston European call price via Lewis (2001) single-integral formula (Simpson quadrature).
#[allow(clippy::too_many_arguments)]
pub fn heston_call_fourier(
    s: f64,
    k: f64,
    r: f64,
    v0: f64,
    kappa: f64,
    theta: f64,
    xi: f64,
    rho: f64,
    tau: f64,
) -> f64 {
    if tau <= 1e-12 {
        return (s - k).max(0.0);
    }
    let log_fk = (s / k).ln() + r * tau;
    let u_max = 100.0_f64;
    let n = 1000_usize; // even for Simpson (h = 0.1)
    let h = u_max / n as f64;
    let mut sum = 0.0_f64;
    for j in 0..=n {
        let u = j as f64 * h;
        let z = Complex::new(u, -0.5);
        let cf = heston_cf_norm(z, v0, kappa, theta, xi, rho, tau);
        let e_iu = Complex::new(0.0, u * log_fk).exp();
        let integ = e_iu.mul(cf).re / (u * u + 0.25);
        let w = if j == 0 || j == n {
            1.0
        } else if j % 2 == 1 {
            4.0
        } else {
            2.0
        };
        sum += w * integ;
    }
    let integral = sum * h / 3.0;
    s - (s * k).sqrt() * (-r * tau).exp() * integral / PI
}

/// Heston Fourier delta via central finite difference in spot `S`.
#[allow(clippy::too_many_arguments)]
pub fn heston_delta_fourier(
    s: f64,
    k: f64,
    r: f64,
    v: f64,
    kappa: f64,
    theta: f64,
    xi: f64,
    rho: f64,
    tau: f64,
) -> f64 {
    let ds = 1e-4 * s;
    let cp = heston_call_fourier(s + ds, k, r, v, kappa, theta, xi, rho, tau);
    let cm = heston_call_fourier(s - ds, k, r, v, kappa, theta, xi, rho, tau);
    (cp - cm) / (2.0 * ds)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bs::bs_call;

    #[test]
    fn test_heston_fourier_matches_bs_in_zero_vol_of_vol_limit() {
        let c_h = heston_call_fourier(1.0, 1.0, 0.0, 0.04, 2.0, 0.04, 1e-4, -0.7, 0.25);
        let c_bs = bs_call(1.0, 1.0, 0.0, 0.2, 0.25);
        assert!((c_h - c_bs).abs() < 1e-6, "c_h={c_h} c_bs={c_bs}");
    }

    #[test]
    fn test_heston_delta_smooth_across_h() {
        let d_a = heston_delta_fourier(1.0, 1.0, 0.0, 0.04, 2.0, 0.04, 0.5, -0.7, 0.25);
        let s = 1.0_f64;
        let ds2 = 1e-5 * s;
        let cp = heston_call_fourier(s + ds2, 1.0, 0.0, 0.04, 2.0, 0.04, 0.5, -0.7, 0.25);
        let cm = heston_call_fourier(s - ds2, 1.0, 0.0, 0.04, 2.0, 0.04, 0.5, -0.7, 0.25);
        let d_b = (cp - cm) / (2.0 * ds2);
        assert!((d_a - d_b).abs() < 1e-7, "d_a={d_a} d_b={d_b}");
    }
}
