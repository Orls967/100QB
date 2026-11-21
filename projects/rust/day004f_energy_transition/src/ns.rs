//! Nelson-Siegel and Nelson-Siegel-Svensson yield curves, small-tau Taylor stability, and GL integral check.

use crate::linalg;

#[derive(Clone, Copy, Debug)]
pub struct NsParams {
    pub l: f64,
    pub s: f64,
    pub c: f64,
    pub lambda: f64,
    pub c2: f64,
    pub lambda2: f64,
}

impl NsParams {
    pub fn default_ns() -> Self {
        Self {
            l: 0.045,
            s: -0.025,
            c: 0.015,
            lambda: 0.6,
            c2: 0.0,
            lambda2: 0.0,
        }
    }

    pub fn default_nss() -> Self {
        let mut p = Self::default_ns();
        p.c2 = 0.01;
        p.lambda2 = 0.15;
        p
    }

    pub fn yield_at(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return self.l + self.s;
        }
        let mut y = self.l + self.s * g1(tau, self.lambda) + self.c * g2(tau, self.lambda);
        if self.c2 != 0.0 {
            y += self.c2 * g2(tau, self.lambda2);
        }
        y
    }

    pub fn yield_naive(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return self.l + self.s;
        }
        let mut y =
            self.l + self.s * g1_naive(tau, self.lambda) + self.c * g2_naive(tau, self.lambda);
        if self.c2 != 0.0 {
            y += self.c2 * g2_naive(tau, self.lambda2);
        }
        y
    }

    pub fn forward(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return self.l + self.s;
        }
        let u = self.lambda * tau;
        let eu = (-u).exp();
        let mut f = self.l + self.s * eu + self.c * u * eu;
        if self.c2 != 0.0 {
            let u2 = self.lambda2 * tau;
            f += self.c2 * u2 * (-u2).exp();
        }
        f
    }

    pub fn discount(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return 1.0;
        }
        (-tau * self.yield_at(tau)).exp()
    }

    pub fn forward_integral_avg(&self, tau: f64, quad_order: usize) -> f64 {
        if tau <= 0.0 {
            return self.l + self.s;
        }
        let nodes = linalg::gauss_legendre(quad_order, 0.0, tau);
        let sum: f64 = nodes.iter().map(|nd| nd.w * self.forward(nd.x)).sum();
        sum / tau
    }
}

#[inline]
pub fn g1(tau: f64, lambda: f64) -> f64 {
    if tau <= 0.0 || lambda == 0.0 {
        return 1.0;
    }
    let u = lambda * tau;
    if u.abs() < 1e-3 {
        1.0 + u
            * (-0.5
                + u * (1.0 / 6.0
                    + u * (-1.0 / 24.0
                        + u * (1.0 / 120.0
                            + u * (-1.0 / 720.0 + u * (1.0 / 5040.0 - u / 40320.0))))))
    } else {
        -((-u).exp_m1()) / u
    }
}

#[inline]
pub fn g2(tau: f64, lambda: f64) -> f64 {
    if tau <= 0.0 || lambda == 0.0 {
        return 0.0;
    }
    let u = lambda * tau;
    if u.abs() < 1e-3 {
        u * (0.5
            + u * (-1.0 / 3.0
                + u * (1.0 / 8.0
                    + u * (-1.0 / 30.0 + u * (1.0 / 144.0 + u * (-1.0 / 840.0 + u / 5760.0))))))
    } else {
        g1(tau, lambda) - (-u).exp()
    }
}

#[inline]
pub fn g1_naive(tau: f64, lambda: f64) -> f64 {
    let u = lambda * tau;
    if u == 0.0 {
        1.0
    } else {
        (1.0 - (-u).exp()) / u
    }
}

#[inline]
pub fn g2_naive(tau: f64, lambda: f64) -> f64 {
    let u = lambda * tau;
    if u == 0.0 {
        0.0
    } else {
        (1.0 - (-u).exp()) / u - (-u).exp()
    }
}

#[derive(Clone, Debug)]
pub struct SmallTauRow {
    pub tau: f64,
    pub g1_stable: f64,
    pub g1_naive: f64,
    pub g1_rel_err: f64,
    pub g2_stable: f64,
    pub g2_naive: f64,
    pub g2_rel_err: f64,
    pub yield_stable: f64,
    pub yield_naive: f64,
    pub yield_abs_err: f64,
}

pub fn compare_small_tau_stability(p: &NsParams, taus: &[f64]) -> Vec<SmallTauRow> {
    taus.iter()
        .map(|&tau| {
            let g1s = g1(tau, p.lambda);
            let g1n = g1_naive(tau, p.lambda);
            let g2s = g2(tau, p.lambda);
            let g2n = g2_naive(tau, p.lambda);
            let ys = p.yield_at(tau);
            let yn = p.yield_naive(tau);
            SmallTauRow {
                tau,
                g1_stable: g1s,
                g1_naive: g1n,
                g1_rel_err: (g1n - g1s).abs() / g1s.abs().max(1e-30),
                g2_stable: g2s,
                g2_naive: g2n,
                g2_rel_err: (g2n - g2s).abs() / g2s.abs().max(1e-30),
                yield_stable: ys,
                yield_naive: yn,
                yield_abs_err: (yn - ys).abs(),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ns_forward_integral_and_limits() {
        let ns = NsParams::default_ns();
        let nss = NsParams::default_nss();
        for &tau in &[0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0] {
            let y_ns = ns.yield_at(tau);
            let avg_ns = ns.forward_integral_avg(tau, 32);
            assert!((y_ns - avg_ns).abs() < 1e-13);
            let y_nss = nss.yield_at(tau);
            let avg_nss = nss.forward_integral_avg(tau, 32);
            assert!((y_nss - avg_nss).abs() < 1e-13);
        }
        assert!((ns.yield_at(1e-10) - (ns.l + ns.s)).abs() < 1e-8);
        assert!((ns.yield_at(1e6) - ns.l).abs() < 1e-6);
    }
}
