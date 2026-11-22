use crate::parallel::reduce_tuple4_fixed_blocks;
use crate::rng::{derive_seed, Rng};

#[derive(Clone, Copy, Debug)]
pub struct ConcRange {
    pub pa: f64,
    pub pb: f64,
    pub l: f64,
}

impl ConcRange {
    /// Creates a concentrated range `[pa, pb]` normalized so initial value `V(p0 = 1.0) = 1.0`.
    pub fn normalized(pa: f64, pb: f64) -> Self {
        let unit_v0 = 2.0 - pa.sqrt() - 1.0 / pb.sqrt();
        Self {
            pa,
            pb,
            l: 1.0 / unit_v0,
        }
    }

    pub fn full_range() -> Self {
        Self {
            pa: 0.0,
            pb: f64::INFINITY,
            l: 0.5,
        }
    }

    #[inline]
    pub fn capital_efficiency(&self) -> f64 {
        self.l / 0.5
    }

    #[inline]
    pub fn reserves(&self, p: f64) -> (f64, f64) {
        if p <= self.pa {
            let x = if self.pa > 0.0 {
                self.l * (1.0 / self.pa.sqrt() - 1.0 / self.pb.sqrt())
            } else {
                f64::INFINITY
            };
            (x, 0.0)
        } else if p >= self.pb {
            let y = self.l * (self.pb.sqrt() - self.pa.sqrt());
            (0.0, y)
        } else {
            let sp = p.sqrt();
            (
                self.l * (1.0 / sp - 1.0 / self.pb.sqrt()),
                self.l * (sp - self.pa.sqrt()),
            )
        }
    }

    #[inline]
    pub fn value(&self, p: f64) -> f64 {
        let (x, y) = self.reserves(p);
        x * p + y
    }

    #[inline]
    pub fn gamma(&self, p: f64) -> f64 {
        if p > self.pa && p < self.pb {
            -0.5 * self.l * p.powf(-1.5)
        } else {
            0.0
        }
    }

    #[inline]
    pub fn lvr_rate(&self, p: f64, sigma: f64) -> f64 {
        if p > self.pa && p < self.pb {
            0.25 * sigma * sigma * self.l * p.sqrt()
        } else {
            0.0
        }
    }
}

/// Closed-form expected first exit time from `x0 = 0` to `{a = ln(pa), b = ln(pb)}`
/// for `dX = nu dt + sigma dW` (`nu = -0.5 * sigma^2`), derived from `0.5*sigma^2 u'' + nu u' = -1`.
pub fn expected_exit_time_gbm(pa: f64, pb: f64, sigma: f64) -> f64 {
    if pa <= 0.0 || !pb.is_finite() {
        return f64::INFINITY;
    }
    let a = pa.ln();
    let b = pb.ln();
    let nu = -0.5 * sigma * sigma;
    let theta = -2.0 * nu / (sigma * sigma); // = +1.0 for mu=0!
                                             // Homogeneous solutions: 1 and exp(theta * x). Particular: -x / nu.
                                             // u(x) = C1 + C2 * exp(theta * x) - x / nu, with u(a)=0, u(b)=0
    let ea = (theta * a).exp();
    let eb = (theta * b).exp();
    let c2 = ((b - a) / nu) / (eb - ea);
    let c1 = a / nu - c2 * ea;
    c1 + c2
}

/// MC verification of concentrated LVR rate, time in range, and exit time.
pub fn mc_conc_metrics(
    range: ConcRange,
    sigma: f64,
    t_days: f64,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let t_yr = t_days / 365.0;
    let n_steps = (t_days * 48.0).round().max(100.0) as usize; // 30-min steps
    let dt = t_yr / n_steps as f64;
    let drift = -0.5 * sigma * sigma * dt;
    let vol_dt = sigma * dt.sqrt();

    let sums = reduce_tuple4_fixed_blocks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(
            seed,
            0x4001 ^ ((range.pa * 100.0) as u64),
            r as u64,
        ));
        let mut p = 1.0_f64;
        let mut v_rebal = range.value(1.0);
        let mut lvr_int = 0.0_f64;
        let mut in_range_cnt = 0_usize;
        for _ in 0..n_steps {
            if p > range.pa && p < range.pb {
                in_range_cnt += 1;
            }
            let d = range.reserves(p).0;
            lvr_int += range.lvr_rate(p, sigma) * dt;
            let p_next = p * (drift + vol_dt * rg.normal()).exp();
            v_rebal += d * (p_next - p);
            p = p_next;
        }
        let lvr_rep = (v_rebal - range.value(p)) / t_yr;
        let lvr_a = lvr_int / t_yr;
        let frac_in = in_range_cnt as f64 / n_steps as f64;
        [lvr_rep, lvr_a, frac_in, lvr_rep * lvr_rep]
    });
    let nf = n_paths as f64;
    let m_rep = sums[0] / nf;
    let m_a = sums[1] / nf;
    let m_in = sums[2] / nf;
    let se_rep = ((sums[3] / nf - m_rep * m_rep).max(0.0) / nf).sqrt();
    (m_rep, se_rep, m_a, m_in)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conc_gamma_and_exit_time() {
        let r10 = ConcRange::normalized(0.9, 1.1);
        assert!((r10.value(1.0) - 1.0).abs() < 1e-14);
        assert!((r10.capital_efficiency() - r10.l / 0.5).abs() < 1e-14);
        let et = expected_exit_time_gbm(0.9, 1.1, 0.8);
        assert!(et > 0.0 && et < 0.1, "et={et}");
    }
}
