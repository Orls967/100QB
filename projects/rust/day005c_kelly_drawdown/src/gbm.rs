/// Continuous-time GBM market parameters and Kelly log-growth formulas.
#[derive(Clone, Copy, Debug)]
pub struct GbmParams {
    pub mu: f64,
    pub sigma: f64,
    pub r: f64,
}

impl Default for GbmParams {
    fn default() -> Self {
        Self {
            mu: 0.08,
            sigma: 0.20,
            r: 0.03,
        }
    }
}

impl GbmParams {
    #[inline]
    pub fn kelly_fraction(&self) -> f64 {
        (self.mu - self.r) / (self.sigma * self.sigma)
    }

    #[inline]
    pub fn kelly_growth(&self) -> f64 {
        self.r + 0.5 * (self.mu - self.r).powi(2) / (self.sigma * self.sigma)
    }

    #[inline]
    pub fn drift_log(&self, f: f64) -> f64 {
        self.r + f * (self.mu - self.r) - 0.5 * f * f * self.sigma * self.sigma
    }

    #[inline]
    pub fn vol_log(&self, f: f64) -> f64 {
        f.abs() * self.sigma
    }
}
