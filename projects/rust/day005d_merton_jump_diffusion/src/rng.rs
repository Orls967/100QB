/// Deterministic PRNG: SplitMix64 seeder + xoshiro256++ + Box-Muller + Poisson + Gamma + ncx2.
pub struct Rng {
    s: [u64; 4],
    spare: Option<f64>,
}

#[inline]
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[inline]
pub fn derive_seed(main_seed: u64, stream_tag: u64, idx: u64) -> u64 {
    let x = main_seed
        ^ stream_tag.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ idx.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    splitmix64(x)
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let s0 = splitmix64(seed);
        let s1 = splitmix64(s0);
        let s2 = splitmix64(s1);
        let s3 = splitmix64(s2);
        Self {
            s: [s0, s1, s2, s3],
            spare: None,
        }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = (self.s[0].wrapping_add(self.s[3]))
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    #[inline]
    pub fn uniform(&mut self) -> f64 {
        let u = (self.next_u64() >> 11) as f64 * (1.0 / ((1u64 << 53) as f64));
        if u <= 0.0 {
            1e-16
        } else if u >= 1.0 {
            1.0 - 1e-16
        } else {
            u
        }
    }

    #[inline]
    pub fn normal_pair(&mut self) -> (f64, f64) {
        let u1 = self.uniform();
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        (r * theta.cos(), r * theta.sin())
    }

    #[inline]
    pub fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let (z0, z1) = self.normal_pair();
        self.spare = Some(z1);
        z0
    }

    /// Exact Poisson(lambda) via Knuth's product method (chunked for large lambda).
    pub fn poisson(&mut self, lambda: f64) -> usize {
        if lambda <= 0.0 {
            return 0;
        }
        if lambda <= 30.0 {
            let l = (-lambda).exp();
            let mut k = 0_usize;
            let mut p = 1.0_f64;
            loop {
                p *= self.uniform();
                if p <= l {
                    return k;
                }
                k += 1;
            }
        } else {
            let n_sub = (lambda / 15.0).ceil() as usize;
            let sub_lam = lambda / n_sub as f64;
            let mut tot = 0_usize;
            for _ in 0..n_sub {
                tot += self.poisson(sub_lam);
            }
            tot
        }
    }

    /// Exact Gamma(shape, scale) via Marsaglia & Tsang (2000).
    pub fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        if shape <= 0.0 || scale <= 0.0 {
            return 0.0;
        }
        if shape < 1.0 {
            let g1 = self.gamma(shape + 1.0, 1.0);
            let u = self.uniform();
            return scale * g1 * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v_base = 1.0 + c * x;
            if v_base <= 0.0 {
                continue;
            }
            let v = v_base * v_base * v_base;
            let u = self.uniform();
            let x2 = x * x;
            if u < 1.0 - 0.0331 * x2 * x2 {
                return scale * d * v;
            }
            if u.ln() < 0.5 * x2 + d * (1.0 - v + v.ln()) {
                return scale * d * v;
            }
        }
    }

    /// Exact non-central chi-squared `chi'^2(df, ncp)` via Poisson-Gamma mixture (`df >= 0`).
    pub fn ncx2(&mut self, df: f64, ncp: f64) -> f64 {
        let k = self.poisson(0.5 * ncp.max(0.0));
        let shape = 0.5 * df.max(0.0) + k as f64;
        if shape <= 0.0 {
            0.0
        } else {
            self.gamma(shape, 2.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamma_and_ncx2_moments() {
        let mut rg = Rng::new(20261123);
        let n = 100_000_usize;
        let mut sum_g = 0.0;
        let mut sum_x = 0.0;
        let df = 2.5;
        let ncp = 4.0;
        for _ in 0..n {
            sum_g += rg.gamma(0.4, 2.0);
            sum_x += rg.ncx2(df, ncp);
        }
        let mg = sum_g / n as f64;
        let mx = sum_x / n as f64;
        assert!((mg - 0.8).abs() < 0.02, "gamma mean={mg}");
        assert!((mx - (df + ncp)).abs() < 0.05, "ncx2 mean={mx}");
    }
}
