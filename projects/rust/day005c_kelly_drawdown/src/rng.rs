/// Deterministic PRNG: SplitMix64 seeder + xoshiro256++ + normal, Bernoulli, Gamma & Student-t.
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
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p
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

    /// Gamma(shape, scale) via Marsaglia & Tsang (2000).
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
            let v = 1.0 + c * x;
            if v <= 0.0 {
                continue;
            }
            let v3 = v * v * v;
            let u = self.uniform();
            if u < 1.0 - 0.0331 * x * x * x * x {
                return scale * d * v3;
            }
            if u.ln() < 0.5 * x * x + d * (1.0 - v3 + v3.ln()) {
                return scale * d * v3;
            }
        }
    }

    /// Student-t(nu) via Z / sqrt(V / nu), V ~ chi2(nu) = Gamma(nu/2, 2).
    #[inline]
    pub fn student_t(&mut self, nu: f64) -> f64 {
        let z = self.normal();
        let v = self.gamma(0.5 * nu, 2.0);
        z / (v / nu).sqrt()
    }
}
