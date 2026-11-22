/// Deterministic PRNG: SplitMix64 seeder + xoshiro256++ + Box-Muller + Poisson.
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

    /// Exact Poisson(lambda) via Knuth's product method (for lambda <= 30) or normal-approx/Ahrens-Dieter fallback.
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
            // Sum of smaller chunks for exactness
            let n_sub = (lambda / 15.0).ceil() as usize;
            let sub_lam = lambda / n_sub as f64;
            let mut tot = 0_usize;
            for _ in 0..n_sub {
                tot += self.poisson(sub_lam);
            }
            tot
        }
    }
}
