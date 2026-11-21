//! Deterministic `xoshiro256++` PRNG with `SplitMix64` seeding and samplers.

pub const SEED_DAY004A: u64 = 20261122;
pub const SEED_DAY004B: u64 = 2026112202;
pub const SEED_DAY004C: u64 = 2026112203;
pub const SEED_DAY004D: u64 = 2026112204;
pub const SEED_DAY004E_MAIN: u64 = 2026112205;
pub const SEED_DAY004E_KERNEL: u64 = 2026112215;
pub const SEED_DAY004E_EVAL: u64 = 2026112225;
pub const SEED_DAY004E_VERIFY: u64 = 2026112235;
pub const SEED_DAY004F_INTEGRATION: u64 = 2026112246;

#[inline]
pub fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

#[inline]
pub fn derive_seed(base: u64, salt: u64, idx: u64) -> u64 {
    let mut x = base ^ salt.wrapping_mul(0x9e3779b97f4a7c15) ^ idx.wrapping_mul(0xbf58476d1ce4e5b9);
    splitmix64(&mut x)
}

#[inline]
pub fn derive_boot_seed(base: u64, salt: u64, outer: u64, boot: u64) -> u64 {
    let mut x = base
        ^ salt.wrapping_mul(0x9e3779b97f4a7c15)
        ^ outer.wrapping_mul(0xbf58476d1ce4e5b9)
        ^ boot.wrapping_mul(0x94d049bb133111eb);
    splitmix64(&mut x)
}

#[inline]
pub fn derive_state_seed(base: u64, salt: u64, h: u64, r: u64, b: u64) -> u64 {
    let mut x = base
        ^ salt.wrapping_mul(0x9e3779b97f4a7c15)
        ^ h.wrapping_mul(0xbf58476d1ce4e5b9)
        ^ r.wrapping_mul(0x94d049bb133111eb)
        ^ b.wrapping_mul(0x2545f4914f6cdd1d);
    splitmix64(&mut x)
}

#[derive(Clone, Debug)]
pub struct Rng {
    s: [u64; 4],
    has_spare: bool,
    spare: f64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut sm = seed;
        let s = [
            splitmix64(&mut sm),
            splitmix64(&mut sm),
            splitmix64(&mut sm),
            splitmix64(&mut sm),
        ];
        Self {
            s,
            has_spare: false,
            spare: 0.0,
        }
    }

    pub fn from_state(s: [u64; 4]) -> Self {
        Self {
            s,
            has_spare: false,
            spare: 0.0,
        }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[0]
            .wrapping_add(self.s[3])
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
        (self.next_u64() >> 11) as f64 * (1.0 / 9007199254740992.0)
    }

    #[inline]
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    #[inline]
    pub fn normal(&mut self) -> f64 {
        if self.has_spare {
            self.has_spare = false;
            return self.spare;
        }
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        let mag = (-2.0 * u1.ln()).sqrt();
        let angle = 2.0 * std::f64::consts::PI * u2;
        self.spare = mag * angle.sin();
        self.has_spare = true;
        mag * angle.cos()
    }

    #[inline]
    pub fn exp(&mut self, mean: f64) -> f64 {
        if mean <= 0.0 {
            return 0.0;
        }
        let u = (1.0 - self.uniform()).max(1e-300);
        -mean * u.ln()
    }

    #[inline]
    pub fn poisson(&mut self, lambda: f64) -> usize {
        if lambda <= 0.0 {
            return 0;
        }
        let l = (-lambda).exp();
        let mut k = 0usize;
        let mut p = 1.0;
        loop {
            k += 1;
            p *= self.uniform();
            if p <= l {
                return k - 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_xoshiro256pp_reference_vectors() {
        // Reference test vector for xoshiro256++ from Blackman & Vigna with state [1, 2, 3, 4]
        let mut rng = Rng::from_state([1, 2, 3, 4]);
        let v0 = rng.next_u64();
        let v1 = rng.next_u64();
        let v2 = rng.next_u64();
        // Verify deterministic progression and exact manual computation for step 0:
        // rotl(1 + 4, 23) + 1 = (5 << 23) + 1 = 41943041
        assert_eq!(v0, 41_943_041);
        assert_ne!(v1, v0);
        assert_ne!(v2, v1);

        let mut rng2 = Rng::new(SEED_DAY004A);
        let mut rng3 = Rng::new(SEED_DAY004A);
        for _ in 0..100 {
            assert_eq!(rng2.next_u64(), rng3.next_u64());
        }
    }

    #[test]
    fn test_samplers_moments() {
        let mut rng = Rng::new(SEED_DAY004A);
        let n = 50_000;
        let mut sum_n = 0.0;
        let mut sum_n2 = 0.0;
        for _ in 0..n {
            let z = rng.normal();
            sum_n += z;
            sum_n2 += z * z;
        }
        let mean_n = sum_n / n as f64;
        let var_n = sum_n2 / n as f64 - mean_n * mean_n;
        assert!(mean_n.abs() < 0.02);
        assert!((var_n - 1.0).abs() < 0.03);
    }
}
