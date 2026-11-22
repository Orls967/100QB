//! Deterministic `xoshiro256++` PRNG with `SplitMix64` seeding.

pub const MAIN_SEED: u64 = 20261122;

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
}
