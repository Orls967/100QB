//! Pembangkit bilangan acak xoshiro256++ dengan seeding SplitMix64 dan Box-Muller normal.

/// Fungsi hash SplitMix64 untuk menurunkan seed per chunk secara deterministik.
pub fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// RNG xoshiro256++ + cache Box-Muller.
#[derive(Clone, Debug)]
pub struct Xoshiro256PlusPlus {
    s: [u64; 4],
    spare: Option<f64>,
}

impl Xoshiro256PlusPlus {
    /// Bangun RNG dari seed 64-bit.
    pub fn from_seed(seed: u64) -> Self {
        let mut x = seed;
        let mut s = [0u64; 4];
        for slot in &mut s {
            x = splitmix64(x);
            *slot = x;
        }
        Self { s, spare: None }
    }

    /// Ambil bilangan acak ``u64``.
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

    /// Ambil ``f64`` seragam di ``[0, 1)``.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    /// Ambil sampel normal standar ``N(0, 1)`` via Box-Muller.
    #[inline]
    pub fn next_normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        let mut u1 = self.next_f64();
        while u1 <= 1e-15 {
            u1 = self.next_f64();
        }
        let u2 = self.next_f64();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }
}
