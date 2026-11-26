pub const BP: i64 = 10000;

pub fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

pub fn fnv1a_bytes(buf: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in buf { h = (h ^ b as u64).wrapping_mul(0x100000001b3); }
    h
}

pub fn fnv1a_u32pair(a: u32, b: u32) -> u64 {
    let mut buf = [0u8; 8];
    buf[0..4].copy_from_slice(&a.to_le_bytes());
    buf[4..8].copy_from_slice(&b.to_le_bytes());
    fnv1a_bytes(&buf)
}

pub fn fnv1a_i64(mut h: u64, v: i64) -> u64 {
    for b in (v as u64).to_le_bytes() { h = (h ^ b as u64).wrapping_mul(0x100000001b3); }
    h
}

#[derive(Clone, Default)]
pub struct Xoshiro256 { pub s: [u64; 4] }

impl Xoshiro256 {
    pub fn seed(sd: u64) -> Self {
        let mut x = sd; let mut s = [0u64; 4];
        for i in 0..4 { s[i] = splitmix64(&mut x); }
        Self { s }
    }
    pub fn next_u64(&mut self) -> u64 {
        let res = (self.s[0].wrapping_add(self.s[3])).rotate_left(23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1]; self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].rotate_left(45);
        res
    }
    pub fn rand_int(&mut self, lo: i64, hi: i64) -> i64 {
        if hi <= lo { return lo; }
        lo + ((self.next_u64() >> 11) % ((hi - lo + 1) as u64)) as i64
    }
    pub fn bernoulli(&mut self, p_bp: i64) -> i32 {
        if p_bp <= 0 { return 0; }
        if p_bp >= BP { return 1; }
        (((self.next_u64() >> 11) % (BP as u64)) < (p_bp as u64)) as i32
    }
}
