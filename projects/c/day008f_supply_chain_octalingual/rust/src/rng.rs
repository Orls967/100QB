use std::cmp::Ordering;

#[derive(Clone, Copy)]
pub struct Xoshiro256 { pub s: [u64; 4] }

pub fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

pub fn fnv1a_bytes(buf: &[u8]) -> u64 {
    let mut h = 0xCBF29CE484222325u64;
    for &b in buf { h ^= b as u64; h = h.wrapping_mul(0x100000001B3); }
    h
}

pub fn fnv1a_u32pair(a: u32, b: u32) -> u64 {
    let mut buf = [0u8; 8];
    buf[0..4].copy_from_slice(&a.to_le_bytes());
    buf[4..8].copy_from_slice(&b.to_le_bytes());
    fnv1a_bytes(&buf)
}

pub fn fnv1a_i64(mut h: u64, v: i64) -> u64 {
    for &b in &(v as u64).to_le_bytes() { h ^= b as u64; h = h.wrapping_mul(0x100000001B3); }
    h
}

impl Xoshiro256 {
    pub fn seed(seed: u64) -> Self {
        let mut x = seed;
        let mut s = [0u64; 4];
        for i in 0..4 { s[i] = splitmix64(&mut x); }
        Self { s }
    }
    pub fn next(&mut self) -> u64 {
        let res = (self.s[0].wrapping_add(self.s[3])).rotate_left(23).wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0]; self.s[3] ^= self.s[1]; self.s[1] ^= self.s[2]; self.s[0] ^= self.s[3];
        self.s[2] ^= t; self.s[3] = self.s[3].rotate_left(45);
        res
    }
    pub fn rand_int(&mut self, lo: i64, hi: i64) -> i64 {
        if hi <= lo { return lo; }
        let span = (hi - lo + 1) as u64;
        lo + ((self.next() >> 11) % span) as i64
    }
}

#[derive(Clone, Copy, Eq, PartialEq)]
pub struct SimEvent {
    pub time: i64, pub seq: u64, pub etype: i32, pub src: i32, pub dst: i32, pub sku: i32, pub qty: i32, pub aux: i32,
}

impl Ord for SimEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        other.time.cmp(&self.time).then_with(|| other.seq.cmp(&self.seq))
    }
}
impl PartialOrd for SimEvent { fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) } }
