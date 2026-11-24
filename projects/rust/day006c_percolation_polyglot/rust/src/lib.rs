#![forbid(unsafe_code)]

pub struct Xoshiro256Pp {
    s: [u64; 4],
}

fn splitmix64(x: &mut u64) -> u64 {
    *x = x.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *x;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

impl Xoshiro256Pp {
    pub fn from_seed(seed: u64) -> Self {
        let mut sm = seed;
        Self {
            s: [
                splitmix64(&mut sm),
                splitmix64(&mut sm),
                splitmix64(&mut sm),
                splitmix64(&mut sm),
            ],
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let res = (self.s[0].wrapping_add(self.s[3]))
            .rotate_left(23)
            .wrapping_add(self.s[0]);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        res
    }

    pub fn bounded(&mut self, bound: u32) -> u32 {
        let b = bound as u64;
        let mut x = self.next_u64();
        let mut m = (x as u128) * (b as u128);
        let mut l = m as u64;
        if l < b {
            let t = b.wrapping_neg() % b;
            while l < t {
                x = self.next_u64();
                m = (x as u128) * (b as u128);
                l = m as u64;
            }
        }
        (m >> 64) as u32
    }

    pub fn next_f64(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }

    pub fn shuffle<T>(&mut self, slice: &mut [T]) {
        let mut i = slice.len();
        while i > 1 {
            let j = self.bounded(i as u32) as usize;
            i -= 1;
            slice.swap(i, j);
        }
    }
}

pub struct Uf {
    parent: Vec<u32>,
    sz: Vec<u32>,
}

impl Uf {
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            sz: vec![1; n],
        }
    }

    pub fn reset(&mut self) {
        for (i, p) in self.parent.iter_mut().enumerate() {
            *p = i as u32;
        }
        self.sz.fill(1);
    }

    pub fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let gp = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = gp;
            x = gp;
        }
        x
    }

    pub fn union(&mut self, a: u32, b: u32) -> u32 {
        let mut ra = self.find(a);
        let mut rb = self.find(b);
        if ra == rb {
            return self.sz[ra as usize];
        }
        if self.sz[ra as usize] < self.sz[rb as usize] {
            std::mem::swap(&mut ra, &mut rb);
        }
        self.parent[rb as usize] = ra;
        self.sz[ra as usize] += self.sz[rb as usize];
        self.sz[ra as usize]
    }
}

pub fn build_rect_bonds(rows: usize, cols: usize) -> Vec<(u32, u32)> {
    let mut bonds = Vec::with_capacity(rows * (cols - 1) + (rows - 1) * cols);
    for r in 0..rows {
        for c in 0..(cols - 1) {
            let u = (r * cols + c) as u32;
            bonds.push((u, u + 1));
        }
    }
    for r in 0..(rows - 1) {
        for c in 0..cols {
            let u = (r * cols + c) as u32;
            bonds.push((u, u + cols as u32));
        }
    }
    bonds
}

pub struct JobResult {
    pub hist: Vec<u64>,
    pub sum_largest: u64,
    pub sum_largest_sq: u64,
}

pub fn run_newman_ziff(l: usize, seed: u64, samples: usize) -> JobResult {
    let (rows, cols) = (l, l + 1);
    let nv = rows * cols;
    let bonds = build_rect_bonds(rows, cols);
    let nb = bonds.len();
    let n_half = nb / 2;
    let (v_left, v_right) = (nv as u32, (nv + 1) as u32);

    let mut rng = Xoshiro256Pp::from_seed(seed);
    let mut perm: Vec<u32> = (0..nb as u32).collect();
    let mut uf_cross = Uf::new(nv + 2);
    let mut uf_phys = Uf::new(nv);
    let mut hist = vec![0u64; nb + 1];
    let mut sum_largest = 0u64;
    let mut sum_largest_sq = 0u64;

    for _ in 0..samples {
        rng.shuffle(&mut perm);
        uf_cross.reset();
        uf_phys.reset();
        for r in 0..rows {
            uf_cross.union(v_left, (r * cols) as u32);
            uf_cross.union(v_right, (r * cols + cols - 1) as u32);
        }
        let mut max_sz = 1u32;
        let mut n_star = nb;
        let mut crossed = false;
        for (idx, &bi) in perm.iter().enumerate() {
            let step = idx + 1;
            let (u, v) = bonds[bi as usize];
            if step <= n_half {
                let s = uf_phys.union(u, v);
                if s > max_sz {
                    max_sz = s;
                }
            }
            if !crossed {
                uf_cross.union(u, v);
                if uf_cross.find(v_left) == uf_cross.find(v_right) {
                    crossed = true;
                    n_star = step;
                }
            }
            if step == n_half {
                let m = max_sz as u64;
                sum_largest += m;
                sum_largest_sq += m * m;
            }
            if crossed && step >= n_half {
                break;
            }
        }
        hist[n_star] += 1;
    }
    JobResult {
        hist,
        sum_largest,
        sum_largest_sq,
    }
}

pub fn sample_bonds(l: usize, p: f64, seed: u64) -> Vec<(u32, u32)> {
    let bonds = build_rect_bonds(l, l + 1);
    let mut rng = Xoshiro256Pp::from_seed(seed);
    bonds.into_iter().filter(|_| rng.next_f64() < p).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_union_find() {
        let mut uf = Uf::new(6);
        assert_eq!(uf.union(0, 1), 2);
        assert_eq!(uf.union(2, 3), 2);
        assert_eq!(uf.union(1, 2), 4);
        assert_eq!(uf.find(0), uf.find(3));
        assert_ne!(uf.find(0), uf.find(4));
        assert_eq!(uf.union(0, 3), 4);
    }

    #[test]
    fn test_permutation_uniform() {
        let mut rng = Xoshiro256Pp::from_seed(42);
        let mut counts = [0i64; 6];
        for _ in 0..60_000 {
            let mut a = [0u8, 1, 2];
            rng.shuffle(&mut a);
            let idx = match a {
                [0, 1, 2] => 0,
                [0, 2, 1] => 1,
                [1, 0, 2] => 2,
                [1, 2, 0] => 3,
                [2, 0, 1] => 4,
                _ => 5,
            };
            counts[idx] += 1;
        }
        for c in counts {
            assert!((c - 10_000).abs() < 450, "non-uniform count {c}");
        }
    }

    #[test]
    fn test_crossing_2x3_matches_exact() {
        let bonds = build_rect_bonds(2, 3);
        assert_eq!(bonds.len(), 7);
        let mut exact_a = [0u32; 8];
        for mask in 0..(1u32 << 7) {
            let mut uf = Uf::new(8);
            uf.union(6, 0);
            uf.union(6, 3);
            uf.union(7, 2);
            uf.union(7, 5);
            for (i, &(u, v)) in bonds.iter().enumerate() {
                if (mask >> i) & 1 == 1 {
                    uf.union(u, v);
                }
            }
            if uf.find(6) == uf.find(7) {
                exact_a[mask.count_ones() as usize] += 1;
            }
        }
        let total_cross: u32 = exact_a.iter().sum();
        assert_eq!(total_cross * 2, 128);
        let res = run_newman_ziff(2, 20261124, 40_000);
        let binom = [1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0];
        let mut cum = 0u64;
        let mut p_half = 0.0;
        for n in 0..=7 {
            cum += res.hist[n];
            let qn = (cum as f64) / 40_000.0;
            p_half += binom[n] * qn / 128.0;
        }
        assert!((p_half - 0.5).abs() < 0.01, "p_half={p_half}");
    }

    #[test]
    fn test_seed_determinism() {
        let a = run_newman_ziff(8, 20261124, 500);
        let b = run_newman_ziff(8, 20261124, 500);
        let c = run_newman_ziff(8, 20261125, 500);
        assert_eq!(a.hist, b.hist);
        assert_eq!(a.sum_largest, b.sum_largest);
        assert_ne!(a.hist, c.hist);
    }
}
