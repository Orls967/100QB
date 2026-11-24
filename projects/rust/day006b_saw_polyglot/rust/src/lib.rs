#![forbid(unsafe_code)]

const OFF: i32 = 32;
const DIRS: [(i32, i32, u8); 4] = [(1, 0, b'E'), (0, 1, b'N'), (-1, 0, b'W'), (0, -1, b'S')];
pub type Walk = Vec<(i32, i32)>;

pub struct Xoshiro256Pp {
    s: [u64; 4],
}

impl Xoshiro256Pp {
    pub fn from_seed(seed: u64) -> Self {
        let mut z = seed;
        let mut sm = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Self {
            s: [sm(), sm(), sm(), sm()],
        }
    }

    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let res = self.s[0]
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
        res
    }

    #[inline]
    pub fn bounded(&mut self, bound: u64) -> u64 {
        let mut m = (self.next_u64() as u128) * (bound as u128);
        let mut l = m as u64;
        if l < bound {
            let thr = bound.wrapping_neg() % bound;
            while l < thr {
                m = (self.next_u64() as u128) * (bound as u128);
                l = m as u64;
            }
        }
        (m >> 64) as u64
    }
}

#[derive(Clone)]
struct Bitmap([u64; 64]);

impl Bitmap {
    fn new() -> Self {
        Self([0; 64])
    }
    #[inline]
    fn has(&self, x: i32, y: i32) -> bool {
        (self.0[(y + OFF) as usize] & (1u64 << (x + OFF))) != 0
    }
    #[inline]
    fn set(&mut self, x: i32, y: i32) {
        self.0[(y + OFF) as usize] |= 1u64 << (x + OFF);
    }
    #[inline]
    fn clr(&mut self, x: i32, y: i32) {
        self.0[(y + OFF) as usize] &= !(1u64 << (x + OFF));
    }
}

#[inline]
pub fn apply_d4(sym: u8, dx: i32, dy: i32) -> (i32, i32) {
    match sym & 7 {
        0 => (dx, dy),
        1 => (-dy, dx),
        2 => (-dx, -dy),
        3 => (dy, -dx),
        4 => (dx, -dy),
        5 => (-dx, dy),
        6 => (dy, dx),
        _ => (-dy, -dx),
    }
}

struct SymState {
    nmax: usize,
    hm_n: Option<usize>,
    vis: Bitmap,
    counts: Vec<u64>,
    hm: Option<Vec<Vec<u64>>>,
}

impl SymState {
    fn dfs(&mut self, x: i32, y: i32, step: usize, turned: bool) {
        self.counts[step] += if turned { 8 } else { 4 };
        if self.hm_n == Some(step) {
            if let Some(g) = self.hm.as_mut() {
                let h = step as i32;
                if !turned {
                    for &(sx, sy) in &[(x, 0), (-x, 0), (0, x), (0, -x)] {
                        g[(sy + h) as usize][(sx + h) as usize] += 1;
                    }
                } else {
                    for s in 0..8u8 {
                        let (sx, sy) = apply_d4(s, x, y);
                        g[(sy + h) as usize][(sx + h) as usize] += 1;
                    }
                }
            }
        }
        if step == self.nmax {
            return;
        }
        if !turned {
            self.vis.set(x + 1, 0);
            self.dfs(x + 1, 0, step + 1, false);
            self.vis.clr(x + 1, 0);
            self.vis.set(x, 1);
            self.dfs(x, 1, step + 1, true);
            self.vis.clr(x, 1);
            return;
        }
        if step + 1 == self.nmax && self.hm_n != Some(self.nmax) {
            let free = DIRS
                .iter()
                .filter(|&&(dx, dy, _)| !self.vis.has(x + dx, y + dy))
                .count() as u64;
            self.counts[self.nmax] += 8 * free;
            return;
        }
        for &(dx, dy, _) in &DIRS {
            let (nx, ny) = (x + dx, y + dy);
            if !self.vis.has(nx, ny) {
                self.vis.set(nx, ny);
                self.dfs(nx, ny, step + 1, true);
                self.vis.clr(nx, ny);
            }
        }
    }
}

pub fn enumerate_symmetry(nmax: usize, hm_n: Option<usize>) -> (Vec<u64>, Option<Vec<Vec<u64>>>) {
    let mut st = SymState {
        nmax,
        hm_n,
        vis: Bitmap::new(),
        counts: vec![0u64; nmax + 1],
        hm: hm_n.map(|h| vec![vec![0u64; 2 * h + 1]; 2 * h + 1]),
    };
    st.counts[0] = 1;
    if nmax > 0 {
        st.vis.set(0, 0);
        st.vis.set(1, 0);
        st.dfs(1, 0, 1, false);
    }
    (st.counts, st.hm)
}

pub fn enumerate_from_prefix(steps: &str, weight: u64, nmax: usize) -> Result<Vec<u64>, String> {
    let mut vis = Bitmap::new();
    let (mut x, mut y) = (0i32, 0i32);
    vis.set(0, 0);
    for b in steps.bytes() {
        let (dx, dy) = match b {
            b'E' => (1, 0),
            b'N' => (0, 1),
            b'W' => (-1, 0),
            b'S' => (0, -1),
            _ => return Err("bad step".into()),
        };
        x += dx;
        y += dy;
        if vis.has(x, y) {
            return Err("collision".into());
        }
        vis.set(x, y);
    }
    let mut counts = vec![0u64; nmax + 1];
    fn dfs(vis: &mut Bitmap, counts: &mut [u64], x: i32, y: i32, step: usize, w: u64, nmax: usize) {
        counts[step] += w;
        if step == nmax {
            return;
        }
        if step + 1 == nmax {
            let free = DIRS
                .iter()
                .filter(|&&(dx, dy, _)| !vis.has(x + dx, y + dy))
                .count() as u64;
            counts[nmax] += w * free;
            return;
        }
        for &(dx, dy, _) in &DIRS {
            let (nx, ny) = (x + dx, y + dy);
            if !vis.has(nx, ny) {
                vis.set(nx, ny);
                dfs(vis, counts, nx, ny, step + 1, w, nmax);
                vis.clr(nx, ny);
            }
        }
    }
    if steps.len() <= nmax {
        dfs(&mut vis, &mut counts, x, y, steps.len(), weight, nmax);
    }
    Ok(counts)
}

struct StampSet {
    mask: usize,
    stamp: u32,
    stamps: Vec<u32>,
    keys: Vec<u64>,
}

impl StampSet {
    fn new(cap: usize) -> Self {
        let c = cap.next_power_of_two().max(64);
        Self {
            mask: c - 1,
            stamp: 1,
            stamps: vec![0; c],
            keys: vec![0; c],
        }
    }
    #[inline]
    fn clear(&mut self) {
        self.stamp = self.stamp.wrapping_add(1);
        if self.stamp == 0 {
            self.stamps.fill(0);
            self.stamp = 1;
        }
    }
    #[inline]
    fn pack(x: i32, y: i32) -> u64 {
        ((x as u32 as u64) << 32) | (y as u32 as u64)
    }
    #[inline]
    fn ins(&mut self, x: i32, y: i32) -> bool {
        let k = Self::pack(x, y);
        let mut i = ((k.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) as usize) & self.mask;
        loop {
            if self.stamps[i] != self.stamp {
                self.stamps[i] = self.stamp;
                self.keys[i] = k;
                return true;
            }
            if self.keys[i] == k {
                return false;
            }
            i = (i + 1) & self.mask;
        }
    }
    #[inline]
    fn has(&self, x: i32, y: i32) -> bool {
        let k = Self::pack(x, y);
        let mut i = ((k.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 32) as usize) & self.mask;
        loop {
            if self.stamps[i] != self.stamp {
                return false;
            }
            if self.keys[i] == k {
                return true;
            }
            i = (i + 1) & self.mask;
        }
    }
}

pub struct PivotCfg {
    pub n: usize,
    pub seed: u64,
    pub burn_in: usize,
    pub batches: usize,
    pub batch_size: usize,
    pub samples: usize,
    pub stride: usize,
}

fn pivot_step(
    pts: &mut [(i32, i32)],
    cand: &mut [(i32, i32)],
    set: &mut StampSet,
    rng: &mut Xoshiro256Pp,
) -> bool {
    let n = pts.len() - 1;
    let k = rng.bounded(n as u64) as usize;
    let sym = (rng.bounded(7) as u8) + 1;
    if k == 0 {
        for p in pts.iter_mut().skip(1) {
            *p = apply_d4(sym, p.0, p.1);
        }
        return true;
    }
    let (xk, yk) = pts[k];
    set.clear();
    for &(px, py) in &pts[..k] {
        set.ins(px, py);
    }
    for i in (k + 1)..=n {
        let (rx, ry) = apply_d4(sym, pts[i].0 - xk, pts[i].1 - yk);
        let (nx, ny) = (xk + rx, yk + ry);
        if set.has(nx, ny) {
            return false;
        }
        cand[i] = (nx, ny);
    }
    pts[(k + 1)..=n].copy_from_slice(&cand[(k + 1)..=n]);
    true
}

pub fn run_pivot(cfg: &PivotCfg) -> (f64, f64, f64, Vec<Walk>) {
    let mut rng = Xoshiro256Pp::from_seed(cfg.seed);
    let (mut x, mut y) = (0i32, 0i32);
    let mut pts = vec![(0, 0)];
    for i in 0..cfg.n {
        if i.is_multiple_of(2) {
            x += 1;
        } else {
            y += 1;
        }
        pts.push((x, y));
    }
    let mut cand = vec![(0, 0); cfg.n + 1];
    let mut set = StampSet::new((cfg.n + 1) * 4);
    for _ in 0..cfg.burn_in {
        pivot_step(&mut pts, &mut cand, &mut set, &mut rng);
    }
    let (mut acc, mut cnt) = (0usize, 0usize);
    let mut b_means = Vec::with_capacity(cfg.batches);
    let mut walks = Vec::with_capacity(cfg.samples);
    let stride = cfg.stride.max(1);
    for _ in 0..cfg.batches {
        let mut s_r2 = 0.0;
        for _ in 0..cfg.batch_size {
            if pivot_step(&mut pts, &mut cand, &mut set, &mut rng) {
                acc += 1;
            }
            let (ex, ey) = pts[cfg.n];
            s_r2 += (ex as f64).powi(2) + (ey as f64).powi(2);
            cnt += 1;
            if walks.len() < cfg.samples && cnt.is_multiple_of(stride) {
                walks.push(pts.clone());
            }
        }
        b_means.push(s_r2 / (cfg.batch_size as f64));
    }
    while walks.len() < cfg.samples {
        pivot_step(&mut pts, &mut cand, &mut set, &mut rng);
        cnt += 1;
        if cnt.is_multiple_of(stride) {
            walks.push(pts.clone());
        }
    }
    let bf = cfg.batches as f64;
    let mean = b_means.iter().sum::<f64>() / bf;
    let var = b_means.iter().map(|&m| (m - mean).powi(2)).sum::<f64>() / (bf - 1.0);
    let rate = (acc as f64) / ((cfg.batches * cfg.batch_size) as f64);
    (mean, (var / bf).sqrt(), rate, walks)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_cn_and_symmetry_proof() {
        let (sym, _) = enumerate_symmetry(9, None);
        assert_eq!(&sym[1..=9], &[4, 12, 36, 100, 284, 780, 2172, 5916, 16268]);
        let mut full = [0u64; 9];
        for r in ["E", "N", "W", "S"] {
            let c = enumerate_from_prefix(r, 1, 8).expect("valid");
            for n in 1..=8 {
                full[n] += c[n];
            }
        }
        assert_eq!(&sym[1..=8], &full[1..=8]);
    }

    #[test]
    fn test_prefix_and_heatmap() {
        assert!(enumerate_from_prefix("EW", 4, 6).is_err());
        let (c, hm) = enumerate_symmetry(8, Some(6));
        assert_eq!(
            hm.expect("hm").iter().flat_map(|r| r.iter()).sum::<u64>(),
            c[6]
        );
    }

    #[test]
    fn test_pivot_validity() {
        let cfg = PivotCfg {
            n: 40,
            seed: 20261124,
            burn_in: 300,
            batches: 4,
            batch_size: 200,
            samples: 8,
            stride: 25,
        };
        let (_, _, _, walks) = run_pivot(&cfg);
        for w in &walks {
            let mut set = StampSet::new(256);
            assert_eq!(w.len(), 41);
            for &(x, y) in w {
                assert!(set.ins(x, y));
            }
        }
    }

    #[test]
    fn test_pivot_seed_determinism() {
        let cfg = PivotCfg {
            n: 30,
            seed: 77,
            burn_in: 200,
            batches: 3,
            batch_size: 100,
            samples: 2,
            stride: 20,
        };
        assert_eq!(run_pivot(&cfg).0.to_bits(), run_pivot(&cfg).0.to_bits());
    }
}
