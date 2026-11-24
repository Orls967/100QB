#![forbid(unsafe_code)]

//! Self-Avoiding Walk (SAW) core on the 2D square lattice:
//! - Exact DFS enumeration with D4 (x8) symmetry reduction and prefix-based worker mode.
//! - Madras-Sokal pivot MCMC sampler using custom xoshiro256++ RNG and batch means.

const BOARD_OFFSET: i32 = 32;
const BOARD_SIDE: usize = 64;
const DIRS: [(i32, i32, u8); 4] = [(1, 0, b'E'), (0, 1, b'N'), (-1, 0, b'W'), (0, -1, b'S')];

/// SplitMix64 seed expander used to initialize xoshiro256++ state.
#[derive(Clone, Debug)]
pub struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Custom xoshiro256++ 1.0 PRNG (Blackman & Vigna).
#[derive(Clone, Debug)]
pub struct Xoshiro256PlusPlus {
    s: [u64; 4],
}

impl Xoshiro256PlusPlus {
    pub fn from_seed(seed: u64) -> Self {
        let mut sm = SplitMix64::new(seed);
        Self {
            s: [sm.next_u64(), sm.next_u64(), sm.next_u64(), sm.next_u64()],
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

    /// Unbiased bounded integer in `0..bound` via Lemire's rejection method.
    #[inline]
    pub fn next_bounded(&mut self, bound: u64) -> u64 {
        debug_assert!(bound > 0);
        let mut x = self.next_u64();
        let mut m = (x as u128) * (bound as u128);
        let mut l = m as u64;
        if l < bound {
            let threshold = bound.wrapping_neg() % bound;
            while l < threshold {
                x = self.next_u64();
                m = (x as u128) * (bound as u128);
                l = m as u64;
            }
        }
        (m >> 64) as u64
    }
}

/// Compact 64x64 bitmap for fast DFS visited tracking around origin.
#[derive(Clone)]
struct VisitedBitmap {
    rows: [u64; BOARD_SIDE],
}

impl VisitedBitmap {
    fn new() -> Self {
        Self {
            rows: [0u64; BOARD_SIDE],
        }
    }

    #[inline]
    fn is_visited(&self, x: i32, y: i32) -> bool {
        let ux = (x + BOARD_OFFSET) as usize;
        let uy = (y + BOARD_OFFSET) as usize;
        (self.rows[uy] & (1u64 << ux)) != 0
    }

    #[inline]
    fn set(&mut self, x: i32, y: i32) {
        let ux = (x + BOARD_OFFSET) as usize;
        let uy = (y + BOARD_OFFSET) as usize;
        self.rows[uy] |= 1u64 << ux;
    }

    #[inline]
    fn unset(&mut self, x: i32, y: i32) {
        let ux = (x + BOARD_OFFSET) as usize;
        let uy = (y + BOARD_OFFSET) as usize;
        self.rows[uy] &= !(1u64 << ux);
    }
}

/// Result of exact SAW enumeration.
#[derive(Clone, Debug)]
pub struct EnumResult {
    pub nmax: usize,
    pub counts: Vec<u64>,
    pub heatmap_n: Option<usize>,
    pub heatmap_grid: Option<Vec<Vec<u64>>>,
}

/// Apply one of the 8 symmetries of the dihedral group D4 on Z^2.
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

fn record_heatmap(grid: &mut [Vec<u64>], hn: i32, x: i32, y: i32, has_turned: bool) {
    if !has_turned {
        for &(sx, sy) in &[(x, 0), (-x, 0), (0, x), (0, -x)] {
            grid[(sy + hn) as usize][(sx + hn) as usize] += 1;
        }
    } else {
        for sym in 0..8u8 {
            let (sx, sy) = apply_d4(sym, x, y);
            grid[(sy + hn) as usize][(sx + hn) as usize] += 1;
        }
    }
}

struct SymDfsCtx<'a> {
    nmax: usize,
    heatmap_n: Option<usize>,
    visited: VisitedBitmap,
    counts: &'a mut [u64],
    heatmap: Option<&'a mut [Vec<u64>]>,
}

fn dfs_sym(ctx: &mut SymDfsCtx<'_>, x: i32, y: i32, step: usize, has_turned: bool) {
    let w = if has_turned { 8u64 } else { 4u64 };
    ctx.counts[step] += w;

    if let Some(hn) = ctx.heatmap_n {
        if step == hn {
            if let Some(grid) = ctx.heatmap.as_deref_mut() {
                record_heatmap(grid, hn as i32, x, y, has_turned);
            }
        }
    }

    if step == ctx.nmax {
        return;
    }

    if !has_turned {
        // Straight East continuation: weight remains 4.
        let nx_east = x + 1;
        ctx.visited.set(nx_east, 0);
        dfs_sym(ctx, nx_east, 0, step + 1, false);
        ctx.visited.unset(nx_east, 0);

        // First non-straight turn forced North (+y): breaks reflection symmetry, weight becomes 8.
        ctx.visited.set(x, 1);
        dfs_sym(ctx, x, 1, step + 1, true);
        ctx.visited.unset(x, 1);
        return;
    }

    // Leaf-counting shortcut when on penultimate step and heatmap is not requested at nmax.
    if step + 1 == ctx.nmax && ctx.heatmap_n != Some(ctx.nmax) {
        let mut free = 0u64;
        for &(dx, dy, _) in &DIRS {
            if !ctx.visited.is_visited(x + dx, y + dy) {
                free += 1;
            }
        }
        ctx.counts[ctx.nmax] += 8 * free;
        return;
    }

    for &(dx, dy, _) in &DIRS {
        let nx = x + dx;
        let ny = y + dy;
        if !ctx.visited.is_visited(nx, ny) {
            ctx.visited.set(nx, ny);
            dfs_sym(ctx, nx, ny, step + 1, true);
            ctx.visited.unset(nx, ny);
        }
    }
}

/// Exact enumeration of c_n for n = 1..=nmax using x8 D4 symmetry reduction.
pub fn enumerate_symmetry(nmax: usize, heatmap_n: Option<usize>) -> EnumResult {
    assert!(nmax <= 30, "nmax must be <= 30 for 64x64 bitmap");
    let mut counts = vec![0u64; nmax + 1];
    counts[0] = 1;

    let mut heatmap_grid = heatmap_n.map(|hn| {
        let side = 2 * hn + 1;
        vec![vec![0u64; side]; side]
    });

    if nmax >= 1 {
        let mut visited = VisitedBitmap::new();
        visited.set(0, 0);
        visited.set(1, 0);
        let mut ctx = SymDfsCtx {
            nmax,
            heatmap_n,
            visited,
            counts: &mut counts,
            heatmap: heatmap_grid.as_deref_mut(),
        };
        dfs_sym(&mut ctx, 1, 0, 1, false);
    }

    EnumResult {
        nmax,
        counts,
        heatmap_n,
        heatmap_grid,
    }
}

/// Naive 4-direction brute-force enumeration without symmetry (used to prove x8 reduction in tests).
pub fn enumerate_naive(nmax: usize) -> Vec<u64> {
    assert!(nmax <= 16, "naive enumeration is only intended for small n");
    let mut counts = vec![0u64; nmax + 1];
    counts[0] = 1;
    let mut visited = VisitedBitmap::new();
    visited.set(0, 0);

    fn dfs(
        visited: &mut VisitedBitmap,
        counts: &mut [u64],
        x: i32,
        y: i32,
        step: usize,
        nmax: usize,
    ) {
        if step > 0 {
            counts[step] += 1;
        }
        if step == nmax {
            return;
        }
        for &(dx, dy, _) in &DIRS {
            let nx = x + dx;
            let ny = y + dy;
            if !visited.is_visited(nx, ny) {
                visited.set(nx, ny);
                dfs(visited, counts, nx, ny, step + 1, nmax);
                visited.unset(nx, ny);
            }
        }
    }

    dfs(&mut visited, &mut counts, 0, 0, 0, nmax);
    counts
}

fn dfs_prefix(
    visited: &mut VisitedBitmap,
    counts: &mut [u64],
    x: i32,
    y: i32,
    step: usize,
    weight: u64,
    nmax: usize,
) {
    counts[step] += weight;
    if step == nmax {
        return;
    }
    if step + 1 == nmax {
        let mut free = 0u64;
        for &(dx, dy, _) in &DIRS {
            if !visited.is_visited(x + dx, y + dy) {
                free += 1;
            }
        }
        counts[nmax] += weight * free;
        return;
    }
    for &(dx, dy, _) in &DIRS {
        let nx = x + dx;
        let ny = y + dy;
        if !visited.is_visited(nx, ny) {
            visited.set(nx, ny);
            dfs_prefix(visited, counts, nx, ny, step + 1, weight, nmax);
            visited.unset(nx, ny);
        }
    }
}

/// Enumerate all extensions up to `nmax` from an explicit step-string prefix (`E`, `N`, `W`, `S`).
pub fn enumerate_from_prefix(steps: &str, weight: u64, nmax: usize) -> Result<Vec<u64>, String> {
    if nmax > 30 {
        return Err("nmax exceeds bitmap capacity 30".to_string());
    }
    let mut visited = VisitedBitmap::new();
    let mut x = 0i32;
    let mut y = 0i32;
    visited.set(0, 0);

    for ch in steps.bytes() {
        let (dx, dy) = match ch {
            b'E' => (1, 0),
            b'N' => (0, 1),
            b'W' => (-1, 0),
            b'S' => (0, -1),
            _ => return Err(format!("invalid step byte: {}", ch as char)),
        };
        x += dx;
        y += dy;
        if !(-30..=30).contains(&x) || !(-30..=30).contains(&y) {
            return Err("prefix steps out of range".to_string());
        }
        if visited.is_visited(x, y) {
            return Err("prefix is not self-avoiding".to_string());
        }
        visited.set(x, y);
    }

    let k = steps.len();
    let mut counts = vec![0u64; nmax + 1];
    if k <= nmax {
        dfs_prefix(&mut visited, &mut counts, x, y, k, weight, nmax);
    }
    Ok(counts)
}

/// Generate all valid SAW prefixes of length `k` with first step forced East (`'E'`), weight 4.
pub fn generate_east_prefixes(k: usize) -> Vec<String> {
    if k == 0 {
        return vec![String::new()];
    }
    let mut out = Vec::new();
    let mut visited = VisitedBitmap::new();
    visited.set(0, 0);
    visited.set(1, 0);
    let mut buf = vec![b'E'];

    fn build(
        visited: &mut VisitedBitmap,
        buf: &mut Vec<u8>,
        out: &mut Vec<String>,
        x: i32,
        y: i32,
        target: usize,
    ) {
        if buf.len() == target {
            out.push(String::from_utf8(buf.clone()).expect("ascii steps"));
            return;
        }
        for &(dx, dy, code) in &DIRS {
            let nx = x + dx;
            let ny = y + dy;
            if !visited.is_visited(nx, ny) {
                visited.set(nx, ny);
                buf.push(code);
                build(visited, buf, out, nx, ny, target);
                buf.pop();
                visited.unset(nx, ny);
            }
        }
    }

    build(&mut visited, &mut buf, &mut out, 1, 0, k);
    out
}

/// Check that a walk of `n + 1` lattice coordinates starts at `(0,0)`, takes unit steps, and has 0 collisions.
pub fn is_valid_saw(pts: &[(i32, i32)], expected_n: usize) -> bool {
    if pts.len() != expected_n + 1 || pts[0] != (0, 0) {
        return false;
    }
    for w in pts.windows(2) {
        let manhattan = (w[1].0 - w[0].0).abs() + (w[1].1 - w[0].1).abs();
        if manhattan != 1 {
            return false;
        }
    }
    let mut table = StampHashSet::with_capacity(pts.len() * 4);
    table.clear();
    for &(x, y) in pts {
        if !table.insert(x, y) {
            return false;
        }
    }
    true
}

/// Allocation-free stamped open-addressing hash set for O(N) pivot collision detection.
struct StampHashSet {
    mask: usize,
    current_stamp: u32,
    stamps: Vec<u32>,
    keys: Vec<u64>,
}

impl StampHashSet {
    fn with_capacity(min_cap: usize) -> Self {
        let cap = min_cap.next_power_of_two().max(64);
        Self {
            mask: cap - 1,
            current_stamp: 1,
            stamps: vec![0; cap],
            keys: vec![0; cap],
        }
    }

    #[inline]
    fn clear(&mut self) {
        self.current_stamp = self.current_stamp.wrapping_add(1);
        if self.current_stamp == 0 {
            self.stamps.fill(0);
            self.current_stamp = 1;
        }
    }

    #[inline]
    fn pack(x: i32, y: i32) -> u64 {
        ((x as u32 as u64) << 32) | (y as u32 as u64)
    }

    #[inline]
    fn bucket(&self, key: u64) -> usize {
        let h = key.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        ((h ^ (h >> 32)) as usize) & self.mask
    }

    /// Insert `(x, y)`; returns `false` if already present in the current stamp.
    #[inline]
    fn insert(&mut self, x: i32, y: i32) -> bool {
        let key = Self::pack(x, y);
        let mut idx = self.bucket(key);
        loop {
            if self.stamps[idx] != self.current_stamp {
                self.stamps[idx] = self.current_stamp;
                self.keys[idx] = key;
                return true;
            }
            if self.keys[idx] == key {
                return false;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    #[inline]
    fn contains(&self, x: i32, y: i32) -> bool {
        let key = Self::pack(x, y);
        let mut idx = self.bucket(key);
        loop {
            if self.stamps[idx] != self.current_stamp {
                return false;
            }
            if self.keys[idx] == key {
                return true;
            }
            idx = (idx + 1) & self.mask;
        }
    }
}

/// Configuration for Madras-Sokal pivot MCMC simulation.
#[derive(Clone, Debug)]
pub struct PivotConfig {
    pub n: usize,
    pub seed: u64,
    pub burn_in: usize,
    pub n_batches: usize,
    pub batch_size: usize,
    pub n_samples: usize,
    pub sample_stride: usize,
}

/// Output of Madras-Sokal pivot simulation.
#[derive(Clone, Debug)]
pub struct PivotResult {
    pub n: usize,
    pub seed: u64,
    pub r2_mean: f64,
    pub r2_se: f64,
    pub acceptance_rate: f64,
    pub total_steps: usize,
    pub batch_means: Vec<f64>,
    pub samples: Vec<Vec<(i32, i32)>>,
}

fn pivot_step(
    pts: &mut [(i32, i32)],
    cand: &mut [(i32, i32)],
    hash: &mut StampHashSet,
    rng: &mut Xoshiro256PlusPlus,
) -> bool {
    let n = pts.len() - 1;
    let k = rng.next_bounded(n as u64) as usize;
    // Choose uniformly among the 7 non-identity D4 symmetries (indices 1..=7).
    let sym = (rng.next_bounded(7) as u8) + 1;
    let (xk, yk) = pts[k];

    if k == 0 {
        for pt in pts.iter_mut().skip(1) {
            let (rx, ry) = apply_d4(sym, pt.0, pt.1);
            *pt = (rx, ry);
        }
        return true;
    }

    hash.clear();
    for &(px, py) in &pts[..k] {
        hash.insert(px, py);
    }

    for i in (k + 1)..=n {
        let (dx, dy) = (pts[i].0 - xk, pts[i].1 - yk);
        let (rx, ry) = apply_d4(sym, dx, dy);
        let nx = xk + rx;
        let ny = yk + ry;
        if hash.contains(nx, ny) {
            return false;
        }
        cand[i] = (nx, ny);
    }

    pts[(k + 1)..=n].copy_from_slice(&cand[(k + 1)..=n]);
    true
}

/// Run the Madras-Sokal pivot sampler on `N`-step SAWs in Z^2.
pub fn run_pivot(cfg: &PivotConfig) -> PivotResult {
    assert!(cfg.n >= 1, "walk length N must be >= 1");
    assert!(cfg.n_batches >= 2, "need at least 2 batches for SE");
    assert!(cfg.batch_size >= 1, "batch_size must be >= 1");

    let mut rng = Xoshiro256PlusPlus::from_seed(cfg.seed);
    // Initialize as a valid staircase SAW to avoid artificial straight-rod transient.
    let mut pts = Vec::with_capacity(cfg.n + 1);
    let mut x = 0i32;
    let mut y = 0i32;
    pts.push((0, 0));
    for i in 0..cfg.n {
        if i.is_multiple_of(2) {
            x += 1;
        } else {
            y += 1;
        }
        pts.push((x, y));
    }

    let mut cand = vec![(0i32, 0i32); cfg.n + 1];
    let mut hash = StampHashSet::with_capacity((cfg.n + 1) * 4);

    for _ in 0..cfg.burn_in {
        pivot_step(&mut pts, &mut cand, &mut hash, &mut rng);
    }

    let mut accepted = 0usize;
    let total_steps = cfg.n_batches * cfg.batch_size;
    let mut batch_means = Vec::with_capacity(cfg.n_batches);
    let mut samples = Vec::with_capacity(cfg.n_samples);
    let stride = cfg.sample_stride.max(1);
    let mut step_counter = 0usize;

    for _ in 0..cfg.n_batches {
        let mut sum_r2 = 0.0f64;
        for _ in 0..cfg.batch_size {
            if pivot_step(&mut pts, &mut cand, &mut hash, &mut rng) {
                accepted += 1;
            }
            let (ex, ey) = pts[cfg.n];
            let r2 = (ex as f64) * (ex as f64) + (ey as f64) * (ey as f64);
            sum_r2 += r2;
            step_counter += 1;
            if samples.len() < cfg.n_samples && step_counter.is_multiple_of(stride) {
                samples.push(pts.clone());
            }
        }
        batch_means.push(sum_r2 / (cfg.batch_size as f64));
    }

    while samples.len() < cfg.n_samples {
        pivot_step(&mut pts, &mut cand, &mut hash, &mut rng);
        step_counter += 1;
        if step_counter.is_multiple_of(stride) {
            samples.push(pts.clone());
        }
    }

    let b_f = cfg.n_batches as f64;
    let r2_mean = batch_means.iter().sum::<f64>() / b_f;
    let var_b = batch_means
        .iter()
        .map(|&m| {
            let d = m - r2_mean;
            d * d
        })
        .sum::<f64>()
        / (b_f - 1.0);
    let r2_se = (var_b / b_f).sqrt();

    PivotResult {
        n: cfg.n,
        seed: cfg.seed,
        r2_mean,
        r2_se,
        acceptance_rate: (accepted as f64) / (total_steps as f64),
        total_steps,
        batch_means,
        samples,
    }
}

/// Minimal parser for worker JSON lines: `{"id":12,"nmax":20,"weight":4,"steps":"EENW..."}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrefixTask {
    pub id: u64,
    pub nmax: usize,
    pub weight: u64,
    pub steps: String,
}

fn extract_json_u64(line: &str, key: &str) -> Result<u64, String> {
    let pat = format!("\"{}\":", key);
    let start = line
        .find(&pat)
        .ok_or_else(|| format!("missing key {}", key))?
        + pat.len();
    let rest = line[start..].trim_start();
    let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        return Err(format!("invalid numeric value for {}", key));
    }
    digits
        .parse::<u64>()
        .map_err(|e| format!("parse error for {}: {}", key, e))
}

fn extract_json_str(line: &str, key: &str) -> Result<String, String> {
    let pat = format!("\"{}\":", key);
    let start = line
        .find(&pat)
        .ok_or_else(|| format!("missing key {}", key))?
        + pat.len();
    let rest = line[start..].trim_start();
    if !rest.starts_with('"') {
        return Err(format!("expected string for {}", key));
    }
    let after_quote = &rest[1..];
    let end = after_quote
        .find('"')
        .ok_or_else(|| format!("unterminated string for {}", key))?;
    Ok(after_quote[..end].to_string())
}

pub fn parse_prefix_task(line: &str) -> Result<PrefixTask, String> {
    let id = extract_json_u64(line, "id")?;
    let nmax = extract_json_u64(line, "nmax")? as usize;
    let weight = extract_json_u64(line, "weight")?;
    let steps = extract_json_str(line, "steps")?;
    Ok(PrefixTask {
        id,
        nmax,
        weight,
        steps,
    })
}

pub fn format_u64_slice(vals: &[u64]) -> String {
    let mut out = String::from("[");
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&v.to_string());
    }
    out.push(']');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_small_cn_known_values() {
        let res = enumerate_symmetry(12, None);
        let expected: [u64; 13] = [
            1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932,
        ];
        assert_eq!(res.counts, expected);
    }

    #[test]
    fn test_symmetry_x8_matches_naive_bruteforce() {
        let naive = enumerate_naive(10);
        let sym = enumerate_symmetry(10, None);
        assert_eq!(sym.counts, naive);
    }

    #[test]
    fn test_heatmap_sum_and_d4_symmetry() {
        let res = enumerate_symmetry(10, Some(8));
        let grid = res.heatmap_grid.expect("heatmap present");
        let total: u64 = grid.iter().flat_map(|r| r.iter()).sum();
        assert_eq!(total, res.counts[8]);
        let hn = 8usize;
        for dy in -8i32..=8 {
            for dx in -8i32..=8 {
                let base = grid[(dy + hn as i32) as usize][(dx + hn as i32) as usize];
                for sym in 1..8u8 {
                    let (sx, sy) = apply_d4(sym, dx, dy);
                    assert_eq!(
                        grid[(sy + hn as i32) as usize][(sx + hn as i32) as usize],
                        base
                    );
                }
            }
        }
    }

    #[test]
    fn test_prefix_round_trip_matches_single_process() {
        let nmax = 12;
        let k = 5;
        let full = enumerate_symmetry(nmax, None).counts;
        let prefixes = generate_east_prefixes(k);
        assert_eq!(prefixes.len() as u64 * 4, full[k]);

        let mut summed = vec![0u64; nmax + 1];
        for p in &prefixes {
            let line = format!(
                "{{\"id\":0,\"nmax\":{},\"weight\":4,\"steps\":\"{}\"}}",
                nmax, p
            );
            let task = parse_prefix_task(&line).expect("valid task JSON");
            let ext = enumerate_from_prefix(&task.steps, task.weight, task.nmax)
                .expect("valid prefix extension");
            for n in k..=nmax {
                summed[n] += ext[n];
            }
        }
        assert_eq!(&summed[k..=nmax], &full[k..=nmax]);
    }

    #[test]
    fn test_pivot_walks_are_valid_and_collision_free() {
        let cfg = PivotConfig {
            n: 80,
            seed: 20261124,
            burn_in: 500,
            n_batches: 4,
            batch_size: 250,
            n_samples: 25,
            sample_stride: 40,
        };
        let res = run_pivot(&cfg);
        assert_eq!(res.samples.len(), 25);
        for walk in &res.samples {
            assert!(is_valid_saw(walk, 80));
        }
        assert!(res.r2_mean > 20.0 && res.r2_se > 0.0);
    }

    #[test]
    fn test_pivot_seed_determinism() {
        let cfg = PivotConfig {
            n: 60,
            seed: 424242,
            burn_in: 300,
            n_batches: 5,
            batch_size: 200,
            n_samples: 5,
            sample_stride: 50,
        };
        let r1 = run_pivot(&cfg);
        let r2 = run_pivot(&cfg);
        assert_eq!(r1.r2_mean.to_bits(), r2.r2_mean.to_bits());
        assert_eq!(r1.r2_se.to_bits(), r2.r2_se.to_bits());
        assert_eq!(r1.samples, r2.samples);
    }

    #[test]
    fn test_pivot_small_n_matches_exact_r2() {
        // Exact <R^2> for N=4: compute directly via brute force over all 100 SAWs.
        let mut visited = VisitedBitmap::new();
        visited.set(0, 0);
        fn sum_r2(visited: &mut VisitedBitmap, x: i32, y: i32, step: usize) -> (u64, u64) {
            if step == 4 {
                return (1, ((x * x + y * y) as u64));
            }
            let (mut c, mut s) = (0, 0);
            for &(dx, dy, _) in &DIRS {
                let (nx, ny) = (x + dx, y + dy);
                if !visited.is_visited(nx, ny) {
                    visited.set(nx, ny);
                    let (dc, ds) = sum_r2(visited, nx, ny, step + 1);
                    c += dc;
                    s += ds;
                    visited.unset(nx, ny);
                }
            }
            (c, s)
        }
        let (cnt, sr2) = sum_r2(&mut visited, 0, 0, 0);
        assert_eq!(cnt, 100);
        let exact_r2 = (sr2 as f64) / (cnt as f64);

        let cfg = PivotConfig {
            n: 4,
            seed: 20261124,
            burn_in: 1000,
            n_batches: 20,
            batch_size: 2000,
            n_samples: 0,
            sample_stride: 100,
        };
        let mc = run_pivot(&cfg);
        assert!((mc.r2_mean - exact_r2).abs() < 0.15);
    }
}
