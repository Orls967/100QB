use std::f64::consts::PI;

pub const SEED: u64 = 20261126;

pub fn fnv1a(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h = (h ^ b as u64).wrapping_mul(0x100_0000_01b3);
    }
    h
}

pub fn rng_vec(n: usize, seed: u64) -> Vec<f64> {
    let (mut x, mut s) = (seed, [0u64; 4]);
    for slot in &mut s {
        x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        *slot = z ^ (z >> 31);
    }
    let mut v: Vec<f64> = (0..n)
        .map(|_| {
            let r = (s[0].wrapping_add(s[3])).rotate_left(23).wrapping_add(s[0]);
            let t = s[1] << 17;
            s[2] ^= s[0];
            s[3] ^= s[1];
            s[1] ^= s[2];
            s[0] ^= s[3];
            s[2] ^= t;
            s[3] = s[3].rotate_left(45);
            ((r >> 11) as f64) * (2.0 / ((1u64 << 53) as f64)) - 1.0
        })
        .collect();
    let nrm = v.iter().map(|y| y * y).sum::<f64>().sqrt();
    v.iter_mut().for_each(|y| *y /= nrm);
    v
}

pub struct Sector {
    pub n: usize,
    pub m: usize,
    pub states: Vec<u32>,
    binom: [[usize; 26]; 26],
}

impl Sector {
    pub fn new(n: usize, m: usize) -> Self {
        let mut binom = [[0usize; 26]; 26];
        for i in 0..=25 {
            binom[i][0] = 1;
            for j in 1..=i {
                binom[i][j] = binom[i - 1][j - 1] + binom[i - 1][j];
            }
        }
        let dim = binom[n][m];
        let mut states = Vec::with_capacity(dim);
        let mut s = if m == 0 { 0u32 } else { (1u32 << m) - 1 };
        for _ in 0..dim {
            states.push(s);
            if m == 0 || m == n {
                break;
            }
            let c = s & s.wrapping_neg();
            let r = s.wrapping_add(c);
            s = (((r ^ s) >> 2) / c) | r;
        }
        Self {
            n,
            m,
            states,
            binom,
        }
    }

    #[inline]
    pub fn rank(&self, mut s: u32) -> usize {
        let (mut j, mut idx) = (1usize, 0usize);
        while s != 0 {
            idx += self.binom[s.trailing_zeros() as usize][j];
            j += 1;
            s &= s - 1;
        }
        idx
    }
}

pub fn apply_h(sec: &Sector, delta: f64, v: &[f64], w: &mut [f64]) {
    let n = sec.n;
    let mask = if n == 32 { u32::MAX } else { (1u32 << n) - 1 };
    for (idx, &s) in sec.states.iter().enumerate() {
        let mut om = s ^ (((s >> 1) | ((s & 1) << (n - 1))) & mask);
        let mut acc = 0.25 * delta * ((n as i32 - 2 * om.count_ones() as i32) as f64) * v[idx];
        while om != 0 {
            let i = om.trailing_zeros() as usize;
            let j = if i + 1 == n { 0 } else { i + 1 };
            acc += 0.5 * v[sec.rank(s ^ (1u32 << i) ^ (1u32 << j))];
            om &= om - 1;
        }
        w[idx] = acc;
    }
}

pub fn build_dense_h(sec: &Sector, delta: f64) -> Vec<f64> {
    let (d, n) = (sec.states.len(), sec.n);
    let mut h = vec![0.0; d * d];
    for (row, &s) in sec.states.iter().enumerate() {
        let mut diag = 0.0;
        for i in 0..n {
            let j = (i + 1) % n;
            if ((s >> i) & 1) == ((s >> j) & 1) {
                diag += 0.25 * delta;
            } else {
                diag -= 0.25 * delta;
                h[row * d + sec.rank(s ^ (1u32 << i) ^ (1u32 << j))] += 0.5;
            }
        }
        h[row * d + row] = diag;
    }
    h
}

pub fn jacobi_eigs(a: &mut [f64], n: usize, need_vecs: bool) -> (Vec<f64>, Vec<f64>) {
    let mut v = if need_vecs { vec![0.0; n * n] } else { vec![] };
    for i in 0..n {
        if need_vecs {
            v[i * n + i] = 1.0;
        }
    }
    for _ in 0..80 {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                off += apq * apq;
                if apq.abs() <= 1e-15 {
                    continue;
                }
                let (app, aqq) = (a[p * n + p], a[q * n + q]);
                let tau = (aqq - app) / (2.0 * apq);
                let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
                let (c, s) = (1.0 / (1.0 + t * t).sqrt(), t / (1.0 + t * t).sqrt());
                (a[p * n + p], a[q * n + q], a[p * n + q], a[q * n + p]) =
                    (app - t * apq, aqq + t * apq, 0.0, 0.0);
                for k in 0..n {
                    if k != p && k != q {
                        let (np, nq) = (
                            c * a[k * n + p] - s * a[k * n + q],
                            s * a[k * n + p] + c * a[k * n + q],
                        );
                        (a[k * n + p], a[p * n + k], a[k * n + q], a[q * n + k]) = (np, np, nq, nq);
                    }
                    if need_vecs {
                        let (vp, vq) = (v[k * n + p], v[k * n + q]);
                        (v[k * n + p], v[k * n + q]) = (c * vp - s * vq, s * vp + c * vq);
                    }
                }
            }
        }
        if off.sqrt() < 1e-13 {
            break;
        }
    }
    let mut pairs: Vec<(f64, usize)> = (0..n).map(|i| (a[i * n + i], i)).collect();
    pairs.sort_by(|x, y| x.0.total_cmp(&y.0));
    let vals = pairs.iter().map(|p| p.0).collect();
    let mut vecs = if need_vecs { vec![0.0; n * n] } else { vec![] };
    if need_vecs {
        for (nc, &(_, oc)) in pairs.iter().enumerate() {
            (0..n).for_each(|r| vecs[r * n + nc] = v[r * n + oc]);
        }
    }
    (vals, vecs)
}

pub fn lanczos(
    sec: &Sector,
    delta: f64,
    need_vec: bool,
) -> (f64, usize, Vec<f64>, Vec<f64>, usize) {
    let d = sec.states.len();
    let q0 = rng_vec(d, SEED ^ ((sec.n as u64) << 8) ^ sec.m as u64);
    let (mut qs, mut q_prev, mut w) = (vec![q0], vec![0.0; d], vec![0.0; d]);
    let (mut alpha, mut beta, mut res_hist) = (Vec::new(), Vec::new(), Vec::new());
    let (mut e0, mut y0) = (0.0, vec![1.0]);
    for k in 0..300 {
        let qk = qs.last().unwrap();
        apply_h(sec, delta, qk, &mut w);
        if k > 0 {
            let b = beta[k - 1];
            (0..d).for_each(|i| w[i] -= b * q_prev[i]);
        }
        let a: f64 = qk.iter().zip(&w).map(|(x, y)| x * y).sum();
        alpha.push(a);
        (0..d).for_each(|i| w[i] -= a * qk[i]);
        if qs.len() <= 120 {
            for qj in &qs {
                let dot: f64 = qj.iter().zip(&w).map(|(x, y)| x * y).sum();
                (0..d).for_each(|i| w[i] -= dot * qj[i]);
            }
        }
        let b = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        beta.push(b);
        let m = alpha.len();
        let mut tm = vec![0.0; m * m];
        for i in 0..m {
            tm[i * m + i] = alpha[i];
            if i + 1 < m {
                (tm[i * m + i + 1], tm[(i + 1) * m + i]) = (beta[i], beta[i]);
            }
        }
        let (vals, vecs) = jacobi_eigs(&mut tm, m, true);
        e0 = vals[0];
        y0 = (0..m).map(|r| vecs[r * m]).collect();
        let ritz = (b * y0[m - 1]).abs();
        res_hist.push(ritz);
        if (m >= 4 && ritz < 1e-10) || b < 1e-13 || m >= d {
            break;
        }
        q_prev.copy_from_slice(qs.last().unwrap());
        let q_next: Vec<f64> = w.iter().map(|&x| x / b).collect();
        if qs.len() < 120 {
            qs.push(q_next);
        } else {
            *qs.last_mut().unwrap() = q_next;
        }
    }
    let mut psi = vec![];
    if need_vec && y0.len() <= qs.len() {
        psi = vec![0.0; d];
        for (j, &c) in y0.iter().enumerate() {
            (0..d).for_each(|i| psi[i] += c * qs[j][i]);
        }
        let pn = psi.iter().map(|x| x * x).sum::<f64>().sqrt();
        psi.iter_mut().for_each(|x| *x /= pn);
    }
    (e0, alpha.len(), psi, res_hist, qs.len())
}

pub fn jw_xx_ground(n: usize) -> f64 {
    let (m, phi) = (n / 2, if (n / 2) % 2 == 1 { 0.0 } else { 0.5 });
    let mut eps: Vec<f64> = (0..n)
        .map(|k| (2.0 * PI * (k as f64 + phi) / n as f64).cos())
        .collect();
    eps.sort_by(|a, b| a.total_cmp(b));
    eps.into_iter().take(m).sum()
}

pub fn exp_s2_and_corr(sec: &Sector, psi: &[f64]) -> (f64, Vec<f64>) {
    let (n, mut corr) = (sec.n, vec![0.0; sec.n]);
    let sz_tot = (sec.m as f64) - 0.5 * (n as f64);
    let mut s2 = 0.5 * (n as f64) + sz_tot * sz_tot;
    corr[0] = 0.75;
    for (idx, &s) in sec.states.iter().enumerate() {
        let (p, mut ones) = (psi[idx], s);
        let zmask = (!s) & ((1u32 << n) - 1);
        while ones != 0 {
            let i = ones.trailing_zeros() as usize;
            let mut zeros = zmask;
            while zeros != 0 {
                let j = zeros.trailing_zeros() as usize;
                s2 += p * psi[sec.rank(s ^ (1u32 << i) ^ (1u32 << j))];
                zeros &= zeros - 1;
            }
            ones &= ones - 1;
        }
        for r in 1..n {
            corr[r] += if (s & 1) == ((s >> r) & 1) {
                0.25 * p * p
            } else {
                -0.25 * p * p + 0.5 * p * psi[sec.rank(s ^ 1 ^ (1u32 << r))]
            };
        }
    }
    (s2, corr)
}
