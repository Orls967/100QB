#[derive(Clone, Default)]
pub struct DtmcRes { pub p_block: f64, pub l: f64, pub l_q: f64, pub u: f64, pub w: f64, pub w_q: f64, pub lam_eff: f64 }

fn binom_probs(m: usize, mu: f64) -> [f64; 8] {
    let mut out = [0.0f64; 8]; out[0] = 1.0;
    for _ in 0..m {
        for k in (1..=m).rev() { out[k] = out[k] * (1.0 - mu) + out[k - 1] * mu; }
        out[0] *= 1.0 - mu;
    }
    out
}

pub fn solve_dtmc_single(c: usize, k_cap: usize, lam: f64, mu: f64) -> DtmcRes {
    let mut p = [[0.0f64; 64]; 64];
    for n in 0..=k_cap {
        let m = n.min(c);
        let bp = binom_probs(m, mu);
        for d in 0..=m {
            if n < k_cap {
                p[n][n - d] += (1.0 - lam) * bp[d];
                p[n][n + 1 - d] += lam * bp[d];
            } else {
                p[n][n - d] += bp[d];
            }
        }
    }
    let mut pi = [0.0f64; 64];
    for i in 0..=k_cap { pi[i] = 1.0 / (k_cap + 1) as f64; }
    for _ in 0..25000 {
        let (mut maxd, mut sum) = (0.0f64, 0.0f64);
        for j in 0..=k_cap {
            let mut nv = 0.0;
            for i in 0..=k_cap { nv += pi[i] * p[i][j]; }
            maxd = maxd.max((nv - pi[j]).abs());
            pi[j] = nv; sum += nv;
        }
        for j in 0..=k_cap { pi[j] /= sum; }
        if maxd < 1e-14 { break; }
    }
    let mut r = DtmcRes::default();
    r.p_block = pi[k_cap];
    r.lam_eff = lam * (1.0 - pi[k_cap]);
    for n in 0..=k_cap {
        let m = n.min(c);
        r.l += n as f64 * pi[n];
        r.l_q += (n - m) as f64 * pi[n];
        r.u += (m as f64 / c as f64) * pi[n];
    }
    r.w = if r.lam_eff > 0.0 { r.l / r.lam_eff } else { 0.0 };
    r.w_q = if r.lam_eff > 0.0 { r.l_q / r.lam_eff } else { 0.0 };
    r
}

pub fn solve_geo_closed_form(lam: f64, mu: f64) -> DtmcRes {
    let mut r = DtmcRes::default();
    r.lam_eff = lam; r.u = lam / mu;
    r.l = (lam * (1.0 - lam)) / (mu - lam);
    r.w = (1.0 - lam) / (mu - lam);
    r.w_q = r.w - 1.0 / mu; r.l_q = lam * r.w_q;
    r
}

pub fn solve_dtmc_tandem(k1: usize, k2: usize, lam: f64, mu1: f64, mu2: f64) -> (DtmcRes, DtmcRes, f64) {
    let s_tot = (k1 + 1) * (k2 + 1);
    let mut p = [[0.0f64; 64]; 64];
    for n1 in 0..=k1 {
        for n2 in 0..=k2 {
            let s_idx = n1 * (k2 + 1) + n2;
            for a in 0..=1 {
                let pa = if a != 0 { lam } else { 1.0 - lam };
                for d1 in 0..=(if n1 > 0 { 1 } else { 0 }) {
                    let pd1 = if n1 == 0 { 1.0 } else if d1 != 0 { mu1 } else { 1.0 - mu1 };
                    for d2 in 0..=(if n2 > 0 { 1 } else { 0 }) {
                        let pd2 = if n2 == 0 { 1.0 } else if d2 != 0 { mu2 } else { 1.0 - mu2 };
                        let nn1 = if n1 < k1 { n1 + a - d1 } else { n1 - d1 };
                        let n2_mid = n2 - d2;
                        let nn2 = if n2_mid < k2 { n2_mid + d1 } else { n2_mid };
                        p[s_idx][nn1 * (k2 + 1) + nn2] += pa * pd1 * pd2;
                    }
                }
            }
        }
    }
    let mut pi = [0.0f64; 64];
    for i in 0..s_tot { pi[i] = 1.0 / s_tot as f64; }
    for _ in 0..25000 {
        let (mut maxd, mut sum) = (0.0f64, 0.0f64);
        for j in 0..s_tot {
            let mut nv = 0.0;
            for i in 0..s_tot { nv += pi[i] * p[i][j]; }
            maxd = maxd.max((nv - pi[j]).abs());
            pi[j] = nv; sum += nv;
        }
        for j in 0..s_tot { pi[j] /= sum; }
        if maxd < 1e-14 { break; }
    }
    let (mut r1, mut r2) = (DtmcRes::default(), DtmcRes::default());
    let (mut blk1, mut blk2) = (0.0, 0.0);
    for n1 in 0..=k1 {
        for n2 in 0..=k2 {
            let pv = pi[n1 * (k2 + 1) + n2];
            r1.l += n1 as f64 * pv; r2.l += n2 as f64 * pv;
            if n1 > 0 { r1.u += pv; r1.l_q += (n1 - 1) as f64 * pv; }
            if n2 > 0 { r2.u += pv; r2.l_q += (n2 - 1) as f64 * pv; }
            if n1 == k1 { blk1 += pv * lam; }
            if n1 > 0 && n2 == k2 { blk2 += pv * mu1 * (1.0 - mu2); }
        }
    }
    r1.p_block = blk1; r1.lam_eff = lam - blk1;
    r1.w = if r1.lam_eff > 0.0 { r1.l / r1.lam_eff } else { 0.0 };
    r1.w_q = if r1.lam_eff > 0.0 { r1.l_q / r1.lam_eff } else { 0.0 };
    r2.p_block = blk2; r2.lam_eff = r1.lam_eff - blk2;
    r2.w = if r2.lam_eff > 0.0 { r2.l / r2.lam_eff } else { 0.0 };
    r2.w_q = if r2.lam_eff > 0.0 { r2.l_q / r2.lam_eff } else { 0.0 };
    let pbt = blk2 / if r1.lam_eff > 0.0 { r1.lam_eff } else { 1.0 };
    (r1, r2, pbt)
}
