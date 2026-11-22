#![allow(clippy::too_many_arguments, clippy::needless_range_loop)]

use crate::cev::cev_call_exact;
use crate::parallel::default_workers;
use crate::rng::{derive_seed, Rng};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EulerScheme {
    Absorbing,
    FullTruncation,
}

#[derive(Clone, Debug)]
pub struct SabrMcResult {
    pub strikes: Vec<f64>,
    pub prices: Vec<f64>,
    pub ses: Vec<f64>,
    pub p_absorb: f64,
    pub mean_ft: f64,
    pub se_ft: f64,
}

/// Deterministic fixed-block parallel reduction over `Vec<f64>` (`BLOCK = 64`).
fn reduce_vec_fixed(
    n_items: usize,
    vec_len: usize,
    n_workers: usize,
    eval: impl Fn(usize, &mut [f64]) + Sync,
) -> Vec<f64> {
    const BLOCK: usize = 64;
    if n_items == 0 {
        return vec![0.0; vec_len];
    }
    let n_blocks = n_items.div_ceil(BLOCK);
    let w = n_workers.max(1).min(n_blocks);
    let block_sums: Vec<Vec<f64>> = if w == 1 {
        (0..n_blocks)
            .map(|b| {
                let start = b * BLOCK;
                let end = ((b + 1) * BLOCK).min(n_items);
                let mut acc = vec![0.0_f64; vec_len];
                let mut tmp = vec![0.0_f64; vec_len];
                for idx in start..end {
                    tmp.fill(0.0);
                    eval(idx, &mut tmp);
                    for j in 0..vec_len {
                        acc[j] += tmp[j];
                    }
                }
                acc
            })
            .collect()
    } else {
        let bpw = n_blocks.div_ceil(w);
        let mut outer: Vec<Vec<Vec<f64>>> = Vec::with_capacity(w);
        std::thread::scope(|s| {
            let mut handles = Vec::with_capacity(w);
            for tid in 0..w {
                let bs = tid * bpw;
                if bs >= n_blocks {
                    break;
                }
                let be = ((tid + 1) * bpw).min(n_blocks);
                let eval_ref = &eval;
                handles.push(s.spawn(move || {
                    let mut out = Vec::with_capacity(be - bs);
                    let mut tmp = vec![0.0_f64; vec_len];
                    for b in bs..be {
                        let start = b * BLOCK;
                        let end = ((b + 1) * BLOCK).min(n_items);
                        let mut acc = vec![0.0_f64; vec_len];
                        for idx in start..end {
                            tmp.fill(0.0);
                            eval_ref(idx, &mut tmp);
                            for j in 0..vec_len {
                                acc[j] += tmp[j];
                            }
                        }
                        out.push(acc);
                    }
                    out
                }));
            }
            for h in handles {
                outer.push(h.join().expect("worker panicked"));
            }
        });
        let mut flat = Vec::with_capacity(n_blocks);
        for mut part in outer {
            flat.append(&mut part);
        }
        flat
    };
    let mut tot = vec![0.0_f64; vec_len];
    for bs in block_sums {
        for j in 0..vec_len {
            tot[j] += bs[j];
        }
    }
    tot
}

/// Conditional CEV Monte Carlo for `rho = 0`:
/// simulates `alpha_t = alpha0 * exp(nu * W2_t - 0.5 * nu^2 * t)` on an even grid
/// and integrates `I = int_0^T alpha_s^2 ds` via composite Simpson's rule,
/// then averages `C_CEV(F0, K; alpha=1, T=I)`.
pub fn sabr_cond_rho0_mc(
    f0: f64,
    ks: &[f64],
    t: f64,
    alpha0: f64,
    beta: f64,
    nu: f64,
    n_steps_simpson: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> SabrMcResult {
    let n_st = if n_steps_simpson.is_multiple_of(2) {
        n_steps_simpson.max(2)
    } else {
        n_steps_simpson + 1
    };
    let dt = t / n_st as f64;
    let drift = -0.5 * nu * nu * dt;
    let vol_dt = nu * dt.sqrt();
    let nk = ks.len();
    let vlen = 2 * nk;

    let sums = reduce_vec_fixed(n_paths, vlen, n_workers, |r, out| {
        let mut rg = Rng::new(derive_seed(seed, 0x3001, r as u64));
        let mut a = alpha0;
        let mut int_a2 = a * a; // endpoint w=1
        for i in 1..=n_st {
            a *= (drift + vol_dt * rg.normal()).exp();
            let w = if i == n_st {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            int_a2 += w * a * a;
        }
        let i_var = int_a2 * (dt / 3.0);
        for (j, &k) in ks.iter().enumerate() {
            let c = cev_call_exact(f0, k, 1.0, beta, i_var);
            out[2 * j] = c;
            out[2 * j + 1] = c * c;
        }
    });

    let nf = n_paths as f64;
    let mut prices = Vec::with_capacity(nk);
    let mut ses = Vec::with_capacity(nk);
    for j in 0..nk {
        let m = sums[2 * j] / nf;
        let var = (sums[2 * j + 1] / nf - m * m).max(0.0);
        prices.push(m);
        ses.push((var / nf).sqrt());
    }
    SabrMcResult {
        strikes: ks.to_vec(),
        prices,
        ses,
        p_absorb: 0.0,
        mean_ft: f0,
        se_ft: 0.0,
    }
}

/// Full 2D SABR Monte Carlo with log-Euler for `alpha` and `Absorbing` or `FullTruncation` for `F`.
pub fn sabr_euler_mc(
    f0: f64,
    ks: &[f64],
    t: f64,
    alpha0: f64,
    beta: f64,
    rho: f64,
    nu: f64,
    n_steps: usize,
    n_paths: usize,
    scheme: EulerScheme,
    n_workers: usize,
    seed: u64,
) -> SabrMcResult {
    let dt = t / n_steps as f64;
    let sqrt_dt = dt.sqrt();
    let drift_a = -0.5 * nu * nu * dt;
    let vol_a = nu * sqrt_dt;
    let one_m_rho2_sqrt = (1.0 - rho * rho).max(0.0).sqrt();
    let is_half = (beta - 0.5).abs() < 1e-14;
    let nk = ks.len();
    let vlen = 2 * nk + 3;
    let tag = 0x3010_u64
        ^ (n_steps as u64)
        ^ match scheme {
            EulerScheme::Absorbing => 0x11,
            EulerScheme::FullTruncation => 0x22,
        };

    let sums = reduce_vec_fixed(n_paths, vlen, n_workers, |r, out| {
        let mut rg = Rng::new(derive_seed(seed, tag, r as u64));
        let mut f = f0;
        let mut a = alpha0;
        for _ in 0..n_steps {
            let (z1, z2_ind) = rg.normal_pair();
            let dw1 = z1;
            let dw2 = rho * z1 + one_m_rho2_sqrt * z2_ind;
            match scheme {
                EulerScheme::Absorbing => {
                    if f > 0.0 {
                        let f_b = if is_half { f.sqrt() } else { f.powf(beta) };
                        f += a * f_b * sqrt_dt * dw1;
                        if f <= 0.0 {
                            f = 0.0;
                        }
                    }
                }
                EulerScheme::FullTruncation => {
                    let f_pos = f.max(0.0);
                    let f_b = if is_half {
                        f_pos.sqrt()
                    } else {
                        f_pos.powf(beta)
                    };
                    f += a * f_b * sqrt_dt * dw1;
                }
            }
            a *= (drift_a + vol_a * dw2).exp();
        }
        let ft = f.max(0.0);
        for (j, &k) in ks.iter().enumerate() {
            let pay = (ft - k).max(0.0);
            out[2 * j] = pay;
            out[2 * j + 1] = pay * pay;
        }
        out[2 * nk] = if ft == 0.0 { 1.0 } else { 0.0 };
        out[2 * nk + 1] = ft;
        out[2 * nk + 2] = ft * ft;
    });

    let nf = n_paths as f64;
    let mut prices = Vec::with_capacity(nk);
    let mut ses = Vec::with_capacity(nk);
    for j in 0..nk {
        let m = sums[2 * j] / nf;
        let var = (sums[2 * j + 1] / nf - m * m).max(0.0);
        prices.push(m);
        ses.push((var / nf).sqrt());
    }
    let m_ft = sums[2 * nk + 1] / nf;
    let var_ft = (sums[2 * nk + 2] / nf - m_ft * m_ft).max(0.0);
    SabrMcResult {
        strikes: ks.to_vec(),
        prices,
        ses,
        p_absorb: sums[2 * nk] / nf,
        mean_ft: m_ft,
        se_ft: (var_ft / nf).sqrt(),
    }
}

pub fn default_sabr_workers() -> usize {
    default_workers().min(8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sabr_rho0_cond_vs_euler() {
        let f0 = 0.04;
        let ks = [0.02, 0.04, 0.06];
        let t = 1.0;
        let alpha = 0.0325;
        let beta = 0.5;
        let nu = 0.4;
        let cond = sabr_cond_rho0_mc(f0, &ks, t, alpha, beta, nu, 200, 100_000, 4, 20261123);
        let eul = sabr_euler_mc(
            f0,
            &ks,
            t,
            alpha,
            beta,
            0.0,
            nu,
            400,
            200_000,
            EulerScheme::Absorbing,
            4,
            20261123,
        );
        for j in 0..ks.len() {
            let se_comb = (cond.ses[j].powi(2) + eul.ses[j].powi(2)).sqrt();
            let z = (cond.prices[j] - eul.prices[j]).abs() / se_comb;
            assert!(z <= 3.0, "k={} z={z}", ks[j]);
        }
    }
}
