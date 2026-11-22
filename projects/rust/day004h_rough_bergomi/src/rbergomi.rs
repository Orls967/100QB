//! Rough Bergomi model simulation and exact pricing estimators:
//! (1) Variance martingale anchor E[exp(eta*W_tilde - 0.5*eta^2*u^{2H})] = 1 & variance swap strike,
//! (2) Price martingale anchor E[S_T] = S_0,
//! (3) Pure MC log-Euler vs Conditional Black-Scholes (Romano-Touzi / McCrickerd-Pakkanen) & VRF,
//! (4) Put-call parity & exact Cholesky reference comparison (n=50).

use crate::bs;
use crate::linalg;
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;
use crate::volterra::{build_joint_cov_matrix, exact_var, HybridWeights};

#[derive(Clone, Copy, Debug)]
pub struct RBergomiParams {
    pub s0: f64,
    pub xi0: f64,
    pub eta: f64,
    pub h: f64,
    pub rho: f64,
    pub r: f64,
}

impl RBergomiParams {
    pub fn default_params() -> Self {
        Self {
            s0: 1.0,
            xi0: 0.04,
            eta: 1.9,
            h: 0.1,
            rho: -0.7,
            r: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MartingaleAnchorSummary {
    pub t_horizon: f64,
    pub mean_exp_w: f64,
    pub se_exp_w: f64,
    pub z_exp_w: f64,
    pub mean_int_v: f64,
    pub se_int_v: f64,
    pub exact_var_swap: f64,
    pub z_var_swap: f64,
    pub mean_s_t: f64,
    pub se_s_t: f64,
    pub z_s_t: f64,
}

pub fn run_martingale_anchors(
    p: RBergomiParams,
    ts: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<MartingaleAnchorSummary> {
    let mut out = Vec::with_capacity(ts.len());
    let rho_perp = (1.0 - p.rho * p.rho).max(0.0).sqrt();
    for &t_h in ts {
        let hw = HybridWeights::new(n_steps, t_h, p.h);
        let dt = hw.dt;
        let half_eta2 = 0.5 * p.eta * p.eta;
        let rows = parallel::run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x2100 ^ (t_h.to_bits()), r as u64));
            let mut u = vec![0.0_f64; n_steps];
            let mut v = vec![0.0_f64; n_steps];
            let mut w_perp = vec![0.0_f64; n_steps];
            let mut wt = vec![0.0_f64; n_steps];
            let mut db = vec![0.0_f64; n_steps];
            for i in 0..n_steps {
                u[i] = rg.normal();
                v[i] = rg.normal();
                w_perp[i] = rg.normal();
            }
            hw.generate_path(&u, &v, &mut wt, &mut db);

            // Evaluate V_u at grid midpoints/left-points (V_0 = xi0, V_i = xi0 * exp(eta*wt[i-1] - 0.5*eta^2*t_i^{2H}))
            let mut int_v = 0.0_f64;
            let mut log_s = p.s0.ln();
            let sqrt_dt = dt.sqrt();
            for i in 0..n_steps {
                let v_i = if i == 0 {
                    p.xi0
                } else {
                    let t_i = i as f64 * dt;
                    p.xi0 * (p.eta * wt[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
                };
                int_v += v_i * dt;
                let dw = p.rho * db[i] + rho_perp * sqrt_dt * w_perp[i];
                log_s += -0.5 * v_i * dt + v_i.sqrt() * dw;
            }
            let exp_w_t = (p.eta * wt[n_steps - 1] - half_eta2 * exact_var(t_h, p.h)).exp();
            (exp_w_t, int_v, log_s.exp())
        });

        let ew: Vec<f64> = rows.iter().map(|x| x.0).collect();
        let iv: Vec<f64> = rows.iter().map(|x| x.1).collect();
        let st: Vec<f64> = rows.iter().map(|x| x.2).collect();

        let m_ew = stats::mean(&ew);
        let se_ew = stats::se(&ew);
        let m_iv = stats::mean(&iv);
        let se_iv = stats::se(&iv);
        let ex_vs = p.xi0 * t_h;
        let m_st = stats::mean(&st);
        let se_st = stats::se(&st);

        out.push(MartingaleAnchorSummary {
            t_horizon: t_h,
            mean_exp_w: m_ew,
            se_exp_w: se_ew,
            z_exp_w: (m_ew - 1.0).abs() / se_ew.max(1e-12),
            mean_int_v: m_iv,
            se_int_v: se_iv,
            exact_var_swap: ex_vs,
            z_var_swap: (m_iv - ex_vs).abs() / se_iv.max(1e-12),
            mean_s_t: m_st,
            se_s_t: se_st,
            z_s_t: (m_st - p.s0).abs() / se_st.max(1e-12),
        });
    }
    out
}

#[derive(Clone, Debug)]
pub struct CholeskyMartingaleCheck {
    pub mean_exp_w: f64,
    pub se_exp_w: f64,
    pub z_exp_w: f64,
    pub boot_q025: f64,
    pub boot_q500: f64,
    pub boot_q975: f64,
}

/// Tests `E[exp(eta * W_tilde^H_1 - 0.5 * eta^2)] = 1` using exact Gaussian `W_tilde^H_1 ~ N(0, 1)`
/// (matching exact Cholesky terminal law) with 1000-resample bootstrap quantiles for heavy-tail diagnostics.
pub fn run_exact_cholesky_martingale_check(
    p: RBergomiParams,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> CholeskyMartingaleCheck {
    let half_eta2 = 0.5 * p.eta * p.eta;
    let vals = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x2999, r as u64));
        let w = rg.normal();
        (p.eta * w - half_eta2).exp()
    });
    let m = stats::mean(&vals);
    let se = stats::se(&vals);
    let z = (m - 1.0).abs() / se.max(1e-12);

    let n_boot = 1000;
    let mut boot_means = Vec::with_capacity(n_boot);
    let mut brg = Rng::new(derive_seed(seed, 0x2AAA, 0));
    let inv_n = 1.0 / n_paths as f64;
    for _ in 0..n_boot {
        let mut s = 0.0_f64;
        for _ in 0..n_paths {
            let idx = (brg.uniform() * n_paths as f64) as usize;
            s += vals[idx.min(n_paths - 1)];
        }
        boot_means.push(s * inv_n);
    }
    boot_means.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    CholeskyMartingaleCheck {
        mean_exp_w: m,
        se_exp_w: se,
        z_exp_w: z,
        boot_q025: boot_means[(0.025 * n_boot as f64) as usize],
        boot_q500: boot_means[n_boot / 2],
        boot_q975: boot_means[(0.975 * n_boot as f64) as usize],
    }
}

#[derive(Clone, Debug)]
pub struct CallEstimatorRow {
    pub t_horizon: f64,
    pub log_k: f64,
    pub strike: f64,
    pub mc_call: f64,
    pub mc_se: f64,
    pub cond_call: f64,
    pub cond_se: f64,
    pub z_diff: f64,
    pub vrf: f64,
    pub put_call_parity_err: f64,
}

#[derive(Clone, Debug)]
pub struct PricingComparisonSummary {
    pub rows: Vec<CallEstimatorRow>,
    pub mc_secs: f64,
    pub cond_secs: f64,
}

/// Compares (a) pure MC log-Euler vs (b) conditional Black-Scholes (Romano-Touzi / McCrickerd-Pakkanen)
/// across strikes `log_ks` and maturities `ts` using CRN.
pub fn run_call_estimators_comparison(
    p: RBergomiParams,
    ts: &[f64],
    log_ks: &[f64],
    n_steps: usize,
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> PricingComparisonSummary {
    let rho_perp = (1.0 - p.rho * p.rho).max(0.0).sqrt();
    let half_eta2 = 0.5 * p.eta * p.eta;
    let one_m_rho2 = (1.0 - p.rho * p.rho).max(0.0);

    // Precompute conditional inputs (S_eff, V_cond) and pure MC S_T per maturity
    let mut rows = Vec::new();
    let mut mc_secs_total = 0.0_f64;
    let mut cond_secs_total = 0.0_f64;

    for &t_h in ts {
        let hw = HybridWeights::new(n_steps, t_h, p.h);
        let dt = hw.dt;
        let sqrt_dt = dt.sqrt();

        // Time conditional BS path generation
        let tc = std::time::Instant::now();
        let cond_paths = parallel::run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x2200 ^ t_h.to_bits(), r as u64));
            let mut u = vec![0.0_f64; n_steps];
            let mut v = vec![0.0_f64; n_steps];
            let mut wt = vec![0.0_f64; n_steps];
            let mut db = vec![0.0_f64; n_steps];
            for i in 0..n_steps {
                u[i] = rg.normal();
                v[i] = rg.normal();
            }
            hw.generate_path(&u, &v, &mut wt, &mut db);

            let mut int_v = 0.0_f64;
            let mut stoch_int = 0.0_f64;
            for i in 0..n_steps {
                let v_i = if i == 0 {
                    p.xi0
                } else {
                    let t_i = i as f64 * dt;
                    p.xi0 * (p.eta * wt[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
                };
                int_v += v_i * dt;
                stoch_int += v_i.sqrt() * db[i];
            }
            let s_eff = p.s0 * (p.rho * stoch_int - 0.5 * p.rho * p.rho * int_v).exp();
            let v_cond = one_m_rho2 * int_v;
            (s_eff, v_cond)
        });
        cond_secs_total += tc.elapsed().as_secs_f64();

        // Time pure MC path generation (same B-driver + independent W_perp)
        let tm = std::time::Instant::now();
        let mc_paths = parallel::run_chunks(n_paths, n_workers, |r| {
            let mut rg = Rng::new(derive_seed(seed, 0x2200 ^ t_h.to_bits(), r as u64));
            let mut u = vec![0.0_f64; n_steps];
            let mut v = vec![0.0_f64; n_steps];
            let mut wt = vec![0.0_f64; n_steps];
            let mut db = vec![0.0_f64; n_steps];
            for i in 0..n_steps {
                u[i] = rg.normal();
                v[i] = rg.normal();
            }
            hw.generate_path(&u, &v, &mut wt, &mut db);

            let mut log_s = p.s0.ln();
            for i in 0..n_steps {
                let v_i = if i == 0 {
                    p.xi0
                } else {
                    let t_i = i as f64 * dt;
                    p.xi0 * (p.eta * wt[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
                };
                let w_perp = rg.normal();
                let dw = p.rho * db[i] + rho_perp * sqrt_dt * w_perp;
                log_s += -0.5 * v_i * dt + v_i.sqrt() * dw;
            }
            log_s.exp()
        });
        mc_secs_total += tm.elapsed().as_secs_f64();

        let mean_st = stats::mean(&mc_paths);
        let n_f = n_paths as f64;
        for &lk in log_ks {
            let strike = p.s0 * lk.exp();
            let mut sum_mc = 0.0_f64;
            let mut sum_mc2 = 0.0_f64;
            let mut sum_put = 0.0_f64;
            let mut sum_cond = 0.0_f64;
            let mut sum_cond2 = 0.0_f64;
            let mut sum_diff = 0.0_f64;
            let mut sum_diff2 = 0.0_f64;

            for r in 0..n_paths {
                let st = mc_paths[r];
                let c_mc = (st - strike).max(0.0);
                let p_mc = (strike - st).max(0.0);
                let (s_eff, v_cond) = cond_paths[r];
                let c_cond = bs::bs_call_total_var(s_eff, strike, v_cond);
                let d = c_mc - c_cond;

                sum_mc += c_mc;
                sum_mc2 += c_mc * c_mc;
                sum_put += p_mc;
                sum_cond += c_cond;
                sum_cond2 += c_cond * c_cond;
                sum_diff += d;
                sum_diff2 += d * d;
            }

            let mean_mc = sum_mc / n_f;
            let var_mc = (sum_mc2 / n_f - mean_mc * mean_mc).max(0.0) * n_f / (n_f - 1.0);
            let se_mc = (var_mc / n_f).sqrt();

            let mean_cond = sum_cond / n_f;
            let var_cond = (sum_cond2 / n_f - mean_cond * mean_cond).max(0.0) * n_f / (n_f - 1.0);
            let se_cond = (var_cond / n_f).sqrt();

            let mean_d = sum_diff / n_f;
            let var_d = (sum_diff2 / n_f - mean_d * mean_d).max(0.0) * n_f / (n_f - 1.0);
            let se_d = (var_d / n_f).sqrt();

            let mean_put = sum_put / n_f;
            let parity_err = ((mean_mc - mean_put) - (mean_st - strike)).abs();
            let vrf = if var_cond > 0.0 {
                var_mc / var_cond
            } else {
                1.0
            };

            rows.push(CallEstimatorRow {
                t_horizon: t_h,
                log_k: lk,
                strike,
                mc_call: mean_mc,
                mc_se: se_mc,
                cond_call: mean_cond,
                cond_se: se_cond,
                z_diff: mean_d.abs() / se_d.max(1e-12),
                vrf,
                put_call_parity_err: parity_err,
            });
        }
    }

    PricingComparisonSummary {
        rows,
        mc_secs: mc_secs_total,
        cond_secs: cond_secs_total,
    }
}

#[derive(Clone, Debug)]
pub struct RefVsHybridPriceRow {
    pub log_k: f64,
    pub ref_call: f64,
    pub ref_se: f64,
    pub hyb_call: f64,
    pub hyb_se: f64,
    pub z_comb: f64,
}

/// Compares hybrid (`n=50`) conditional BS price vs exact `2n x 2n` Cholesky reference (`n=50`, Euler in S) at `T=1`.
pub fn compare_hybrid_vs_cholesky_ref(
    p: RBergomiParams,
    n: usize,
    t_horizon: f64,
    log_ks: &[f64],
    n_paths: usize,
    n_workers: usize,
    seed: u64,
) -> Vec<RefVsHybridPriceRow> {
    let cov = build_joint_cov_matrix(n, t_horizon, p.h);
    let chol = linalg::cholesky(&cov).expect("joint covariance SPD");
    let dim = 2 * n;
    let dt = t_horizon / n as f64;
    let half_eta2 = 0.5 * p.eta * p.eta;
    let one_m_rho2 = (1.0 - p.rho * p.rho).max(0.0);

    let ref_cond = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x2401, r as u64));
        let mut z = vec![0.0_f64; dim];
        let mut x = vec![0.0_f64; dim];
        for zj in z.iter_mut() {
            *zj = rg.normal();
        }
        for i in 0..dim {
            let row = &chol[i];
            let mut s = 0.0_f64;
            for j in 0..=i {
                s += row[j] * z[j];
            }
            x[i] = s;
        }
        let mut int_v = 0.0_f64;
        let mut stoch_int = 0.0_f64;
        for i in 0..n {
            let v_i = if i == 0 {
                p.xi0
            } else {
                let t_i = i as f64 * dt;
                p.xi0 * (p.eta * x[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
            };
            int_v += v_i * dt;
            stoch_int += v_i.sqrt() * x[n + i];
        }
        let s_eff = p.s0 * (p.rho * stoch_int - 0.5 * p.rho * p.rho * int_v).exp();
        (s_eff, one_m_rho2 * int_v)
    });

    let hw = HybridWeights::new(n, t_horizon, p.h);
    let hyb_cond = parallel::run_chunks(n_paths, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x2402, r as u64));
        let mut u = vec![0.0_f64; n];
        let mut v = vec![0.0_f64; n];
        let mut wt = vec![0.0_f64; n];
        let mut db = vec![0.0_f64; n];
        for i in 0..n {
            u[i] = rg.normal();
            v[i] = rg.normal();
        }
        hw.generate_path(&u, &v, &mut wt, &mut db);
        let mut int_v = 0.0_f64;
        let mut stoch_int = 0.0_f64;
        for i in 0..n {
            let v_i = if i == 0 {
                p.xi0
            } else {
                let t_i = i as f64 * dt;
                p.xi0 * (p.eta * wt[i - 1] - half_eta2 * exact_var(t_i, p.h)).exp()
            };
            int_v += v_i * dt;
            stoch_int += v_i.sqrt() * db[i];
        }
        let s_eff = p.s0 * (p.rho * stoch_int - 0.5 * p.rho * p.rho * int_v).exp();
        (s_eff, one_m_rho2 * int_v)
    });

    let mut out = Vec::with_capacity(log_ks.len());
    for &lk in log_ks {
        let strike = p.s0 * lk.exp();
        let ref_vals: Vec<f64> = ref_cond
            .iter()
            .map(|&(seff, vc)| bs::bs_call_total_var(seff, strike, vc))
            .collect();
        let hyb_vals: Vec<f64> = hyb_cond
            .iter()
            .map(|&(seff, vc)| bs::bs_call_total_var(seff, strike, vc))
            .collect();
        let mr = stats::mean(&ref_vals);
        let ser = stats::se(&ref_vals);
        let mh = stats::mean(&hyb_vals);
        let seh = stats::se(&hyb_vals);
        let zc = (mr - mh).abs() / (ser * ser + seh * seh).sqrt().max(1e-12);
        out.push(RefVsHybridPriceRow {
            log_k: lk,
            ref_call: mr,
            ref_se: ser,
            hyb_call: mh,
            hyb_se: seh,
            z_comb: zc,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_martingale_and_conditional_bs_small() {
        let p = RBergomiParams::default_params();
        let m = run_martingale_anchors(p, &[0.5], 50, 8000, 4, 20261122);
        assert!(m[0].z_exp_w < 3.5, "z_exp_w={}", m[0].z_exp_w);
        assert!(m[0].z_s_t < 3.5, "z_s_t={}", m[0].z_s_t);

        let cm = run_exact_cholesky_martingale_check(p, 10_000, 4, 20261122);
        assert!(cm.z_exp_w < 3.5, "chol z_exp_w={}", cm.z_exp_w);
        assert!(cm.boot_q025 < cm.boot_q975);

        let cmp = run_call_estimators_comparison(p, &[0.5], &[0.0], 50, 8000, 4, 20261122);
        assert!(cmp.rows[0].z_diff < 3.5, "z_diff={}", cmp.rows[0].z_diff);
        assert!(cmp.rows[0].vrf > 1.5, "vrf={}", cmp.rows[0].vrf);
        assert!(cmp.rows[0].put_call_parity_err < 1e-12);
    }
}
