//! Arbitrage-Free Nelson-Siegel (AFNS) closed-form loadings, convexity adjustment,
//! Riccati RK4 verification, P-measure exact transition, and 10y term premium.

use crate::linalg;
use crate::ns;
use crate::parallel;
use crate::rng::{derive_seed, Rng};
use crate::stats;

#[derive(Clone, Copy, Debug)]
pub struct AfnsParams {
    pub lambda: f64,
    pub sigma: [f64; 3],
    pub kp: [[f64; 3]; 3],
    pub theta_p: [f64; 3],
}

impl AfnsParams {
    pub fn default_params() -> Self {
        Self {
            lambda: 0.6,
            sigma: [0.0060, 0.0100, 0.0200],
            kp: [[0.05, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.4]],
            theta_p: [0.045, -0.020, 0.010],
        }
    }

    pub fn is_diagonal_kp(&self) -> bool {
        self.kp[1][0] == 0.0
            && self.kp[2][0] == 0.0
            && self.kp[2][1] == 0.0
            && self.kp[0][1] == 0.0
            && self.kp[0][2] == 0.0
            && self.kp[1][2] == 0.0
    }

    pub fn b_vec(&self, tau: f64) -> [f64; 3] {
        if tau <= 0.0 {
            return [0.0, 0.0, 0.0];
        }
        [
            -tau,
            -tau * ns::g1(tau, self.lambda),
            -tau * ns::g2(tau, self.lambda),
        ]
    }

    pub fn yield_loadings(&self, tau: f64) -> [f64; 3] {
        [1.0, ns::g1(tau, self.lambda), ns::g2(tau, self.lambda)]
    }

    pub fn yield_adjustment_components(&self, tau: f64) -> (f64, f64, f64, f64) {
        if tau <= 0.0 {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let s1_sq = self.sigma[0] * self.sigma[0];
        let s2_sq = self.sigma[1] * self.sigma[1];
        let s3_sq = self.sigma[2] * self.sigma[2];
        let adj1 = -s1_sq * tau * tau / 6.0;

        let lam = self.lambda;
        let u = lam * tau;
        if u < 0.05 {
            let nodes = linalg::gauss_legendre(16, 0.0, tau);
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            for nd in nodes {
                let bv = self.b_vec(nd.x);
                sum2 += nd.w * bv[1] * bv[1];
                sum3 += nd.w * bv[2] * bv[2];
            }
            let adj2 = -0.5 * s2_sq * sum2 / tau;
            let adj3 = -0.5 * s3_sq * sum3 / tau;
            return (adj1, adj2, adj3, adj1 + adj2 + adj3);
        }

        let lam2 = lam * lam;
        let lam3 = lam2 * lam;
        let eu = (-u).exp();
        let e2u = eu * eu;
        let adj2 = -s2_sq
            * (1.0 / (2.0 * lam2) - (1.0 - eu) / (lam3 * tau) + (1.0 - e2u) / (4.0 * lam3 * tau));
        let adj3 = -s3_sq
            * (0.5 / lam2 + eu / lam2
                - 0.25 * u * e2u / lam2
                - 0.75 * e2u / lam2
                - 11.0 / (8.0 * lam3 * tau)
                + 2.0 * eu / (lam3 * tau)
                - 5.0 * e2u / (8.0 * lam3 * tau));
        (adj1, adj2, adj3, adj1 + adj2 + adj3)
    }

    pub fn yield_adjustment(&self, tau: f64) -> (f64, f64) {
        let (_, _, _, total) = self.yield_adjustment_components(tau);
        (total, -tau * total)
    }

    pub fn yield_at(&self, x: [f64; 3], tau: f64) -> f64 {
        let (adj, _) = self.yield_adjustment(tau);
        let h = self.yield_loadings(tau);
        adj + h[0] * x[0] + h[1] * x[1] + h[2] * x[2]
    }

    pub fn solve_riccati_rk4(&self, tau: f64, n_steps: usize) -> ([f64; 3], f64, [f64; 3], f64) {
        let lam = self.lambda;
        let s_sq = [
            self.sigma[0] * self.sigma[0],
            self.sigma[1] * self.sigma[1],
            self.sigma[2] * self.sigma[2],
        ];
        let deriv = |_t: f64, y: &[f64], dydt: &mut [f64]| {
            dydt[0] = -1.0;
            dydt[1] = -1.0 - lam * y[1];
            dydt[2] = lam * y[1] - lam * y[2];
            dydt[3] = 0.5 * (s_sq[0] * y[0] * y[0] + s_sq[1] * y[1] * y[1] + s_sq[2] * y[2] * y[2]);
        };
        let y0 = [0.0; 4];
        let sol = linalg::rk4_system(0.0, tau, &y0, n_steps, deriv);
        let sol_r = linalg::rk4_richardson(0.0, tau, &y0, n_steps / 2, deriv);
        (
            [sol[0], sol[1], sol[2]],
            sol[3],
            [sol_r[0], sol_r[1], sol_r[2]],
            sol_r[3],
        )
    }

    pub fn transition_matrix_closed(&self, dt: f64) -> [[f64; 3]; 3] {
        let k11 = self.kp[0][0];
        let k22 = self.kp[1][1];
        let k33 = self.kp[2][2];
        let k21 = self.kp[1][0];
        let k31 = self.kp[2][0];
        let k32 = self.kp[2][1];
        let mut phi = [[0.0; 3]; 3];
        phi[0][0] = (-k11 * dt).exp();
        phi[1][1] = (-k22 * dt).exp();
        phi[2][2] = (-k33 * dt).exp();
        if k21 != 0.0 {
            phi[1][0] = -k21 * phi1(k11, k22, dt);
        }
        if k32 != 0.0 {
            phi[2][1] = -k32 * phi1(k22, k33, dt);
        }
        if k31 != 0.0 {
            phi[2][0] = -k31 * phi1(k11, k33, dt);
        }
        if k21 != 0.0 && k32 != 0.0 {
            let psi2 = if (k11 - k22).abs() < 1e-6 && (k22 - k33).abs() < 1e-6 {
                0.5 * dt * dt * (-k11 * dt).exp()
            } else if (k11 - k22).abs() >= 1e-6 {
                (phi1(k11, k33, dt) - phi1(k22, k33, dt)) / (k22 - k11)
            } else {
                (phi1(k11, k33, dt) - dt * (-k22 * dt).exp()) / (k33 - k22)
            };
            phi[2][0] += k32 * k21 * psi2;
        }
        phi
    }

    pub fn transition_cov_closed(&self, dt: f64) -> [[f64; 3]; 3] {
        let mut q = [[0.0; 3]; 3];
        if dt <= 0.0 {
            return q;
        }
        if self.is_diagonal_kp() {
            for i in 0..3 {
                let k = self.kp[i][i];
                let u = 2.0 * k * dt;
                let fac = if u.abs() < 1e-5 {
                    dt * (1.0 - 0.5 * u + u * u / 6.0 - u * u * u / 24.0)
                } else {
                    -((-u).exp_m1()) / (2.0 * k)
                };
                q[i][i] = self.sigma[i] * self.sigma[i] * fac;
            }
            return q;
        }
        self.transition_cov_gl(dt, 32)
    }

    pub fn transition_cov_gl(&self, dt: f64, quad_order: usize) -> [[f64; 3]; 3] {
        let mut q = [[0.0; 3]; 3];
        if dt <= 0.0 {
            return q;
        }
        let nodes = linalg::gauss_legendre(quad_order, 0.0, dt);
        let s_sq = [
            self.sigma[0] * self.sigma[0],
            self.sigma[1] * self.sigma[1],
            self.sigma[2] * self.sigma[2],
        ];
        for nd in nodes {
            let phi = self.transition_matrix_closed(nd.x);
            for i in 0..3 {
                for j in i..3 {
                    let mut sum = 0.0;
                    for k in 0..3 {
                        sum += phi[i][k] * s_sq[k] * phi[j][k];
                    }
                    q[i][j] += nd.w * sum;
                }
            }
        }
        q[1][0] = q[0][1];
        q[2][0] = q[0][2];
        q[2][1] = q[1][2];
        q
    }

    pub fn transition_cov_van_loan(&self, dt: f64) -> [[f64; 3]; 3] {
        let mut m = linalg::zeros(6, 6);
        for i in 0..3 {
            m[i][3 + i] = self.sigma[i] * self.sigma[i] * dt;
            for j in 0..3 {
                m[i][j] = -self.kp[i][j] * dt;
                m[3 + i][3 + j] = self.kp[j][i] * dt;
            }
        }
        let exp_m = linalg::mat_exp_series(&m, 35);
        let mut q = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let mut sum = 0.0;
                for k in 0..3 {
                    sum += exp_m[i][3 + k] * exp_m[j][k];
                }
                q[i][j] = sum;
            }
        }
        q
    }

    pub fn stationary_cov(&self) -> [[f64; 3]; 3] {
        let mut q = [[0.0; 3]; 3];
        let s_sq = [
            self.sigma[0] * self.sigma[0],
            self.sigma[1] * self.sigma[1],
            self.sigma[2] * self.sigma[2],
        ];
        for i in 0..3 {
            for j in 0..=i {
                let mut rhs = if i == j { s_sq[i] } else { 0.0 };
                for m in 0..i {
                    rhs -= self.kp[i][m] * q[m][j];
                }
                for m in 0..j {
                    rhs -= self.kp[j][m] * q[i][m];
                }
                let val = rhs / (self.kp[i][i] + self.kp[j][j]);
                q[i][j] = val;
                q[j][i] = val;
            }
        }
        q
    }

    pub fn market_price_of_risk(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        let kq = [
            [0.0, 0.0, 0.0],
            [0.0, self.lambda, -self.lambda],
            [0.0, 0.0, self.lambda],
        ];
        let mut lam0 = [0.0; 3];
        let mut lam1 = [[0.0; 3]; 3];
        for i in 0..3 {
            let mut sum = 0.0;
            for j in 0..3 {
                sum += self.kp[i][j] * self.theta_p[j];
                lam1[i][j] = (kq[i][j] - self.kp[i][j]) / self.sigma[i];
            }
            lam0[i] = sum / self.sigma[i];
        }
        (lam0, lam1)
    }

    pub fn long_run_term_premium(&self, tau: f64) -> (f64, f64, f64, f64, f64, f64) {
        let (adj, _) = self.yield_adjustment(tau);
        let h = self.yield_loadings(tau);
        let lr_sr = self.theta_p[0] + self.theta_p[1];
        let lr_y = adj + h[0] * self.theta_p[0] + h[1] * self.theta_p[1] + h[2] * self.theta_p[2];
        let slope_term = (h[1] - 1.0) * self.theta_p[1];
        let curv_term = h[2] * self.theta_p[2];
        (lr_y - lr_sr, lr_y, lr_sr, adj, slope_term, curv_term)
    }
}

fn phi1(a: f64, b: f64, dt: f64) -> f64 {
    let d = (a - b) * dt;
    if d.abs() < 1e-6 {
        let e = (-a * dt).exp();
        dt * e * (1.0 + d * (0.5 + d * (1.0 / 6.0 + d / 24.0)))
    } else {
        ((-a * dt).exp() - (-b * dt).exp()) / (b - a)
    }
}

pub fn cholesky3(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut l = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..=i {
            let mut s = a[i][j];
            for k in 0..j {
                s -= l[i][k] * l[j][k];
            }
            if i == j {
                l[i][j] = s.max(0.0).sqrt();
            } else if l[j][j] > 0.0 {
                l[i][j] = s / l[j][j];
            }
        }
    }
    l
}

pub fn simulate_stationary_path(
    p: &AfnsParams,
    taus: &[f64],
    t_steps: usize,
    dt: f64,
    sigma_eps: f64,
    rg: &mut Rng,
) -> (Vec<[f64; 3]>, Vec<Vec<f64>>) {
    simulate_stationary_path_opt(p, taus, t_steps, dt, sigma_eps, true, rg)
}

pub fn simulate_stationary_path_opt(
    p: &AfnsParams,
    taus: &[f64],
    t_steps: usize,
    dt: f64,
    sigma_eps: f64,
    use_afns_adj: bool,
    rg: &mut Rng,
) -> (Vec<[f64; 3]>, Vec<Vec<f64>>) {
    let m = taus.len();
    let mut a_vec = vec![0.0; m];
    let mut h_mat = Vec::with_capacity(m);
    for (j, &tau) in taus.iter().enumerate() {
        h_mat.push(p.yield_loadings(tau));
        if use_afns_adj {
            a_vec[j] = p.yield_adjustment(tau).0;
        }
    }
    let phi = p.transition_matrix_closed(dt);
    let q_dt = p.transition_cov_closed(dt);
    let q_inf = p.stationary_cov();
    let chol_inf = cholesky3(q_inf);
    let chol_dt = cholesky3(q_dt);

    let mut xs = Vec::with_capacity(t_steps);
    let mut ys = Vec::with_capacity(t_steps);

    let z0 = [rg.normal(), rg.normal(), rg.normal()];
    let mut x = [
        p.theta_p[0] + chol_inf[0][0] * z0[0],
        p.theta_p[1] + chol_inf[1][0] * z0[0] + chol_inf[1][1] * z0[1],
        p.theta_p[2] + chol_inf[2][0] * z0[0] + chol_inf[2][1] * z0[1] + chol_inf[2][2] * z0[2],
    ];

    for t in 0..t_steps {
        if t > 0 {
            let dx = [
                x[0] - p.theta_p[0],
                x[1] - p.theta_p[1],
                x[2] - p.theta_p[2],
            ];
            let z = [rg.normal(), rg.normal(), rg.normal()];
            let mut next_x = [0.0; 3];
            for i in 0..3 {
                let mut e = 0.0;
                for k in 0..=i {
                    e += chol_dt[i][k] * z[k];
                }
                next_x[i] =
                    p.theta_p[i] + phi[i][0] * dx[0] + phi[i][1] * dx[1] + phi[i][2] * dx[2] + e;
            }
            x = next_x;
        }
        xs.push(x);
        let mut yt = Vec::with_capacity(m);
        for j in 0..m {
            let mut val = a_vec[j] + h_mat[j][0] * x[0] + h_mat[j][1] * x[1] + h_mat[j][2] * x[2];
            if sigma_eps > 0.0 {
                val += sigma_eps * rg.normal();
            }
            yt.push(val);
        }
        ys.push(yt);
    }
    (xs, ys)
}

#[derive(Clone, Debug)]
pub struct TransitionMcReport {
    pub max_z_score: f64,
    pub mc_mean: [f64; 3],
    pub theo_mean: [f64; 3],
    pub mc_cov: [[f64; 3]; 3],
    pub theo_cov: [[f64; 3]; 3],
}

pub fn verify_transition_mc(
    p: &AfnsParams,
    x0: [f64; 3],
    dt: f64,
    n_paths: usize,
    n_chunks: usize,
    n_workers: usize,
    seed: u64,
) -> TransitionMcReport {
    let phi = p.transition_matrix_closed(dt);
    let q_dt = p.transition_cov_closed(dt);
    let chol = cholesky3(q_dt);
    let dx = [
        x0[0] - p.theta_p[0],
        x0[1] - p.theta_p[1],
        x0[2] - p.theta_p[2],
    ];
    let mut theo_mean = [0.0; 3];
    for i in 0..3 {
        theo_mean[i] = p.theta_p[i] + phi[i][0] * dx[0] + phi[i][1] * dx[1] + phi[i][2] * dx[2];
    }

    let per = n_paths / n_chunks;
    let rem = n_paths % n_chunks;
    struct Acc {
        sum_x: [f64; 3],
        sum_cov: [[f64; 3]; 3],
        sum_c2: [[f64; 3]; 3],
    }
    let chunks = parallel::run_chunks(n_chunks, n_workers, |c| {
        let cnt = per + if c < rem { 1 } else { 0 };
        let mut rg = Rng::new(derive_seed(seed, 0x2301, c as u64));
        let mut acc = Acc {
            sum_x: [0.0; 3],
            sum_cov: [[0.0; 3]; 3],
            sum_c2: [[0.0; 3]; 3],
        };
        for _ in 0..cnt {
            let z0 = rg.normal();
            let z1 = rg.normal();
            let z2 = rg.normal();
            let e = [
                chol[0][0] * z0,
                chol[1][0] * z0 + chol[1][1] * z1,
                chol[2][0] * z0 + chol[2][1] * z1 + chol[2][2] * z2,
            ];
            for i in 0..3 {
                acc.sum_x[i] += theo_mean[i] + e[i];
                for j in i..3 {
                    let pij = e[i] * e[j];
                    acc.sum_cov[i][j] += pij;
                    acc.sum_c2[i][j] += pij * pij;
                }
            }
        }
        acc
    });

    let mut tot_x = [0.0; 3];
    let mut tot_cov = [[0.0; 3]; 3];
    let mut tot_c2 = [[0.0; 3]; 3];
    for ch in chunks {
        for i in 0..3 {
            tot_x[i] += ch.sum_x[i];
            for j in i..3 {
                tot_cov[i][j] += ch.sum_cov[i][j];
                tot_c2[i][j] += ch.sum_c2[i][j];
            }
        }
    }
    let fn_p = n_paths as f64;
    let mut mc_mean = [0.0; 3];
    let mut mc_cov = [[0.0; 3]; 3];
    let mut max_z = 0.0f64;
    for i in 0..3 {
        mc_mean[i] = tot_x[i] / fn_p;
        let se_m = (q_dt[i][i] / fn_p).sqrt();
        let zm = (mc_mean[i] - theo_mean[i]).abs() / se_m.max(1e-15);
        max_z = max_z.max(zm);
        for j in i..3 {
            mc_cov[i][j] = tot_cov[i][j] / fn_p;
            mc_cov[j][i] = mc_cov[i][j];
            let var_p = (tot_c2[i][j] / fn_p - mc_cov[i][j] * mc_cov[i][j]).max(0.0);
            let se_c = (var_p / fn_p).sqrt();
            let zc = (mc_cov[i][j] - q_dt[i][j]).abs() / se_c.max(1e-18);
            max_z = max_z.max(zc);
        }
    }
    TransitionMcReport {
        max_z_score: max_z,
        mc_mean,
        theo_mean,
        mc_cov,
        theo_cov: q_dt,
    }
}

pub fn verify_long_run_tp_mc(
    p: &AfnsParams,
    tau: f64,
    n_paths: usize,
    n_chunks: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64) {
    let (tp_closed, _, _, _, _, _) = p.long_run_term_premium(tau);
    let q_inf = p.stationary_cov();
    let chol = cholesky3(q_inf);
    let per = n_paths / n_chunks;
    let rem = n_paths % n_chunks;
    let chunks = parallel::run_chunks(n_chunks, n_workers, |c| {
        let cnt = per + if c < rem { 1 } else { 0 };
        let mut rg = Rng::new(derive_seed(seed, 0x2401, c as u64));
        let mut vals = Vec::with_capacity(cnt);
        for _ in 0..cnt {
            let z0 = rg.normal();
            let z1 = rg.normal();
            let z2 = rg.normal();
            let x = [
                p.theta_p[0] + chol[0][0] * z0,
                p.theta_p[1] + chol[1][0] * z0 + chol[1][1] * z1,
                p.theta_p[2] + chol[2][0] * z0 + chol[2][1] * z1 + chol[2][2] * z2,
            ];
            let y = p.yield_at(x, tau);
            let r_t = x[0] + x[1];
            vals.push(y - r_t);
        }
        vals
    });
    let mut all = Vec::with_capacity(n_paths);
    for ch in chunks {
        all.extend(ch);
    }
    let m = stats::mean(&all);
    let se_m = stats::se(&all);
    let z = (m - tp_closed).abs() / se_m.max(1e-15);
    (tp_closed, m, se_m, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_afns_riccati_and_van_loan() {
        let p = AfnsParams::default_params();
        for &tau in &[0.5, 1.0, 5.0, 10.0, 20.0, 30.0] {
            let b_closed = p.b_vec(tau);
            let (adj_closed, a_closed) = p.yield_adjustment(tau);
            let (_, _, b_rich, a_rich) = p.solve_riccati_rk4(tau, 400);
            for i in 0..3 {
                assert!((b_closed[i] - b_rich[i]).abs() < 1e-11);
            }
            assert!((a_closed - a_rich).abs() < 1e-11);
            assert!((adj_closed - (-a_rich / tau)).abs() < 1e-11);
            assert!(adj_closed <= 0.0);
        }

        let q_closed = p.transition_cov_closed(1.0);
        let q_vl = p.transition_cov_van_loan(1.0);
        let q_gl = p.transition_cov_gl(1.0, 32);
        for i in 0..3 {
            for j in 0..3 {
                assert!((q_closed[i][j] - q_vl[i][j]).abs() < 1e-14);
                assert!((q_closed[i][j] - q_gl[i][j]).abs() < 1e-14);
            }
        }
    }
}
