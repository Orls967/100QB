//! Kalman filter (Joseph form + fast steady-state Riccati path), RTS smoother,
//! T=5 joint Gaussian anchor, and 11-parameter BFGS MLE with Hessian SE & diagnostics.

use crate::afns::AfnsParams;
use crate::linalg;
use crate::optim;
use crate::rng::Rng;

pub const STANDARD_TENORS: [f64; 10] = [0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0];
pub const WEEKLY_DT: f64 = 1.0 / 52.0;
pub const DEFAULT_SIGMA_EPS: f64 = 0.0005;
const LOG_2PI: f64 = 1.837_877_066_409_345_6;

pub const PARAM_NAMES: [&str; 11] = [
    "lambda",
    "kp11",
    "kp22",
    "kp33",
    "theta1",
    "theta2",
    "theta3",
    "sigma1",
    "sigma2",
    "sigma3",
    "sigma_eps",
];

pub const IS_POSITIVE_PARAM: [bool; 11] = [
    true, true, true, true, false, false, false, true, true, true, true,
];

pub const LOWER_BOUNDS: [f64; 11] = [
    0.05, 5e-4, 5e-4, 5e-4, -0.5, -0.5, -0.5, 1e-5, 1e-5, 1e-5, 1e-6,
];

pub const UPPER_BOUNDS: [f64; 11] = [5.0, 5.0, 5.0, 5.0, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.05];

#[derive(Clone, Copy, Debug)]
pub struct ModelSpec {
    pub afns: AfnsParams,
    pub sigma_eps: f64,
    pub use_afns_adj: bool,
}

impl ModelSpec {
    pub fn default_spec() -> Self {
        Self {
            afns: AfnsParams::default_params(),
            sigma_eps: DEFAULT_SIGMA_EPS,
            use_afns_adj: true,
        }
    }

    pub fn to_vec(&self) -> [f64; 11] {
        [
            self.afns.lambda,
            self.afns.kp[0][0],
            self.afns.kp[1][1],
            self.afns.kp[2][2],
            self.afns.theta_p[0],
            self.afns.theta_p[1],
            self.afns.theta_p[2],
            self.afns.sigma[0],
            self.afns.sigma[1],
            self.afns.sigma[2],
            self.sigma_eps,
        ]
    }

    pub fn from_vec(v: [f64; 11], use_afns_adj: bool) -> Self {
        Self {
            afns: AfnsParams {
                lambda: v[0],
                kp: [[v[1], 0.0, 0.0], [0.0, v[2], 0.0], [0.0, 0.0, v[3]]],
                theta_p: [v[4], v[5], v[6]],
                sigma: [v[7], v[8], v[9]],
            },
            sigma_eps: v[10],
            use_afns_adj,
        }
    }
}

pub fn pack_unconstrained(v: [f64; 11]) -> [f64; 11] {
    [
        v[0].ln(),
        v[1].ln(),
        v[2].ln(),
        v[3].ln(),
        v[4] / 0.02,
        v[5] / 0.02,
        v[6] / 0.02,
        (v[7] / 0.01).ln(),
        (v[8] / 0.01).ln(),
        (v[9] / 0.01).ln(),
        (v[10] / 0.0005).ln(),
    ]
}

pub fn unpack_unconstrained(u: &[f64]) -> [f64; 11] {
    [
        u[0].exp(),
        u[1].exp(),
        u[2].exp(),
        u[3].exp(),
        u[4] * 0.02,
        u[5] * 0.02,
        u[6] * 0.02,
        0.01 * u[7].exp(),
        0.01 * u[8].exp(),
        0.01 * u[9].exp(),
        0.0005 * u[10].exp(),
    ]
}

fn u_bounds() -> ([f64; 11], [f64; 11]) {
    let mut lo = pack_unconstrained(LOWER_BOUNDS);
    let mut hi = pack_unconstrained(UPPER_BOUNDS);
    for i in 4..=6 {
        lo[i] = LOWER_BOUNDS[i] / 0.02;
        hi[i] = UPPER_BOUNDS[i] / 0.02;
    }
    (lo, hi)
}

fn jacobian_diag(u: &[f64]) -> [f64; 11] {
    let v = unpack_unconstrained(u);
    [
        v[0], v[1], v[2], v[3], 0.02, 0.02, 0.02, v[7], v[8], v[9], v[10],
    ]
}

pub fn perturb_initial_guess(true_spec: ModelSpec, rg: &mut Rng) -> ModelSpec {
    let v = true_spec.to_vec();
    let mut init_v = [0.0; 11];
    for i in 0..11 {
        let mut mult = rg.uniform_range(0.5, 1.8);
        if (mult - 1.0).abs() < 0.15 {
            mult = if mult < 1.0 { 0.72 } else { 1.35 };
        }
        init_v[i] = v[i] * mult;
    }
    ModelSpec::from_vec(init_v, true_spec.use_afns_adj)
}

#[derive(Clone, Debug)]
pub struct FilterOutput {
    pub log_likelihood: f64,
    pub x_pred: Vec<[f64; 3]>,
    pub p_pred: Vec<[[f64; 3]; 3]>,
    pub x_filt: Vec<[f64; 3]>,
    pub p_filt: Vec<[[f64; 3]; 3]>,
    pub x_smooth: Vec<[f64; 3]>,
    pub p_smooth: Vec<[[f64; 3]; 3]>,
}

pub fn fast_log_likelihood(
    spec: &ModelSpec,
    taus: &[f64],
    dt: f64,
    yields: &[Vec<f64>],
) -> Result<f64, String> {
    let m = taus.len();
    let t_steps = yields.len();
    if m == 0 || m > 16 || t_steps == 0 || spec.sigma_eps <= 0.0 {
        return Err("fast_log_likelihood: invalid inputs".into());
    }

    let mut a_vec = [0.0; 16];
    let mut h_mat = [[0.0; 3]; 16];
    for (j, &tau) in taus.iter().enumerate() {
        h_mat[j] = spec.afns.yield_loadings(tau);
        if spec.use_afns_adj {
            a_vec[j] = spec.afns.yield_adjustment(tau).0;
        }
    }

    let phi = spec.afns.transition_matrix_closed(dt);
    let q_dt = spec.afns.transition_cov_closed(dt);
    let q_inf = spec.afns.stationary_cov();
    let theta = spec.afns.theta_p;
    let r_var = spec.sigma_eps * spec.sigma_eps;

    let mut x_pred = theta;
    let mut p_pred = q_inf;
    let mut prev_p_pred = [[0.0; 3]; 3];
    let mut x_filt = [0.0; 3];
    let mut p_filt = [[0.0; 3]; 3];

    let mut f_mat = [[0.0; 16]; 16];
    let mut l_mat = [[0.0; 16]; 16];
    let mut inv_diag_l = [0.0; 16];
    let mut hp = [[0.0; 3]; 16];
    let mut w_mat = [[0.0; 3]; 16];
    let mut v = [0.0; 16];
    let mut z = [0.0; 16];

    let mut total_ll = 0.0;
    let const_half_m = m as f64 * LOG_2PI;
    let mut log_det = 0.0;
    let mut steady = false;

    for t in 0..t_steps {
        if t > 0 {
            let dx0 = x_filt[0] - theta[0];
            let dx1 = x_filt[1] - theta[1];
            let dx2 = x_filt[2] - theta[2];
            for i in 0..3 {
                x_pred[i] = theta[i] + phi[i][0] * dx0 + phi[i][1] * dx1 + phi[i][2] * dx2;
            }
            if !steady {
                prev_p_pred = p_pred;
                let mut tmp = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        tmp[i][j] = phi[i][0] * p_filt[0][j]
                            + phi[i][1] * p_filt[1][j]
                            + phi[i][2] * p_filt[2][j];
                    }
                }
                for i in 0..3 {
                    for j in i..3 {
                        let val = q_dt[i][j]
                            + tmp[i][0] * phi[j][0]
                            + tmp[i][1] * phi[j][1]
                            + tmp[i][2] * phi[j][2];
                        p_pred[i][j] = val;
                        p_pred[j][i] = val;
                    }
                }
            }
        }

        let yt = &yields[t];
        for j in 0..m {
            let val = yt[j];
            if val.is_nan() || val.is_infinite() {
                return Err(format!("non-finite observation at ({t},{j})"));
            }
            let pred_y = a_vec[j]
                + h_mat[j][0] * x_pred[0]
                + h_mat[j][1] * x_pred[1]
                + h_mat[j][2] * x_pred[2];
            v[j] = val - pred_y;
        }

        if !steady {
            for j in 0..m {
                hp[j][0] = h_mat[j][0] * p_pred[0][0]
                    + h_mat[j][1] * p_pred[1][0]
                    + h_mat[j][2] * p_pred[2][0];
                hp[j][1] = h_mat[j][0] * p_pred[0][1]
                    + h_mat[j][1] * p_pred[1][1]
                    + h_mat[j][2] * p_pred[2][1];
                hp[j][2] = h_mat[j][0] * p_pred[0][2]
                    + h_mat[j][1] * p_pred[1][2]
                    + h_mat[j][2] * p_pred[2][2];
            }
            for i in 0..m {
                for j in 0..=i {
                    let mut s =
                        hp[i][0] * h_mat[j][0] + hp[i][1] * h_mat[j][1] + hp[i][2] * h_mat[j][2];
                    if i == j {
                        s += r_var;
                    }
                    f_mat[i][j] = s;
                    f_mat[j][i] = s;
                }
            }
            log_det = 0.0;
            for i in 0..m {
                for j in 0..=i {
                    let mut sum = f_mat[i][j];
                    for k in 0..j {
                        sum -= l_mat[i][k] * l_mat[j][k];
                    }
                    if i == j {
                        if sum <= 0.0 || sum.is_nan() || sum.is_infinite() {
                            return Err(format!("non-PD innovation covariance at t={t}"));
                        }
                        let diag = sum.sqrt();
                        l_mat[i][j] = diag;
                        inv_diag_l[i] = 1.0 / diag;
                        log_det += diag.ln();
                    } else {
                        l_mat[i][j] = sum * inv_diag_l[j];
                    }
                }
            }
            log_det *= 2.0;
        }

        let mut quad = 0.0;
        for i in 0..m {
            let mut sum = v[i];
            for j in 0..i {
                sum -= l_mat[i][j] * z[j];
            }
            let zi = sum * inv_diag_l[i];
            z[i] = zi;
            quad += zi * zi;
        }
        total_ll += -0.5 * (const_half_m + log_det + quad);

        if t + 1 == t_steps {
            break;
        }

        if !steady {
            for col in 0..3 {
                for i in 0..m {
                    let mut sum = hp[i][col];
                    for j in 0..i {
                        sum -= l_mat[i][j] * w_mat[j][col];
                    }
                    w_mat[i][col] = sum * inv_diag_l[i];
                }
                for i in (0..m).rev() {
                    let mut sum = w_mat[i][col];
                    for j in (i + 1)..m {
                        sum -= l_mat[j][i] * w_mat[j][col];
                    }
                    w_mat[i][col] = sum * inv_diag_l[i];
                }
            }

            let mut ikh = [[0.0; 3]; 3];
            for i in 0..3 {
                for k in 0..3 {
                    let mut kh = 0.0;
                    for j in 0..m {
                        kh += w_mat[j][i] * h_mat[j][k];
                    }
                    ikh[i][k] = if i == k { 1.0 - kh } else { -kh };
                }
            }
            let mut tmp_p = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    tmp_p[i][j] = ikh[i][0] * p_pred[0][j]
                        + ikh[i][1] * p_pred[1][j]
                        + ikh[i][2] * p_pred[2][j];
                }
            }
            for i in 0..3 {
                for j in i..3 {
                    let mut val =
                        tmp_p[i][0] * ikh[j][0] + tmp_p[i][1] * ikh[j][1] + tmp_p[i][2] * ikh[j][2];
                    let mut kk = 0.0;
                    for r in 0..m {
                        kk += w_mat[r][i] * w_mat[r][j];
                    }
                    val += r_var * kk;
                    p_filt[i][j] = val;
                    p_filt[j][i] = val;
                }
            }

            if t >= 5 {
                let mut max_diff = 0.0f64;
                let mut max_scale = 0.0f64;
                for i in 0..3 {
                    max_scale = max_scale.max(p_pred[i][i].abs());
                    for j in 0..3 {
                        max_diff = max_diff.max((p_pred[i][j] - prev_p_pred[i][j]).abs());
                    }
                }
                if max_diff <= 1e-15 * max_scale.max(1e-18) {
                    steady = true;
                }
            }
        }

        for i in 0..3 {
            let mut upd = 0.0;
            for j in 0..m {
                upd += w_mat[j][i] * v[j];
            }
            x_filt[i] = x_pred[i] + upd;
        }
    }

    Ok(total_ll)
}

pub fn filter_and_smooth(
    spec: &ModelSpec,
    taus: &[f64],
    dt: f64,
    yields: &[Vec<f64>],
) -> Result<FilterOutput, String> {
    let m = taus.len();
    let t_steps = yields.len();
    if m == 0 || t_steps == 0 || spec.sigma_eps <= 0.0 {
        return Err("filter_and_smooth: invalid inputs".into());
    }

    let mut a_vec = vec![0.0; m];
    let mut h_mat = Vec::with_capacity(m);
    for (j, &tau) in taus.iter().enumerate() {
        h_mat.push(spec.afns.yield_loadings(tau));
        if spec.use_afns_adj {
            a_vec[j] = spec.afns.yield_adjustment(tau).0;
        }
    }

    let phi = spec.afns.transition_matrix_closed(dt);
    let q_dt = spec.afns.transition_cov_closed(dt);
    let q_inf = spec.afns.stationary_cov();
    let theta = spec.afns.theta_p;
    let r_var = spec.sigma_eps * spec.sigma_eps;

    let mut x_pred = vec![[0.0; 3]; t_steps];
    let mut p_pred = vec![[[0.0; 3]; 3]; t_steps];
    let mut x_filt = vec![[0.0; 3]; t_steps];
    let mut p_filt = vec![[[0.0; 3]; 3]; t_steps];

    let mut total_ll = 0.0;
    for t in 0..t_steps {
        if t == 0 {
            x_pred[0] = theta;
            p_pred[0] = q_inf;
        } else {
            let dx = [
                x_filt[t - 1][0] - theta[0],
                x_filt[t - 1][1] - theta[1],
                x_filt[t - 1][2] - theta[2],
            ];
            for i in 0..3 {
                x_pred[t][i] = theta[i] + phi[i][0] * dx[0] + phi[i][1] * dx[1] + phi[i][2] * dx[2];
            }
            let mut tmp = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    tmp[i][j] = phi[i][0] * p_filt[t - 1][0][j]
                        + phi[i][1] * p_filt[t - 1][1][j]
                        + phi[i][2] * p_filt[t - 1][2][j];
                }
            }
            for i in 0..3 {
                for j in i..3 {
                    let val = q_dt[i][j]
                        + tmp[i][0] * phi[j][0]
                        + tmp[i][1] * phi[j][1]
                        + tmp[i][2] * phi[j][2];
                    p_pred[t][i][j] = val;
                    p_pred[t][j][i] = val;
                }
            }
        }

        let mut v = vec![0.0; m];
        for j in 0..m {
            let pred_y = a_vec[j]
                + h_mat[j][0] * x_pred[t][0]
                + h_mat[j][1] * x_pred[t][1]
                + h_mat[j][2] * x_pred[t][2];
            v[j] = yields[t][j] - pred_y;
        }

        let mut hp = vec![[0.0; 3]; m];
        for j in 0..m {
            for k in 0..3 {
                hp[j][k] = h_mat[j][0] * p_pred[t][0][k]
                    + h_mat[j][1] * p_pred[t][1][k]
                    + h_mat[j][2] * p_pred[t][2][k];
            }
        }
        let mut f_mat = linalg::zeros(m, m);
        for i in 0..m {
            for j in 0..=i {
                let mut s =
                    hp[i][0] * h_mat[j][0] + hp[i][1] * h_mat[j][1] + hp[i][2] * h_mat[j][2];
                if i == j {
                    s += r_var;
                }
                f_mat[i][j] = s;
                f_mat[j][i] = s;
            }
        }
        let l_chol = linalg::cholesky(&f_mat)?;
        let mut log_det = 0.0;
        for i in 0..m {
            log_det += l_chol[i][i].ln();
        }
        log_det *= 2.0;
        let z = {
            let mut z_vec = vec![0.0; m];
            for i in 0..m {
                let mut s = v[i];
                for k in 0..i {
                    s -= l_chol[i][k] * z_vec[k];
                }
                z_vec[i] = s / l_chol[i][i];
            }
            z_vec
        };
        let quad: f64 = z.iter().map(|&zi| zi * zi).sum();
        total_ll += -0.5 * (m as f64 * LOG_2PI + log_det + quad);

        let mut k_t = vec![[0.0; 3]; m];
        for col in 0..3 {
            let rhs: Vec<f64> = (0..m).map(|i| hp[i][col]).collect();
            let sol = linalg::solve_cholesky(&l_chol, &rhs);
            for i in 0..m {
                k_t[i][col] = sol[i];
            }
        }
        for i in 0..3 {
            let upd: f64 = (0..m).map(|j| k_t[j][i] * v[j]).sum();
            x_filt[t][i] = x_pred[t][i] + upd;
        }
        let mut ikh = [[0.0; 3]; 3];
        for i in 0..3 {
            for k in 0..3 {
                let kh: f64 = (0..m).map(|j| k_t[j][i] * h_mat[j][k]).sum();
                ikh[i][k] = if i == k { 1.0 - kh } else { -kh };
            }
        }
        let mut tmp_p = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                tmp_p[i][j] = ikh[i][0] * p_pred[t][0][j]
                    + ikh[i][1] * p_pred[t][1][j]
                    + ikh[i][2] * p_pred[t][2][j];
            }
        }
        for i in 0..3 {
            for j in i..3 {
                let mut val =
                    tmp_p[i][0] * ikh[j][0] + tmp_p[i][1] * ikh[j][1] + tmp_p[i][2] * ikh[j][2];
                let kk: f64 = (0..m).map(|r| k_t[r][i] * k_t[r][j]).sum();
                val += r_var * kk;
                p_filt[t][i][j] = val;
                p_filt[t][j][i] = val;
            }
        }
    }

    let mut x_smooth = vec![[0.0; 3]; t_steps];
    let mut p_smooth = vec![[[0.0; 3]; 3]; t_steps];
    x_smooth[t_steps - 1] = x_filt[t_steps - 1];
    p_smooth[t_steps - 1] = p_filt[t_steps - 1];

    for t in (0..(t_steps - 1)).rev() {
        let p_next_mat: Vec<Vec<f64>> = (0..3).map(|i| p_pred[t + 1][i].to_vec()).collect();
        let l_p = linalg::cholesky(&p_next_mat)?;
        let mut pf_phi_t = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                pf_phi_t[i][j] = p_filt[t][i][0] * phi[j][0]
                    + p_filt[t][i][1] * phi[j][1]
                    + p_filt[t][i][2] * phi[j][2];
            }
        }
        let mut j_t = [[0.0; 3]; 3];
        for row in 0..3 {
            let rhs = vec![pf_phi_t[row][0], pf_phi_t[row][1], pf_phi_t[row][2]];
            let sol = linalg::solve_cholesky(&l_p, &rhs);
            j_t[row].copy_from_slice(&sol);
        }
        let dx_s = [
            x_smooth[t + 1][0] - x_pred[t + 1][0],
            x_smooth[t + 1][1] - x_pred[t + 1][1],
            x_smooth[t + 1][2] - x_pred[t + 1][2],
        ];
        for i in 0..3 {
            x_smooth[t][i] =
                x_filt[t][i] + j_t[i][0] * dx_s[0] + j_t[i][1] * dx_s[1] + j_t[i][2] * dx_s[2];
        }
        let mut dp = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                dp[i][j] = p_smooth[t + 1][i][j] - p_pred[t + 1][i][j];
            }
        }
        let mut tmp = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                tmp[i][j] = j_t[i][0] * dp[0][j] + j_t[i][1] * dp[1][j] + j_t[i][2] * dp[2][j];
            }
        }
        for i in 0..3 {
            for j in i..3 {
                let val = p_filt[t][i][j]
                    + tmp[i][0] * j_t[j][0]
                    + tmp[i][1] * j_t[j][1]
                    + tmp[i][2] * j_t[j][2];
                p_smooth[t][i][j] = val;
                p_smooth[t][j][i] = val;
            }
        }
    }

    Ok(FilterOutput {
        log_likelihood: total_ll,
        x_pred,
        p_pred,
        x_filt,
        p_filt,
        x_smooth,
        p_smooth,
    })
}

#[derive(Clone, Debug)]
pub struct JointGaussianAnchor {
    pub t_steps: usize,
    pub num_tenors: usize,
    pub kalman_log_lik: f64,
    pub fast_log_lik: f64,
    pub joint_log_lik: f64,
    pub log_lik_abs_diff: f64,
    pub fast_log_lik_abs_diff: f64,
    pub max_smooth_mean_abs_diff: f64,
    pub max_smooth_cov_abs_diff: f64,
}

pub fn verify_joint_gaussian_anchor(
    spec: &ModelSpec,
    taus: &[f64],
    dt: f64,
    yields: &[Vec<f64>],
) -> Result<JointGaussianAnchor, String> {
    let t_steps = yields.len();
    let m = taus.len();
    let nx = 3 * t_steps;
    let ny = m * t_steps;

    let phi = spec.afns.transition_matrix_closed(dt);
    let q_dt = spec.afns.transition_cov_closed(dt);
    let q_inf = spec.afns.stationary_cov();
    let theta = spec.afns.theta_p;
    let r_var = spec.sigma_eps * spec.sigma_eps;

    let mut sigma_x = linalg::zeros(nx, nx);
    for t in 0..t_steps {
        for s in t..t_steps {
            let mut blk = if t == 0 { q_inf } else { [[0.0; 3]; 3] };
            if t > 0 {
                // Sigma_{t,t} = Phi Sigma_{t-1,t-1} Phi^T + Q_dt (stationary so = Q_inf, keep exact recursion)
                blk = q_inf;
            }
            if s > t {
                // Cov(X_s, X_t) = Phi^(s-t) Q_inf
                let mut p_pow = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
                for _ in 0..(s - t) {
                    let mut next = [[0.0; 3]; 3];
                    for i in 0..3 {
                        for j in 0..3 {
                            next[i][j] = phi[i][0] * p_pow[0][j]
                                + phi[i][1] * p_pow[1][j]
                                + phi[i][2] * p_pow[2][j];
                        }
                    }
                    p_pow = next;
                }
                for i in 0..3 {
                    for j in 0..3 {
                        let v = p_pow[i][0] * q_inf[0][j]
                            + p_pow[i][1] * q_inf[1][j]
                            + p_pow[i][2] * q_inf[2][j];
                        sigma_x[3 * s + i][3 * t + j] = v;
                        sigma_x[3 * t + j][3 * s + i] = v;
                    }
                }
            } else {
                for i in 0..3 {
                    for j in 0..3 {
                        sigma_x[3 * t + i][3 * t + j] = blk[i][j];
                    }
                }
            }
        }
    }

    let mut h_big = linalg::zeros(ny, nx);
    let mut a_big = vec![0.0; ny];
    let mut y_big = vec![0.0; ny];
    for t in 0..t_steps {
        for j in 0..m {
            let row = t * m + j;
            let h = spec.afns.yield_loadings(taus[j]);
            for k in 0..3 {
                h_big[row][3 * t + k] = h[k];
            }
            if spec.use_afns_adj {
                a_big[row] = spec.afns.yield_adjustment(taus[j]).0;
            }
            y_big[row] = yields[t][j];
        }
    }

    let mut mu_x = vec![0.0; nx];
    for t in 0..t_steps {
        for k in 0..3 {
            mu_x[3 * t + k] = theta[k];
        }
    }
    let h_mu = linalg::mat_vec(&h_big, &mu_x);
    let mut diff_y = vec![0.0; ny];
    for i in 0..ny {
        diff_y[i] = y_big[i] - (a_big[i] + h_mu[i]);
    }

    let sigma_xy = linalg::mat_mul(&sigma_x, &linalg::transpose(&h_big));
    let mut sigma_y = linalg::mat_mul(&h_big, &sigma_xy);
    for i in 0..ny {
        sigma_y[i][i] += r_var;
    }

    let l_y = linalg::cholesky(&sigma_y)?;
    let mut log_det_y = 0.0;
    for i in 0..ny {
        log_det_y += l_y[i][i].ln();
    }
    log_det_y *= 2.0;
    let mut z_y = vec![0.0; ny];
    for i in 0..ny {
        let mut s = diff_y[i];
        for k in 0..i {
            s -= l_y[i][k] * z_y[k];
        }
        z_y[i] = s / l_y[i][i];
    }
    let quad_y: f64 = z_y.iter().map(|&x| x * x).sum();
    let joint_ll = -0.5 * (ny as f64 * LOG_2PI + log_det_y + quad_y);

    let alpha_y = linalg::solve_cholesky(&l_y, &diff_y);
    let corr_mean = linalg::mat_vec(&sigma_xy, &alpha_y);
    let cond_mean: Vec<f64> = (0..nx).map(|i| mu_x[i] + corr_mean[i]).collect();

    let mut k_mat = linalg::zeros(ny, nx);
    for col in 0..nx {
        let rhs: Vec<f64> = (0..ny).map(|r| sigma_xy[col][r]).collect();
        let sol = linalg::solve_cholesky(&l_y, &rhs);
        for r in 0..ny {
            k_mat[r][col] = sol[r];
        }
    }
    let gain_cov = linalg::mat_mul(&sigma_xy, &k_mat);
    let mut cov_post = linalg::zeros(nx, nx);
    for i in 0..nx {
        for j in 0..nx {
            cov_post[i][j] = sigma_x[i][j] - gain_cov[i][j];
        }
    }

    let filt_out = filter_and_smooth(spec, taus, dt, yields)?;
    let fast_ll = fast_log_likelihood(spec, taus, dt, yields)?;

    let mut max_mean_diff = 0.0f64;
    let mut max_cov_diff = 0.0f64;
    for t in 0..t_steps {
        for i in 0..3 {
            max_mean_diff =
                max_mean_diff.max((filt_out.x_smooth[t][i] - cond_mean[3 * t + i]).abs());
            for j in 0..3 {
                max_cov_diff = max_cov_diff
                    .max((filt_out.p_smooth[t][i][j] - cov_post[3 * t + i][3 * t + j]).abs());
            }
        }
    }

    let _ = q_dt;
    Ok(JointGaussianAnchor {
        t_steps,
        num_tenors: m,
        kalman_log_lik: filt_out.log_likelihood,
        fast_log_lik: fast_ll,
        joint_log_lik: joint_ll,
        log_lik_abs_diff: (filt_out.log_likelihood - joint_ll).abs(),
        fast_log_lik_abs_diff: (fast_ll - joint_ll).abs(),
        max_smooth_mean_abs_diff: max_mean_diff,
        max_smooth_cov_abs_diff: max_cov_diff,
    })
}

#[derive(Clone, Debug)]
pub struct StageDiagnostic {
    pub stage: usize,
    pub iterations: usize,
    pub func_evals: usize,
    pub log_lik: f64,
    pub raw_grad_norm: f64,
    pub grad_norm: f64,
    pub scaled_grad_norm: f64,
    pub converged: bool,
    pub hit_bound: bool,
    pub termination_reason: String,
    pub min_alpha_tested: f64,
    pub last_step_norm: f64,
    pub candidate_bit_identical: bool,
    pub last_obj_change: f64,
    pub obj_bit_identical: bool,
    pub line_search_failure_cause: String,
}

#[derive(Clone, Debug)]
pub struct MleResult {
    pub finite: bool,
    pub init_params: [f64; 11],
    pub est_params: [f64; 11],
    pub true_params: [f64; 11],
    pub init_log_lik: f64,
    pub true_log_lik: f64,
    pub est_log_lik: f64,
    pub grad_rel_err_init: f64,
    pub grad_rel_err_true: f64,
    pub grad_abs_diff_opt: f64,
    pub raw_grad_norm_opt: f64,
    pub grad_norm_opt: f64,
    pub scaled_grad_norm_opt: f64,
    pub hessian_se: [f64; 11],
    pub hessian_min_eig: f64,
    pub hessian_regularized: bool,
    pub iterations: usize,
    pub best_stage_iterations: usize,
    pub cumulative_iterations: usize,
    pub best_stage: usize,
    pub stage1: StageDiagnostic,
    pub stage2: StageDiagnostic,
    pub converged: bool,
    pub termination_reason: String,
    pub hit_bound: bool,
    pub hit_bound_by_param: [bool; 11],
    pub min_alpha_tested: f64,
    pub last_step_norm: f64,
    pub candidate_bit_identical: bool,
    pub last_obj_change: f64,
    pub obj_bit_identical: bool,
    pub line_search_failure_cause: String,
}

fn make_stage_diag(stage: usize, res: &optim::BfgsResult) -> StageDiagnostic {
    StageDiagnostic {
        stage,
        iterations: res.iterations,
        func_evals: res.func_evals,
        log_lik: -res.fun,
        raw_grad_norm: res.raw_grad_norm,
        grad_norm: res.grad_norm,
        scaled_grad_norm: res.scaled_grad_norm,
        converged: res.converged,
        hit_bound: res.hit_bound,
        termination_reason: res.termination_reason.clone(),
        min_alpha_tested: res.min_alpha_tested,
        last_step_norm: res.last_step_norm,
        candidate_bit_identical: res.candidate_bit_identical,
        last_obj_change: res.last_obj_change,
        obj_bit_identical: res.obj_bit_identical,
        line_search_failure_cause: res.line_search_failure_cause.clone(),
    }
}

pub fn estimate_mle(
    init_spec: ModelSpec,
    true_spec: ModelSpec,
    taus: &[f64],
    dt: f64,
    yields: &[Vec<f64>],
    compute_hessian: bool,
) -> MleResult {
    let use_afns = init_spec.use_afns_adj;
    let mut neg_ll = |u: &[f64]| -> f64 {
        let v = unpack_unconstrained(u);
        let sp = ModelSpec::from_vec(v, use_afns);
        match fast_log_likelihood(&sp, taus, dt, yields) {
            Ok(ll) if ll.is_finite() => -ll,
            _ => 1e30,
        }
    };

    let u_init = pack_unconstrained(init_spec.to_vec());
    let u_true = pack_unconstrained(true_spec.to_vec());
    let (lo, hi) = u_bounds();

    let (mut grad_err_init, mut grad_err_true) = (0.0, 0.0);
    if compute_hessian {
        grad_err_init = optim::verify_grad_vs_richardson(&mut neg_ll, &u_init, None).2;
        grad_err_true = optim::verify_grad_vs_richardson(&mut neg_ll, &u_true, None).2;
    }

    let init_ll = fast_log_likelihood(&init_spec, taus, dt, yields).unwrap_or(f64::NAN);
    let true_ll = fast_log_likelihood(&true_spec, taus, dt, yields).unwrap_or(f64::NAN);

    let res1 = optim::minimize_bfgs(&mut neg_ll, &u_init, &lo, &hi, 200, 1e-6);
    let res2 = optim::minimize_bfgs(&mut neg_ll, &res1.x, &lo, &hi, 150, 1e-6);
    let (best, best_stage) = if res2.fun < res1.fun {
        (&res2, 2usize)
    } else {
        (&res1, 1usize)
    };
    let cumulative_iterations = res1.iterations + res2.iterations;
    let best_stage_iterations = best.iterations;

    let est_v = unpack_unconstrained(&best.x);
    let est_ll = -best.fun;
    let finite = best.termination_reason != optim::TERM_INVALID_OBJECTIVE
        && best.fun < 1e29
        && est_ll.is_finite()
        && est_v.iter().all(|x| x.is_finite());

    let mut grad_abs_diff_opt = 0.0;
    let mut ses = [0.0; 11];
    let mut min_eig = 0.0;
    let mut reg = false;

    if finite && compute_hessian {
        let (cd_opt, rich_opt, _) = optim::verify_grad_vs_richardson(&mut neg_ll, &best.x, None);
        grad_abs_diff_opt = cd_opt
            .iter()
            .zip(rich_opt.iter())
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt();

        let hess_u = optim::central_diff_hessian(&mut neg_ll, &best.x, None);
        if let Ok((vals, vecs, _)) = linalg::symmetric_eigen_jacobi(&hess_u, 80, 1e-13) {
            if vals.len() == 11 {
                min_eig = vals[10];
                let mut inv_hess_u = [[0.0; 11]; 11];
                const EIG_FLOOR: f64 = 1e-6;
                for k in 0..11 {
                    let mut ev = vals[k];
                    if ev < EIG_FLOOR {
                        ev = EIG_FLOOR;
                        reg = true;
                    }
                    let inv_ev = 1.0 / ev;
                    for i in 0..11 {
                        for j in 0..11 {
                            inv_hess_u[i][j] += inv_ev * vecs[k][i] * vecs[k][j];
                        }
                    }
                }
                let jac = jacobian_diag(&best.x);
                for i in 0..11 {
                    let var_psi = jac[i] * jac[i] * inv_hess_u[i][i];
                    if var_psi > 0.0 {
                        ses[i] = var_psi.sqrt();
                    }
                }
            }
        }
    }

    let mut hit_by_param = [false; 11];
    for (i, &h) in best.hit_bound_by_param.iter().take(11).enumerate() {
        hit_by_param[i] = h;
    }

    MleResult {
        finite,
        init_params: init_spec.to_vec(),
        est_params: est_v,
        true_params: true_spec.to_vec(),
        init_log_lik: init_ll,
        true_log_lik: true_ll,
        est_log_lik: est_ll,
        grad_rel_err_init: grad_err_init,
        grad_rel_err_true: grad_err_true,
        grad_abs_diff_opt,
        raw_grad_norm_opt: best.raw_grad_norm,
        grad_norm_opt: best.grad_norm,
        scaled_grad_norm_opt: best.scaled_grad_norm,
        hessian_se: ses,
        hessian_min_eig: min_eig,
        hessian_regularized: reg,
        iterations: best_stage_iterations,
        best_stage_iterations,
        cumulative_iterations,
        best_stage,
        stage1: make_stage_diag(1, &res1),
        stage2: make_stage_diag(2, &res2),
        converged: best.converged,
        termination_reason: best.termination_reason.clone(),
        hit_bound: best.hit_bound,
        hit_bound_by_param: hit_by_param,
        min_alpha_tested: best.min_alpha_tested,
        last_step_norm: best.last_step_norm,
        candidate_bit_identical: best.candidate_bit_identical,
        last_obj_change: best.last_obj_change,
        obj_bit_identical: best.obj_bit_identical,
        line_search_failure_cause: best.line_search_failure_cause.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::afns::simulate_stationary_path;
    use crate::rng::SEED_DAY004A;

    #[test]
    fn test_joint_gaussian_anchor_t5() {
        let spec = ModelSpec::default_spec();
        let taus = [0.5, 2.0, 10.0];
        let mut rg = Rng::new(SEED_DAY004A);
        let (_, y) =
            simulate_stationary_path(&spec.afns, &taus, 5, WEEKLY_DT, spec.sigma_eps, &mut rg);
        let anc = verify_joint_gaussian_anchor(&spec, &taus, WEEKLY_DT, &y).unwrap();
        assert!(anc.log_lik_abs_diff < 1e-9);
        assert!(anc.fast_log_lik_abs_diff < 1e-9);
        assert!(anc.max_smooth_mean_abs_diff < 1e-9);
        assert!(anc.max_smooth_cov_abs_diff < 1e-9);
    }
}
