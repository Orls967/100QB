//! Finite-difference & Richardson gradients, Hessian, and box-constrained BFGS with line-search diagnostics.

#[derive(Clone, Debug)]
pub struct BfgsResult {
    pub x: Vec<f64>,
    pub fun: f64,
    pub raw_grad_norm: f64,
    pub grad_norm: f64,
    pub scaled_grad_norm: f64,
    pub iterations: usize,
    pub func_evals: usize,
    pub converged: bool,
    pub hit_bound: bool,
    pub hit_bound_by_param: Vec<bool>,
    pub termination_reason: String,
    pub min_alpha_tested: f64,
    pub last_step_norm: f64,
    pub candidate_bit_identical: bool,
    pub last_obj_change: f64,
    pub obj_bit_identical: bool,
    pub line_search_failure_cause: String,
}

pub const TERM_GRADIENT_TOL: &str = "gradient_tolerance";
pub const TERM_BOUNDARY_STATIONARY: &str = "boundary_stationary";
pub const TERM_LINE_SEARCH_FAILED: &str = "line_search_failed";
pub const TERM_MAX_ITERATIONS: &str = "max_iterations";
pub const TERM_INVALID_OBJECTIVE: &str = "invalid_objective_or_gradient";

pub const LS_CAUSE_NONE: &str = "none";
pub const LS_CAUSE_ARMIJO_NOT_MET: &str = "armijo_not_met";
pub const LS_CAUSE_STEP_UNDERFLOW: &str = "step_underflow_plateau";
pub const LS_CAUSE_INVALID_OBJECTIVE: &str = "invalid_objective_or_gradient";

pub fn central_diff_grad<F>(f: &mut F, x: &[f64], steps: Option<&[f64]>) -> Vec<f64>
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x.len();
    let mut g = vec![0.0; n];
    let mut xp = x.to_vec();
    for i in 0..n {
        let h = if let Some(st) = steps {
            if st.len() == n && st[i] > 0.0 {
                st[i]
            } else {
                1e-5 * x[i].abs().max(1.0)
            }
        } else {
            1e-5 * x[i].abs().max(1.0)
        };
        let orig = xp[i];
        xp[i] = orig + h;
        let fp = f(&xp);
        xp[i] = orig - h;
        let fm = f(&xp);
        xp[i] = orig;
        g[i] = (fp - fm) / (2.0 * h);
    }
    g
}

pub fn richardson_grad<F>(f: &mut F, x: &[f64], steps: Option<&[f64]>) -> Vec<f64>
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x.len();
    let mut h_coarse = vec![0.0; n];
    let mut h_fine = vec![0.0; n];
    for i in 0..n {
        let h = if let Some(st) = steps {
            if st.len() == n && st[i] > 0.0 {
                st[i]
            } else {
                2e-5 * x[i].abs().max(1.0)
            }
        } else {
            2e-5 * x[i].abs().max(1.0)
        };
        h_coarse[i] = h;
        h_fine[i] = 0.5 * h;
    }
    let g1 = central_diff_grad(f, x, Some(&h_coarse));
    let g2 = central_diff_grad(f, x, Some(&h_fine));
    (0..n).map(|i| (4.0 * g2[i] - g1[i]) / 3.0).collect()
}

pub fn verify_grad_vs_richardson<F>(
    f: &mut F,
    x: &[f64],
    steps: Option<&[f64]>,
) -> (Vec<f64>, Vec<f64>, f64)
where
    F: FnMut(&[f64]) -> f64,
{
    let cd = central_diff_grad(f, x, steps);
    let rich = richardson_grad(f, x, steps);
    let mut diff_sq = 0.0;
    let mut norm_sq = 0.0;
    for i in 0..x.len() {
        let d = cd[i] - rich[i];
        diff_sq += d * d;
        norm_sq += rich[i] * rich[i];
    }
    let rel_err = diff_sq.sqrt() / norm_sq.sqrt().max(1e-12);
    (cd, rich, rel_err)
}

pub fn central_diff_hessian<F>(f: &mut F, x: &[f64], steps: Option<&[f64]>) -> Vec<Vec<f64>>
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x.len();
    let mut h = vec![0.0; n];
    for i in 0..n {
        h[i] = if let Some(st) = steps {
            if st.len() == n && st[i] > 0.0 {
                st[i]
            } else {
                1e-4 * x[i].abs().max(1.0)
            }
        } else {
            1e-4 * x[i].abs().max(1.0)
        };
    }
    let f0 = f(x);
    let mut hess = vec![vec![0.0; n]; n];
    let mut xp = x.to_vec();
    for i in 0..n {
        let orig = xp[i];
        xp[i] = orig + h[i];
        let fp = f(&xp);
        xp[i] = orig - h[i];
        let fm = f(&xp);
        xp[i] = orig;
        hess[i][i] = (fp - 2.0 * f0 + fm) / (h[i] * h[i]);
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let oi = xp[i];
            let oj = xp[j];
            xp[i] = oi + h[i];
            xp[j] = oj + h[j];
            let fpp = f(&xp);
            xp[j] = oj - h[j];
            let fpm = f(&xp);
            xp[i] = oi - h[i];
            xp[j] = oj + h[j];
            let fmp = f(&xp);
            xp[j] = oj - h[j];
            let fmm = f(&xp);
            xp[i] = oi;
            xp[j] = oj;
            let val = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j]);
            hess[i][j] = val;
            hess[j][i] = val;
        }
    }
    hess
}

fn clamp_vec(x: &mut [f64], lo: &[f64], hi: &[f64]) {
    if lo.len() == x.len() && hi.len() == x.len() {
        for i in 0..x.len() {
            if x[i] < lo[i] {
                x[i] = lo[i];
            } else if x[i] > hi[i] {
                x[i] = hi[i];
            }
        }
    }
}

fn is_invalid_state(val: f64, g: &[f64]) -> bool {
    if val.is_nan() || val.is_infinite() || val >= 1e29 {
        return true;
    }
    g.iter().any(|&gi| gi.is_nan() || gi.is_infinite())
}

pub fn minimize_bfgs<F>(
    mut f: F,
    x0: &[f64],
    lo: &[f64],
    hi: &[f64],
    max_iter: usize,
    grad_tol: f64,
) -> BfgsResult
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x0.len();
    let mut x = x0.to_vec();
    clamp_vec(&mut x, lo, hi);

    let mut evals = 0usize;
    let mut eval_f = |v: &[f64]| -> f64 {
        evals += 1;
        f(v)
    };

    let raw_inf_norm = |v: &[f64]| -> f64 {
        let mut m = 0.0;
        for &a in v {
            if a.is_nan() {
                return f64::NAN;
            }
            if a.abs() > m {
                m = a.abs();
            }
        }
        m
    };

    let check_bounds = |vec: &[f64]| -> (bool, Vec<bool>) {
        let mut hit = false;
        let mut by_param = vec![false; n];
        if lo.len() == n && hi.len() == n {
            for i in 0..n {
                if vec[i] <= lo[i] + 1e-9 || vec[i] >= hi[i] - 1e-9 {
                    hit = true;
                    by_param[i] = true;
                }
            }
        }
        (hit, by_param)
    };

    let inf_norm = |v: &[f64], cur_x: &[f64]| -> f64 {
        let mut m = 0.0;
        for (i, &a) in v.iter().enumerate() {
            if a.is_nan() {
                return f64::NAN;
            }
            if lo.len() == v.len() && hi.len() == v.len() {
                if cur_x[i] <= lo[i] + 1e-12 && a > 0.0 {
                    continue;
                }
                if cur_x[i] >= hi[i] - 1e-12 && a < 0.0 {
                    continue;
                }
            }
            if a.abs() > m {
                m = a.abs();
            }
        }
        m
    };

    let mut fx = eval_f(&x);
    let mut g = central_diff_grad(&mut eval_f, &x, None);

    let mut h_mat = vec![vec![0.0; n]; n];
    for i in 0..n {
        h_mat[i][i] = 1.0;
    }

    let mut term_reason = String::new();
    let mut ls_cause = LS_CAUSE_NONE.to_string();
    let mut min_alpha_tested = 0.0;
    let mut last_step_norm = 0.0;
    let mut cand_bit_ident = false;
    let mut last_obj_change = 0.0;
    let mut obj_bit_ident = false;
    let mut iters = 0usize;

    let mut x_cand = vec![0.0; n];
    let mut d = vec![0.0; n];
    let mut s = vec![0.0; n];
    let mut y = vec![0.0; n];
    let mut hy = vec![0.0; n];

    if is_invalid_state(fx, &g) {
        term_reason = TERM_INVALID_OBJECTIVE.to_string();
        ls_cause = LS_CAUSE_INVALID_OBJECTIVE.to_string();
    }

    for iter in 0..max_iter {
        if !term_reason.is_empty() {
            break;
        }
        let gn = inf_norm(&g, &x);
        if gn <= grad_tol {
            let (hb, _) = check_bounds(&x);
            term_reason = if hb {
                TERM_BOUNDARY_STATIONARY.to_string()
            } else {
                TERM_GRADIENT_TOL.to_string()
            };
            break;
        }
        iters = iter + 1;

        let mut slope = 0.0;
        for i in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                sum -= h_mat[i][j] * g[j];
            }
            d[i] = sum;
            slope += g[i] * d[i];
        }
        if slope >= 0.0 || slope.is_nan() {
            for i in 0..n {
                for j in 0..n {
                    h_mat[i][j] = 0.0;
                }
                h_mat[i][i] = 1.0;
                d[i] = -g[i];
            }
        }

        let mut alpha = 1.0;
        let mut accepted = false;
        let mut f_cand = fx;
        let mut alpha_last = 1.0;
        let mut any_finite_cand = false;
        const C1: f64 = 1e-4;
        for _ in 0..35 {
            alpha_last = alpha;
            for i in 0..n {
                x_cand[i] = x[i] + alpha * d[i];
            }
            clamp_vec(&mut x_cand, lo, hi);
            let mut dir_deriv = 0.0;
            for i in 0..n {
                dir_deriv += g[i] * (x_cand[i] - x[i]);
            }
            f_cand = eval_f(&x_cand);
            if !f_cand.is_nan() && !f_cand.is_infinite() && f_cand < 1e29 {
                any_finite_cand = true;
                if f_cand <= fx + C1 * dir_deriv {
                    accepted = true;
                    break;
                }
            }
            alpha *= 0.5;
        }
        min_alpha_tested = alpha_last;
        if !accepted {
            term_reason = TERM_LINE_SEARCH_FAILED.to_string();
            let mut step_sq = 0.0;
            let mut bit_eq = true;
            for i in 0..n {
                let diff = x_cand[i] - x[i];
                step_sq += diff * diff;
                if x_cand[i].to_bits() != x[i].to_bits() {
                    bit_eq = false;
                }
            }
            last_step_norm = step_sq.sqrt();
            cand_bit_ident = bit_eq;
            if !any_finite_cand || f_cand.is_nan() || f_cand.is_infinite() || f_cand >= 1e29 {
                last_obj_change = 0.0;
                obj_bit_ident = false;
                ls_cause = LS_CAUSE_INVALID_OBJECTIVE.to_string();
            } else {
                last_obj_change = f_cand - fx;
                obj_bit_ident = f_cand.to_bits() == fx.to_bits();
                ls_cause = LS_CAUSE_ARMIJO_NOT_MET.to_string();
            }
            break;
        }

        let mut step_norm_sq = 0.0;
        let mut bit_eq = true;
        for i in 0..n {
            s[i] = x_cand[i] - x[i];
            if x_cand[i].to_bits() != x[i].to_bits() {
                bit_eq = false;
            }
            x[i] = x_cand[i];
            step_norm_sq += s[i] * s[i];
        }
        let step_norm = step_norm_sq.sqrt();
        last_step_norm = step_norm;
        cand_bit_ident = bit_eq;
        let f_prev = fx;
        fx = f_cand;
        last_obj_change = fx - f_prev;
        obj_bit_ident = fx.to_bits() == f_prev.to_bits();

        let g_new = central_diff_grad(&mut eval_f, &x, None);
        if is_invalid_state(fx, &g_new) {
            g = g_new;
            term_reason = TERM_INVALID_OBJECTIVE.to_string();
            ls_cause = LS_CAUSE_INVALID_OBJECTIVE.to_string();
            break;
        }

        let mut sy = 0.0;
        let mut s_norm_sq = 0.0;
        let mut y_norm_sq = 0.0;
        for i in 0..n {
            y[i] = g_new[i] - g[i];
            sy += s[i] * y[i];
            s_norm_sq += s[i] * s[i];
            y_norm_sq += y[i] * y[i];
        }
        g = g_new;

        if step_norm < 1e-13 && (fx - f_prev).abs() < 1e-12 * fx.abs().max(1.0) {
            if inf_norm(&g, &x) <= grad_tol {
                let (hb, _) = check_bounds(&x);
                term_reason = if hb {
                    TERM_BOUNDARY_STATIONARY.to_string()
                } else {
                    TERM_GRADIENT_TOL.to_string()
                };
            } else {
                term_reason = TERM_LINE_SEARCH_FAILED.to_string();
                ls_cause = LS_CAUSE_STEP_UNDERFLOW.to_string();
            }
            break;
        }

        if sy > 1e-12 * (s_norm_sq * y_norm_sq).sqrt() {
            let rho = 1.0 / sy;
            if iter == 0 {
                let scale = sy / y_norm_sq;
                for i in 0..n {
                    h_mat[i][i] = scale;
                }
            }
            let mut y_hy = 0.0;
            for i in 0..n {
                let mut sum = 0.0;
                for j in 0..n {
                    sum += h_mat[i][j] * y[j];
                }
                hy[i] = sum;
                y_hy += y[i] * sum;
            }
            let coef = (1.0 + rho * y_hy) * rho;
            for i in 0..n {
                for j in i..n {
                    let val =
                        h_mat[i][j] + coef * s[i] * s[j] - rho * (hy[i] * s[j] + s[i] * hy[j]);
                    h_mat[i][j] = val;
                    h_mat[j][i] = val;
                }
            }
        }
    }

    let (hit_bound, hit_by_param) = check_bounds(&x);
    let gn_final = inf_norm(&g, &x);
    let raw_final = raw_inf_norm(&g);
    if term_reason.is_empty() {
        if gn_final <= grad_tol {
            term_reason = if hit_bound {
                TERM_BOUNDARY_STATIONARY.to_string()
            } else {
                TERM_GRADIENT_TOL.to_string()
            };
        } else {
            term_reason = TERM_MAX_ITERATIONS.to_string();
        }
    }
    let scaled_final = gn_final / fx.abs().max(1.0);
    BfgsResult {
        x,
        fun: fx,
        raw_grad_norm: raw_final,
        grad_norm: gn_final,
        scaled_grad_norm: scaled_final,
        iterations: iters,
        func_evals: evals,
        converged: term_reason == TERM_GRADIENT_TOL,
        hit_bound,
        hit_bound_by_param: hit_by_param,
        termination_reason: term_reason,
        min_alpha_tested,
        last_step_norm,
        candidate_bit_identical: cand_bit_ident,
        last_obj_change,
        obj_bit_identical: obj_bit_ident,
        line_search_failure_cause: ls_cause,
    }
}
