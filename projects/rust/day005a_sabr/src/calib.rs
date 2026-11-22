use crate::hagan::hagan_implied_vol;

/// Nelder-Mead simplex optimizer for `n`-dimensional minimization.
pub fn nelder_mead<F>(
    x0: &[f64],
    steps: &[f64],
    max_iter: usize,
    tol: f64,
    mut f: F,
) -> (Vec<f64>, f64, usize)
where
    F: FnMut(&[f64]) -> f64,
{
    let n = x0.len();
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    simplex.push((x0.to_vec(), f(x0)));
    let mut evals = 1_usize;
    for i in 0..n {
        let mut xi = x0.to_vec();
        xi[i] += steps[i];
        let fi = f(&xi);
        evals += 1;
        simplex.push((xi, fi));
    }

    for _ in 0..max_iter {
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let f_best = simplex[0].1;
        let f_worst = simplex[n].1;
        if (f_worst - f_best).abs() < tol {
            break;
        }
        let mut max_diam = 0.0_f64;
        for i in 1..=n {
            for j in 0..n {
                max_diam = max_diam.max((simplex[i].0[j] - simplex[0].0[j]).abs());
            }
        }
        if max_diam < 1e-9 {
            break;
        }

        let mut centroid = vec![0.0_f64; n];
        for (xi, _) in simplex.iter().take(n) {
            for j in 0..n {
                centroid[j] += xi[j];
            }
        }
        for c in &mut centroid {
            *c /= n as f64;
        }

        let mut xr = vec![0.0_f64; n];
        for j in 0..n {
            xr[j] = centroid[j] + 1.0 * (centroid[j] - simplex[n].0[j]);
        }
        let fr = f(&xr);
        evals += 1;

        if fr >= simplex[0].1 && fr < simplex[n - 1].1 {
            simplex[n] = (xr, fr);
            continue;
        }
        if fr < simplex[0].1 {
            let mut xe = vec![0.0_f64; n];
            for j in 0..n {
                xe[j] = centroid[j] + 2.0 * (xr[j] - centroid[j]);
            }
            let fe = f(&xe);
            evals += 1;
            if fe < fr {
                simplex[n] = (xe, fe);
            } else {
                simplex[n] = (xr, fr);
            }
            continue;
        }
        let mut xc = vec![0.0_f64; n];
        if fr < simplex[n].1 {
            for j in 0..n {
                xc[j] = centroid[j] + 0.5 * (xr[j] - centroid[j]);
            }
        } else {
            for j in 0..n {
                xc[j] = centroid[j] - 0.5 * (centroid[j] - simplex[n].0[j]);
            }
        }
        let fc = f(&xc);
        evals += 1;
        if fc < simplex[n].1 {
            simplex[n] = (xc, fc);
            continue;
        }
        let x0_best = simplex[0].0.clone();
        for (xi, fi) in simplex.iter_mut().skip(1) {
            for j in 0..n {
                xi[j] = x0_best[j] + 0.5 * (xi[j] - x0_best[j]);
            }
            *fi = f(xi);
            evals += 1;
        }
    }
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    (simplex[0].0.clone(), simplex[0].1, evals)
}

#[derive(Clone, Debug)]
pub struct CalibFit {
    pub alpha: f64,
    pub rho: f64,
    pub nu: f64,
    pub beta: f64,
    pub rmse_iv_bps: f64,
    pub evals: usize,
    pub corr_matrix: [[f64; 3]; 3],
}

/// Calibrates `(alpha, rho, nu)` for fixed `beta` to a target `(strikes, target_ivs)` smile at `t`.
/// Computes 3x3 numerical Hessian at the optimum and parameter correlation matrix.
pub fn calibrate_sabr_smile(
    f0: f64,
    t: f64,
    beta: f64,
    strikes: &[f64],
    target_ivs: &[f64],
    init_guess: [f64; 3],
) -> CalibFit {
    let sse = |p: &[f64]| -> f64 {
        let (a, r, nu) = (p[0], p[1], p[2]);
        if a <= 1e-5 || a > 2.0 || r.abs() >= 0.999 || nu <= 1e-5 || nu > 5.0 {
            return 1e6;
        }
        let mut s = 0.0_f64;
        for (&k, &iv_tgt) in strikes.iter().zip(target_ivs.iter()) {
            let iv_mod = hagan_implied_vol(f0, k, t, a, beta, r, nu);
            let d = iv_mod - iv_tgt;
            s += d * d;
        }
        s
    };

    let steps = [0.005, 0.1, 0.1];
    let (opt, min_sse, evals) = nelder_mead(&init_guess, &steps, 800, 1e-16, sse);
    let rmse_iv_bps = (min_sse / strikes.len() as f64).sqrt() * 1e4;

    // Numerical 3x3 Hessian of SSE -> invert -> correlation matrix
    let h_eps = [1e-4, 1e-3, 1e-3];
    let mut hess = [[0.0_f64; 3]; 3];
    let f0_val = sse(&opt);
    for i in 0..3 {
        let mut xp = opt.clone();
        let mut xm = opt.clone();
        xp[i] += h_eps[i];
        xm[i] -= h_eps[i];
        hess[i][i] = (sse(&xp) - 2.0 * f0_val + sse(&xm)) / (h_eps[i] * h_eps[i]);
        for j in (i + 1)..3 {
            let mut xpp = opt.clone();
            let mut xpm = opt.clone();
            let mut xmp = opt.clone();
            let mut xmm = opt.clone();
            xpp[i] += h_eps[i];
            xpp[j] += h_eps[j];
            xpm[i] += h_eps[i];
            xpm[j] -= h_eps[j];
            xmp[i] -= h_eps[i];
            xmp[j] += h_eps[j];
            xmm[i] -= h_eps[i];
            xmm[j] -= h_eps[j];
            let hij = (sse(&xpp) - sse(&xpm) - sse(&xmp) + sse(&xmm)) / (4.0 * h_eps[i] * h_eps[j]);
            hess[i][j] = hij;
            hess[j][i] = hij;
        }
    }
    let corr_matrix = inv_corr_3x3(&hess);

    CalibFit {
        alpha: opt[0],
        rho: opt[1],
        nu: opt[2],
        beta,
        rmse_iv_bps,
        evals,
        corr_matrix,
    }
}

fn inv_corr_3x3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut a = *m;
    let mut inv = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for i in 0..3 {
        let mut piv = a[i][i];
        if piv.abs() < 1e-14 {
            piv = 1e-14;
        }
        for j in 0..3 {
            a[i][j] /= piv;
            inv[i][j] /= piv;
        }
        for r in 0..3 {
            if r != i {
                let fac = a[r][i];
                for j in 0..3 {
                    a[r][j] -= fac * a[i][j];
                    inv[r][j] -= fac * inv[i][j];
                }
            }
        }
    }
    let mut corr = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let d = (inv[i][i].abs() * inv[j][j].abs()).sqrt().max(1e-15);
            corr[i][j] = (inv[i][j] / d).clamp(-1.0, 1.0);
        }
    }
    corr
}
