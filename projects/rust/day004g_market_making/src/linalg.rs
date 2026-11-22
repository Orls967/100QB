//! Pure Rust linear algebra: symmetric Jacobi eigensolver and scaling-and-squaring matrix exponential.

pub fn zeros(r: usize, c: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; c]; r]
}

pub fn identity(n: usize) -> Vec<Vec<f64>> {
    let mut m = zeros(n, n);
    for i in 0..n {
        m[i][i] = 1.0;
    }
    m
}

pub fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let k_len = b.len();
    let c = b[0].len();
    let mut out = zeros(r, c);
    for i in 0..r {
        for k in 0..k_len {
            let aik = a[i][k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..c {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(&a, &b)| a * b).sum())
        .collect()
}

pub fn mat_inf_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
}

/// Matrix exponential `exp(A)` via scaling-and-squaring with Taylor order 18.
pub fn mat_exp_scaling_squaring(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let norm = mat_inf_norm(a);
    let s = if norm > 0.5 {
        (norm / 0.5).log2().ceil().max(0.0) as u32
    } else {
        0
    };
    let scale = 2.0_f64.powi(-(s as i32));
    let mut a_scaled = zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            a_scaled[i][j] = a[i][j] * scale;
        }
    }

    let mut res = identity(n);
    let mut term = identity(n);
    for k in 1..=18 {
        let next = mat_mul(&term, &a_scaled);
        let inv_k = 1.0 / k as f64;
        let mut max_t = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                term[i][j] = next[i][j] * inv_k;
                res[i][j] += term[i][j];
                max_t = max_t.max(term[i][j].abs());
            }
        }
        if max_t < 1e-17 {
            break;
        }
    }

    for _ in 0..s {
        res = mat_mul(&res, &res);
    }
    res
}

pub type EigenResult = (Vec<f64>, Vec<Vec<f64>>, usize);

/// Symmetric Jacobi eigensolver returning `(eigenvalues_asc, eigenvectors_by_row, sweeps)`.
pub fn symmetric_eigen_jacobi_asc(
    mat: &[Vec<f64>],
    max_sweeps: usize,
    tol: f64,
) -> Result<EigenResult, String> {
    let n = mat.len();
    let mut a = zeros(n, n);
    let mut v = identity(n);
    for i in 0..n {
        for j in 0..n {
            a[i][j] = 0.5 * (mat[i][j] + mat[j][i]);
        }
    }
    let mut sweeps = 0usize;
    for sw in 0..max_sweeps {
        sweeps = sw + 1;
        let mut any_rotated = false;
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.abs() <= tol {
                    continue;
                }
                any_rotated = true;
                let app = a[p][p];
                let aqq = a[q][q];
                let tau = (aqq - app) / (2.0 * apq);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;

                a[p][p] = app - t * apq;
                a[q][q] = aqq + t * apq;
                a[p][q] = 0.0;
                a[q][p] = 0.0;

                for r in 0..n {
                    if r != p && r != q {
                        let arp = a[r][p];
                        let arq = a[r][q];
                        a[r][p] = c * arp - s * arq;
                        a[p][r] = a[r][p];
                        a[r][q] = s * arp + c * arq;
                        a[q][r] = a[r][q];
                    }
                }
                for r in 0..n {
                    let vrp = v[r][p];
                    let vrq = v[r][q];
                    v[r][p] = c * vrp - s * vrq;
                    v[r][q] = s * vrp + c * vrq;
                }
            }
        }
        if !any_rotated {
            break;
        }
    }

    let mut pairs: Vec<(f64, Vec<f64>)> = (0..n)
        .map(|j| (a[j][j], (0..n).map(|i| v[i][j]).collect()))
        .collect();
    pairs.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let vals = pairs.iter().map(|p| p.0).collect();
    let vecs = pairs.into_iter().map(|p| p.1).collect();
    Ok((vals, vecs, sweeps))
}
