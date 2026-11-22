//! Linear algebra: LU solve, Cholesky decomposition, symmetric Jacobi eigensolver, and tridiagonal Thomas solver.

#![allow(clippy::needless_range_loop)]

pub fn solve_tridiagonal(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    if n == 0 {
        return Vec::new();
    }
    let mut c_p = vec![0.0_f64; n];
    let mut d_p = vec![0.0_f64; n];
    c_p[0] = c[0] / b[0];
    d_p[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * c_p[i - 1];
        if i + 1 < n {
            c_p[i] = c[i] / denom;
        }
        d_p[i] = (d[i] - a[i] * d_p[i - 1]) / denom;
    }
    let mut x = vec![0.0_f64; n];
    x[n - 1] = d_p[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_p[i] - c_p[i] * x[i + 1];
    }
    x
}

pub fn cholesky_2x2(a: &[[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let l00 = a[0][0].max(1e-18).sqrt();
    let l10 = a[1][0] / l00;
    let l11 = (a[1][1] - l10 * l10).max(1e-18).sqrt();
    [[l00, 0.0], [l10, l11]]
}

pub fn inv_2x2(a: &[[f64; 2]; 2]) -> [[f64; 2]; 2] {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let inv_det = 1.0 / det;
    [
        [a[1][1] * inv_det, -a[0][1] * inv_det],
        [-a[1][0] * inv_det, a[0][0] * inv_det],
    ]
}

/// Eigenvalues (descending) of a 2x2 symmetric matrix.
pub fn eig_sym_2x2_desc(a: &[[f64; 2]; 2]) -> [f64; 2] {
    let tr = a[0][0] + a[1][1];
    let diff = a[0][0] - a[1][1];
    let disc = (diff * diff + 4.0 * a[0][1] * a[1][0]).max(0.0).sqrt();
    [0.5 * (tr + disc), 0.5 * (tr - disc)]
}

/// General LU linear solve `A x = b` with partial pivoting.
pub fn solve_lu(a_in: &[Vec<f64>], b_in: &[f64]) -> Result<Vec<f64>, String> {
    let n = a_in.len();
    let mut a = a_in.to_vec();
    let mut b = b_in.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut best = a[col][col].abs();
        for r in (col + 1)..n {
            if a[r][col].abs() > best {
                best = a[r][col].abs();
                piv = r;
            }
        }
        if best < 1e-15 {
            return Err("singular matrix in solve_lu".to_string());
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let diag = a[col][col];
        for r in (col + 1)..n {
            let f = a[r][col] / diag;
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Ok(x)
}

/// Cyclic Jacobi eigensolver for general symmetric matrix, returning eigenvalues descending.
pub fn jacobi_sym_desc(a_in: &[Vec<f64>]) -> Vec<f64> {
    let n = a_in.len();
    let mut a = a_in.to_vec();
    for _ in 0..80 {
        let mut off = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                off = off.max(a[i][j].abs());
            }
        }
        if off < 1e-14 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p][q];
                if apq.abs() < 1e-15 {
                    continue;
                }
                let tau = (a[q][q] - a[p][p]) / (2.0 * apq);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    -1.0 / (-tau + (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                let app = a[p][p];
                let aqq = a[q][q];
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
            }
        }
    }
    let mut vals: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    vals.sort_by(|x, y| y.partial_cmp(x).unwrap_or(std::cmp::Ordering::Equal));
    vals
}
