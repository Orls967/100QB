//! Dense linear algebra: Cholesky decomposition, linear solve, and inverse for small symmetric matrices.

#[allow(clippy::needless_range_loop)]
pub fn cholesky(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let mut l = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum <= 0.0 || !sum.is_finite() {
                    return Err(format!("non-PSD pivot at index {i}: {sum}"));
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Solves symmetric positive-definite `A x = b` via Cholesky.
pub fn solve_spd(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, String> {
    let n = a.len();
    let l = cholesky(a)?;
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i][j] * y[j];
        }
        y[i] = s / l[i][i];
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for j in (i + 1)..n {
            s -= l[j][i] * x[j];
        }
        x[i] = s / l[i][i];
    }
    Ok(x)
}

/// Inverts a small symmetric matrix (e.g. 3x3 Hessian) via Gaussian elimination with partial pivoting.
#[allow(clippy::needless_range_loop)]
pub fn invert_dense(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let mut aug = vec![vec![0.0_f64; 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    for col in 0..n {
        let mut piv = col;
        let mut best = aug[col][col].abs();
        for (r, row) in aug.iter().enumerate().skip(col + 1) {
            if row[col].abs() > best {
                best = row[col].abs();
                piv = r;
            }
        }
        if best < 1e-14 || !best.is_finite() {
            return Err(format!("singular matrix at col {col}"));
        }
        aug.swap(col, piv);
        let div = aug[col][col];
        for c in 0..(2 * n) {
            aug[col][c] /= div;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = aug[r][col];
            if f != 0.0 {
                for c in 0..(2 * n) {
                    aug[r][c] -= f * aug[col][c];
                }
            }
        }
    }
    let mut inv = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            inv[i][j] = aug[i][n + j];
        }
    }
    Ok(inv)
}
