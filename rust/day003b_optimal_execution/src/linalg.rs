//! Aljabar linear murni: Cholesky, Thomas tridiagonal, dan Jacobi eigenvalue simetris.

/// Produk titik dua vektor.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Perkalian matriks-vektor ``A * x``.
pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

/// Bentuk kuadrat ``x^T * A * x``.
pub fn quad_form(a: &[Vec<f64>], x: &[f64]) -> f64 {
    let ax = mat_vec(a, x);
    dot(x, &ax)
}

/// Dekomposisi Cholesky ``A = L * L^T`` untuk matriks simetris.
/// Mengembalikan ``Err(idx)`` jika pivot pada baris ``idx`` <= 0 (bukan positif-definit).
pub fn cholesky(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, usize> {
    let n = a.len();
    let mut l = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let sum_prod: f64 = l[i][..j]
                .iter()
                .zip(l[j][..j].iter())
                .map(|(a, b)| a * b)
                .sum();
            let sum = a[i][j] - sum_prod;
            if i == j {
                if sum <= 0.0 || !sum.is_finite() {
                    return Err(i);
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Selesaikan ``L * L^T * x = b`` dari faktor Cholesky bawah ``L``.
pub fn solve_cholesky(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    assert_eq!(b.len(), n);
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = y[i];
        for k in (i + 1)..n {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    x
}

/// Selesaikan sistem SPD ``A * x = b`` lewat Cholesky.
pub fn solve_spd(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, usize> {
    let l = cholesky(a)?;
    Ok(solve_cholesky(&l, b))
}

/// Algoritma Thomas O(N) untuk sistem tridiagonal.
/// ``lower[i]`` untuk i=1..n-1 (lower[0] diabaikan), ``diag[0..n]``, ``upper[i]`` untuk i=0..n-2.
pub fn solve_tridiagonal(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64]) -> Vec<f64> {
    let n = diag.len();
    assert_eq!(lower.len(), n);
    assert_eq!(upper.len(), n);
    assert_eq!(rhs.len(), n);
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![rhs[0] / diag[0]];
    }

    let mut c_prime = vec![0.0; n];
    let mut d_prime = vec![0.0; n];

    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = rhs[0] / diag[0];

    for i in 1..n {
        let denom = diag[i] - lower[i] * c_prime[i - 1];
        c_prime[i] = if i < n - 1 { upper[i] / denom } else { 0.0 };
        d_prime[i] = (rhs[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    let mut x = vec![0.0; n];
    x[n - 1] = d_prime[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    x
}

/// Eigenvalue minimum dan maksimum matriks simetris (untuk bilangan kondisi).
pub fn jacobi_extreme_eigenvalues(a: &[Vec<f64>]) -> (f64, f64) {
    let (min_val, _) = jacobi_min_eigenvalue(a);
    // eigenvalue maksimum = -(eigenvalue minimum dari -A)
    let neg: Vec<Vec<f64>> = a.iter().map(|r| r.iter().map(|v| -v).collect()).collect();
    let (neg_min, _) = jacobi_min_eigenvalue(&neg);
    (min_val, -neg_min)
}

/// Cari eigenvalue minimum dan eigenvector-nya untuk matriks simetris via rotasi Jacobi.
pub fn jacobi_min_eigenvalue(a: &[Vec<f64>]) -> (f64, Vec<f64>) {
    let n = a.len();
    assert!(n > 0);
    let mut m = a.to_vec();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    let max_sweeps = 60 * n * n;
    for _ in 0..max_sweeps {
        let mut p = 0;
        let mut q = 1;
        let mut max_off = 0.0;
        for (i, row) in m.iter().enumerate() {
            for (j, &val) in row.iter().enumerate().skip(i + 1) {
                let av = val.abs();
                if av > max_off {
                    max_off = av;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-13 {
            break;
        }

        let app = m[p][p];
        let aqq = m[q][q];
        let apq = m[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        m[p][p] = app - t * apq;
        m[q][q] = aqq + t * apq;
        m[p][q] = 0.0;
        m[q][p] = 0.0;

        for (r, row) in m.iter_mut().enumerate() {
            if r != p && r != q {
                let mrp = row[p];
                let mrq = row[q];
                row[p] = c * mrp - s * mrq;
                row[q] = s * mrp + c * mrq;
            }
        }
        let col_p: Vec<f64> = m.iter().map(|row| row[p]).collect();
        let col_q: Vec<f64> = m.iter().map(|row| row[q]).collect();
        m[p][..n].copy_from_slice(&col_p[..n]);
        m[q][..n].copy_from_slice(&col_q[..n]);
        for row in &mut v {
            let vrp = row[p];
            let vrq = row[q];
            row[p] = c * vrp - s * vrq;
            row[q] = s * vrp + c * vrq;
        }
    }

    let mut min_idx = 0;
    let mut min_val = m[0][0];
    for (i, row) in m.iter().enumerate().skip(1) {
        if row[i] < min_val {
            min_val = row[i];
            min_idx = i;
        }
    }
    let eigvec = (0..n).map(|r| v[r][min_idx]).collect();
    (min_val, eigvec)
}
