//! Pure Rust linear algebra, Gauss-Legendre quadrature, and RK4 ODE integrator.

#[derive(Clone, Copy, Debug)]
pub struct QuadNode {
    pub x: f64,
    pub w: f64,
}

pub fn gauss_legendre(n: usize, a: f64, b: f64) -> Vec<QuadNode> {
    let mut nodes = vec![QuadNode { x: 0.0, w: 0.0 }; n];
    let m = n.div_ceil(2);
    let xm = 0.5 * (b + a);
    let xl = 0.5 * (b - a);
    for i in 0..m {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut pp = 0.0;
        for _ in 0..80 {
            let mut p1 = 1.0;
            let mut p2 = 0.0;
            for j in 1..=n {
                let p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
            }
            pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
            let z1 = z;
            z = z1 - p1 / pp;
            if (z - z1).abs() <= 1e-15 {
                break;
            }
        }
        nodes[i].x = xm - xl * z;
        nodes[n - 1 - i].x = xm + xl * z;
        let w = 2.0 * xl / ((1.0 - z * z) * pp * pp);
        nodes[i].w = w;
        nodes[n - 1 - i].w = w;
    }
    nodes
}

pub fn rk4_system<F>(t0: f64, t1: f64, y0: &[f64], n_steps: usize, deriv: F) -> Vec<f64>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let dim = y0.len();
    let mut y = y0.to_vec();
    if n_steps == 0 || t1 == t0 {
        return y;
    }
    let h = (t1 - t0) / n_steps as f64;
    let mut k1 = vec![0.0; dim];
    let mut k2 = vec![0.0; dim];
    let mut k3 = vec![0.0; dim];
    let mut k4 = vec![0.0; dim];
    let mut tmp = vec![0.0; dim];

    for step in 0..n_steps {
        let t = t0 + step as f64 * h;
        deriv(t, &y, &mut k1);
        for i in 0..dim {
            tmp[i] = y[i] + 0.5 * h * k1[i];
        }
        deriv(t + 0.5 * h, &tmp, &mut k2);
        for i in 0..dim {
            tmp[i] = y[i] + 0.5 * h * k2[i];
        }
        deriv(t + 0.5 * h, &tmp, &mut k3);
        for i in 0..dim {
            tmp[i] = y[i] + h * k3[i];
        }
        deriv(t + h, &tmp, &mut k4);
        for i in 0..dim {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }
    y
}

pub fn rk4_richardson<F>(t0: f64, t1: f64, y0: &[f64], n_coarse: usize, deriv: F) -> Vec<f64>
where
    F: Fn(f64, &[f64], &mut [f64]),
{
    let y1 = rk4_system(t0, t1, y0, n_coarse, &deriv);
    let y2 = rk4_system(t0, t1, y0, n_coarse * 2, &deriv);
    y1.iter()
        .zip(y2.iter())
        .map(|(&a, &b)| (16.0 * b - a) / 15.0)
        .collect()
}

pub fn zeros(r: usize, c: usize) -> Vec<Vec<f64>> {
    vec![vec![0.0; c]; r]
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

pub fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let c = a[0].len();
    let mut t = zeros(c, r);
    for i in 0..r {
        for j in 0..c {
            t[j][i] = a[i][j];
        }
    }
    t
}

pub fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let k_len = b.len();
    let c = b[0].len();
    let mut out = zeros(r, c);
    for i in 0..r {
        for k in 0..k_len {
            let aik = a[i][k];
            for j in 0..c {
                out[i][j] += aik * b[k][j];
            }
        }
    }
    out
}

pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter().map(|row| dot(row, x)).collect()
}

pub fn mat_exp_series(a: &[Vec<f64>], order: usize) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut res = zeros(n, n);
    let mut term = zeros(n, n);
    for i in 0..n {
        res[i][i] = 1.0;
        term[i][i] = 1.0;
    }
    for k in 1..=order {
        let next = mat_mul(&term, a);
        let inv_k = 1.0 / k as f64;
        for i in 0..n {
            for j in 0..n {
                term[i][j] = next[i][j] * inv_k;
                res[i][j] += term[i][j];
            }
        }
    }
    res
}

pub fn solve_lu(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, String> {
    let n = a.len();
    if n == 0 || b.len() != n {
        return Err("solve_lu: dimension mismatch".into());
    }
    let mut lu = a.to_vec();
    let mut rhs = b.to_vec();

    for k in 0..n {
        let mut piv = k;
        let mut max_v = lu[k][k].abs();
        for i in (k + 1)..n {
            let v = lu[i][k].abs();
            if v > max_v {
                max_v = v;
                piv = i;
            }
        }
        if max_v < 1e-15 {
            return Err("solve_lu: singular matrix".into());
        }
        if piv != k {
            lu.swap(k, piv);
            rhs.swap(k, piv);
        }
        let inv_diag = 1.0 / lu[k][k];
        for i in (k + 1)..n {
            let fac = lu[i][k] * inv_diag;
            lu[i][k] = fac;
            for j in (k + 1)..n {
                lu[i][j] -= fac * lu[k][j];
            }
            rhs[i] -= fac * rhs[k];
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= lu[i][j] * x[j];
        }
        x[i] = s / lu[i][i];
    }
    Ok(x)
}

pub fn solve_tridiagonal(
    sub: &[f64],
    diag: &[f64],
    sup: &[f64],
    rhs: &[f64],
) -> Result<Vec<f64>, String> {
    let n = diag.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let mut c_prime = vec![0.0; n];
    let mut d_prime = vec![0.0; n];
    if diag[0].abs() < 1e-15 {
        return Err("solve_tridiagonal: zero pivot".into());
    }
    if n > 1 {
        c_prime[0] = sup[0] / diag[0];
    }
    d_prime[0] = rhs[0] / diag[0];

    for i in 1..n {
        let denom = diag[i] - sub[i - 1] * c_prime[i - 1];
        if denom.abs() < 1e-15 {
            return Err("solve_tridiagonal: zero pivot".into());
        }
        if i + 1 < n {
            c_prime[i] = sup[i] / denom;
        }
        d_prime[i] = (rhs[i] - sub[i - 1] * d_prime[i - 1]) / denom;
    }

    let mut x = vec![0.0; n];
    x[n - 1] = d_prime[n - 1];
    for i in (0..(n - 1)).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    Ok(x)
}

pub fn cholesky(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let mut l = zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for k in 0..j {
                sum -= l[i][k] * l[j][k];
            }
            if i == j {
                if sum <= 0.0 || sum.is_nan() || sum.is_infinite() {
                    return Err(format!("cholesky: non-PD at {i}"));
                }
                l[i][j] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}

pub fn solve_cholesky(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut y = vec![0.0; n];
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i][k] * y[k];
        }
        y[i] = s / l[i][i];
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for k in (i + 1)..n {
            s -= l[k][i] * x[k];
        }
        x[i] = s / l[i][i];
    }
    x
}

pub type EigenResult = (Vec<f64>, Vec<Vec<f64>>, usize);

/// Symmetric Jacobi eigensolver returning `(eigenvalues_desc, eigenvectors_by_row, sweeps)`.
pub fn symmetric_eigen_jacobi(
    mat: &[Vec<f64>],
    max_sweeps: usize,
    tol: f64,
) -> Result<EigenResult, String> {
    let n = mat.len();
    let mut a = zeros(n, n);
    let mut v = zeros(n, n);
    for i in 0..n {
        v[i][i] = 1.0;
        for j in 0..n {
            a[i][j] = 0.5 * (mat[i][j] + mat[j][i]);
        }
    }
    let mut sweeps = 0usize;
    for sw in 0..max_sweeps {
        sweeps = sw + 1;
        let mut max_off = 0.0;
        let mut p = 0usize;
        let mut q = 1usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let val = a[i][j].abs();
                if val > max_off {
                    max_off = val;
                    p = i;
                    q = j;
                }
            }
        }
        if max_off <= tol {
            break;
        }
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
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

    let mut pairs: Vec<(f64, Vec<f64>)> = (0..n)
        .map(|j| {
            let vec_j: Vec<f64> = (0..n).map(|i| v[i][j]).collect();
            (a[j][j], vec_j)
        })
        .collect();
    pairs.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));
    let vals = pairs.iter().map(|p| p.0).collect();
    let vecs = pairs.into_iter().map(|p| p.1).collect();
    Ok((vals, vecs, sweeps))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lu_cholesky_jacobi() {
        let a = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        let b = vec![1.0, 2.0];
        let x = solve_lu(&a, &b).unwrap();
        assert!((x[0] - 1.0 / 11.0).abs() < 1e-14);
        assert!((x[1] - 7.0 / 11.0).abs() < 1e-14);

        let l = cholesky(&a).unwrap();
        let x2 = solve_cholesky(&l, &b);
        assert!((x2[0] - x[0]).abs() < 1e-14);
        assert!((x2[1] - x[1]).abs() < 1e-14);

        let (vals, vecs, _) = symmetric_eigen_jacobi(&a, 50, 1e-14).unwrap();
        assert!((vals[0] + vals[1] - 7.0).abs() < 1e-12);
        assert!((vals[0] * vals[1] - 11.0).abs() < 1e-12);
        assert!(dot(&vecs[0], &vecs[1]).abs() < 1e-12);
    }
}
