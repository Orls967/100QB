/// Tridiagonal linear system solver (Thomas algorithm) in O(N).
/// Solves `a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]` for `i = 0..n`.
pub fn solve_tridiagonal(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    assert_eq!(a.len(), n);
    assert_eq!(b.len(), n);
    assert_eq!(c.len(), n);
    if n == 0 {
        return Vec::new();
    }
    let mut c_prime = vec![0.0_f64; n];
    let mut d_prime = vec![0.0_f64; n];

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for i in 1..n {
        let denom = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = if i + 1 < n { c[i] / denom } else { 0.0 };
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    let mut x = vec![0.0_f64; n];
    x[n - 1] = d_prime[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    x
}
