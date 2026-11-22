//! Gauss-Legendre quadrature on `[-1, 1]` and `[a, b]` computed from Legendre polynomial roots.

use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct GaussLegendre {
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
}

impl GaussLegendre {
    /// Constructs an `n`-point Gauss-Legendre rule on `[-1, 1]` via Newton-Raphson on `P_n(x)`.
    pub fn new(n: usize) -> Self {
        assert!(n >= 1, "Gauss-Legendre order must be >= 1");
        let mut nodes = vec![0.0_f64; n];
        let mut weights = vec![0.0_f64; n];
        let m = n.div_ceil(2);
        for i in 0..m {
            let mut z = (PI * ((i as f64) + 0.75) / ((n as f64) + 0.5)).cos();
            let mut pp = 0.0_f64;
            for _ in 0..60 {
                let mut p1 = 1.0_f64;
                let mut p2 = 0.0_f64;
                for j in 1..=n {
                    let p3 = p2;
                    p2 = p1;
                    p1 = ((2 * j - 1) as f64 * z * p2 - (j - 1) as f64 * p3) / j as f64;
                }
                pp = n as f64 * (z * p1 - p2) / (z * z - 1.0);
                let z1 = z;
                z = z1 - p1 / pp;
                if (z - z1).abs() < 1e-15 {
                    break;
                }
            }
            nodes[i] = -z;
            nodes[n - 1 - i] = z;
            let w = 2.0 / ((1.0 - z * z) * pp * pp);
            weights[i] = w;
            weights[n - 1 - i] = w;
        }
        Self { nodes, weights }
    }

    pub fn integrate<F>(&self, a: f64, b: f64, mut f: F) -> f64
    where
        F: FnMut(f64) -> f64,
    {
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        let mut sum = 0.0_f64;
        for (&xi, &wi) in self.nodes.iter().zip(self.weights.iter()) {
            sum += wi * f(mid + half * xi);
        }
        half * sum
    }
}
