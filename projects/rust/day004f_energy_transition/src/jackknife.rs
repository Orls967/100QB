//! Temporal two-half jackknife (unclamped) for AFNS mean-reversion bias correction.

use crate::kalman::{IS_POSITIVE_PARAM, LOWER_BOUNDS, UPPER_BOUNDS};

pub const TARGET_PARAM_INDICES: [usize; 4] = [1, 2, 3, 4];

#[inline]
pub fn combine_split(full: f64, half1: f64, half2: f64) -> f64 {
    2.0 * full - 0.5 * (half1 + half2)
}

#[derive(Clone, Debug)]
pub struct JackknifeEvalResult {
    pub est_params: [f64; 11],
    pub out_of_bounds: [bool; 11],
    pub non_physical: [bool; 11],
}

pub fn combine_vec(full: &[f64; 11], half1: &[f64; 11], half2: &[f64; 11]) -> JackknifeEvalResult {
    let mut est_params = [0.0; 11];
    let mut out_of_bounds = [false; 11];
    let mut non_physical = [false; 11];
    for i in 0..11 {
        let val = combine_split(full[i], half1[i], half2[i]);
        est_params[i] = val;
        if val < LOWER_BOUNDS[i] || val > UPPER_BOUNDS[i] {
            out_of_bounds[i] = true;
        }
        if IS_POSITIVE_PARAM[i] && val <= 0.0 {
            non_physical[i] = true;
        }
    }
    JackknifeEvalResult {
        est_params,
        out_of_bounds,
        non_physical,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jackknife_c_over_h_algebraic_exactness() {
        let psi_true = 0.05;
        let c = 4.0;
        for &h in &[20.0_f64, 80.0_f64] {
            let full = psi_true + c / h;
            let half = psi_true + c / (0.5 * h);
            let jk = combine_split(full, half, half);
            assert!((jk - psi_true).abs() < 1e-15);
        }
    }
}
