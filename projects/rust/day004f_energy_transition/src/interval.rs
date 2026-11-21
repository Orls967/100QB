//! Five confidence interval methods (Wald, Percentile, Basic, Studentized, Log-Basic)
//! evaluated without clamping and with explicit validity / non-physical flags.

use crate::kalman::{IS_POSITIVE_PARAM, LOWER_BOUNDS, UPPER_BOUNDS};
use crate::stats;

pub const METHOD_WALD: &str = "wald";
pub const METHOD_PERCENTILE: &str = "percentile";
pub const METHOD_BASIC: &str = "basic";
pub const METHOD_STUDENTIZED: &str = "studentized";
pub const METHOD_LOG_BASIC: &str = "log_basic";
pub const MIN_VALID_FRACTION: f64 = 0.95;

pub fn min_required_valid(total_boot: usize) -> usize {
    if total_boot == 0 {
        return 1;
    }
    ((MIN_VALID_FRACTION * total_boot as f64).ceil() as usize).max(2)
}

#[derive(Clone, Debug)]
pub struct IntervalResult {
    pub method: &'static str,
    pub lower: f64,
    pub upper: f64,
    pub width: f64,
    pub valid: bool,
    pub invalid_reason: &'static str,
    pub covered: bool,
    pub non_physical: bool,
    pub out_of_bounds: bool,
    pub valid_draws_count: usize,
}

fn invalid_interval(
    method: &'static str,
    reason: &'static str,
    valid_draws: usize,
) -> IntervalResult {
    IntervalResult {
        method,
        lower: f64::NAN,
        upper: f64::NAN,
        width: f64::NAN,
        valid: false,
        invalid_reason: reason,
        covered: false,
        non_physical: false,
        out_of_bounds: false,
        valid_draws_count: valid_draws,
    }
}

fn finalize_interval(
    method: &'static str,
    param_idx: usize,
    lo: f64,
    hi: f64,
    true_val: f64,
    valid_draws: usize,
) -> IntervalResult {
    if !lo.is_finite() || !hi.is_finite() {
        return invalid_interval(method, "outer_non_finite", valid_draws);
    }
    let non_phys = IS_POSITIVE_PARAM[param_idx] && lo <= 0.0;
    let oob = lo < LOWER_BOUNDS[param_idx] || hi > UPPER_BOUNDS[param_idx];
    IntervalResult {
        method,
        lower: lo,
        upper: hi,
        width: hi - lo,
        valid: true,
        invalid_reason: "none",
        covered: true_val >= lo && true_val <= hi,
        non_physical: non_phys,
        out_of_bounds: oob,
        valid_draws_count: valid_draws,
    }
}

pub fn eval_wald(
    param_idx: usize,
    psi_hat: f64,
    se_hat: f64,
    true_val: f64,
    outer_finite: bool,
) -> IntervalResult {
    if !outer_finite || !psi_hat.is_finite() {
        return invalid_interval(METHOD_WALD, "outer_non_finite", 0);
    }
    if !se_hat.is_finite() || se_hat <= 0.0 {
        return invalid_interval(METHOD_WALD, "outer_se_non_positive", 0);
    }
    finalize_interval(
        METHOD_WALD,
        param_idx,
        psi_hat - stats::Z95 * se_hat,
        psi_hat + stats::Z95 * se_hat,
        true_val,
        1,
    )
}

pub fn eval_percentile(
    param_idx: usize,
    true_val: f64,
    outer_finite: bool,
    boot_ests: &[f64],
    total_boot: usize,
) -> IntervalResult {
    if !outer_finite {
        return invalid_interval(METHOD_PERCENTILE, "outer_non_finite", 0);
    }
    let finite: Vec<f64> = boot_ests
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .collect();
    if finite.len() < min_required_valid(total_boot) {
        return invalid_interval(METHOD_PERCENTILE, "insufficient_finite_draws", finite.len());
    }
    let lo = stats::quantile(&finite, 0.025);
    let hi = stats::quantile(&finite, 0.975);
    finalize_interval(METHOD_PERCENTILE, param_idx, lo, hi, true_val, finite.len())
}

pub fn eval_basic(
    param_idx: usize,
    psi_hat: f64,
    true_val: f64,
    outer_finite: bool,
    boot_ests: &[f64],
    total_boot: usize,
) -> IntervalResult {
    if !outer_finite || !psi_hat.is_finite() {
        return invalid_interval(METHOD_BASIC, "outer_non_finite", 0);
    }
    let finite: Vec<f64> = boot_ests
        .iter()
        .copied()
        .filter(|x| x.is_finite())
        .collect();
    if finite.len() < min_required_valid(total_boot) {
        return invalid_interval(METHOD_BASIC, "insufficient_finite_draws", finite.len());
    }
    let q_lo = stats::quantile(&finite, 0.025);
    let q_hi = stats::quantile(&finite, 0.975);
    finalize_interval(
        METHOD_BASIC,
        param_idx,
        2.0 * psi_hat - q_hi,
        2.0 * psi_hat - q_lo,
        true_val,
        finite.len(),
    )
}

pub fn eval_studentized(
    param_idx: usize,
    psi_hat: f64,
    se_hat: f64,
    true_val: f64,
    outer_finite: bool,
    boot_ests: &[f64],
    boot_ses: &[f64],
    total_boot: usize,
) -> IntervalResult {
    if !outer_finite || !psi_hat.is_finite() {
        return invalid_interval(METHOD_STUDENTIZED, "outer_non_finite", 0);
    }
    if !se_hat.is_finite() || se_hat <= 0.0 {
        return invalid_interval(METHOD_STUDENTIZED, "outer_se_non_positive", 0);
    }
    let n = boot_ests.len().min(boot_ses.len());
    let mut t_stars = Vec::with_capacity(n);
    for b in 0..n {
        if boot_ests[b].is_finite() && boot_ses[b].is_finite() && boot_ses[b] > 0.0 {
            t_stars.push((boot_ests[b] - psi_hat) / boot_ses[b]);
        }
    }
    if t_stars.len() < min_required_valid(total_boot) {
        return invalid_interval(
            METHOD_STUDENTIZED,
            "insufficient_valid_se_draws",
            t_stars.len(),
        );
    }
    let qt_lo = stats::quantile(&t_stars, 0.025);
    let qt_hi = stats::quantile(&t_stars, 0.975);
    finalize_interval(
        METHOD_STUDENTIZED,
        param_idx,
        psi_hat - qt_hi * se_hat,
        psi_hat - qt_lo * se_hat,
        true_val,
        t_stars.len(),
    )
}

pub fn eval_log_basic(
    param_idx: usize,
    psi_hat: f64,
    true_val: f64,
    outer_finite: bool,
    boot_ests: &[f64],
    total_boot: usize,
) -> IntervalResult {
    if !outer_finite || !psi_hat.is_finite() {
        return invalid_interval(METHOD_LOG_BASIC, "outer_non_finite", 0);
    }
    if psi_hat <= 0.0 {
        return invalid_interval(METHOD_LOG_BASIC, "outer_non_positive_for_log", 0);
    }
    let log_draws: Vec<f64> = boot_ests
        .iter()
        .copied()
        .filter(|&x| x.is_finite() && x > 0.0)
        .map(f64::ln)
        .collect();
    if log_draws.len() < min_required_valid(total_boot) {
        return invalid_interval(
            METHOD_LOG_BASIC,
            "insufficient_positive_draws",
            log_draws.len(),
        );
    }
    let eta_hat = psi_hat.ln();
    let q_lo = stats::quantile(&log_draws, 0.025);
    let q_hi = stats::quantile(&log_draws, 0.975);
    finalize_interval(
        METHOD_LOG_BASIC,
        param_idx,
        (2.0 * eta_hat - q_hi).exp(),
        (2.0 * eta_hat - q_lo).exp(),
        true_val,
        log_draws.len(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_identities_and_non_physical_flags() {
        // Basic & Percentile must have identical width on finite draws
        let draws: Vec<f64> = (1..=100).map(|i| 0.04 + 0.002 * i as f64).collect();
        let pct = eval_percentile(1, 0.05, true, &draws, 100);
        let bas = eval_basic(1, 0.25, 0.05, true, &draws, 100);
        assert!(pct.valid && bas.valid);
        assert!((pct.width - bas.width).abs() < 1e-12);
        // Basic with large psi_hat - q_hi gives negative lower endpoint -> non_physical = true for kp11 (idx 1)
        let bas_neg = eval_basic(1, 0.06, 0.05, true, &draws, 100);
        assert!(bas_neg.non_physical);
        // Log-basic is always strictly positive
        let lb = eval_log_basic(1, 0.06, 0.05, true, &draws, 100);
        assert!(lb.valid && !lb.non_physical && lb.lower > 0.0);
    }
}
