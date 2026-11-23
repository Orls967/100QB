use crate::gbm::GbmParams;
use crate::hitting::prob_hit_infinite;
use crate::special::{inv_norm_cdf, norm_cdf};

#[derive(Clone, Debug)]
pub struct FractionalCurveRow {
    pub c: f64,
    pub f: f64,
    pub g: f64,
    pub s: f64,
    pub excess_ratio_measured: f64,
    pub excess_ratio_formula: f64,
    pub prob_hit_half_inf: f64,
}

#[derive(Clone, Debug)]
pub struct FractionalMapRow {
    pub c: f64,
    pub t: f64,
    pub median_wealth: f64,
    pub p05_wealth: f64,
    pub sd_log_wealth: f64,
    pub prob_loss: f64,
}

#[derive(Clone, Debug)]
pub struct FractionalRoots {
    pub c_zero_excess_closed: f64,
    pub c_zero_excess_numeric: f64,
    pub c_hit_half_eq_half_closed: f64,
    pub c_hit_half_eq_half_numeric: f64,
    pub c_hit_half_eq_quarter_closed: f64,
    pub c_hit_half_eq_quarter_numeric: f64,
}

pub fn build_fractional_curve(gbm: &GbmParams, cs: &[f64]) -> Vec<FractionalCurveRow> {
    let f_star = gbm.kelly_fraction();
    let excess_star = gbm.kelly_growth() - gbm.r;
    cs.iter()
        .map(|&c| {
            let f = c * f_star;
            let g = gbm.drift_log(f);
            let s = gbm.vol_log(f);
            let excess_ratio_measured = (g - gbm.r) / excess_star;
            let excess_ratio_formula = c * (2.0 - c);
            let prob_hit_half_inf = prob_hit_infinite(g, s, 0.5);
            FractionalCurveRow {
                c,
                f,
                g,
                s,
                excess_ratio_measured,
                excess_ratio_formula,
                prob_hit_half_inf,
            }
        })
        .collect()
}

pub fn build_fractional_map(gbm: &GbmParams, cs: &[f64], ts: &[f64]) -> Vec<FractionalMapRow> {
    let z05 = inv_norm_cdf(0.05);
    let f_star = gbm.kelly_fraction();
    let mut out = Vec::with_capacity(cs.len() * ts.len());
    for &c in cs {
        let f = c * f_star;
        let g = gbm.drift_log(f);
        let s = gbm.vol_log(f);
        for &t in ts {
            let sd_log_wealth = s * t.sqrt();
            let median_wealth = (g * t).exp();
            let p05_wealth = (g * t + z05 * sd_log_wealth).exp();
            let prob_loss = norm_cdf(-(g * t) / sd_log_wealth);
            out.push(FractionalMapRow {
                c,
                t,
                median_wealth,
                p05_wealth,
                sd_log_wealth,
                prob_loss,
            });
        }
    }
    out
}

/// Closed-form and numerical roots for `g(c*f*) = r` and `P_inf(hit 0.5) = p_target`.
pub fn compute_fractional_roots(gbm: &GbmParams) -> FractionalRoots {
    let f_star = gbm.kelly_fraction();
    let c_zero_excess_closed = 2.0;

    // Bisection for g(c*f*) - r = 0 on [0.5, 3.0]
    let mut lo = 0.5_f64;
    let mut hi = 3.0_f64;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if gbm.drift_log(mid * f_star) - gbm.r > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let c_zero_excess_numeric = 0.5 * (lo + hi);

    let solve_hit_target = |p_target: f64| -> (f64, f64) {
        // 0.5^(2g/s^2) = p_target => 2g/s^2 = -ln(p_target)/ln(2) = k
        let k = -p_target.ln() / 2.0_f64.ln();
        let a_coef = (k + 1.0) * f_star * f_star * gbm.sigma * gbm.sigma;
        let b_coef = -2.0 * f_star * (gbm.mu - gbm.r);
        let c_coef = -2.0 * gbm.r;
        let disc = b_coef * b_coef - 4.0 * a_coef * c_coef;
        let closed = (-b_coef + disc.sqrt()) / (2.0 * a_coef);

        let mut blo = 0.05_f64;
        let mut bhi = 4.0_f64;
        for _ in 0..100 {
            let mid = 0.5 * (blo + bhi);
            let f = mid * f_star;
            let p_mid = prob_hit_infinite(gbm.drift_log(f), gbm.vol_log(f), 0.5);
            if p_mid < p_target {
                blo = mid;
            } else {
                bhi = mid;
            }
        }
        (closed, 0.5 * (blo + bhi))
    };

    let (c_hit_half_eq_half_closed, c_hit_half_eq_half_numeric) = solve_hit_target(0.5);
    let (c_hit_half_eq_quarter_closed, c_hit_half_eq_quarter_numeric) = solve_hit_target(0.25);

    FractionalRoots {
        c_zero_excess_closed,
        c_zero_excess_numeric,
        c_hit_half_eq_half_closed,
        c_hit_half_eq_half_numeric,
        c_hit_half_eq_quarter_closed,
        c_hit_half_eq_quarter_numeric,
    }
}
