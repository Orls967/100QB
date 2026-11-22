//! Market making quote strategies: GLT optimal, Avellaneda-Stoikov classic,
//! symmetric fixed, tight fixed, and adverse-selection adjusted.

use crate::glt::{GltParams, QuoteTable};

#[derive(Clone)]
pub enum Strategy<'a> {
    GltOptimal { table: &'a QuoteTable },
    AvellanedaStoikov { p: GltParams },
    SymmetricFixed { delta: f64, q_max: isize },
    GltConstantShift { table: &'a QuoteTable, shift: f64 },
    GltAsExact { table: &'a QuoteTable },
}

impl<'a> Strategy<'a> {
    #[inline]
    pub fn quotes(&self, t: f64, q: isize) -> (f64, f64) {
        match self {
            Strategy::GltOptimal { table } => table.interp_quotes(t, q),
            Strategy::AvellanedaStoikov { p } => {
                let tau = (p.t_horizon - t).max(0.0);
                let half_spread =
                    0.5 * p.gamma * p.sigma * p.sigma * tau + p.half_spread_zero_inv();
                let skew = (q as f64) * p.gamma * p.sigma * p.sigma * tau;
                let db = if q >= p.q_max {
                    f64::INFINITY
                } else {
                    half_spread + skew
                };
                let da = if q <= -p.q_max {
                    f64::INFINITY
                } else {
                    half_spread - skew
                };
                (db, da)
            }
            Strategy::SymmetricFixed { delta, q_max } => {
                let db = if q >= *q_max { f64::INFINITY } else { *delta };
                let da = if q <= -*q_max { f64::INFINITY } else { *delta };
                (db, da)
            }
            Strategy::GltConstantShift { table, shift } => {
                let (db, da) = table.interp_quotes(t, q);
                (db + shift, da + shift)
            }
            Strategy::GltAsExact { table } => table.interp_quotes(t, q),
        }
    }
}
