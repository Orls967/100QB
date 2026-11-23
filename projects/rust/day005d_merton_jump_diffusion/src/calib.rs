use crate::bs::implied_vol_call;
use crate::merton::{merton_call_series, MertonParams};

#[derive(Clone, Debug)]
pub struct SmileRow {
    pub t: f64,
    pub k: f64,
    pub call_price: f64,
    pub implied_vol: f64,
}

#[derive(Clone, Debug)]
pub struct CalibResult {
    pub true_params: MertonParams,
    pub init_params: MertonParams,
    pub fitted_params: MertonParams,
    pub price_rmse: f64,
    pub iv_rmse_bps: f64,
    pub iters: usize,
}

pub fn build_smile_table(p: &MertonParams, ts: &[f64], ks: &[f64]) -> Vec<SmileRow> {
    let mut out = Vec::with_capacity(ts.len() * ks.len());
    for &t in ts {
        for &k in ks {
            let call_price = merton_call_series(p, k, t);
            let implied_vol = implied_vol_call(call_price, p.s0, k, p.r, p.q, t);
            out.push(SmileRow {
                t,
                k,
                call_price,
                implied_vol,
            });
        }
    }
    out
}

/// Nelder-Mead simplex optimization in 4D `(sigma, lambda, mu_j, delta_j)`
/// minimizing SSE over synthetic market call prices.
pub fn calibrate_merton_nelder_mead(
    true_params: &MertonParams,
    init_params: &MertonParams,
    ts: &[f64],
    ks: &[f64],
    max_iters: usize,
) -> CalibResult {
    let mut market = Vec::with_capacity(ts.len() * ks.len());
    for &t in ts {
        for &k in ks {
            let px = merton_call_series(true_params, k, t);
            let iv = implied_vol_call(px, true_params.s0, k, true_params.r, true_params.q, t);
            market.push((t, k, px, iv));
        }
    }

    let obj = |x: &[f64; 4]| -> f64 {
        let sigma = x[0];
        let lambda = x[1];
        let mu_j = x[2];
        let delta_j = x[3];
        if sigma <= 0.01
            || sigma > 1.5
            || !(0.0..=10.0).contains(&lambda)
            || delta_j <= 0.01
            || delta_j > 1.5
            || !(-1.0..=1.0).contains(&mu_j)
        {
            return 1e12;
        }
        let cand = MertonParams {
            sigma,
            lambda,
            mu_j,
            delta_j,
            ..*true_params
        };
        let mut sse = 0.0_f64;
        for &(t, k, px, _) in &market {
            let diff = merton_call_series(&cand, k, t) - px;
            sse += diff * diff;
        }
        sse
    };

    let x0 = [
        init_params.sigma,
        init_params.lambda,
        init_params.mu_j,
        init_params.delta_j,
    ];
    let step = [0.03, 0.25, 0.05, 0.05];
    let mut simplex: Vec<([f64; 4], f64)> = Vec::with_capacity(5);
    simplex.push((x0, obj(&x0)));
    for i in 0..4 {
        let mut xi = x0;
        xi[i] += step[i];
        simplex.push((xi, obj(&xi)));
    }

    let mut iters = 0_usize;
    for it in 0..max_iters {
        iters = it + 1;
        simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
        let spread: f64 = (0..4)
            .map(|j| (simplex[4].0[j] - simplex[0].0[j]).abs())
            .fold(0.0_f64, f64::max);
        if spread < 1e-8 && (simplex[4].1 - simplex[0].1).abs() < 1e-12 {
            break;
        }

        let mut centroid = [0.0_f64; 4];
        for v in simplex.iter().take(4) {
            for (c, &vj) in centroid.iter_mut().zip(v.0.iter()) {
                *c += 0.25 * vj;
            }
        }

        let mut xr = [0.0_f64; 4];
        for j in 0..4 {
            xr[j] = centroid[j] + 1.0 * (centroid[j] - simplex[4].0[j]);
        }
        let fr = obj(&xr);

        if fr < simplex[0].1 {
            let mut xe = [0.0_f64; 4];
            for j in 0..4 {
                xe[j] = centroid[j] + 2.0 * (xr[j] - centroid[j]);
            }
            let fe = obj(&xe);
            simplex[4] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < simplex[3].1 {
            simplex[4] = (xr, fr);
        } else {
            let use_out = fr < simplex[4].1;
            let mut xc = [0.0_f64; 4];
            for j in 0..4 {
                xc[j] = if use_out {
                    centroid[j] + 0.5 * (xr[j] - centroid[j])
                } else {
                    centroid[j] - 0.5 * (centroid[j] - simplex[4].0[j])
                };
            }
            let fc = obj(&xc);
            if fc < if use_out { fr } else { simplex[4].1 } {
                simplex[4] = (xc, fc);
            } else {
                let best = simplex[0].0;
                for v in simplex.iter_mut().skip(1) {
                    for (vj, &bj) in v.0.iter_mut().zip(best.iter()) {
                        *vj = bj + 0.5 * (*vj - bj);
                    }
                    v.1 = obj(&v.0);
                }
            }
        }
    }

    simplex.sort_by(|a, b| a.1.total_cmp(&b.1));
    let best = simplex[0].0;
    let fitted_params = MertonParams {
        sigma: best[0],
        lambda: best[1],
        mu_j: best[2],
        delta_j: best[3],
        ..*true_params
    };

    let mut sse_px = 0.0_f64;
    let mut sse_iv = 0.0_f64;
    for &(t, k, px, iv) in &market {
        let f_px = merton_call_series(&fitted_params, k, t);
        let f_iv = implied_vol_call(
            f_px,
            fitted_params.s0,
            k,
            fitted_params.r,
            fitted_params.q,
            t,
        );
        sse_px += (f_px - px).powi(2);
        sse_iv += ((f_iv - iv) * 10_000.0).powi(2);
    }
    let m_cnt = market.len() as f64;

    CalibResult {
        true_params: *true_params,
        init_params: *init_params,
        fitted_params,
        price_rmse: (sse_px / m_cnt).sqrt(),
        iv_rmse_bps: (sse_iv / m_cnt).sqrt(),
        iters,
    }
}
