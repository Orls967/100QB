//! Bond duration/convexity (closed-form vs bump-and-reprice), Key-Rate Durations,
//! 10y swap PV01, two-bond KRD hedge, and generic cashflow curve risk.

use crate::linalg;
use crate::ns::NsParams;

pub const STANDARD_KRD_NODES: [f64; 8] = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0];

pub fn bond_cashflows(maturity: f64, coupon_rate: f64, freq: f64) -> (Vec<f64>, Vec<f64>) {
    let n = (maturity * freq).round() as usize;
    let dt = 1.0 / freq;
    let mut times = Vec::with_capacity(n);
    let mut cfs = Vec::with_capacity(n);
    for i in 0..n {
        times.push((i + 1) as f64 * dt);
        let mut cf = coupon_rate * dt;
        if i + 1 == n {
            cf += 1.0;
        }
        cfs.push(cf);
    }
    (times, cfs)
}

pub fn price_cashflows_with_shift<F>(
    yield_fn: &dyn Fn(f64) -> f64,
    times: &[f64],
    cfs: &[f64],
    shift_fn: Option<&F>,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let mut pv = 0.0;
    for (&t, &cf) in times.iter().zip(cfs.iter()) {
        let mut y = yield_fn(t);
        if let Some(sf) = shift_fn {
            y += sf(t);
        }
        pv += cf * (-t * y).exp();
    }
    pv
}

pub fn solve_ytm(price: f64, times: &[f64], cfs: &[f64], freq: f64) -> f64 {
    let mut y = 0.04;
    for _ in 0..60 {
        let base = 1.0 + y / freq;
        let mut pv = 0.0;
        let mut dpv = 0.0;
        for (&t, &cf) in times.iter().zip(cfs.iter()) {
            let df = base.powf(-t * freq);
            pv += cf * df;
            dpv -= cf * t * df / base;
        }
        let diff = pv - price;
        if diff.abs() < 1e-14 {
            break;
        }
        y -= diff / dpv;
    }
    y
}

#[derive(Clone, Debug)]
pub struct BondRiskMetrics {
    pub price: f64,
    pub ytm: f64,
    pub mac_duration_ytm: f64,
    pub mod_duration_ytm_closed: f64,
    pub mod_duration_ytm_fd: f64,
    pub mod_duration_ytm_rel_err: f64,
    pub convexity_ytm_closed: f64,
    pub convexity_ytm_fd: f64,
    pub convexity_ytm_rel_err: f64,
    pub eff_duration_zero_closed: f64,
    pub eff_duration_zero_fd: f64,
    pub eff_duration_zero_rel_err: f64,
    pub eff_convexity_zero_closed: f64,
    pub eff_convexity_zero_fd: f64,
    pub eff_convexity_zero_rel_err: f64,
}

pub fn analyze_bond_risk(
    curve: &NsParams,
    maturity: f64,
    coupon_rate: f64,
    freq: f64,
    bump_dur: f64,
    bump_conv: f64,
) -> BondRiskMetrics {
    let (times, cfs) = bond_cashflows(maturity, coupon_rate, freq);
    let yfn = |t: f64| curve.yield_at(t);
    let price = price_cashflows_with_shift::<fn(f64) -> f64>(&yfn, &times, &cfs, None);

    let mut sum_t = 0.0;
    let mut sum_t2 = 0.0;
    for (&t, &cf) in times.iter().zip(cfs.iter()) {
        let pv_cf = cf * (-t * curve.yield_at(t)).exp();
        sum_t += t * pv_cf;
        sum_t2 += t * t * pv_cf;
    }
    let eff_dur_closed = sum_t / price;
    let eff_conv_closed = sum_t2 / price;

    let up_d = |_: f64| bump_dur;
    let dn_d = |_: f64| -bump_dur;
    let p_up = price_cashflows_with_shift(&yfn, &times, &cfs, Some(&up_d));
    let p_dn = price_cashflows_with_shift(&yfn, &times, &cfs, Some(&dn_d));
    let eff_dur_fd = -(p_up - p_dn) / (2.0 * bump_dur * price);

    let up_c = |_: f64| bump_conv;
    let dn_c = |_: f64| -bump_conv;
    let p_upc = price_cashflows_with_shift(&yfn, &times, &cfs, Some(&up_c));
    let p_dnc = price_cashflows_with_shift(&yfn, &times, &cfs, Some(&dn_c));
    let eff_conv_fd = (p_upc - 2.0 * price + p_dnc) / (bump_conv * bump_conv * price);

    let ytm = solve_ytm(price, &times, &cfs, freq);
    let ytm_price = |y: f64| -> f64 {
        let base = 1.0 + y / freq;
        times
            .iter()
            .zip(cfs.iter())
            .map(|(&t, &cf)| cf * base.powf(-t * freq))
            .sum()
    };
    let base = 1.0 + ytm / freq;
    let mut mac_sum = 0.0;
    let mut conv_sum = 0.0;
    for (&t, &cf) in times.iter().zip(cfs.iter()) {
        let pv_cf = cf * base.powf(-t * freq);
        mac_sum += t * pv_cf;
        conv_sum += t * (t + 1.0 / freq) * pv_cf;
    }
    let mac_dur = mac_sum / price;
    let mod_dur_closed = mac_dur / base;
    let conv_ytm_closed = conv_sum / (price * base * base);
    let mod_dur_fd =
        -(ytm_price(ytm + bump_dur) - ytm_price(ytm - bump_dur)) / (2.0 * bump_dur * price);
    let conv_ytm_fd = (ytm_price(ytm + bump_conv) - 2.0 * price + ytm_price(ytm - bump_conv))
        / (bump_conv * bump_conv * price);

    BondRiskMetrics {
        price,
        ytm,
        mac_duration_ytm: mac_dur,
        mod_duration_ytm_closed: mod_dur_closed,
        mod_duration_ytm_fd: mod_dur_fd,
        mod_duration_ytm_rel_err: (mod_dur_fd - mod_dur_closed).abs() / mod_dur_closed.abs(),
        convexity_ytm_closed: conv_ytm_closed,
        convexity_ytm_fd: conv_ytm_fd,
        convexity_ytm_rel_err: (conv_ytm_fd - conv_ytm_closed).abs() / conv_ytm_closed.abs(),
        eff_duration_zero_closed: eff_dur_closed,
        eff_duration_zero_fd: eff_dur_fd,
        eff_duration_zero_rel_err: (eff_dur_fd - eff_dur_closed).abs() / eff_dur_closed.abs(),
        eff_convexity_zero_closed: eff_conv_closed,
        eff_convexity_zero_fd: eff_conv_fd,
        eff_convexity_zero_rel_err: (eff_conv_fd - eff_conv_closed).abs() / eff_conv_closed.abs(),
    }
}

pub fn tent_weight(nodes: &[f64], idx: usize, tau: f64) -> f64 {
    let m = nodes.len();
    if idx == 0 {
        if tau <= nodes[0] {
            return 1.0;
        }
        if m > 1 && tau < nodes[1] {
            return (nodes[1] - tau) / (nodes[1] - nodes[0]);
        }
        return 0.0;
    }
    if idx + 1 == m {
        if tau >= nodes[m - 1] {
            return 1.0;
        }
        if tau > nodes[m - 2] {
            return (tau - nodes[m - 2]) / (nodes[m - 1] - nodes[m - 2]);
        }
        return 0.0;
    }
    if tau <= nodes[idx - 1] || tau >= nodes[idx + 1] {
        return 0.0;
    }
    if tau <= nodes[idx] {
        (tau - nodes[idx - 1]) / (nodes[idx] - nodes[idx - 1])
    } else {
        (nodes[idx + 1] - tau) / (nodes[idx + 1] - nodes[idx])
    }
}

#[derive(Clone, Debug)]
pub struct KrdReport {
    pub nodes: Vec<f64>,
    pub krd_fd: Vec<f64>,
    pub krd_closed: Vec<f64>,
    pub dollar_krd_fd: Vec<f64>,
    pub dollar_krd_closed: Vec<f64>,
    pub sum_krd_fd: f64,
    pub sum_krd_closed: f64,
    pub eff_parallel_fd: f64,
    pub eff_parallel_closed: f64,
    pub sum_vs_par_fd_err: f64,
    pub sum_vs_par_closed_err: f64,
}

pub fn compute_krd_generic(
    yield_fn: &dyn Fn(f64) -> f64,
    times: &[f64],
    cfs: &[f64],
    nodes: &[f64],
    bump: f64,
) -> KrdReport {
    let price = price_cashflows_with_shift::<fn(f64) -> f64>(yield_fn, times, cfs, None);
    let mut krd_fd = Vec::with_capacity(nodes.len());
    let mut krd_cl = Vec::with_capacity(nodes.len());
    let mut dkrd_fd = Vec::with_capacity(nodes.len());
    let mut dkrd_cl = Vec::with_capacity(nodes.len());
    let mut sum_fd = 0.0;
    let mut sum_cl = 0.0;

    for idx in 0..nodes.len() {
        let up = |t: f64| bump * tent_weight(nodes, idx, t);
        let dn = |t: f64| -bump * tent_weight(nodes, idx, t);
        let p_up = price_cashflows_with_shift(yield_fn, times, cfs, Some(&up));
        let p_dn = price_cashflows_with_shift(yield_fn, times, cfs, Some(&dn));
        let k_f = -(p_up - p_dn) / (2.0 * bump * price);
        krd_fd.push(k_f);
        dkrd_fd.push(k_f * price);

        let num: f64 = times
            .iter()
            .zip(cfs.iter())
            .map(|(&t, &cf)| t * tent_weight(nodes, idx, t) * cf * (-t * yield_fn(t)).exp())
            .sum();
        let k_c = num / price;
        krd_cl.push(k_c);
        dkrd_cl.push(num);
        sum_fd += k_f;
        sum_cl += k_c;
    }

    let par_up = |_: f64| bump;
    let par_dn = |_: f64| -bump;
    let p_par_up = price_cashflows_with_shift(yield_fn, times, cfs, Some(&par_up));
    let p_par_dn = price_cashflows_with_shift(yield_fn, times, cfs, Some(&par_dn));
    let eff_par_fd = -(p_par_up - p_par_dn) / (2.0 * bump * price);
    let sum_t: f64 = times
        .iter()
        .zip(cfs.iter())
        .map(|(&t, &cf)| t * cf * (-t * yield_fn(t)).exp())
        .sum();
    let eff_par_cl = sum_t / price;

    KrdReport {
        nodes: nodes.to_vec(),
        krd_fd,
        krd_closed: krd_cl,
        dollar_krd_fd: dkrd_fd,
        dollar_krd_closed: dkrd_cl,
        sum_krd_fd: sum_fd,
        sum_krd_closed: sum_cl,
        eff_parallel_fd: eff_par_fd,
        eff_parallel_closed: eff_par_cl,
        sum_vs_par_fd_err: (sum_fd - eff_par_fd).abs(),
        sum_vs_par_closed_err: (sum_cl - eff_par_cl).abs(),
    }
}

#[derive(Clone, Debug)]
pub struct SwapPv01Report {
    pub maturity: f64,
    pub par_swap_rate: f64,
    pub annuity: f64,
    pub annuity_pv01_per_1: f64,
    pub annuity_pv01_per_1m: f64,
    pub payer_curve_dv01_exact: f64,
    pub payer_curve_dv01_fd: f64,
    pub payer_curve_dv01_1m: f64,
}

pub fn analyze_swap_pv01(curve: &NsParams, maturity: f64, freq: f64) -> SwapPv01Report {
    let n = (maturity * freq).round() as usize;
    let dt = 1.0 / freq;
    let mut annuity = 0.0;
    let mut d_ann_dy = 0.0;
    for i in 1..=n {
        let t = i as f64 * dt;
        let df = curve.discount(t);
        annuity += dt * df;
        d_ann_dy -= dt * t * df;
    }
    let p_mat = curve.discount(maturity);
    let par_rate = (1.0 - p_mat) / annuity;
    let d_pmat_dy = -maturity * p_mat;
    let d_v_dy = -d_pmat_dy - par_rate * d_ann_dy;

    let eval_payer = |shift: f64| -> f64 {
        let mut ann = 0.0;
        for i in 1..=n {
            let t = i as f64 * dt;
            ann += dt * (-t * (curve.yield_at(t) + shift)).exp();
        }
        let p_t = (-maturity * (curve.yield_at(maturity) + shift)).exp();
        1.0 - p_t - par_rate * ann
    };
    let dv01_fd = eval_payer(0.5e-4) - eval_payer(-0.5e-4);

    SwapPv01Report {
        maturity,
        par_swap_rate: par_rate,
        annuity,
        annuity_pv01_per_1: annuity * 1e-4,
        annuity_pv01_per_1m: annuity * 1e-4 * 1_000_000.0,
        payer_curve_dv01_exact: d_v_dy * 1e-4,
        payer_curve_dv01_fd: dv01_fd,
        payer_curve_dv01_1m: d_v_dy * 1e-4 * 1_000_000.0,
    }
}

#[derive(Clone, Debug)]
pub struct TwoBondHedgeReport {
    pub target_name: String,
    pub hedge1_name: String,
    pub hedge2_name: String,
    pub target_node1: f64,
    pub target_node2: f64,
    pub weight_hedge1: f64,
    pub weight_hedge2: f64,
    pub target_dollar_krds: Vec<f64>,
    pub hedged_dollar_krds: Vec<f64>,
    pub residual_node1: f64,
    pub residual_node2: f64,
    pub residual_total_par_d: f64,
    pub target_total_par_d: f64,
}

pub fn hedge_two_krds(
    curve: &NsParams,
    nodes: &[f64],
    node_idx1: usize,
    node_idx2: usize,
    bump: f64,
) -> Result<TwoBondHedgeReport, String> {
    let (t_tgt, cf_tgt) = bond_cashflows(10.0, 0.04, 2.0);
    let (t_h1, cf_h1) = bond_cashflows(5.0, 0.04, 2.0);
    let (t_h2, cf_h2) = bond_cashflows(10.0, 0.02, 2.0);
    let yfn = |t: f64| curve.yield_at(t);
    let k_tgt = compute_krd_generic(&yfn, &t_tgt, &cf_tgt, nodes, bump);
    let k_h1 = compute_krd_generic(&yfn, &t_h1, &cf_h1, nodes, bump);
    let k_h2 = compute_krd_generic(&yfn, &t_h2, &cf_h2, nodes, bump);

    let mat = vec![
        vec![
            k_h1.dollar_krd_closed[node_idx1],
            k_h2.dollar_krd_closed[node_idx1],
        ],
        vec![
            k_h1.dollar_krd_closed[node_idx2],
            k_h2.dollar_krd_closed[node_idx2],
        ],
    ];
    let rhs = vec![
        -k_tgt.dollar_krd_closed[node_idx1],
        -k_tgt.dollar_krd_closed[node_idx2],
    ];
    let w = linalg::solve_lu(&mat, &rhs)?;
    let mut hedged = Vec::with_capacity(nodes.len());
    let mut sum_tgt = 0.0;
    let mut sum_hedged = 0.0;
    for i in 0..nodes.len() {
        let h = k_tgt.dollar_krd_closed[i]
            + w[0] * k_h1.dollar_krd_closed[i]
            + w[1] * k_h2.dollar_krd_closed[i];
        hedged.push(h);
        sum_tgt += k_tgt.dollar_krd_closed[i];
        sum_hedged += h;
    }

    Ok(TwoBondHedgeReport {
        target_name: "10y_4pct_bond".into(),
        hedge1_name: "5y_4pct_bond".into(),
        hedge2_name: "10y_2pct_bond".into(),
        target_node1: nodes[node_idx1],
        target_node2: nodes[node_idx2],
        weight_hedge1: w[0],
        weight_hedge2: w[1],
        target_dollar_krds: k_tgt.dollar_krd_closed,
        hedged_dollar_krds: hedged.clone(),
        residual_node1: hedged[node_idx1],
        residual_node2: hedged[node_idx2],
        residual_total_par_d: sum_hedged,
        target_total_par_d: sum_tgt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_risk_krd_and_hedge() {
        let ns = NsParams::default_ns();
        let rm = analyze_bond_risk(&ns, 10.0, 0.04, 2.0, 1e-6, 1e-4);
        assert!(rm.mod_duration_ytm_rel_err < 1e-6);
        assert!(rm.convexity_ytm_rel_err < 1e-6);
        assert!(rm.eff_duration_zero_rel_err < 1e-6);
        assert!(rm.eff_convexity_zero_rel_err < 1e-6);

        let (times, cfs) = bond_cashflows(10.0, 0.04, 2.0);
        let yfn = |t: f64| ns.yield_at(t);
        let krd = compute_krd_generic(&yfn, &times, &cfs, &STANDARD_KRD_NODES, 1e-6);
        assert!(krd.sum_vs_par_closed_err < 1e-12);
        assert!(krd.sum_vs_par_fd_err < 1e-8);

        let hedge = hedge_two_krds(&ns, &STANDARD_KRD_NODES, 3, 5, 1e-6).unwrap();
        assert!(hedge.residual_node1.abs() < 1e-12);
        assert!(hedge.residual_node2.abs() < 1e-12);
    }
}
