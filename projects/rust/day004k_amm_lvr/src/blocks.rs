use crate::parallel::reduce_tuple4_fixed_blocks;
use crate::rng::{derive_seed, Rng};
use crate::special::{gauss_legendre, norm_pdf};

/// Exact single-block arbitrage outcome for log-price shock `z = ln(P' / p)` and pool fee `phi` (paid on input).
/// Returns `(arb_profit, lp_fee_income, arb_volume_numeraire, new_pool_price)` normalized to `L = 0.5, p = 1.0` (so `V(p) = 1.0`).
#[inline]
pub fn block_arb_outcome(z: f64, phi: f64, l: f64, p: f64) -> (f64, f64, f64, f64) {
    let one_m_phi = 1.0 - phi;
    let thresh = -one_m_phi.ln(); // -ln(1 - phi) >= 0
    let scale = l * p.sqrt();
    if z > thresh {
        // P' > p / (1 - phi): arbitrageur pays dy in numeraire to buy dx
        let u = (one_m_phi * z.exp()).sqrt();
        let dy_in = (scale / one_m_phi) * (u - 1.0);
        let profit = (scale / one_m_phi) * (u - 1.0) * (u - 1.0);
        let fee = phi * dy_in;
        let p_new = one_m_phi * p * z.exp();
        (profit, fee, dy_in, p_new)
    } else if z < -thresh {
        // P' < p * (1 - phi): arbitrageur sells dx into pool for dy
        let u = (one_m_phi * (-z).exp()).sqrt();
        let dx_in = (l / (one_m_phi * p.sqrt())) * (u - 1.0);
        let profit = scale * (u - 1.0) * (u - 1.0) / (u * u);
        let p_ext = p * z.exp();
        let fee = phi * dx_in * p_ext;
        let vol_num = dx_in * p_ext;
        let p_new = p_ext / one_m_phi;
        (profit, fee, vol_num, p_new)
    } else {
        (0.0, 0.0, 0.0, p)
    }
}

/// Exact expected `(arb_profit, lp_fee_income, arb_volume)` per block via Gauss-Legendre quadrature
/// over `z ~ N(-0.5 * sigma^2 * dt, sigma^2 * dt)`.
pub fn exact_block_expectations_gl(sigma: f64, dt: f64, phi: f64, l: f64) -> (f64, f64, f64) {
    let mean_z = -0.5 * sigma * sigma * dt;
    let sd_z = sigma * dt.sqrt();
    let thresh = -(1.0 - phi).ln();
    let z_max = (mean_z.abs() + 9.0 * sd_z).max(thresh + 9.0 * sd_z);

    let (nodes, weights) = gauss_legendre(64);
    let integrate_interval = |a: f64, b: f64| -> (f64, f64, f64) {
        if b <= a {
            return (0.0, 0.0, 0.0);
        }
        let half = 0.5 * (b - a);
        let mid = 0.5 * (a + b);
        let mut s_p = 0.0_f64;
        let mut s_f = 0.0_f64;
        let mut s_v = 0.0_f64;
        for (&xi, &wi) in nodes.iter().zip(weights.iter()) {
            let z = mid + half * xi;
            let std_z = (z - mean_z) / sd_z;
            let dens = norm_pdf(std_z) / sd_z;
            let (pr, fe, vo, _) = block_arb_outcome(z, phi, l, 1.0);
            let w = wi * half * dens;
            s_p += w * pr;
            s_f += w * fe;
            s_v += w * vo;
        }
        (s_p, s_f, s_v)
    };

    let (p1, f1, v1) = integrate_interval(-z_max, -thresh);
    let (p2, f2, v2) = integrate_interval(thresh, z_max);
    (p1 + p2, f1 + f2, v1 + v2)
}

/// Monte Carlo per-block `(arb_profit, lp_fee_income, arb_volume)` with SEs (annualized by dividing by `dt`).
pub fn mc_block_annualized(
    sigma: f64,
    dt: f64,
    phi: f64,
    n_blocks: usize,
    n_workers: usize,
    seed: u64,
) -> (f64, f64, f64, f64, f64, f64) {
    let l = 0.5_f64;
    let mean_z = -0.5 * sigma * sigma * dt;
    let sd_z = sigma * dt.sqrt();
    let inv_dt = 1.0 / dt;
    let sums = reduce_tuple4_fixed_blocks(n_blocks, n_workers, |r| {
        let mut rg = Rng::new(derive_seed(seed, 0x2001, r as u64));
        let z = mean_z + sd_z * rg.normal();
        let (pr, fe, vo, _) = block_arb_outcome(z, phi, l, 1.0);
        let pr_a = pr * inv_dt;
        let fe_a = fe * inv_dt;
        let vo_a = vo * inv_dt;
        [pr_a, fe_a, vo_a, pr_a * pr_a]
    });
    let nf = n_blocks as f64;
    let m_pr = sums[0] / nf;
    let m_fe = sums[1] / nf;
    let m_vo = sums[2] / nf;
    let var_pr = (sums[3] / nf - m_pr * m_pr).max(0.0);
    let se_pr = (var_pr / nf).sqrt();
    (m_pr, se_pr, m_fe, m_vo, sums[1] / nf, sums[2] / nf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_block_gl_matches_mc_and_zero_fee_limit() {
        let sec_yr = 365.0 * 86_400.0;
        let dt = 12.0 / sec_yr;
        let (ep, _ef, _ev) = exact_block_expectations_gl(0.8, dt, 0.0, 0.5);
        let ann_lvr = ep / dt;
        assert!((ann_lvr - 0.08).abs() < 1e-4, "ann_lvr={ann_lvr}");

        let (gp, gf, gv) = exact_block_expectations_gl(0.8, 600.0 / sec_yr, 0.003, 0.5);
        let (mp, se_p, mf, mv, _, _) =
            mc_block_annualized(0.8, 600.0 / sec_yr, 0.003, 400_000, 4, 20261122);
        let scale = 600.0 / sec_yr;
        assert!(
            (gp / scale - mp).abs() <= 3.0 * se_p,
            "gl={} mc={} se={se_p}",
            gp / scale,
            mp
        );
        assert!((gf / scale - mf).abs() < 0.005);
        assert!((gv / scale - mv).abs() < 0.5);
    }
}
