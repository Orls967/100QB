use day005d_merton_jump_diffusion::calib::calibrate_merton_nelder_mead;
use day005d_merton_jump_diffusion::merton::MertonParams;
use day005d_merton_jump_diffusion::pide::pide_grid_convergence;

fn main() {
    let p = MertonParams::default();
    let conv = pide_grid_convergence(&p, &[80.0, 100.0, 120.0], 1.0, &[(401, 200), (801, 400)]);
    for row in conv {
        println!(
            "PIDE {}x{}: max_abs_err={:.6}, atm_abs_err={:.6}",
            row.n_x, row.n_t, row.max_abs_err, row.atm_abs_err
        );
    }

    let starts = [
        ("low_vol", 0.10, 1.10, -0.15, 0.30),
        ("high_vol", 0.22, 0.40, -0.02, 0.15),
    ];
    for (name, s, l, m, d) in starts {
        let init = MertonParams {
            sigma: s,
            lambda: l,
            mu_j: m,
            delta_j: d,
            ..p
        };
        let cal = calibrate_merton_nelder_mead(
            &p,
            &init,
            &[0.25, 1.0],
            &[80.0, 90.0, 100.0, 110.0, 120.0],
            1200,
        );
        println!(
            "Calib [{name}]: sigma={:.4}, lambda={:.4}, mu_j={:.4}, delta_j={:.4}, price_rmse={:.2e}, iv_rmse={:.3} bps",
            cal.fitted_params.sigma,
            cal.fitted_params.lambda,
            cal.fitted_params.mu_j,
            cal.fitted_params.delta_j,
            cal.price_rmse,
            cal.iv_rmse_bps
        );
    }
}
