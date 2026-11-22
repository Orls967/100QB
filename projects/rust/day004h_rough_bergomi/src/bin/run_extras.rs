use std::path::PathBuf;
use std::time::Instant;

use day004h_rough_bergomi::calib::run_calibration_a6_extras;
use day004h_rough_bergomi::jsonout::{write_json_file, JsonVal};
use day004h_rough_bergomi::parallel;
use day004h_rough_bergomi::rbergomi::RBergomiParams;
use day004h_rough_bergomi::rng::MAIN_SEED;

fn find_repo_root() -> PathBuf {
    let mut dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
    for _ in 0..8 {
        if dir.join("outputs").is_dir() && dir.join("projects").join("rust").is_dir() {
            return dir;
        }
        if !dir.pop() {
            break;
        }
    }
    PathBuf::from(".")
}

fn rep_json(r: &day004h_rough_bergomi::calib::ReplicationSummary) -> JsonVal {
    JsonVal::Obj(vec![
        ("n_rep", r.n_rep.into()),
        ("bias_h", r.bias_h.into()),
        ("sd_h", r.sd_h.into()),
        ("bias_eta", r.bias_eta.into()),
        ("sd_eta", r.sd_eta.into()),
        ("bias_rho", r.bias_rho.into()),
        ("sd_rho", r.sd_rho.into()),
    ])
}

fn fit_json(f: &day004h_rough_bergomi::calib::CalibResult) -> JsonVal {
    JsonVal::Obj(vec![
        ("h_hat", f.h_hat.into()),
        ("eta_hat", f.eta_hat.into()),
        ("rho_hat", f.rho_hat.into()),
        ("rmse_opt", f.rmse_opt.into()),
    ])
}

fn main() {
    let t0 = Instant::now();
    let n_workers = parallel::default_workers();
    let root = find_repo_root();
    let out_json = root.join("outputs").join("day004h").join("extras.json");
    let p = RBergomiParams::default_params();

    let ex = run_calibration_a6_extras(p, n_workers, MAIN_SEED);
    let total_secs = t0.elapsed().as_secs_f64();

    let doc = JsonVal::Obj(vec![
        ("total_seconds", total_secs.into()),
        ("workers_used", n_workers.into()),
        ("self_consistency_fit", fit_json(&ex.self_consistency_fit)),
        ("fit_r80k", fit_json(&ex.fit_r80k)),
        ("fit_vega_price", fit_json(&ex.fit_vega_price)),
        ("fit_long_only_t_ge_05", fit_json(&ex.fit_long_only)),
        ("fit_short_only_t_le_01", fit_json(&ex.fit_short_only)),
        ("rep_old_diff_n", rep_json(&ex.rep_old_diff_n)),
        ("rep_blp_diff_n", rep_json(&ex.rep_blp_diff_n)),
        ("rep_diff_crn", rep_json(&ex.rep_diff_crn)),
        ("rep_cond_bs_fixed_crn", rep_json(&ex.rep_cond_bs_fixed_crn)),
        ("rep_cond_bs_diff_crn", rep_json(&ex.rep_cond_bs_diff_crn)),
        (
            "true_rmse_grid_10seed_stats",
            ex.true_rmse_grid_stats
                .iter()
                .map(|s| {
                    JsonVal::Obj(vec![
                        ("n_steps", s.n_steps.into()),
                        ("n_paths", s.n_paths.into()),
                        ("n_seeds", s.n_seeds.into()),
                        ("median_rmse", s.median_rmse.into()),
                        ("iqr_rmse", s.iqr_rmse.into()),
                        ("min_rmse", s.min_rmse.into()),
                        ("max_rmse", s.max_rmse.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
        (
            "hill_tail_estimates_t1",
            ex.hill_estimates
                .iter()
                .map(|h| {
                    JsonVal::Obj(vec![
                        ("k", h.k.into()),
                        ("alpha_st", h.alpha_st.into()),
                        ("alpha_payoff2", h.alpha_payoff2.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
    ]);

    write_json_file(&out_json, &doc).expect("write day004h extras.json");
    println!("day004h run_extras selesai dalam {:.2}s", total_secs);
}
