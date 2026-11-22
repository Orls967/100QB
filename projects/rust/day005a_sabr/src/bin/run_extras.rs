use std::path::Path;
use std::time::Instant;

use day005a_sabr::cev::{cev_call_exact, cev_call_pde_cn};
use day005a_sabr::density::scan_hagan_negative_density;
use day005a_sabr::jsonout::{obj, write_json_file};
use day005a_sabr::sabr_mc::default_sabr_workers;

fn repo_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

fn main() {
    let t0 = Instant::now();
    let n_workers = default_sabr_workers();
    let root = repo_root();
    let out_json = root.join("outputs/day005a/extras.json");

    let f0 = 0.04_f64;
    let alpha0 = 0.0325_f64;
    let mut cev_full_grid = Vec::new();
    for &b in &[0.3_f64, 0.5, 0.7] {
        for &t_h in &[1.0_f64, 5.0, 10.0] {
            for &k in &[0.01_f64, 0.02, 0.04, 0.06, 0.08] {
                let c_ex = cev_call_exact(f0, k, alpha0, b, t_h);
                let p1 = cev_call_pde_cn(f0, k, alpha0, b, t_h, 800, 400);
                let p2 = cev_call_pde_cn(f0, k, alpha0, b, t_h, 1600, 800);
                let c_pde = (4.0 * p2 - p1) / 3.0;
                cev_full_grid.push((
                    b,
                    t_h,
                    k,
                    c_ex * 1e4,
                    c_pde * 1e4,
                    (c_ex - c_pde).abs() / c_ex.max(1e-12),
                ));
            }
        }
    }

    let mut density_nu_sweep = Vec::new();
    for &nu_v in &[0.2_f64, 0.4, 0.8] {
        let sc = scan_hagan_negative_density(f0, 10.0, alpha0, 0.5, -0.3, nu_v, 0.0005, 0.08, 400);
        density_nu_sweep.push((nu_v, sc.min_q, sc.total_neg_mass));
    }

    let total_secs = t0.elapsed().as_secs_f64();
    let doc = obj([
        ("total_seconds", total_secs.into()),
        ("workers_used", n_workers.into()),
        (
            "cev_full_pde_grid",
            cev_full_grid
                .iter()
                .map(|r| {
                    obj([
                        ("beta", r.0.into()),
                        ("t", r.1.into()),
                        ("k", r.2.into()),
                        ("exact_bps", r.3.into()),
                        ("pde_bps", r.4.into()),
                        ("rel_err", r.5.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
        (
            "density_nu_sweep_t10",
            density_nu_sweep
                .iter()
                .map(|r| {
                    obj([
                        ("nu", r.0.into()),
                        ("min_q", r.1.into()),
                        ("neg_mass", r.2.into()),
                    ])
                })
                .collect::<Vec<_>>()
                .into(),
        ),
    ]);
    write_json_file(&out_json, &doc).expect("write day005a extras.json");
    println!("day005a run_extras selesai dalam {:.2}s", total_secs);
}
