use day001d_asian_options::{
    bs_call, delta_asian_bump, delta_asian_pathwise, fit_loglog_slope, geometric_asian_call,
    price_asian_antithetic, price_asian_control_variate, price_asian_naive, price_european_mc,
    ControlVariateResult, DeltaResult, EstimatorResult, OptionParams,
};
use plotters::prelude::*;
use serde::Serialize;
use std::error::Error;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

#[derive(Serialize)]
struct SanityCheck {
    bs_call: f64,
    mc_call: EstimatorResult,
    diff: f64,
    z_score: f64,
}

#[derive(Serialize)]
struct EstimatorEntry {
    name: String,
    price: f64,
    se: f64,
    ci_low: f64,
    ci_high: f64,
    vrf: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    beta: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    y_exact: Option<f64>,
}

#[derive(Serialize)]
struct DeltaComparison {
    pathwise: DeltaResult,
    bump_1pct: DeltaResult,
    diff: f64,
    z_score: f64,
}

#[derive(Serialize)]
struct ConvergencePoint {
    n: usize,
    se_naive: f64,
    se_antithetic: f64,
    se_control_variate: f64,
}

#[derive(Serialize)]
struct Slopes {
    naive: f64,
    antithetic: f64,
    control_variate: f64,
}

#[derive(Serialize)]
struct TimingsMs {
    sanity_check: f64,
    naive: f64,
    antithetic: f64,
    control_variate: f64,
    delta_pathwise: f64,
    delta_bump: f64,
    convergence: f64,
    total_200k_core: f64,
}

#[derive(Serialize)]
struct FullResults {
    params: OptionParams,
    n_paths: usize,
    seed: u64,
    sanity_check: SanityCheck,
    estimators: Vec<EstimatorEntry>,
    delta: DeltaComparison,
    convergence: Vec<ConvergencePoint>,
    slopes: Slopes,
    timings_ms: TimingsMs,
}

fn find_repo_root() -> Result<PathBuf, Box<dyn Error>> {
    let mut dir = std::env::current_dir()?;
    loop {
        if dir.join("src").join("qb100").exists() {
            return Ok(dir);
        }
        if !dir.pop() {
            break;
        }
    }
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    Ok(manifest
        .parent()
        .and_then(|p| p.parent())
        .map(PathBuf::from)
        .unwrap_or(manifest))
}

fn plot_convergence(
    path: &PathBuf,
    points: &[ConvergencePoint],
    slopes: &Slopes,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(path, (960, 560)).into_drawing_area();
    root.fill(&WHITE)?;

    let min_se = points
        .iter()
        .map(|p| p.se_control_variate)
        .fold(f64::INFINITY, f64::min)
        * 0.7;
    let max_se = points.iter().map(|p| p.se_naive).fold(0.0_f64, f64::max) * 1.4;

    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Konvergensi Standard Error vs Jumlah Path (Log-Log, Rust)",
            ("sans-serif", 22).into_font(),
        )
        .margin(20)
        .x_label_area_size(45)
        .y_label_area_size(70)
        .build_cartesian_2d(
            (900f64..230_000f64).log_scale(),
            (min_se..max_se).log_scale(),
        )?;

    chart
        .configure_mesh()
        .x_desc("Jumlah Path N")
        .y_desc("Standard Error (SE)")
        .draw()?;

    let blue = RGBColor(31, 111, 235);
    let red = RGBColor(209, 73, 91);
    let teal = RGBColor(42, 157, 143);

    chart
        .draw_series(LineSeries::new(
            points.iter().map(|p| (p.n as f64, p.se_naive)),
            blue.stroke_width(2),
        ))?
        .label(format!("Naive MC (slope={:.3})", slopes.naive))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], blue.stroke_width(2)));

    chart
        .draw_series(LineSeries::new(
            points.iter().map(|p| (p.n as f64, p.se_antithetic)),
            red.stroke_width(2),
        ))?
        .label(format!("Antithetic (slope={:.3})", slopes.antithetic))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], red.stroke_width(2)));

    chart
        .draw_series(LineSeries::new(
            points.iter().map(|p| (p.n as f64, p.se_control_variate)),
            teal.stroke_width(2),
        ))?
        .label(format!(
            "Control Variate (slope={:.3})",
            slopes.control_variate
        ))
        .legend(move |(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], teal.stroke_width(2)));

    chart
        .configure_series_labels()
        .background_style(WHITE.mix(0.85))
        .border_style(BLACK)
        .draw()?;

    root.present()?;
    Ok(())
}

fn plot_estimators(path: &PathBuf, estimators: &[EstimatorEntry]) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(path, (960, 480)).into_drawing_area();
    root.fill(&WHITE)?;
    let (left, right) = root.split_horizontally(480);

    // Left: Price & 95% CI
    let min_p = estimators
        .iter()
        .map(|e| e.ci_low)
        .fold(f64::INFINITY, f64::min)
        - 0.01;
    let max_p = estimators
        .iter()
        .map(|e| e.ci_high)
        .fold(f64::NEG_INFINITY, f64::max)
        + 0.01;

    let mut chart_l = ChartBuilder::on(&left)
        .caption("Harga + CI 95% (200k Path)", ("sans-serif", 18).into_font())
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(55)
        .build_cartesian_2d(-0.5f64..2.5f64, min_p..max_p)?;

    chart_l
        .configure_mesh()
        .x_labels(3)
        .x_label_formatter(&|v| match v.round() as i32 {
            0 => "Naive".to_string(),
            1 => "Antithetic".to_string(),
            2 => "CV".to_string(),
            _ => String::new(),
        })
        .y_desc("Harga Opsi Asia")
        .draw()?;

    let colors = [
        RGBColor(31, 111, 235),
        RGBColor(209, 73, 91),
        RGBColor(42, 157, 143),
    ];

    for (i, e) in estimators.iter().enumerate() {
        let x = i as f64;
        let c = colors[i];
        chart_l.draw_series(std::iter::once(PathElement::new(
            vec![(x, e.ci_low), (x, e.ci_high)],
            c.stroke_width(3),
        )))?;
        chart_l.draw_series(std::iter::once(Circle::new((x, e.price), 5, c.filled())))?;
    }

    // Right: Variance Reduction Factor (log scale)
    let max_vrf = estimators.iter().map(|e| e.vrf).fold(1.0_f64, f64::max) * 2.0;

    let mut chart_r = ChartBuilder::on(&right)
        .caption(
            "Variance Reduction Factor (Log)",
            ("sans-serif", 18).into_font(),
        )
        .margin(15)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (0usize..3usize).into_segmented(),
            (0.8f64..max_vrf).log_scale(),
        )?;

    chart_r
        .configure_mesh()
        .x_label_formatter(&|v| match v {
            SegmentValue::CenterOf(0) => "Naive".to_string(),
            SegmentValue::CenterOf(1) => "Antithetic".to_string(),
            SegmentValue::CenterOf(2) => "CV".to_string(),
            _ => String::new(),
        })
        .y_desc("Var(Naive) / Var(Metode)")
        .draw()?;

    chart_r.draw_series(estimators.iter().enumerate().map(|(i, e)| {
        Rectangle::new(
            [
                (SegmentValue::Exact(i), 0.8),
                (SegmentValue::Exact(i + 1), e.vrf),
            ],
            colors[i].filled(),
        )
    }))?;

    root.present()?;
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let params = OptionParams {
        s0: 100.0,
        k: 100.0,
        r: 0.05,
        sigma: 0.20,
        t: 1.0,
        n_steps: 252,
    };
    let n_paths = 200_000usize;
    let seed = 42u64;

    let t_total_start = Instant::now();

    let t0 = Instant::now();
    let bs_euro = bs_call(&params);
    let mc_euro = price_european_mc(&params, n_paths, seed)?;
    let t_sanity = t0.elapsed().as_secs_f64() * 1000.0;
    let diff_euro = mc_euro.price - bs_euro;
    let z_euro = diff_euro.abs() / mc_euro.se;

    let t0 = Instant::now();
    let res_naive = price_asian_naive(&params, n_paths, seed)?;
    let t_naive = t0.elapsed().as_secs_f64() * 1000.0;

    let t0 = Instant::now();
    let res_anti = price_asian_antithetic(&params, n_paths, seed)?;
    let t_anti = t0.elapsed().as_secs_f64() * 1000.0;

    let t0 = Instant::now();
    let res_cv: ControlVariateResult = price_asian_control_variate(&params, n_paths, seed)?;
    let t_cv = t0.elapsed().as_secs_f64() * 1000.0;

    let t0 = Instant::now();
    let delta_pw = delta_asian_pathwise(&params, n_paths, seed)?;
    let t_pw = t0.elapsed().as_secs_f64() * 1000.0;

    let t0 = Instant::now();
    let delta_bmp = delta_asian_bump(&params, n_paths, 0.01, seed)?;
    let t_bmp = t0.elapsed().as_secs_f64() * 1000.0;

    let core_ms = t_total_start.elapsed().as_secs_f64() * 1000.0;

    let delta_diff = delta_pw.delta - delta_bmp.delta;
    let delta_z = delta_diff.abs() / delta_pw.se.hypot(delta_bmp.se);

    let var_naive = res_naive.se * res_naive.se;
    let estimators = vec![
        EstimatorEntry {
            name: "Naive MC (200k)".to_string(),
            price: res_naive.price,
            se: res_naive.se,
            ci_low: res_naive.ci_low,
            ci_high: res_naive.ci_high,
            vrf: 1.0,
            beta: None,
            y_exact: None,
        },
        EstimatorEntry {
            name: "Antithetic (100k pasang)".to_string(),
            price: res_anti.price,
            se: res_anti.se,
            ci_low: res_anti.ci_low,
            ci_high: res_anti.ci_high,
            vrf: var_naive / (res_anti.se * res_anti.se),
            beta: None,
            y_exact: None,
        },
        EstimatorEntry {
            name: "Control Variate (Geom. Asia)".to_string(),
            price: res_cv.price,
            se: res_cv.se,
            ci_low: res_cv.ci_low,
            ci_high: res_cv.ci_high,
            vrf: var_naive / (res_cv.se * res_cv.se),
            beta: Some(res_cv.beta),
            y_exact: Some(res_cv.y_exact),
        },
    ];

    let t0 = Instant::now();
    let n_grid = [
        1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000,
    ];
    let mut conv_points = Vec::with_capacity(n_grid.len());
    let mut se_n_vec = Vec::with_capacity(n_grid.len());
    let mut se_a_vec = Vec::with_capacity(n_grid.len());
    let mut se_c_vec = Vec::with_capacity(n_grid.len());

    for &n in &n_grid {
        let rn = price_asian_naive(&params, n, seed)?;
        let ra = price_asian_antithetic(&params, n, seed)?;
        let rc = price_asian_control_variate(&params, n, seed)?;
        se_n_vec.push(rn.se);
        se_a_vec.push(ra.se);
        se_c_vec.push(rc.se);
        conv_points.push(ConvergencePoint {
            n,
            se_naive: rn.se,
            se_antithetic: ra.se,
            se_control_variate: rc.se,
        });
    }
    let t_conv = t0.elapsed().as_secs_f64() * 1000.0;

    let slopes = Slopes {
        naive: fit_loglog_slope(&n_grid, &se_n_vec)?,
        antithetic: fit_loglog_slope(&n_grid, &se_a_vec)?,
        control_variate: fit_loglog_slope(&n_grid, &se_c_vec)?,
    };

    let results = FullResults {
        params,
        n_paths,
        seed,
        sanity_check: SanityCheck {
            bs_call: bs_euro,
            mc_call: mc_euro,
            diff: diff_euro,
            z_score: z_euro,
        },
        estimators,
        delta: DeltaComparison {
            pathwise: delta_pw,
            bump_1pct: delta_bmp,
            diff: delta_diff,
            z_score: delta_z,
        },
        convergence: conv_points,
        slopes,
        timings_ms: TimingsMs {
            sanity_check: t_sanity,
            naive: t_naive,
            antithetic: t_anti,
            control_variate: t_cv,
            delta_pathwise: t_pw,
            delta_bump: t_bmp,
            convergence: t_conv,
            total_200k_core: core_ms,
        },
    };

    let root = find_repo_root()?;
    let out_dir = root.join("outputs").join("day001d");
    let fig_dir = root.join("outputs").join("figures");
    fs::create_dir_all(&out_dir)?;
    fs::create_dir_all(&fig_dir)?;

    let json_path = out_dir.join("results.json");
    fs::write(&json_path, serde_json::to_string_pretty(&results)?)?;

    plot_convergence(
        &fig_dir.join("day001d_asian-convergence-loglog.png"),
        &results.convergence,
        &results.slopes,
    )?;
    plot_estimators(
        &fig_dir.join("day001d_asian-estimators-ci.png"),
        &results.estimators,
    )?;

    println!(
        "BS Eropa: {:.6} | MC Eropa: {:.6} ± {:.6} ({:.2} SE)",
        bs_euro, mc_euro.price, mc_euro.se, z_euro
    );
    println!(
        "Geom. Asian exact: {:.6} | beta: {:.6}",
        geometric_asian_call(&params)?,
        res_cv.beta
    );
    for e in &results.estimators {
        println!(
            "{:28}: {:.6} ± {:.6} (CI [{:.6}, {:.6}]) | VRF = {:.2}x",
            e.name, e.price, e.se, e.ci_low, e.ci_high, e.vrf
        );
    }
    println!(
        "Delta PW: {:.6} ± {:.6} | Bump: {:.6} ± {:.6} | diff = {:+.6} ({:.2} SE)",
        delta_pw.delta, delta_pw.se, delta_bmp.delta, delta_bmp.se, delta_diff, delta_z
    );
    println!(
        "Slopes -> Naive: {:.4}, Antithetic: {:.4}, CV: {:.4}",
        results.slopes.naive, results.slopes.antithetic, results.slopes.control_variate
    );
    println!("Total 200k core time: {:.1} ms", core_ms);

    Ok(())
}
