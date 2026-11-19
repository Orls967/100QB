use std::path::Path;

use plotters::prelude::*;

use crate::grid::GridPosterior;

pub fn plot_pf_unbiased_and_var(
    out_path: &Path,
    ll_diffs: &[f64],
    np_grid: &[usize],
    var_ll: &[f64],
    slope: f64,
) -> Result<(), String> {
    let root = BitMapBackend::new(out_path, (1080, 460)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| e.to_string())?;
    let areas = root.split_evenly((1, 2));

    // Left: histogram ll_hat - ll_exact
    let min_d = ll_diffs.iter().copied().fold(f64::INFINITY, f64::min);
    let max_d = ll_diffs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let n_bins = 35usize;
    let bw = (max_d - min_d).max(1e-6) / (n_bins as f64);
    let mut counts = vec![0usize; n_bins];
    for &d in ll_diffs {
        let b = (((d - min_d) / bw).floor() as usize).min(n_bins - 1);
        counts[b] += 1;
    }
    let max_c = (*counts.iter().max().unwrap_or(&1) as f64) * 1.1;

    let mut ch1 = ChartBuilder::on(&areas[0])
        .caption(
            "PF Error Distribution (ll_hat - ll_exact, Np=500)",
            ("sans-serif", 18),
        )
        .margin(18)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(min_d..max_d, 0.0..max_c)
        .map_err(|e| e.to_string())?;
    ch1.configure_mesh()
        .x_desc("ll_hat - ll_exact")
        .y_desc("Count")
        .draw()
        .map_err(|e| e.to_string())?;
    ch1.draw_series((0..n_bins).map(|b| {
        let x0 = min_d + (b as f64) * bw;
        let x1 = x0 + bw;
        Rectangle::new([(x0, 0.0), (x1, counts[b] as f64)], BLUE.mix(0.6).filled())
    }))
    .map_err(|e| e.to_string())?;

    // Right: log10(Var) vs log10(Np)
    let xs: Vec<f64> = np_grid.iter().map(|&n| (n as f64).log10()).collect();
    let ys: Vec<f64> = var_ll.iter().map(|&v| v.max(1e-12).log10()).collect();
    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min) - 0.1;
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max) + 0.1;
    let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min) - 0.2;
    let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max) + 0.2;

    let mut ch2 = ChartBuilder::on(&areas[1])
        .caption(
            format!("Var(ll_hat) vs Np (log-log, slope = {slope:.3})"),
            ("sans-serif", 18),
        )
        .margin(18)
        .x_label_area_size(40)
        .y_label_area_size(55)
        .build_cartesian_2d(x_min..x_max, y_min..y_max)
        .map_err(|e| e.to_string())?;
    ch2.configure_mesh()
        .x_desc("log10(Np)")
        .y_desc("log10(Var(ll_hat))")
        .draw()
        .map_err(|e| e.to_string())?;
    ch2.draw_series(LineSeries::new(
        xs.iter().copied().zip(ys.iter().copied()),
        &RED,
    ))
    .map_err(|e| e.to_string())?;
    ch2.draw_series(
        xs.iter()
            .copied()
            .zip(ys.iter().copied())
            .map(|pt| Circle::new(pt, 4, RED.filled())),
    )
    .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())
}

fn kde_hist(samples: &[f64], grid: &[f64]) -> Vec<f64> {
    let min_x = grid[0];
    let max_x = *grid.last().unwrap_or(&0.0);
    let step = (max_x - min_x) / ((grid.len() - 1) as f64);
    let mut counts = vec![0.0; grid.len()];
    for &s in samples {
        let idx = ((s - min_x) / step).round() as isize;
        if idx >= 0 && (idx as usize) < grid.len() {
            counts[idx as usize] += 1.0;
        }
    }
    let norm = counts.iter().sum::<f64>() * step;
    if norm > 0.0 {
        for c in &mut counts {
            *c /= norm;
        }
    }
    counts
}

pub fn plot_model_a_posteriors(
    out_path: &Path,
    grid: &GridPosterior,
    mh_a: &[f64],
    mh_b: &[f64],
    pmmh_a: &[f64],
    pmmh_b: &[f64],
) -> Result<(), String> {
    let root = BitMapBackend::new(out_path, (1080, 480)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| e.to_string())?;
    let areas = root.split_evenly((1, 2));

    let draw_panel = |area: &DrawingArea<BitMapBackend, plotters::coord::Shift>,
                      xs: &[f64],
                      g_pdf: &[f64],
                      mh_s: &[f64],
                      pm_s: &[f64],
                      title: &str,
                      xlab: &str|
     -> Result<(), String> {
        let mh_d = kde_hist(mh_s, xs);
        let pm_d = kde_hist(pm_s, xs);
        let max_y = g_pdf.iter().copied().fold(0.0_f64, f64::max) * 1.25;
        let x0 = grid_slice_min(xs);
        let x1 = grid_slice_max(xs);

        let mut ch = ChartBuilder::on(area)
            .caption(title, ("sans-serif", 18))
            .margin(18)
            .x_label_area_size(40)
            .y_label_area_size(50)
            .build_cartesian_2d(x0..x1, 0.0..max_y)
            .map_err(|e| e.to_string())?;
        ch.configure_mesh()
            .x_desc(xlab)
            .y_desc("Density")
            .draw()
            .map_err(|e| e.to_string())?;
        ch.draw_series(LineSeries::new(
            xs.iter().copied().zip(g_pdf.iter().copied()),
            &BLACK,
        ))
        .map_err(|e| e.to_string())?
        .label("Grid Eksak")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], BLACK));
        ch.draw_series(LineSeries::new(xs.iter().copied().zip(mh_d), &BLUE))
            .map_err(|e| e.to_string())?
            .label("MH Eksak")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], BLUE));
        ch.draw_series(LineSeries::new(xs.iter().copied().zip(pm_d), &RED))
            .map_err(|e| e.to_string())?
            .label("PMMH")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], RED));
        ch.configure_series_labels()
            .background_style(WHITE.mix(0.8))
            .border_style(BLACK)
            .draw()
            .map_err(|e| e.to_string())?;
        Ok(())
    };

    draw_panel(
        &areas[0],
        &grid.a_grid,
        &grid.marginal_a,
        mh_a,
        pmmh_a,
        "Model A Posterior: a = ln(sigma_eta)",
        "a",
    )?;
    draw_panel(
        &areas[1],
        &grid.b_grid,
        &grid.marginal_b,
        mh_b,
        pmmh_b,
        "Model A Posterior: b = ln(sigma_eps)",
        "b",
    )?;
    root.present().map_err(|e| e.to_string())
}

fn grid_slice_min(xs: &[f64]) -> f64 {
    xs.first().copied().unwrap_or(0.0)
}

fn grid_slice_max(xs: &[f64]) -> f64 {
    xs.last().copied().unwrap_or(1.0)
}

pub fn plot_model_b_sv(
    out_path: &Path,
    true_h: &[f64],
    filt_mean: &[f64],
    filt_lo: &[f64],
    filt_hi: &[f64],
    phi_samples: &[f64],
    sigma_samples: &[f64],
) -> Result<(), String> {
    let root = BitMapBackend::new(out_path, (1080, 480)).into_drawing_area();
    root.fill(&WHITE).map_err(|e| e.to_string())?;
    let areas = root.split_evenly((1, 2));

    let min_h = filt_lo
        .iter()
        .chain(true_h.iter())
        .copied()
        .fold(f64::INFINITY, f64::min)
        - 0.2;
    let max_h = filt_hi
        .iter()
        .chain(true_h.iter())
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        + 0.2;
    let mut ch1 = ChartBuilder::on(&areas[0])
        .caption("Model B SV: Filtered h_t vs True h_t", ("sans-serif", 18))
        .margin(18)
        .x_label_area_size(40)
        .y_label_area_size(50)
        .build_cartesian_2d(0.0..(true_h.len() as f64), min_h..max_h)
        .map_err(|e| e.to_string())?;
    ch1.configure_mesh()
        .x_desc("t")
        .y_desc("h_t")
        .draw()
        .map_err(|e| e.to_string())?;
    ch1.draw_series(AreaSeries::new(
        (0..true_h.len()).map(|t| (t as f64, filt_hi[t])),
        min_h,
        BLUE.mix(0.15),
    ))
    .map_err(|e| e.to_string())?;
    ch1.draw_series(LineSeries::new(
        (0..true_h.len()).map(|t| (t as f64, true_h[t])),
        &BLACK,
    ))
    .map_err(|e| e.to_string())?
    .label("h_t sejati")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], BLACK));
    ch1.draw_series(LineSeries::new(
        (0..true_h.len()).map(|t| (t as f64, filt_mean[t])),
        &BLUE,
    ))
    .map_err(|e| e.to_string())?
    .label("Filtered mean")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 18, y)], BLUE));
    ch1.configure_series_labels()
        .background_style(WHITE.mix(0.8))
        .border_style(BLACK)
        .draw()
        .map_err(|e| e.to_string())?;

    let mut ch2 = ChartBuilder::on(&areas[1])
        .caption(
            "Model B Posterior Scatter: (phi, sigma)",
            ("sans-serif", 18),
        )
        .margin(18)
        .x_label_area_size(40)
        .y_label_area_size(55)
        .build_cartesian_2d(0.80..0.995, 0.05..0.55)
        .map_err(|e| e.to_string())?;
    ch2.configure_mesh()
        .x_desc("phi")
        .y_desc("sigma")
        .draw()
        .map_err(|e| e.to_string())?;
    let step = (phi_samples.len() / 2500).max(1);
    ch2.draw_series(
        phi_samples
            .iter()
            .zip(sigma_samples.iter())
            .step_by(step)
            .map(|(&p, &s)| Circle::new((p, s), 2, BLUE.mix(0.35).filled())),
    )
    .map_err(|e| e.to_string())?;
    ch2.draw_series(std::iter::once(Circle::new((0.95, 0.25), 6, RED.filled())))
        .map_err(|e| e.to_string())?;

    root.present().map_err(|e| e.to_string())
}
