use std::fmt::Write as _;
use std::fs;
use std::path::Path;

pub struct Series<'a> {
    pub label: &'a str,
    pub x: &'a [f64],
    pub y: &'a [f64],
    pub color: &'a str,
    pub dashed: bool,
}

pub fn write_line_chart(
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[Series<'_>],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let w = 820.0_f64;
    let h = 480.0_f64;
    let ml = 82.0_f64;
    let mr = 30.0_f64;
    let mt = 56.0_f64;
    let mb = 70.0_f64;
    let pw = w - ml - mr;
    let ph = h - mt - mb;

    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for s in series {
        for (&xv, &yv) in s.x.iter().zip(s.y.iter()) {
            if xv.is_finite() && yv.is_finite() {
                xmin = xmin.min(xv);
                xmax = xmax.max(xv);
                ymin = ymin.min(yv);
                ymax = ymax.max(yv);
            }
        }
    }
    if (xmax - xmin).abs() < 1e-12 {
        xmin -= 1.0;
        xmax += 1.0;
    }
    if (ymax - ymin).abs() < 1e-12 {
        ymin -= 1.0;
        ymax += 1.0;
    }
    let ypad = 0.08 * (ymax - ymin);
    ymin -= ypad;
    ymax += ypad;

    let sx = |x: f64| ml + (x - xmin) / (xmax - xmin) * pw;
    let sy = |y: f64| mt + ph - (y - ymin) / (ymax - ymin) * ph;

    let mut svg = String::with_capacity(8192);
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"##
    );
    let _ = writeln!(svg, r##"<rect width="{w}" height="{h}" fill="#ffffff"/>"##);
    let _ = writeln!(
        svg,
        r##"<text x="{}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold" fill="#111827">{title}</text>"##,
        w * 0.5
    );

    for k in 0..=5 {
        let frac = k as f64 / 5.0;
        let yv = ymin + frac * (ymax - ymin);
        let yp = sy(yv);
        let _ = writeln!(
            svg,
            r##"<line x1="{ml}" y1="{yp:.1}" x2="{:.1}" y2="{yp:.1}" stroke="#e5e7eb" stroke-width="1"/>"##,
            ml + pw
        );
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#4b5563">{yv:.4}</text>"##,
            ml - 8.0,
            yp + 4.0
        );
        let xv = xmin + frac * (xmax - xmin);
        let xp = sx(xv);
        let _ = writeln!(
            svg,
            r##"<line x1="{xp:.1}" y1="{mt}" x2="{xp:.1}" y2="{:.1}" stroke="#e5e7eb" stroke-width="1"/>"##,
            mt + ph
        );
        let _ = writeln!(
            svg,
            r##"<text x="{xp:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#4b5563">{xv:.3}</text>"##,
            mt + ph + 18.0
        );
    }

    let _ = writeln!(
        svg,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#374151" stroke-width="1.2"/>"##
    );
    let _ = writeln!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#1f2937">{x_label}</text>"##,
        ml + pw * 0.5,
        h - 16.0
    );
    let _ = writeln!(
        svg,
        r##"<text x="20" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#1f2937" transform="rotate(-90 20 {:.1})">{y_label}</text>"##,
        mt + ph * 0.5,
        mt + ph * 0.5
    );

    for s in series {
        let mut pts = String::new();
        for (&xv, &yv) in s.x.iter().zip(s.y.iter()) {
            if xv.is_finite() && yv.is_finite() {
                let _ = write!(pts, "{:.2},{:.2} ", sx(xv), sy(yv));
            }
        }
        let dash = if s.dashed {
            r##" stroke-dasharray="6,4""##
        } else {
            ""
        };
        let _ = writeln!(
            svg,
            r##"<polyline fill="none" stroke="{}" stroke-width="2.2"{dash} points="{pts}"/>"##,
            s.color
        );
    }

    let mut lx = ml + 14.0;
    let ly = mt + 18.0;
    for s in series {
        let dash = if s.dashed {
            r##" stroke-dasharray="6,4""##
        } else {
            ""
        };
        let _ = writeln!(
            svg,
            r##"<line x1="{lx:.1}" y1="{ly:.1}" x2="{:.1}" y2="{ly:.1}" stroke="{}" stroke-width="2.5"{dash}/>"##,
            lx + 22.0,
            s.color
        );
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-family="sans-serif" font-size="12" fill="#111827">{}</text>"##,
            lx + 28.0,
            ly + 4.0,
            s.label
        );
        lx += 38.0 + (s.label.len() as f64) * 6.8;
    }

    let _ = writeln!(svg, "</svg>");
    fs::write(path, svg)
}

pub fn write_bar_chart(
    path: &Path,
    title: &str,
    y_label: &str,
    labels: &[&str],
    values: &[f64],
    colors: &[&str],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let w = 820.0_f64;
    let h = 460.0_f64;
    let ml = 82.0_f64;
    let mr = 30.0_f64;
    let mt = 56.0_f64;
    let mb = 85.0_f64;
    let pw = w - ml - mr;
    let ph = h - mt - mb;

    let ymin = values.iter().copied().fold(0.0_f64, f64::min).min(0.0);
    let ymax = values.iter().copied().fold(0.0_f64, f64::max).max(1e-6) * 1.15;
    let span = (ymax - ymin).max(1e-6);

    let mut svg = String::with_capacity(4096);
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"##
    );
    let _ = writeln!(svg, r##"<rect width="{w}" height="{h}" fill="#ffffff"/>"##);
    let _ = writeln!(
        svg,
        r##"<text x="{}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold" fill="#111827">{title}</text>"##,
        w * 0.5
    );
    for k in 0..=5 {
        let v = ymin + span * k as f64 / 5.0;
        let yp = mt + ph - ((v - ymin) / span) * ph;
        let _ = writeln!(
            svg,
            r##"<line x1="{ml}" y1="{yp:.1}" x2="{:.1}" y2="{yp:.1}" stroke="#e5e7eb" stroke-width="1"/>"##,
            ml + pw
        );
        let _ = writeln!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#4b5563">{v:.4}</text>"##,
            ml - 8.0,
            yp + 4.0
        );
    }
    let _ = writeln!(
        svg,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#374151" stroke-width="1.2"/>"##
    );
    let _ = writeln!(
        svg,
        r##"<text x="20" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#1f2937" transform="rotate(-90 20 {:.1})">{y_label}</text>"##,
        mt + ph * 0.5,
        mt + ph * 0.5
    );

    let y_zero = mt + ph - ((0.0 - ymin) / span) * ph;
    let n = labels.len().max(1);
    let slot = pw / n as f64;
    let bw = (slot * 0.58).min(80.0);
    for (i, ((&lbl, &val), &col)) in labels
        .iter()
        .zip(values.iter())
        .zip(colors.iter())
        .enumerate()
    {
        let xc = ml + (i as f64 + 0.5) * slot;
        let yv = mt + ph - ((val - ymin) / span) * ph;
        let y0 = yv.min(y_zero);
        let bh = (yv - y_zero).abs().max(1.0);
        let _ = writeln!(
            svg,
            r##"<rect x="{:.1}" y="{y0:.1}" width="{bw:.1}" height="{bh:.1}" fill="{col}"/>"##,
            xc - bw * 0.5
        );
        let _ = writeln!(
            svg,
            r##"<text x="{xc:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#111827">{val:.4}</text>"##,
            y0 - 6.0
        );
        let _ = writeln!(
            svg,
            r##"<text x="{xc:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#1f2937">{lbl}</text>"##,
            mt + ph + 20.0
        );
    }
    let _ = writeln!(svg, "</svg>");
    fs::write(path, svg)
}
