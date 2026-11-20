//! Penulis grafik SVG buatan sendiri (line chart & bar chart) tanpa dependensi eksternal.

use std::fmt::Write;
use std::fs;
use std::path::Path;

pub const PALETTE: [&str; 6] = [
    "##1f6feb", "##d1495b", "##2a9d8f", "##e9c46a", "##6a4c93", "##f4a261",
];

/// Satu seri garis (x, y).
pub struct LineSeries {
    pub label: String,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub color: String,
    pub dashed: bool,
}

/// Gambar grafik garis 2D ke SVG dan simpan ke ``path``.
pub fn write_line_chart_svg(
    path: &Path,
    title: &str,
    x_label: &str,
    y_label: &str,
    series: &[LineSeries],
    log_x: bool,
    log_y: bool,
) -> std::io::Result<()> {
    let width = 820.0;
    let height = 460.0;
    let ml = 82.0;
    let mr = 30.0;
    let mt = 48.0;
    let mb = 64.0;
    let pw = width - ml - mr;
    let ph = height - mt - mb;

    let tx = |x: f64| if log_x { x.max(1e-30).log10() } else { x };
    let ty = |y: f64| if log_y { y.max(1e-30).log10() } else { y };

    let mut xmin = f64::INFINITY;
    let mut xmax = f64::NEG_INFINITY;
    let mut ymin = f64::INFINITY;
    let mut ymax = f64::NEG_INFINITY;
    for s in series {
        for (&x, &y) in s.xs.iter().zip(s.ys.iter()) {
            let vx = tx(x);
            let vy = ty(y);
            if vx.is_finite() && vy.is_finite() {
                xmin = xmin.min(vx);
                xmax = xmax.max(vx);
                ymin = ymin.min(vy);
                ymax = ymax.max(vy);
            }
        }
    }
    if (xmax - xmin).abs() < 1e-15 {
        xmin -= 1.0;
        xmax += 1.0;
    }
    if (ymax - ymin).abs() < 1e-15 {
        ymin -= 1.0;
        ymax += 1.0;
    }
    let ypad = 0.06 * (ymax - ymin);
    ymin -= ypad;
    ymax += ypad;

    let px = |x: f64| ml + (tx(x) - xmin) / (xmax - xmin) * pw;
    let py = |y: f64| mt + (1.0 - (ty(y) - ymin) / (ymax - ymin)) * ph;

    let mut out = String::new();
    let _ = writeln!(
        out,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" style="background:#ffffff;font-family:system-ui,sans-serif">"##
    );
    let _ = writeln!(
        out,
        r##"<text x="{}" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">{}</text>"##,
        width / 2.0,
        escape_xml(title)
    );

    // Grid & ticks
    for k in 0..=4 {
        let frac = (k as f64) / 4.0;
        let gx = ml + frac * pw;
        let gy = mt + (1.0 - frac) * ph;
        let xv = xmin + frac * (xmax - xmin);
        let yv = ymin + frac * (ymax - ymin);
        let _ = writeln!(
            out,
            r##"<line x1="{gx:.1}" y1="{mt}" x2="{gx:.1}" y2="{:.1}" stroke="#e5e7eb" stroke-dasharray="3,3"/>"##,
            mt + ph
        );
        let _ = writeln!(
            out,
            r##"<line x1="{ml}" y1="{gy:.1}" x2="{:.1}" y2="{gy:.1}" stroke="#e5e7eb" stroke-dasharray="3,3"/>"##,
            ml + pw
        );
        let _ = writeln!(
            out,
            r##"<text x="{gx:.1}" y="{:.1}" text-anchor="middle" font-size="11" fill="#444">{:.3}</text>"##,
            mt + ph + 18.0,
            if log_x { 10f64.powf(xv) } else { xv }
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="end" font-size="11" fill="#444">{:.3}</text>"##,
            ml - 8.0,
            gy + 4.0,
            if log_y { 10f64.powf(yv) } else { yv }
        );
    }

    // Axes
    let _ = writeln!(
        out,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#333" stroke-width="1"/>"##
    );
    let _ = writeln!(
        out,
        r##"<text x="{}" y="{:.1}" text-anchor="middle" font-size="12" fill="#222">{}</text>"##,
        ml + pw / 2.0,
        height - 14.0,
        escape_xml(x_label)
    );
    let _ = writeln!(
        out,
        r##"<text x="18" y="{}" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 {})">{}</text>"##,
        mt + ph / 2.0,
        mt + ph / 2.0,
        escape_xml(y_label)
    );

    // Series
    for s in series {
        let mut pts = String::new();
        for (&x, &y) in s.xs.iter().zip(s.ys.iter()) {
            let _ = write!(pts, "{:.2},{:.2} ", px(x), py(y));
        }
        let dash = if s.dashed {
            r##" stroke-dasharray="6,4""##
        } else {
            ""
        };
        let _ = writeln!(
            out,
            r##"<polyline fill="none" stroke="{}" stroke-width="2.2"{dash} points="{pts}"/>"##,
            s.color
        );
    }

    // Legend
    for (idx, s) in series.iter().enumerate() {
        let lx = ml + 14.0;
        let ly = mt + 18.0 + (idx as f64) * 18.0;
        let _ = writeln!(
            out,
            r##"<line x1="{lx}" y1="{ly}" x2="{:.1}" y2="{ly}" stroke="{}" stroke-width="2.5"/>"##,
            lx + 20.0,
            s.color
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{:.1}" font-size="11" fill="#222">{}</text>"##,
            lx + 26.0,
            ly + 4.0,
            escape_xml(&s.label)
        );
    }

    out.push_str("</svg>\n");
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, out)
}

/// Gambar grafik batang sederhana ke SVG.
pub fn write_bar_chart_svg(
    path: &Path,
    title: &str,
    y_label: &str,
    labels: &[String],
    values: &[f64],
) -> std::io::Result<()> {
    let width = 760.0;
    let height = 420.0;
    let ml = 80.0;
    let mr = 30.0;
    let mt = 48.0;
    let mb = 70.0;
    let pw = width - ml - mr;
    let ph = height - mt - mb;

    let ymax = values.iter().copied().fold(0.0f64, f64::max).max(1e-12) * 1.15;
    let n = labels.len().max(1) as f64;
    let slot = pw / n;
    let bw = slot * 0.56;

    let mut out = String::new();
    let _ = writeln!(
        out,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" style="background:#ffffff;font-family:system-ui,sans-serif">"##
    );
    let _ = writeln!(
        out,
        r##"<text x="{}" y="28" text-anchor="middle" font-size="16" font-weight="bold" fill="#111">{}</text>"##,
        width / 2.0,
        escape_xml(title)
    );
    let _ = writeln!(
        out,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#333" stroke-width="1"/>"##
    );

    for (i, (lab, &val)) in labels.iter().zip(values.iter()).enumerate() {
        let bx = ml + (i as f64 + 0.5) * slot - 0.5 * bw;
        let bh = (val.max(0.0) / ymax) * ph;
        let by = mt + ph - bh;
        let col = PALETTE[i % PALETTE.len()];
        let _ = writeln!(
            out,
            r##"<rect x="{bx:.1}" y="{by:.1}" width="{bw:.1}" height="{bh:.1}" fill="{col}"/>"##
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="11" fill="#222">{:.4}</text>"##,
            bx + 0.5 * bw,
            by - 5.0,
            val
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="11" fill="#333">{}</text>"##,
            bx + 0.5 * bw,
            mt + ph + 18.0,
            escape_xml(lab)
        );
    }

    let _ = writeln!(
        out,
        r##"<text x="18" y="{}" text-anchor="middle" font-size="12" fill="#222" transform="rotate(-90 18 {})">{}</text>"##,
        mt + ph / 2.0,
        mt + ph / 2.0,
        escape_xml(y_label)
    );
    out.push_str("</svg>\n");
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, out)
}

fn escape_xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
