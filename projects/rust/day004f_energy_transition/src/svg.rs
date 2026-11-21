//! Pure Rust SVG chart generator for the 7 day004f figures.

use std::fmt::Write;
use std::fs;
use std::path::Path;

pub const PALETTE: [&str; 6] = [
    "##1f77b4", "##d62728", "##2ca02c", "##ff7f0e", "##9467bd", "##8c564b",
];

pub struct LineSeries<'a> {
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
    series: &[LineSeries<'_>],
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let w = 900.0;
    let h = 520.0;
    let ml = 80.0;
    let mr = 30.0;
    let mt = 55.0;
    let mb = 75.0;
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
        xmax = xmin + 1.0;
    }
    if (ymax - ymin).abs() < 1e-12 {
        ymax = ymin + 1.0;
    }
    let ypad = 0.08 * (ymax - ymin);
    ymin -= ypad;
    ymax += ypad;

    let sx = |xv: f64| ml + (xv - xmin) / (xmax - xmin) * pw;
    let sy = |yv: f64| mt + ph - (yv - ymin) / (ymax - ymin) * ph;

    let mut out = String::with_capacity(8192);
    let _ = writeln!(
        out,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"##
    );
    let _ = writeln!(
        out,
        r##"<rect width="100%" height="100%" fill="#ffffff"/>"##
    );
    let _ = writeln!(
        out,
        r##"<text x="{}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold" fill="#222">{}</text>"##,
        w / 2.0,
        title
    );

    for k in 0..=4 {
        let frac = k as f64 / 4.0;
        let yv = ymin + frac * (ymax - ymin);
        let yp = sy(yv);
        let _ = writeln!(
            out,
            r##"<line x1="{ml}" y1="{yp:.1}" x2="{}" y2="{yp:.1}" stroke="#e6e6e6" stroke-width="1"/>"##,
            ml + pw
        );
        let _ = writeln!(
            out,
            r##"<text x="{}" y="{:.1}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#555">{yv:.2}</text>"##,
            ml - 8.0,
            yp + 4.0
        );
        let xv = xmin + frac * (xmax - xmin);
        let xp = sx(xv);
        let _ = writeln!(
            out,
            r##"<text x="{xp:.1}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#555">{xv:.2}</text>"##,
            mt + ph + 18.0
        );
    }

    let _ = writeln!(
        out,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#444" stroke-width="1.2"/>"##
    );

    if ymin < 0.0 && ymax > 0.0 {
        let y0 = sy(0.0);
        let _ = writeln!(
            out,
            r##"<line x1="{ml}" y1="{y0:.1}" x2="{}" y2="{y0:.1}" stroke="#888" stroke-dasharray="4,4" stroke-width="1"/>"##,
            ml + pw
        );
    }

    for s in series {
        if s.x.is_empty() {
            continue;
        }
        let mut pts = String::new();
        for (&xv, &yv) in s.x.iter().zip(s.y.iter()) {
            if xv.is_finite() && yv.is_finite() {
                let _ = write!(pts, "{:.1},{:.1} ", sx(xv), sy(yv));
            }
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

    for (idx, s) in series.iter().enumerate() {
        let lx = ml + 14.0 + (idx % 3) as f64 * 260.0;
        let ly = h - 22.0 + (idx / 3) as f64 * 16.0;
        let _ = writeln!(
            out,
            r##"<line x1="{lx}" y1="{}" x2="{}" y2="{}" stroke="{}" stroke-width="2.5"/>"##,
            ly - 4.0,
            lx + 20.0,
            ly - 4.0,
            s.color
        );
        let _ = writeln!(
            out,
            r##"<text x="{}" y="{ly}" font-family="sans-serif" font-size="11" fill="#222">{}</text>"##,
            lx + 26.0,
            s.label
        );
    }

    let _ = writeln!(
        out,
        r##"<text x="{}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#333">{x_label}</text>"##,
        ml + pw / 2.0,
        mt + ph + 38.0
    );
    let _ = writeln!(
        out,
        r##"<text x="18" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#333" transform="rotate(-90 18 {})">{y_label}</text>"##,
        mt + ph / 2.0,
        mt + ph / 2.0
    );
    out.push_str("</svg>\n");
    fs::write(path, out)
}

pub fn write_bar_chart(
    path: &Path,
    title: &str,
    y_label: &str,
    labels: &[&str],
    values: &[f64],
    color: &str,
) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let w = 900.0;
    let h = 480.0;
    let ml = 80.0;
    let mr = 30.0;
    let mt = 55.0;
    let mb = 80.0;
    let pw = w - ml - mr;
    let ph = h - mt - mb;

    let mut ymin = values.iter().copied().fold(0.0f64, f64::min);
    let mut ymax = values.iter().copied().fold(0.0f64, f64::max);
    if (ymax - ymin).abs() < 1e-12 {
        ymax = ymin + 1.0;
    }
    let pad = 0.1 * (ymax - ymin);
    if ymin < 0.0 {
        ymin -= pad;
    }
    ymax += pad;
    let sy = |yv: f64| mt + ph - (yv - ymin) / (ymax - ymin) * ph;
    let y0 = sy(0.0);

    let mut out = String::with_capacity(4096);
    let _ = writeln!(
        out,
        r##"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">"##
    );
    let _ = writeln!(
        out,
        r##"<rect width="100%" height="100%" fill="#ffffff"/>"##
    );
    let _ = writeln!(
        out,
        r##"<text x="{}" y="30" text-anchor="middle" font-family="sans-serif" font-size="16" font-weight="bold" fill="#222">{title}</text>"##,
        w / 2.0
    );
    let _ = writeln!(
        out,
        r##"<rect x="{ml}" y="{mt}" width="{pw}" height="{ph}" fill="none" stroke="#444" stroke-width="1.2"/>"##
    );
    let _ = writeln!(
        out,
        r##"<line x1="{ml}" y1="{y0:.1}" x2="{}" y2="{y0:.1}" stroke="#666" stroke-width="1"/>"##,
        ml + pw
    );

    let n = labels.len().max(1);
    let slot = pw / n as f64;
    let bw = slot * 0.6;
    for (i, (&lbl, &val)) in labels.iter().zip(values.iter()).enumerate() {
        let bx = ml + i as f64 * slot + (slot - bw) * 0.5;
        let yv = sy(val);
        let (ry, rh) = if val >= 0.0 {
            (yv, y0 - yv)
        } else {
            (y0, yv - y0)
        };
        let _ = writeln!(
            out,
            r##"<rect x="{bx:.1}" y="{ry:.1}" width="{bw:.1}" height="{:.1}" fill="{color}" opacity="0.85"/>"##,
            rh.max(1.0)
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#333">{lbl}</text>"##,
            bx + bw / 2.0,
            mt + ph + 18.0
        );
        let _ = writeln!(
            out,
            r##"<text x="{:.1}" y="{:.1}" text-anchor="middle" font-family="sans-serif" font-size="10" fill="#222">{val:.2}</text>"##,
            bx + bw / 2.0,
            ry - 4.0
        );
    }

    let _ = writeln!(
        out,
        r##"<text x="18" y="{}" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#333" transform="rotate(-90 18 {})">{y_label}</text>"##,
        mt + ph / 2.0,
        mt + ph / 2.0
    );
    out.push_str("</svg>\n");
    fs::write(path, out)
}
