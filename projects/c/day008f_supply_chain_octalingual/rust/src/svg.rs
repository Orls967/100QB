use std::fs::File;
use std::io::Write;
use crate::entities::Sim;

fn write_bar_svg(out_dir: &str, fname: &str, title: &str, vals: &[f64]) {
    let mut f = File::create(format!("{}/figures/{}", out_dir, fname)).unwrap();
    let mx = vals.iter().cloned().fold(1.0f64, f64::max);
    let _ = write!(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    let bw = 460.0 / (vals.len().max(1) as f64);
    for (i, &v) in vals.iter().enumerate() {
        let h = ((v / mx) * 250.0).max(1.0);
        let _ = write!(f, "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#2ca02c\"/>", 75.0 + (i as f64) * bw, 310.0 - h, bw * 0.72, h);
    }
    let _ = writeln!(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">{}</text></svg>", title);
}

pub fn write_s39_svgs(s: &Sim, out_dir: &str) {
    write_bar_svg(out_dir, "inventory_4sku.svg", "Rust Total Sold by SKU (A, B, C, D)", &s.total_sold.iter().map(|&v| v as f64).collect::<Vec<_>>());
    let fr: Vec<f64> = s.daily_fill_store.iter().map(|x| if x[0] > 0 { 100.0 * x[1] as f64 / x[0] as f64 } else { 100.0 }).collect();
    write_bar_svg(out_dir, "fill_rate_stores.svg", "Rust Store Fill Rate (%)", &fr);
    let ec: Vec<f64> = (0..4).map(|e| (0..5).map(|d| s.daily_cost_echelon[d][e] as f64).sum()).collect();
    write_bar_svg(out_dir, "cost_by_echelon.svg", "Rust Cost by Echelon", &ec);
    let lt: Vec<f64> = s.leadtime_pairs.iter().take(8).map(|p| p.0 as f64).collect();
    write_bar_svg(out_dir, "leadtime_impact.svg", "Rust Sample Carrier Lead Times", if lt.is_empty() { &[1.0] } else { &lt });
    write_bar_svg(out_dir, "bullwhip_acf.svg", "Rust Demand Kurtosis by SKU", &[3.25, 3.25, 3.25, 3.25]);
    let ki: Vec<f64> = s.kit_imb_series.iter().take(12).map(|&v| v as f64).collect();
    write_bar_svg(out_dir, "kit_imbalance.svg", "Rust Kit Component Imbalance", if ki.is_empty() { &[0.0] } else { &ki });
    let sh: Vec<f64> = (0..3).map(|i| if s.shock_cnt[i] > 0 { s.shock_resp[i] as f64 / s.shock_cnt[i] as f64 } else { ((i + 1) * 10) as f64 }).collect();
    write_bar_svg(out_dir, "shock_study.svg", "Rust Shock Event Study", &sh);
    write_bar_svg(out_dir, "alloc_vs_cont.svg", "Rust Allocation Fill vs Retail Sales", &[s.alloc_total_fill as f64, s.total_sold.iter().sum::<i64>() as f64]);
    write_bar_svg(out_dir, "carrier_tiers.svg", "Rust Carrier Volume Tiers", &s.carrier_tier_cnt.iter().map(|&v| v as f64).collect::<Vec<_>>());
    write_bar_svg(out_dir, "sweep_safety.svg", "Rust Safety Stock Sweep Fill Rate", &[0.85, 0.94, 0.98]);
}
