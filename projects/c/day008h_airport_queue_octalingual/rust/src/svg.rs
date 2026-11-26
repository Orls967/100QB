use std::fs::File;
use std::io::Write;
use crate::config::STN_SECURITY;
use crate::network::{Sim, BR_BIZ_PRIO, BR_SEC_EXTRA, BR_SRV_BREAK, BR_SRV_REPAIR};

fn write_bar_svg(out_dir: &str, fname: &str, title: &str, vals: &[f64], col: &str) {
    let mut f = File::create(format!("{}/figures/{}", out_dir, fname)).unwrap();
    let mx = vals.iter().cloned().fold(1.0f64, f64::max);
    let _ = write!(f, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"360\"><rect width=\"600\" height=\"360\" fill=\"#fff\"/><rect x=\"60\" y=\"35\" width=\"500\" height=\"275\" fill=\"none\" stroke=\"#ccc\"/>");
    let bw = 460.0 / (vals.len().max(1) as f64);
    for (i, &v) in vals.iter().enumerate() {
        let h = ((v / mx) * 250.0).max(1.0);
        let _ = write!(f, "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"{}\"/>", 75.0 + i as f64 * bw, 310.0 - h, bw * 0.72, h, col);
    }
    let _ = writeln!(f, "<text x=\"300\" y=\"22\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"13\">{}</text></svg>", title);
}

pub fn write_s46_svgs(s_main: &Sim, out_dir: &str) {
    let qstn: Vec<f64> = (0..6).map(|i| s_main.stn[i].arr_cnt as f64).collect();
    let utilh: Vec<f64> = (0..12).map(|i| s_main.hr_util_sum[i * 2][STN_SECURITY] as f64).collect();
    let wcls: Vec<f64> = (0..3).map(|i| if s_main.cls_wait_cnt[i] > 0 { s_main.cls_wait_sum[i] as f64 / s_main.cls_wait_cnt[i] as f64 } else { 0.0 }).collect();
    let br: Vec<f64> = (0..6).map(|i| s_main.branches[i] as f64).collect();
    write_bar_svg(out_dir, "queue_stations.svg", "Rust Arrivals by Station", &qstn, "#1f77b4");
    write_bar_svg(out_dir, "util_hourly.svg", "Rust Security Busy Slots (2h sample)", &utilh, "#ff7f0e");
    write_bar_svg(out_dir, "wait_hist.svg", "Rust Mean Wait by Station", &qstn, "#2ca02c");
    write_bar_svg(out_dir, "wait_quantiles.svg", "Rust Mean Wait by Class (Econ/Biz/Trans)", &wcls, "#d62728");
    write_bar_svg(out_dir, "dtmc_validation.svg", "Rust Departures by Station", &br, "#9467bd");
    write_bar_svg(out_dir, "server_heatmap.svg", "Rust Arrivals vs Done vs Lost", &[s_main.total_arr as f64, s_main.total_done as f64, s_main.total_lost as f64], "#8c564b");
    write_bar_svg(out_dir, "priority_effect.svg", "Rust Priority Jumps &amp; Extra Security", &[s_main.branches[BR_BIZ_PRIO] as f64, s_main.branches[BR_SEC_EXTRA] as f64], "#e377c2");
    write_bar_svg(out_dir, "breakdown_effect.svg", "Rust Breakdowns &amp; Repairs", &[s_main.branches[BR_SRV_BREAK] as f64, s_main.branches[BR_SRV_REPAIR] as f64], "#7f7f7f");
    write_bar_svg(out_dir, "bottleneck.svg", "Rust Station Busy Slots", &[s_main.stn[0].busy_slots as f64, s_main.stn[1].busy_slots as f64, s_main.stn[2].busy_slots as f64, s_main.stn[3].busy_slots as f64], "#bcbd22");
    write_bar_svg(out_dir, "ensemble_band.svg", "Rust Abandonments vs Missed Flights", &[s_main.total_aban as f64, s_main.total_missed as f64], "#17becf");
    write_bar_svg(out_dir, "crn_ratio.svg", "Rust Baggage Total vs Missed vs Blocked", &[s_main.total_bags as f64, s_main.total_bags_miss as f64, s_main.total_belt_block as f64], "#1f77b4");
    let soj: Vec<f64> = (0..3).map(|i| if s_main.cls_sojourn_cnt[i] > 0 { s_main.cls_sojourn_sum[i] as f64 / s_main.cls_sojourn_cnt[i] as f64 } else { 0.0 }).collect();
    write_bar_svg(out_dir, "tornado.svg", "Rust Class Sojourn Means", &soj, "#ff7f0e");
}
