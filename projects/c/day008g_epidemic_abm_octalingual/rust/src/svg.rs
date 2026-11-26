use std::fs::File;
use std::io::Write;
use crate::population::*;

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

pub fn write_s38_svgs(s_no: &Sim, s_main: &Sim, out_dir: &str) {
    let d_cnt = s_main.rp.days as usize; let step = (d_cnt / 12).max(1);
    let mut epi = [0.0; 12]; let mut inc = [0.0; 12]; let mut rt = [0.0; 12]; let mut hosp = [0.0; 12]; let mut ens = [0.0; 12];
    for i in 0..12 {
        let d = (i * step).min(d_cnt - 1);
        epi[i] = (s_main.daily_cnt[d][ST_IA as usize] + s_main.daily_cnt[d][ST_IS as usize]) as f64;
        inc[i] = s_main.daily_inc_age[d].iter().sum::<i64>() as f64;
        rt[i] = (s_no.daily_cnt[d][ST_IA as usize] + s_no.daily_cnt[d][ST_IS as usize]) as f64;
        hosp[i] = s_main.daily_hosp[d][0] as f64;
        ens[i] = s_main.daily_var[d][1] as f64;
    }
    write_bar_svg(out_dir, "epi_curves.svg", "Rust Active Infectious (Interv)", &epi, "#1f77b4");
    write_bar_svg(out_dir, "incidence_age.svg", "Rust Daily Incidence Sample", &inc, "#ff7f0e");
    write_bar_svg(out_dir, "rt_series.svg", "Rust Active Infectious (No Interv)", &rt, "#d62728");
    write_bar_svg(out_dir, "hospital_cap.svg", "Rust Hospital Occupancy Sample", &hosp, "#9467bd");
    write_bar_svg(out_dir, "layer_attack.svg", "Rust Infections by Layer (HH/SW/Comm)", &[s_main.layer_inf[0] as f64, s_main.layer_inf[1] as f64, s_main.layer_inf[2] as f64], "#2ca02c");
    write_bar_svg(out_dir, "offspring_dist.svg", "Rust Infections by Age Group", &[s_main.age_inf[0] as f64, s_main.age_inf[1] as f64, s_main.age_inf[2] as f64], "#8c564b");
    write_bar_svg(out_dir, "sweep_heatmap.svg", "Rust No-Interv vs Interv Total Infections", &[s_no.total_trans as f64, s_main.total_trans as f64], "#e377c2");
    write_bar_svg(out_dir, "vaccine_effect.svg", "Rust Vaccine Doses & Blocked Transmissions", &[s_main.branches[18] as f64, s_main.branches[3] as f64], "#17becf");
    write_bar_svg(out_dir, "variant_share.svg", "Rust Variant 1 Daily Incidence", &ens, "#bcbd22");
    write_bar_svg(out_dir, "ensemble_band.svg", "Rust Final Recovered vs Deaths", &[s_main.daily_cnt[d_cnt - 1][ST_R as usize] as f64, s_main.daily_cnt[d_cnt - 1][ST_D as usize] as f64], "#7f7f7f");
}
