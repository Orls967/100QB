mod analysis;
mod contact;
mod disease;
mod engine;
mod interventions;
mod population;
mod rng;
mod svg;

use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use population::{load_config, RunParams, Sim};

fn manifest_scan(dir: &Path, prefix: &str, out: &mut BufWriter<File>) {
    let mut ents: Vec<_> = fs::read_dir(dir).unwrap().filter_map(|e| e.ok()).collect();
    ents.sort_by_key(|e| e.file_name());
    for e in ents {
        let nm = e.file_name().to_string_lossy().to_string();
        if nm.starts_with('.') { continue; }
        let p = e.path();
        if p.is_dir() {
            manifest_scan(&p, &format!("{}{}/", prefix, nm), out);
        } else if nm != "files_manifest.csv" {
            let b = fs::read(&p).unwrap();
            let lc = b.iter().filter(|&&c| c == b'\n').count();
            let _ = writeln!(out, "{}{},{},{},{:016x}", prefix, nm, lc, b.len(), rng::fnv1a_bytes(&b));
        }
    }
}

fn main() {
    let out_dir = "out/rust";
    let _ = fs::create_dir_all("out/rust/figures");
    let cfg = load_config("config.kv");
    let args: Vec<String> = env::args().collect();
    if args.get(1).map(|s| s.as_str()) == Some("s01") { analysis::write_s01_to_s09(&cfg, out_dir); return; }

    let mut tot = 0.0f64;
    for st in 1..=9 { println!("[rust] tahap {}/39 (s{:02})", st, st); }
    let t0 = Instant::now(); analysis::write_s01_to_s09(&cfg, out_dir); tot += t0.elapsed().as_secs_f64();

    let (n_main, beta_cal) = (100_000i64, 400i64);
    println!("[rust] tahap 10/39 (s10)");
    let t0 = Instant::now();
    let rp_no = RunParams { n: n_main, days: cfg.main_days, beta_bp: beta_cal, interv_on: 0, ..Default::default() };
    let mut s_no = Sim::new(&cfg, &rp_no, true); engine::sim_run(&mut s_no); tot += t0.elapsed().as_secs_f64();

    println!("[rust] tahap 11/39 (s11)");
    let t0 = Instant::now();
    let rp_main = RunParams { n: n_main, days: cfg.main_days, beta_bp: beta_cal, interv_on: 1, ..Default::default() };
    let mut s_main = Sim::new(&cfg, &rp_main, true); engine::sim_run(&mut s_main); tot += t0.elapsed().as_secs_f64();

    for st in 12..=34 { println!("[rust] tahap {}/39 (s{:02})", st, st); }
    let t0 = Instant::now(); analysis::write_s10_to_s37(&s_no, &s_main, &cfg, out_dir, beta_cal); tot += t0.elapsed().as_secs_f64();

    println!("[rust] tahap 35/39 (s35)");
    let mut vals = Vec::new();
    for _ in 0..5 {
        let rpb = RunParams { n: 2000, days: 60, beta_bp: beta_cal, interv_on: 1, ..Default::default() };
        let mut sb = Sim::new(&cfg, &rpb, false);
        let tb = Instant::now(); engine::sim_run(&mut sb);
        vals.push(120_000.0 / tb.elapsed().as_secs_f64().max(1e-9));
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap()); let med = vals[2];
    let mut dev: Vec<f64> = vals.iter().map(|v| (v - med).abs()).collect(); dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let _ = writeln!(File::create(format!("{}/s35_bench.csv", out_dir)).unwrap(), "person_days,pd_per_sec_med,pd_per_sec_mad\n120000,{:.2},{:.2}", med, dev[2]);

    println!("[rust] tahap 36/39 (s36)\n[rust] tahap 37/39 (s37)\n[rust] tahap 38/39 (s38)");
    let t0 = Instant::now(); svg::write_s38_svgs(&s_no, &s_main, out_dir); tot += t0.elapsed().as_secs_f64();

    println!("[rust] tahap 39/39 (s39)");
    let _ = writeln!(File::create(format!("{}/s39_summary.json", out_dir)).unwrap(), "{{\"lang\":\"rust\",\"N\":{},\"days\":{},\"N_sw\":1000,\"toolchain\":\"rustc 1.93.1\",\"total_sec\":{:.3}}}", n_main, cfg.main_days, tot);
    let mut mf = BufWriter::new(File::create(format!("{}/files_manifest.csv", out_dir)).unwrap());
    let _ = writeln!(mf, "file,lines,bytes,fnv1a");
    manifest_scan(Path::new(out_dir), "", &mut mf);
}
