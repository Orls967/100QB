mod analysis;
mod config;
mod dtmc;
mod engine;
mod network;
mod passenger;
mod rng;
mod station;
mod stats;
mod svg;

use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;
use config::{default_rp, load_config};
use network::Sim;

fn manifest_scan(dir: &Path, prefix: &str, out: &mut BufWriter<File>) {
    let mut ents: Vec<_> = fs::read_dir(dir).unwrap().filter_map(|e| e.ok()).collect();
    ents.sort_by_key(|e| e.file_name());
    for e in ents {
        let nm = e.file_name().to_string_lossy().to_string();
        if nm.starts_with('.') { continue; }
        let p = e.path();
        if p.is_dir() { manifest_scan(&p, &format!("{}{}/", prefix, nm), out); }
        else if nm != "files_manifest.csv" {
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
    if args.get(1).map(|s| s.as_str()) == Some("s01") { stats::write_s01_to_s13(&cfg, out_dir); return; }

    let mut tot = 0.0f64;
    for st in 1..=13 { println!("[rust] tahap {}/48 (s{:02})", st, st); }
    let t0 = Instant::now(); stats::write_s01_to_s13(&cfg, out_dir); tot += t0.elapsed().as_secs_f64();

    let main_days = 365i64;
    println!("[rust] tahap 14/48 (s14)");
    let t0 = Instant::now();
    let mut rpm = default_rp(&cfg, main_days); rpm.audit_on = 1;
    let mut s_main = Sim::new(&cfg, &rpm); engine::sim_run(&mut s_main); tot += t0.elapsed().as_secs_f64();

    for st in 15..=41 { println!("[rust] tahap {}/48 (s{:02})", st, st); }
    let t0 = Instant::now(); analysis::write_s14_to_s45(&cfg, out_dir, &mut s_main); tot += t0.elapsed().as_secs_f64();

    println!("[rust] tahap 42/48 (s42)");
    let mut vals = Vec::new();
    for _ in 0..5 {
        let mut sb = Sim::new(&cfg, &default_rp(&cfg, 3));
        let tb = Instant::now(); engine::sim_run(&mut sb);
        vals.push(sb.sys_sum_l as f64 / tb.elapsed().as_secs_f64().max(1e-9));
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap()); let med = vals[2];
    let mut dev: Vec<f64> = vals.iter().map(|v| (v - med).abs()).collect(); dev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let _ = writeln!(File::create(format!("{}/s42_bench.csv", out_dir)).unwrap(), "slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,{:.2},{:.2}", med, dev[2]);

    println!("[rust] tahap 43/48 (s43)\n[rust] tahap 44/48 (s44)\n[rust] tahap 45/48 (s45)\n[rust] tahap 46/48 (s46)");
    let t0 = Instant::now(); svg::write_s46_svgs(&s_main, out_dir); tot += t0.elapsed().as_secs_f64();

    println!("[rust] tahap 47/48 (s47)\n[rust] tahap 48/48 (s48)");
    let _ = writeln!(File::create(format!("{}/s48_summary.json", out_dir)).unwrap(), "{{\"lang\":\"rust\",\"days\":{},\"sw_days\":3,\"toolchain\":\"rustc 1.93.1\",\"total_sec\":{:.3}}}", main_days, tot);
    let mut mf = BufWriter::new(File::create(format!("{}/files_manifest.csv", out_dir)).unwrap());
    let _ = writeln!(mf, "file,lines,bytes,fnv1a");
    manifest_scan(Path::new(out_dir), "", &mut mf);
}
