mod rng; mod inventory; mod allocation; mod entities; mod engine; mod clearing; mod analysis; mod svg;

use std::env;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

fn manifest_scan(dir: &Path, prefix: &str, mf: &mut BufWriter<File>) {
    let mut entries: Vec<_> = fs::read_dir(dir).unwrap().filter_map(|e| e.ok()).collect();
    entries.sort_by_key(|e| e.file_name());
    for e in entries {
        let name = e.file_name().to_string_lossy().to_string();
        if name.starts_with('.') { continue; }
        let p = e.path();
        if p.is_dir() { manifest_scan(&p, &format!("{}{}/", prefix, name), mf); }
        else if name != "files_manifest.csv" {
            let bytes = fs::read(&p).unwrap();
            let lines = bytes.iter().filter(|&&b| b == b'\n').count();
            let h = rng::fnv1a_bytes(&bytes);
            let _ = writeln!(mf, "{}{},{},{},{:016x}", prefix, name, lines, bytes.len(), h);
        }
    }
}

fn main() {
    let out_dir = "out/rust";
    let _ = fs::create_dir_all("out/rust/figures");
    let cfg = entities::load_config("config.kv");
    let args: Vec<String> = env::args().collect();
    let cmd = args.get(1).map(|s| s.as_str()).unwrap_or("all");
    if cmd == "s01" { analysis::write_s01_to_s09(&cfg, out_dir); return; }
    let mut tot = 0.0f64;
    for st in 1..=9 { println!("[rust] tahap {}/40 (s{:02})", st, st); }
    let t0 = Instant::now(); analysis::write_s01_to_s09(&cfg, out_dir); tot += t0.elapsed().as_secs_f64();
    println!("[rust] tahap 10/40 (s10)");
    let t0 = Instant::now(); let mut sim = entities::Sim::new_scaled(&cfg, true, Some(out_dir), 400000);
    engine::sim_run_days(&mut sim, cfg.num_days, 400000, 0, -1, 1); tot += t0.elapsed().as_secs_f64();
    for st in 11..=38 { println!("[rust] tahap {}/40 (s{:02})", st, st); }
    let t0 = Instant::now(); analysis::write_s11_to_s38(&sim, &cfg, out_dir); tot += t0.elapsed().as_secs_f64();
    println!("[rust] tahap 39/40 (s39)");
    let t0 = Instant::now(); svg::write_s39_svgs(&sim, out_dir); tot += t0.elapsed().as_secs_f64();
    println!("[rust] tahap 40/40 (s40)");
    let _ = writeln!(File::create(format!("{}/s40_summary.json", out_dir)).unwrap(), "{{\"lang\":\"rust\",\"days\":5,\"ev_per_day\":400000,\"ev_total\":{},\"toolchain\":\"rustc 1.93.1\",\"total_sec\":{:.3}}}", cfg.num_days * 400000, tot);
    let mut mf = BufWriter::new(File::create(format!("{}/files_manifest.csv", out_dir)).unwrap());
    let _ = writeln!(mf, "file,lines,bytes,fnv1a");
    manifest_scan(Path::new(out_dir), "", &mut mf);
}
