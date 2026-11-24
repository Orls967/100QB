#![forbid(unsafe_code)]

use saw::{enumerate_from_prefix, enumerate_symmetry, run_pivot, PivotCfg};
use std::env;
use std::io::{self, BufRead, Write};

fn arg_u64(args: &[String], flag: &str, def: u64) -> u64 {
    args.windows(2)
        .find(|w| w[0] == flag)
        .map(|w| w[1].parse().expect("u64"))
        .unwrap_or(def)
}

fn fmt_u64s(v: &[u64]) -> String {
    let s: Vec<String> = v.iter().map(|x| x.to_string()).collect();
    format!("[{}]", s.join(","))
}

fn parse_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("\"{}\":", key);
    let rest = line.get(line.find(&pat)? + pat.len()..)?.trim_start();
    if let Some(s) = rest.strip_prefix('"') {
        Some(&s[..s.find('"')?])
    } else {
        let end = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        Some(&rest[..end])
    }
}

fn run_stdin_worker() {
    let stdin = io::stdin();
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in stdin.lock().lines().map_while(Result::ok) {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        let id: u64 = parse_field(t, "id")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let nmax: usize = parse_field(t, "nmax")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let weight: u64 = parse_field(t, "weight")
            .and_then(|s| s.parse().ok())
            .unwrap_or(4);
        let steps = parse_field(t, "steps").unwrap_or("");
        let msg = match enumerate_from_prefix(steps, weight, nmax) {
            Ok(c) => format!(
                "{{\"id\":{},\"ok\":true,\"counts\":{}}}\n",
                id,
                fmt_u64s(&c)
            ),
            Err(e) => format!("{{\"id\":{},\"ok\":false,\"error\":\"{}\"}}\n", id, e),
        };
        out.write_all(msg.as_bytes()).expect("write");
        out.flush().expect("flush");
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("count") => {
            if args.iter().any(|a| a == "--stdin") {
                run_stdin_worker();
                return;
            }
            let nmax = arg_u64(&args[2..], "--nmax", 18) as usize;
            let hm_n = args[2..]
                .windows(2)
                .find(|w| w[0] == "--heatmap-n")
                .map(|w| w[1].parse::<usize>().expect("hm"));
            let (counts, hm) = enumerate_symmetry(nmax, hm_n);
            let mut out = format!("{{\"nmax\":{},\"counts\":{}", nmax, fmt_u64s(&counts));
            if let (Some(hn), Some(grid)) = (hm_n, hm) {
                let rows: Vec<String> = grid.iter().map(|r| fmt_u64s(r)).collect();
                out.push_str(&format!(
                    ",\"heatmap_n\":{},\"heatmap\":[{}]",
                    hn,
                    rows.join(",")
                ));
            }
            out.push('}');
            println!("{}", out);
        }
        Some("pivot") => {
            let n = arg_u64(&args[2..], "--n", 100) as usize;
            let cfg = PivotCfg {
                n,
                seed: arg_u64(&args[2..], "--seed", 20261124),
                burn_in: arg_u64(&args[2..], "--burn-in", (25 * n as u64).max(2000)) as usize,
                batches: arg_u64(&args[2..], "--batches", 50) as usize,
                batch_size: arg_u64(&args[2..], "--batch-size", 3000) as usize,
                samples: arg_u64(&args[2..], "--samples", 0) as usize,
                stride: arg_u64(&args[2..], "--stride", 800) as usize,
            };
            let (mean, se, acc, walks) = run_pivot(&cfg);
            let ws: Vec<String> = walks
                .iter()
                .map(|w| {
                    let ps: Vec<String> = w.iter().map(|(x, y)| format!("[{},{}]", x, y)).collect();
                    format!("[{}]", ps.join(","))
                })
                .collect();
            println!(
                "{{\"n\":{},\"seed\":{},\"r2_mean\":{:.10},\"r2_se\":{:.10},\"acc\":{:.6},\"samples\":[{}]}}",
                cfg.n, cfg.seed, mean, se, acc, ws.join(",")
            );
        }
        _ => {
            eprintln!("usage: saw <count|pivot>");
            std::process::exit(1);
        }
    }
}
