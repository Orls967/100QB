#![forbid(unsafe_code)]

use percolate::{run_newman_ziff, sample_bonds};
use std::io::{self, BufRead, BufWriter, Write};

fn extract_u64(line: &str, key: &str) -> Option<u64> {
    let idx = line.find(&format!("\"{key}\""))? + key.len() + 2;
    let rest = line[idx..].trim_start().strip_prefix(':')?.trim_start();
    rest.chars()
        .take_while(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .ok()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.get(1).map(|s| s.as_str()) == Some("sample") {
        let l: usize = a.get(2).and_then(|s| s.parse().ok()).unwrap_or(96);
        let p: f64 = a.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.5);
        let seed: u64 = a.get(4).and_then(|s| s.parse().ok()).unwrap_or(20261124);
        let pairs: Vec<String> = sample_bonds(l, p, seed)
            .into_iter()
            .map(|(u, v)| format!("[{u},{v}]"))
            .collect();
        println!(
            "{{\"L\":{l},\"rows\":{l},\"cols\":{},\"p\":{p:.6},\"open_bonds\":[{}]}}",
            l + 1,
            pairs.join(",")
        );
        return;
    }
    let mut w = BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines().map_while(Result::ok) {
        if line.trim().is_empty() {
            continue;
        }
        let l = extract_u64(&line, "L").unwrap() as usize;
        let samples = extract_u64(&line, "samples").unwrap() as usize;
        let r = run_newman_ziff(l, extract_u64(&line, "seed").unwrap(), samples);
        let h: Vec<String> = r.hist.iter().map(|v| v.to_string()).collect();
        writeln!(
            w,
            "{{\"L\":{l},\"samples\":{samples},\"sum_largest\":{},\"sum_largest_sq\":{},\"hist\":[{}]}}",
            r.sum_largest,
            r.sum_largest_sq,
            h.join(",")
        )
        .unwrap();
        w.flush().unwrap();
    }
}
