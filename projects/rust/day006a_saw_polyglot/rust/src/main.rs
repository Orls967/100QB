#![forbid(unsafe_code)]

use saw::{
    enumerate_from_prefix, enumerate_symmetry, format_u64_slice, parse_prefix_task, run_pivot,
    PivotConfig,
};
use std::env;
use std::io::{self, BufRead, Write};

fn parse_arg_u64(args: &[String], flag: &str, default: u64) -> u64 {
    for w in args.windows(2) {
        if w[0] == flag {
            return w[1].parse::<u64>().expect("valid u64 flag");
        }
    }
    default
}

fn parse_arg_opt_usize(args: &[String], flag: &str) -> Option<usize> {
    for w in args.windows(2) {
        if w[0] == flag {
            return Some(w[1].parse::<usize>().expect("valid usize flag"));
        }
    }
    None
}

fn run_count_cli(args: &[String]) {
    if args.iter().any(|a| a == "--stdin") {
        run_worker_stdio();
        return;
    }
    let nmax = parse_arg_u64(args, "--nmax", 20) as usize;
    let heatmap_n = parse_arg_opt_usize(args, "--heatmap-n");
    let res = enumerate_symmetry(nmax, heatmap_n);

    let mut json = format!(
        "{{\"nmax\":{},\"counts\":{}",
        res.nmax,
        format_u64_slice(&res.counts)
    );
    if let (Some(hn), Some(grid)) = (res.heatmap_n, res.heatmap_grid) {
        json.push_str(&format!(",\"heatmap_n\":{},\"heatmap\":[", hn));
        for (ri, row) in grid.iter().enumerate() {
            if ri > 0 {
                json.push(',');
            }
            json.push_str(&format_u64_slice(row));
        }
        json.push(']');
    }
    json.push('}');
    println!("{}", json);
}

fn run_worker_stdio() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = io::BufWriter::new(stdout.lock());

    for line_res in stdin.lock().lines() {
        let line = match line_res {
            Ok(l) => l,
            Err(_) => break,
        };
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        match parse_prefix_task(trimmed)
            .and_then(|t| enumerate_from_prefix(&t.steps, t.weight, t.nmax).map(|c| (t.id, c)))
        {
            Ok((id, counts)) => {
                let resp = format!(
                    "{{\"id\":{},\"ok\":true,\"counts\":{}}}\n",
                    id,
                    format_u64_slice(&counts)
                );
                out.write_all(resp.as_bytes()).expect("write worker stdout");
                out.flush().expect("flush worker stdout");
            }
            Err(err) => {
                let safe = err.replace('"', "'");
                let resp = format!("{{\"id\":0,\"ok\":false,\"error\":\"{}\"}}\n", safe);
                out.write_all(resp.as_bytes()).expect("write worker error");
                out.flush().expect("flush worker error");
            }
        }
    }
}

fn run_pivot_cli(args: &[String]) {
    let n = parse_arg_u64(args, "--n", 100) as usize;
    let seed = parse_arg_u64(args, "--seed", 20261124);
    let burn_in = parse_arg_u64(args, "--burn-in", (25 * n as u64).max(2000)) as usize;
    let n_batches = parse_arg_u64(args, "--batches", 50) as usize;
    let batch_size = parse_arg_u64(args, "--batch-size", 4000) as usize;
    let n_samples = parse_arg_u64(args, "--samples", 0) as usize;
    let sample_stride = parse_arg_u64(args, "--sample-stride", 1000) as usize;

    let cfg = PivotConfig {
        n,
        seed,
        burn_in,
        n_batches,
        batch_size,
        n_samples,
        sample_stride,
    };
    let res = run_pivot(&cfg);

    let mut json = format!(
        "{{\"n\":{},\"seed\":{},\"r2_mean\":{:.10},\"r2_se\":{:.10},\"acceptance_rate\":{:.8},\"total_steps\":{},\"samples\":[",
        res.n, res.seed, res.r2_mean, res.r2_se, res.acceptance_rate, res.total_steps
    );
    for (si, walk) in res.samples.iter().enumerate() {
        if si > 0 {
            json.push(',');
        }
        json.push('[');
        for (pi, &(x, y)) in walk.iter().enumerate() {
            if pi > 0 {
                json.push(',');
            }
            json.push_str(&format!("[{},{}]", x, y));
        }
        json.push(']');
    }
    json.push_str("]}");
    println!("{}", json);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: saw <count|worker|pivot> [options]");
        std::process::exit(1);
    }
    match args[1].as_str() {
        "count" => run_count_cli(&args[2..]),
        "worker" => run_worker_stdio(),
        "pivot" => run_pivot_cli(&args[2..]),
        other => {
            eprintln!(
                "Unknown mode '{}'. Expected count, worker, or pivot.",
                other
            );
            std::process::exit(1);
        }
    }
}
