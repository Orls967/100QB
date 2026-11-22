/// Deterministic parallel execution using fixed-order block reduction (`std::thread::scope`).
pub fn default_workers() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1)
}

pub fn run_chunks<T, F>(n_items: usize, n_workers: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync,
{
    if n_items == 0 {
        return Vec::new();
    }
    let w = n_workers.max(1).min(n_items);
    if w == 1 {
        return (0..n_items).map(&f).collect();
    }
    let chunk = n_items.div_ceil(w);
    let mut blocks: Vec<Vec<T>> = Vec::with_capacity(w);
    std::thread::scope(|s| {
        let mut handles = Vec::with_capacity(w);
        for tid in 0..w {
            let start = tid * chunk;
            if start >= n_items {
                break;
            }
            let end = ((tid + 1) * chunk).min(n_items);
            let f_ref = &f;
            handles.push(s.spawn(move || {
                let mut out = Vec::with_capacity(end - start);
                for idx in start..end {
                    out.push(f_ref(idx));
                }
                out
            }));
        }
        for h in handles {
            blocks.push(h.join().expect("worker thread panicked"));
        }
    });
    let mut out = Vec::with_capacity(n_items);
    for mut blk in blocks {
        out.append(&mut blk);
    }
    out
}

/// Bit-identical parallel reduction over fixed block size (`BLOCK = 64`),
/// independent of `n_workers` (`1` vs `8`).
pub fn reduce_vec_fixed_blocks<F>(
    n_items: usize,
    vec_len: usize,
    n_workers: usize,
    eval_item: F,
) -> (f64, Vec<f64>)
where
    F: Fn(usize, &mut Vec<f64>) -> f64 + Sync,
{
    const BLOCK: usize = 64;
    if n_items == 0 {
        return (0.0, vec![0.0; vec_len]);
    }
    let n_blocks = n_items.div_ceil(BLOCK);
    let w = n_workers.max(1).min(n_blocks);
    let block_results: Vec<(f64, Vec<f64>)> = if w == 1 {
        (0..n_blocks)
            .map(|b| {
                let start = b * BLOCK;
                let end = ((b + 1) * BLOCK).min(n_items);
                let mut acc_v = vec![0.0_f64; vec_len];
                let mut local_v = vec![0.0_f64; vec_len];
                let mut acc_s = 0.0_f64;
                for idx in start..end {
                    local_v.fill(0.0);
                    acc_s += eval_item(idx, &mut local_v);
                    for j in 0..vec_len {
                        acc_v[j] += local_v[j];
                    }
                }
                (acc_s, acc_v)
            })
            .collect()
    } else {
        let blocks_per_w = n_blocks.div_ceil(w);
        let mut outer: Vec<Vec<(f64, Vec<f64>)>> = Vec::with_capacity(w);
        std::thread::scope(|s| {
            let mut handles = Vec::with_capacity(w);
            for tid in 0..w {
                let b_start = tid * blocks_per_w;
                if b_start >= n_blocks {
                    break;
                }
                let b_end = ((tid + 1) * blocks_per_w).min(n_blocks);
                let eval_ref = &eval_item;
                handles.push(s.spawn(move || {
                    let mut out = Vec::with_capacity(b_end - b_start);
                    let mut local_v = vec![0.0_f64; vec_len];
                    for b in b_start..b_end {
                        let start = b * BLOCK;
                        let end = ((b + 1) * BLOCK).min(n_items);
                        let mut acc_v = vec![0.0_f64; vec_len];
                        let mut acc_s = 0.0_f64;
                        for idx in start..end {
                            local_v.fill(0.0);
                            acc_s += eval_ref(idx, &mut local_v);
                            for j in 0..vec_len {
                                acc_v[j] += local_v[j];
                            }
                        }
                        out.push((acc_s, acc_v));
                    }
                    out
                }));
            }
            for h in handles {
                outer.push(h.join().expect("worker thread panicked"));
            }
        });
        let mut flat = Vec::with_capacity(n_blocks);
        for mut part in outer {
            flat.append(&mut part);
        }
        flat
    };
    let mut total_s = 0.0_f64;
    let mut total_v = vec![0.0_f64; vec_len];
    for (bs, bv) in block_results {
        total_s += bs;
        for j in 0..vec_len {
            total_v[j] += bv[j];
        }
    }
    (total_s, total_v)
}
