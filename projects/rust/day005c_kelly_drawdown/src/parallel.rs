/// Deterministic parallel execution using fixed-order block reduction (`std::thread::scope`).
pub fn default_workers() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1)
}

/// Bit-identical parallel reduction over fixed block size (`BLOCK = 64`),
/// independent of `n_workers` (`1` vs `8`).
pub fn reduce_vec_fixed_blocks<F>(
    n_items: usize,
    n_workers: usize,
    n_acc: usize,
    eval_item: F,
) -> Vec<f64>
where
    F: Fn(usize, &mut [f64]) + Sync,
{
    const BLOCK: usize = 64;
    if n_items == 0 || n_acc == 0 {
        return vec![0.0; n_acc];
    }
    let n_blocks = n_items.div_ceil(BLOCK);
    let w = n_workers.max(1).min(n_blocks);
    let block_sums: Vec<Vec<f64>> = if w == 1 {
        (0..n_blocks)
            .map(|b| {
                let start = b * BLOCK;
                let end = ((b + 1) * BLOCK).min(n_items);
                let mut acc = vec![0.0_f64; n_acc];
                let mut tmp = vec![0.0_f64; n_acc];
                for idx in start..end {
                    tmp.fill(0.0);
                    eval_item(idx, &mut tmp);
                    for j in 0..n_acc {
                        acc[j] += tmp[j];
                    }
                }
                acc
            })
            .collect()
    } else {
        let blocks_per_w = n_blocks.div_ceil(w);
        let mut outer: Vec<Vec<Vec<f64>>> = Vec::with_capacity(w);
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
                    let mut tmp = vec![0.0_f64; n_acc];
                    for b in b_start..b_end {
                        let start = b * BLOCK;
                        let end = ((b + 1) * BLOCK).min(n_items);
                        let mut acc = vec![0.0_f64; n_acc];
                        for idx in start..end {
                            tmp.fill(0.0);
                            eval_ref(idx, &mut tmp);
                            for j in 0..n_acc {
                                acc[j] += tmp[j];
                            }
                        }
                        out.push(acc);
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
    let mut total = vec![0.0_f64; n_acc];
    for bs in block_sums {
        for j in 0..n_acc {
            total[j] += bs[j];
        }
    }
    total
}

/// Deterministic parallel map over items (preserves index order).
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
