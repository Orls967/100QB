//! Deterministic parallel execution using `std::thread::scope` with fixed-size block reduction.

pub fn default_workers() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get().min(8))
        .unwrap_or(4)
}

pub fn run_chunks<T, F>(n: usize, n_workers: usize, f: F) -> Vec<T>
where
    T: Send,
    F: Fn(usize) -> T + Sync,
{
    if n == 0 {
        return Vec::new();
    }
    let w = n_workers.max(1).min(n);
    if w == 1 {
        return (0..n).map(f).collect();
    }
    let chunk_size = n.div_ceil(w);
    let f_ref = &f;
    let mut parts: Vec<Vec<T>> = Vec::with_capacity(w);

    std::thread::scope(|s| {
        let mut handles = Vec::with_capacity(w);
        for tid in 0..w {
            let start = tid * chunk_size;
            if start >= n {
                break;
            }
            let end = (start + chunk_size).min(n);
            handles.push(s.spawn(move || {
                let mut local = Vec::with_capacity(end - start);
                for idx in start..end {
                    local.push(f_ref(idx));
                }
                local
            }));
        }
        for h in handles {
            parts.push(h.join().expect("worker thread panicked"));
        }
    });

    let mut out = Vec::with_capacity(n);
    for p in parts {
        out.extend(p);
    }
    out
}

pub fn run_chunk_reduce<A, F, R>(n: usize, n_workers: usize, chunk_fn: F, merge_fn: R) -> A
where
    A: Send,
    F: Fn(usize, usize) -> A + Sync,
    R: Fn(Vec<A>) -> A,
{
    const BLOCK: usize = 64;
    let n_blocks = n.max(1).div_ceil(BLOCK);
    let w = n_workers.max(1).min(n_blocks);
    let f_ref = &chunk_fn;
    if w == 1 {
        let mut parts = Vec::with_capacity(n_blocks);
        for b in 0..n_blocks {
            let start = b * BLOCK;
            let end = (start + BLOCK).min(n);
            parts.push(f_ref(start, end));
        }
        return merge_fn(parts);
    }
    let blocks_per_worker = n_blocks.div_ceil(w);
    let mut all_blocks: Vec<A> = Vec::with_capacity(n_blocks);

    std::thread::scope(|s| {
        let mut handles = Vec::with_capacity(w);
        for tid in 0..w {
            let b_start = tid * blocks_per_worker;
            if b_start >= n_blocks {
                break;
            }
            let b_end = (b_start + blocks_per_worker).min(n_blocks);
            handles.push(s.spawn(move || {
                let mut local = Vec::with_capacity(b_end - b_start);
                for b in b_start..b_end {
                    let start = b * BLOCK;
                    let end = (start + BLOCK).min(n);
                    local.push(f_ref(start, end));
                }
                local
            }));
        }
        for h in handles {
            all_blocks.extend(h.join().expect("worker thread panicked"));
        }
    });

    merge_fn(all_blocks)
}
