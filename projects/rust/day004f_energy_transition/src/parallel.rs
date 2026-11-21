//! Deterministic parallel chunk execution using `std::thread::scope`.

pub fn default_workers() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get().min(8))
        .unwrap_or(4)
}

/// Runs `f(idx)` for `0..n` across `n_workers` threads with fixed contiguous chunks,
/// returning results in exact deterministic index order `0..n`.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_run_chunks_deterministic_order() {
        let r1 = run_chunks(103, 1, |i| i * i + 7);
        let r4 = run_chunks(103, 4, |i| i * i + 7);
        let r8 = run_chunks(103, 8, |i| i * i + 7);
        assert_eq!(r1, r4);
        assert_eq!(r1, r8);
    }
}
