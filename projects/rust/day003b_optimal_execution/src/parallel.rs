//! Paralelisme ``std::thread`` deterministik: chunk tetap, seed per chunk, hasil invarian terhadap jumlah thread.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use crate::rng::splitmix64;

/// Jumlah thread default dari ``available_parallelism`` (minimal 1).
pub fn default_num_threads() -> usize {
    thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

/// Jalankan ``total`` unit kerja dalam ``n_chunks`` chunk tetap secara paralel dengan ``n_threads`` thread.
/// Karena jumlah chunk dan seed per chunk hanya bergantung pada indeks chunk (bukan ``n_threads``),
/// hasil reduksi berurutan identik bit-demi-bit untuk 1 maupun banyak thread.
pub fn parallel_fixed_chunks<T, F>(
    total: usize,
    n_chunks: usize,
    n_threads: usize,
    base_seed: u64,
    worker: &F,
) -> Vec<T>
where
    T: Send,
    F: Fn(usize, usize, u64) -> T + Sync,
{
    let chunks = n_chunks.max(1).min(total.max(1));
    let threads = n_threads.max(1).min(chunks);
    let chunk_size = total.div_ceil(chunks);

    let next_idx = AtomicUsize::new(0);
    let mut slots: Vec<Option<T>> = (0..chunks).map(|_| None).collect();

    thread::scope(|s| {
        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let next_ref = &next_idx;
            handles.push(s.spawn(move || {
                let mut local = Vec::new();
                loop {
                    let idx = next_ref.fetch_add(1, Ordering::Relaxed);
                    if idx >= chunks {
                        break;
                    }
                    let start = idx * chunk_size;
                    let end = ((idx + 1) * chunk_size).min(total);
                    let count = end.saturating_sub(start);
                    let seed = splitmix64(base_seed ^ ((idx as u64).wrapping_add(1)));
                    let val = worker(idx, count, seed);
                    local.push((idx, val));
                }
                local
            }));
        }
        for h in handles {
            for (idx, val) in h.join().expect("worker thread panicked") {
                slots[idx] = Some(val);
            }
        }
    });

    slots
        .into_iter()
        .map(|o| o.expect("missing chunk"))
        .collect()
}
