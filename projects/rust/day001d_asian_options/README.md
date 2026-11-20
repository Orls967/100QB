# day001d — Opsi Asia (Rust): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c** dalam bahasa Rust (*crate* Cargo `edition = "2021"`), dengan spesifikasi perhitungan yang identik:

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$.

## Cara Build, Test & Run

```bash
cd projects/rust/day001d_asian_options
cargo test
cargo clippy -- -D warnings
cargo fmt --check
cargo run --release
```

Output dijalankan ke:
- `outputs/day001d/results.json`
- `outputs/figures/day001d_asian-convergence-loglog.png`
- `outputs/figures/day001d_asian-estimators-ci.png`

## Struktur Kode

- `src/lib.rs` — rumus tertutup Black-Scholes & Kemna–Vorst, simulasi GBM paralel deterministik (`rayon` + `ChaCha8Rng` per-chunk), tiga estimator harga, Delta *pathwise* & *bump*, *fitting* slope log-log, serta 5 unit test.
- `src/main.rs` — *runner* seluruh analisis, pencatat waktu eksekusi, penulis `results.json`, dan pembuat grafik PNG (`plotters`).
- `meta.json` — metadata untuk tabel progres di `README.md` utama.

## Ringkasan Hasil (`200.000` path, `seed = 42`)

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.796308` | `0.017942` | `[5.761141, 5.831475]` | `1.00x` |
| Antithetic (100k pasang) | `5.771260` | `0.012407` | `[5.746943, 5.795578]` | `2.09x` |
| Control Variate (Geom. Asia) | `5.782043` | `0.000496` | `[5.781072, 5.783014]` | `1311.15x` |

- **Delta Pathwise vs Bump (`h=1%`)**: `0.591213 ± 0.001179` vs `0.591020 ± 0.001157` (selisih `+0.000193`, `0.12 SE`).
- **Slope konvergensi log-log**: Naive `-0.5010`, Antithetic `-0.4960`, Control Variate `-0.4690`.

## Perbandingan dengan `day001c` (Python / NumPy)

- Rumus tertutup Black-Scholes (`10.450584`) dan Kemna–Vorst (`5.565509`) identik dengan `qb100.options` hingga ketelitian relatif `0.0` (`< 1e-10`).
- Seluruh estimator MC dan Delta konsisten dengan hasil Python `day001c` dalam `< 1.2 SE` gabungan.
- Waktu eksekusi inti `200k` path: **379.8 ms** (Rust `--release`) vs **1381.7 ms** (Python/NumPy, ~3.6× lebih cepat meskipun Rust menjalankan simulasi independen per estimator tanpa menyimpan matriks path di RAM).
