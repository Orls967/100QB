# day001h — Opsi Asia (C17): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c–g** (Python, Rust, Go, TypeScript, C++) dalam **bahasa C murni (C17)** menggunakan CMake, `libm`, dan POSIX threads (`pthread`), tanpa dependensi eksternal:

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$.

## Cara Build, Test & Run

```bash
cd projects/c/day001h_asian_options
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/day001h_run

# Build dengan UndefinedBehaviorSanitizer:
cmake -S . -B build-ubsan -DCMAKE_BUILD_TYPE=RelWithDebInfo -DENABLE_UBSAN=ON
cmake --build build-ubsan -j
ctest --test-dir build-ubsan --output-on-failure
```

Output dijalankan ke:
- `outputs/day001h/results.json`
- `outputs/figures/day001h_asian-convergence-loglog.svg`
- `outputs/figures/day001h_asian-estimators-ci.svg`

## Struktur Kode & Keputusan Desain API / Memori

- `include/asian/rng.h` — PRNG `xoshiro256**` (`SplitMix64` *seed derivation* per-chunk) + *Polar Marsaglia* normal sampler.
- `include/asian/options.h` & `src/options.c` — Black-Scholes & Kemna–Vorst (`erfc`), simulasi GBM paralel deterministik (`pthread` per-chunk), tiga estimator harga, Delta *pathwise* & *bump*, serta *fitting* slope log-log.
- `include/asian/svg.h` & `src/svg.c` — generator grafik SVG murni C (`fprintf`).
- `app/main.c` — *runner* utama dan penulis `results.json`.
- `tests/test_options.c` — *harness* test minimal buatan sendiri (7 test case).

**Manajemen memori & API:**
- Setiap fungsi publik menerima pointer `const` + panjang (`size_t`) dan mengembalikan kode status `asian_status_t` (`ASIAN_OK`, `ASIAN_ERR_INVALID_ARG`, `ASIAN_ERR_ALLOC`, `ASIAN_ERR_ZERO_VAR`, `ASIAN_ERR_THREAD`, `ASIAN_ERR_IO`), tanpa pernah memanggil `exit()` di dalam library.
- Tidak ada alokasi heap per-path di dalam loop panas; buffer chunk (`chunk_acc_t *`, `pthread_t *`, `worker_arg_t *`) dialokasikan sekali dengan `calloc`, dicek `NULL`, dan dibebaskan deterministik oleh pemilik tunggal sebelum fungsi kembali.

## Ringkasan Hasil (`200.000` path, `seed = 42`)

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.772786` | `0.017892` | `[5.737718, 5.807853]` | `1.00x` |
| Antithetic (100k pasang) | `5.784808` | `0.012471` | `[5.760365, 5.809251]` | `2.06x` |
| Control Variate (Geom. Asia) | `5.782097` | `0.000494` | `[5.781128, 5.783066]` | `1310.00x` |

- **Delta Pathwise vs Bump (`h=1%`)**: `0.589589 ± 0.001179` vs `0.589413 ± 0.001157` (selisih `+0.000176`, `0.11 SE`).
- **Slope konvergensi log-log**: Naive `-0.5072`, Antithetic `-0.5072`, Control Variate `-0.4837`.

## Perbandingan 6 Bahasa

- Rumus tertutup Black-Scholes (`10.450583572185565`) dan Kemna–Vorst (`5.5655088313493559`) identik hingga ketelitian relatif `< 3e-15`.
- Seluruh estimator MC dan Delta C17 konsisten dalam `≤ 1.41 SE` terhadap Python, Rust, Go, TypeScript, dan C++.
- Waktu eksekusi inti `200k` path: **379.8 ms** (Rust) · **392.1 ms** (C++20) · **398.2 ms** (C17) · **440.7 ms** (Go) · **1236.6 ms** (TypeScript) · **1381.7 ms** (Python/NumPy).
