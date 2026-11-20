# day001g — Opsi Asia (C++20): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c** (Python), **day001d** (Rust), **day001e** (Go), dan **day001f** (TypeScript) dalam **C++20 modern** (*build* CMake tanpa dependensi eksternal pada kode inti):

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$.

## Cara Build, Test & Run

```bash
cd projects/cpp/day001g_asian_options
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/day001g_run

# Build dengan UndefinedBehaviorSanitizer:
cmake -S . -B build-ubsan -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="-fsanitize=undefined -fno-omit-frame-pointer" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=undefined"
cmake --build build-ubsan -j
ctest --test-dir build-ubsan --output-on-failure
```

Output dijalankan ke:
- `outputs/day001g/results.json`
- `outputs/figures/day001g_asian-convergence-loglog.svg`
- `outputs/figures/day001g_asian-estimators-ci.svg`

## Struktur Kode

- `include/asian/rng.hpp` — PRNG `Xoshiro256**` (`SplitMix64` *seed derivation* per-chunk) + *Polar Marsaglia* normal sampler.
- `include/asian/options.hpp` & `src/options.cpp` — Black-Scholes & Kemna–Vorst (`std::erfc`), simulasi GBM paralel deterministik (`std::thread` per-chunk), tiga estimator harga, Delta *pathwise* & *bump*, serta *fitting* slope log-log.
- `include/asian/svg.hpp` & `src/svg.cpp` — generator grafik SVG murni C++.
- `app/main.cpp` — *runner* utama dan penulis `results.json`.
- `tests/test_options.cpp` — 6 unit test `doctest` (termasuk uji determinisme 1 vs 4 thread).

## Ringkasan Hasil (`200.000` path, `seed = 42`)

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.772786` | `0.017892` | `[5.737718, 5.807853]` | `1.00x` |
| Antithetic (100k pasang) | `5.784808` | `0.012471` | `[5.760365, 5.809251]` | `2.06x` |
| Control Variate (Geom. Asia) | `5.782097` | `0.000494` | `[5.781128, 5.783066]` | `1310.00x` |

- **Delta Pathwise vs Bump (`h=1%`)**: `0.589589 ± 0.001179` vs `0.589413 ± 0.001157` (selisih `+0.000176`, `0.11 SE`).
- **Slope konvergensi log-log**: Naive `-0.5072`, Antithetic `-0.5072`, Control Variate `-0.4837`.

## Perbandingan 5 Bahasa

- Rumus tertutup Black-Scholes dan Kemna–Vorst identik hingga ketelitian relatif `< 4e-15`.
- Seluruh estimator MC dan Delta C++ konsisten dalam `≤ 1.41 SE` terhadap Python, Rust, Go, dan TypeScript.
- Waktu eksekusi inti `200k` path: **392.1 ms** (C++20 Release) vs **379.8 ms** (Rust Release) vs **440.7 ms** (Go) vs **1236.6 ms** (TypeScript) vs **1381.7 ms** (Python/NumPy).
