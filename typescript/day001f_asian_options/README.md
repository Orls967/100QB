# day001f — Opsi Asia (TypeScript): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c** (Python), **day001d** (Rust), dan **day001e** (Go) dalam **TypeScript** (*package* npm ESM dengan `"strict": true`), tanpa dependensi *runtime* eksternal:

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$.

## Cara Install, Test & Run

```bash
cd typescript/day001f_asian_options
npm ci
npm run typecheck
npm run lint
npm run format:check
npm test
npm run build
npm run start
```

Output dijalankan ke:
- `outputs/day001f/results.json`
- `outputs/figures/day001f_asian-convergence-loglog.svg`
- `outputs/figures/day001f_asian-estimators-ci.svg`

## Struktur Kode

- `src/rng.ts` — PRNG *seedable* `Xoshiro128++` (`SplitMix64` *seed derivation* per-chunk) + *Polar Marsaglia* normal sampler.
- `src/stats.ts` — aproksimasi `erf` Cody berpresisi $< 10^{-15}$, `normCdf`, ringkasan momen, dan *fitting* slope log-log.
- `src/core.ts` — rumus tertutup Black-Scholes & Kemna–Vorst serta kernel simulasi GBM bebas alokasi per-path.
- `src/worker.ts` & `src/options.ts` — paralelisme deterministik memakai `worker_threads` (hasil identik 1 worker vs banyak worker).
- `src/svg.ts` & `src/run.ts` — generator grafik SVG murni TypeScript dan *runner* utama.
- `tests/options.test.ts` — 6 unit test Vitest.

## Ringkasan Hasil (`200.000` path, `seed = 42`)

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.776164` | `0.017906` | `[5.741068, 5.811259]` | `1.00x` |
| Antithetic (100k pasang) | `5.787432` | `0.012429` | `[5.763070, 5.811793]` | `2.08x` |
| Control Variate (Geom. Asia) | `5.781629` | `0.000494` | `[5.780661, 5.782597]` | `1314.92x` |

- **Delta Pathwise vs Bump (`h=1%`)**: `0.588719 ± 0.001179` vs `0.588792 ± 0.001157` (selisih `-0.000073`, `0.04 SE`).
- **Slope konvergensi log-log**: Naive `-0.5020`, Antithetic `-0.5048`, Control Variate `-0.4926`.

## Perbandingan dengan Python (`day001c`), Rust (`day001d`) & Go (`day001e`)

- Rumus tertutup Black-Scholes (`10.450584`) dan Kemna–Vorst (`5.565509`) identik dengan ketelitian relatif `0.0` (`< 1e-10`).
- Seluruh estimator MC dan Delta TypeScript konsisten dalam `≤ 1.79 SE` (vs Python), `≤ 1.50 SE` (vs Rust), dan `≤ 2.08 SE` (vs Go).
- Waktu eksekusi inti `200k` path: **1242.9 ms** (TypeScript/Node.js) vs **379.8 ms** (Rust `--release`) vs **440.7 ms** (Go) vs **1381.7 ms** (Python/NumPy).
