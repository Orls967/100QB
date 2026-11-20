# day001e — Opsi Asia (Go): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c** (Python) dan **day001d** (Rust) dalam bahasa **Go** (*module* Go 1.22+), dengan spesifikasi perhitungan yang identik:

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$.

## Cara Build, Test & Run

```bash
cd projects/go/day001e_asian_options
go test ./...
go vet ./...
gofmt -l .
go run ./cmd/run
```

Output dijalankan ke:
- `outputs/day001e/results.json`
- `outputs/figures/day001e_asian-convergence-loglog.png`
- `outputs/figures/day001e_asian-estimators-ci.png`

## Struktur Kode

- `options/options.go` — rumus tertutup Black-Scholes & Kemna–Vorst, simulasi GBM paralel deterministik (*goroutine* + `sync.WaitGroup` + `math/rand/v2` PCG per-chunk), tiga estimator harga, Delta *pathwise* & *bump*, serta *fitting* slope log-log.
- `options/options_test.go` — 5 unit test (termasuk uji determinisme 1 goroutine vs 8 goroutine).
- `cmd/run/main.go` — *runner* seluruh analisis, pencatat waktu eksekusi, penulis `results.json`, dan pembuat grafik PNG (`gonum.org/v1/plot`).
- `meta.json` — metadata untuk tabel progres di `README.md` utama.

## Ringkasan Hasil (`200.000` path, `seed = 42`)

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.791336` | `0.017898` | `[5.756255, 5.826417]` | `1.00x` |
| Antithetic (100k pasang) | `5.799159` | `0.012400` | `[5.774856, 5.823463]` | `2.08x` |
| Control Variate (Geom. Asia) | `5.783085` | `0.000496` | `[5.782114, 5.784056]` | `1304.69x` |

- **Delta Pathwise vs Bump (`h=1%`)**: `0.589879 ± 0.001179` vs `0.589585 ± 0.001157` (selisih `+0.000294`, `0.18 SE`).
- **Slope konvergensi log-log**: Naive `-0.4958`, Antithetic `-0.4996`, Control Variate `-0.4911`.

## Perbandingan dengan Python (`day001c`) & Rust (`day001d`)

- Rumus tertutup Black-Scholes dan Kemna–Vorst identik hingga ketelitian relatif `< 3e-15`.
- Seluruh estimator MC dan Delta Go konsisten dalam `≤ 1.01 SE` terhadap Python dan `≤ 1.59 SE` terhadap Rust.
- Waktu eksekusi inti `200k` path: **440.7 ms** (Go) vs **379.8 ms** (Rust `--release`) vs **1381.7 ms** (Python/NumPy) — ~3.1× lebih cepat dari Python dan mendekati kecepatan Rust.
