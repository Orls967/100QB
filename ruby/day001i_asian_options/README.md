# day001i — Opsi Asia (Ruby): Monte Carlo, Variance Reduction & Greeks

Implementasi ulang **day001c–h** (Python, Rust, Go, TypeScript, C++, C) dalam **Ruby murni (Ruby 3.3+)** tanpa gem numerik eksternal (nol dependensi *runtime*, hanya *standard library* Ruby):

- Model GBM (`S0 = 100`, `K = 100`, `r = 0.05`, `sigma = 0.20`, `T = 1.0`, `252` langkah harian, `seed = 42`).
- *Sanity check* call Eropa Black-Scholes vs Monte Carlo (`200.000` path).
- Tiga estimator *arithmetic Asian call*: **Naive MC**, **Antithetic Variates** (`100.000` pasang), dan **Control Variate** memakai opsi Asia geometrik diskret (Kemna–Vorst, 252 titik).
- **Delta** dengan metode *pathwise derivative* dan *bump-and-revalue* (*central difference*, `h = 1% * S0`, *common random numbers*).
- Analisis konvergensi log-log $\text{SE}$ vs $N \in \{1\text{k}, 2\text{k}, 5\text{k}, 10\text{k}, 20\text{k}, 50\text{k}, 100\text{k}, 200\text{k}\}$ beserta grafik SVG buatan sendiri.

## Cara Install, Test & Run

Butuh **Ruby 3.3+** (diuji pada Ruby 4.0.9 + YJIT) dan Bundler.

```bash
cd ruby/day001i_asian_options
bundle install           # (atau: bundle install --prefer-local bila memakai default gem json)
bundle exec rubocop
bundle exec rake test
RUBY_YJIT_ENABLE=1 bundle exec bin/run
```

> `bin/run` juga mengaktifkan `RubyVM::YJIT.enable` secara otomatis bila tersedia, sehingga `bundle exec bin/run` langsung berjalan dengan YJIT aktif.

Output dihasilkan ke:
- `outputs/day001i/results.json`
- `outputs/figures/day001i_asian-convergence-loglog.svg`
- `outputs/figures/day001i_asian-estimators-ci.svg`

## Struktur Kode

- `lib/asian_options/types.rb` — nilai hasil *immutable* memakai `Data.define` (`Params`, `EstimatorResult`, `ControlVariateResult`, `DeltaResult`, `ConvergencePoint`, `Slopes`, `EstimatorSummary`).
- `lib/asian_options/rng.rb` — `Random.new(seed)` terisolasi (tanpa `rand`/`srand` global), derivasi *seed* per-chunk (*SplitMix64*), dan transformasi normal **Box–Muller polar** buatan sendiri dengan `fill_normals!` *in-place*.
- `lib/asian_options/stats.rb` — CDF normal via `Math.erfc`, agregasi momen Welford/Chan paralel (univariat & bivariat), *fitting* slope log-log, VRF, dan z-score.
- `lib/asian_options/black_scholes.rb` — rumus tertutup Black-Scholes call/put dan Kemna–Vorst *discrete geometric Asian call*.
- `lib/asian_options/parallel.rb` — eksekusi paralel deterministik berbasis `Process.fork` + `IO.pipe` + `Marshal` dengan ukuran chunk tetap (`4.000` item) dan reduksi berurutan.
- `lib/asian_options/path_kernels.rb` & `lib/asian_options/chunk_simulator.rb` — loop simulasi GBM dengan buffer `Array.new(n_steps, 0.0)` pra-alokasi per-chunk.
- `lib/asian_options/estimators.rb` & `lib/asian_options/greeks.rb` — estimator harga (Eropa MC, Naive, Antithetic, Control Variate) dan Greeks (Delta *pathwise* & *bump-and-revalue* CRN).
- `lib/asian_options/svg.rb` & `lib/asian_options/convergence_svg.rb` — generator grafik SVG murni Ruby.
- `lib/asian_options/runner.rb` & `bin/run` — *orchestrator* benchmark, penulis `results.json` dan SVG.
- `test/test_asian_options.rb` & `test/test_stats_and_rng.rb` — 7 unit test Minitest (40 assertions).

## Catatan Performa Ruby (YJIT, `Process.fork` & Alokasi Memori)

Ruby murni dieksekusi di atas VM dinamis sehingga jauh lebih lambat dibanding bahasa terkompilasi (Rust/C/C++/Go) maupun NumPy tervektorisasi—dan itu memang bagian dari studi perbandingan lintas bahasa ini. Optimasi yang diterapkan di dalam batas Ruby murni:

1. **YJIT aktif (`--yjit` / `RubyVM::YJIT.enable`)**: untuk 200k path Naive MC (10 worker), waktu turun dari **1277.2 ms** (YJIT OFF) menjadi **618.0 ms** (YJIT ON) — percepatan **2.07×**.
2. **Paralelisasi multi-proses deterministik (`Process.fork` + `IO.pipe` + `Marshal`)**: menghindari batas GVL Ruby dengan tetap menjamin hasil **bit-for-bit identik** antara 1 proses dan N proses. Untuk 200k path Naive MC (YJIT ON), waktu turun dari **3812.2 ms** (1 proses) menjadi **618.0 ms** (10 proses) — percepatan **6.17×**.
3. **Buffer pra-alokasi**: `normals = Array.new(n_steps, 0.0)` dialokasikan sekali per chunk dan diisi ulang *in-place* oleh `fill_normals!`, meminimalkan tekanan GC di dalam 50,4 juta iterasi langkah.

## Ringkasan Hasil (`200.000` path, `seed = 42`)

- **Black-Scholes Eropa (tertutup)**: `10.450583572185565` vs **MC Eropa**: `10.446569 ± 0.032887` (selisih `-0.004015`, `0.12 SE`).
- **Kemna–Vorst Geometric Asian (tertutup, 252 titik)**: `5.565508831349343` (`beta = 1.035177`).

| Estimator | Harga | SE | CI 95% | VRF |
|-----------|-------|----|--------|-----|
| Naive MC (200k) | `5.776979` | `0.017855` | `[5.741983, 5.811975]` | `1.00x` |
| Antithetic (100k pasang) | `5.782943` | `0.012394` | `[5.758650, 5.807236]` | `2.08x` |
| Control Variate (Geom. Asia) | `5.782149` | `0.000490` | `[5.781189, 5.783109]` | `1329.34x` |

- **Delta Pathwise vs Bump (`h = 1%`)**: `0.589374 ± 0.001179` vs `0.589333 ± 0.001157` (selisih `+0.000041`, `0.02 SE`).
- **Slope konvergensi log-log**: Naive `-0.4981`, Antithetic `-0.4903`, Control Variate `-0.4905`.

## Perbandingan 7 Bahasa

- Rumus tertutup Black-Scholes (`10.450583572185565`, galat relatif `0.0`) dan Kemna–Vorst (`5.565508831349343`, galat relatif `2.39e-15` terhadap anchor dan `0.0` terhadap Python) identik hingga presisi mesin.
- Seluruh 36 pasangan estimator MC & Delta Ruby terhadap Python, Rust, Go, TypeScript, C++, dan C konsisten dengan **maksimum z-score `1.344`** (`< 3.0`).

| Bahasa | Euro MC z | Naive z | Anti z | CV z | Delta PW z | Delta Bump z | Waktu Inti 200k |
|--------|-----------|---------|--------|------|------------|--------------|-----------------|
| Rust | `0.504` | `0.764` | `0.666` | `0.152` | `1.103` | `1.032` | `379.8 ms` |
| C++20 | `0.049` | `0.166` | `0.106` | `0.075` | `0.129` | `0.049` | `392.1 ms` |
| C17 | `0.049` | `0.166` | `0.106` | `0.075` | `0.129` | `0.049` | `398.2 ms` |
| Go | `1.328` | `0.568` | `0.925` | `1.344` | `0.303` | `0.154` | `440.7 ms` |
| TypeScript | `0.545` | `0.032` | `0.256` | `0.748` | `0.393` | `0.331` | `1242.9 ms` |
| Python (NumPy) | `1.246` | `1.098` | `0.506` | `0.329` | `0.152` | `0.156` | `1331.5 ms` |
| **Ruby (YJIT + fork)** | — | — | — | — | — | — | **`4328.4 ms`** |
