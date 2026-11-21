# day004f — Port Lengkap AFNS Term Structure, Kalibrasi Interval, Harga Listrik & Valuasi Proyek Baterai (Rust)

Proyek **Rust stable** mandiri (*zero external crates*, `edition = "2021"`) yang memport seluruh eksperimen, jangkar analitik, diagnostik, ukuran sampel, dan seed dari `day004a`–`day004e` (termasuk audit *expanded-tail* `[-1000, 5000] EUR/MWh`), serta menambahkan **Panel F**: valuasi proyek baterai 10 tahun yang didiskontokan dengan kurva zero-coupon AFNS.

> **Catatan Ruang Lingkup**: Seluruh hasil dalam proyek ini berasal dari **model sintetis** (AFNS + pasar listrik 3-rezim sintetis). Panel F adalah studi valuasi *risk-free-discounted expected cashflow* dengan asumsi bebas korelasi antara harga listrik dan suku bunga, **bukan** harga pasar *arbitrage-free* maupun rekomendasi investasi.

## 1. Inventaris Lengkap Port `day004a`–`day004e` + Panel F

| Sumber | Komponen & Eksperimen yang Diport | Modul Rust (`src/`) | Hasil JSON (`results.json`) | Grafik SVG (`outputs/figures/`) | Tes Unit |
|---|---|---|---|---|---|
| `day004a` | NS/NSS, Taylor `tau` kecil, bootstrap par `0.5-30y`, interpolasi `R=200`, risiko obligasi/KRD/swap/hedge, AFNS Riccati RK4, Van Loan & GL, TP 10y, Kalman/RTS anchor `T=5`, MLE `R=50`, DNS misspec, PCA `T=5000, R=100` | `ns.rs`, `bootstrap.rs`, `interp.rs`, `risk.rs`, `afns.rs`, `kalman.rs`, `optim.rs`, `pca.rs`, `panel_a.rs` | `panel_a_day004a` | `day004f_ns-afns-curves-and-adjustment.svg`, `day004f_kalman-pca-horizon-diagnostics.svg` | `ns::tests`, `bootstrap::tests`, `risk::tests`, `afns::tests`, `kalman::tests` |
| `day004b` | Baseline `1000w` (`R=50`), *long-span* `10/20/40/80y` *shared-prefix* (`R=50`), *infill* `20y` bulanan/mingguan/harian (`R=50`), tren OLS/WLS `1/H`, diagnostik BFGS & *line-search* | `panel_b.rs`, `kalman.rs`, `optim.rs`, `stats.rs` | `panel_b_day004b` | `day004f_kalman-pca-horizon-diagnostics.svg` | `rng::tests`, `parallel::tests`, `linalg::tests` |
| `day004c` | Jackknife temporal `20y` & `80y` (`R=50`), uji aljabar `c/H`, bootstrap parametrik `20y` (`R=50 x B=100 = 5.000` refit), *unclamped* & penanda non-fisik | `jackknife.rs`, `panel_c.rs` | `panel_c_day004c` | `day004f_bias-correction-and-interval-coverage.svg` | `jackknife::tests` |
| `day004d` | Kalibrasi 5 metode interval (Wald, Percentile, Basic, Studentized, Log-Basic) pada `H=20y`, `R=100 x B=100 = 10.000` refit dengan Hessian `SE*`, ringkasan lebar robust & stratifikasi `kp11_hat` | `interval.rs`, `panel_d.rs` | `panel_d_day004d` | `day004f_bias-correction-and-interval-coverage.svg` | `interval::tests` |
| `day004e` | Pasar listrik 3-rezim + *compound Poisson*, momen MC `200k`, baterai `4 MWh`, DP stokastik 17/33 bin (`10.000` draw/state), audit *expanded-tail* `[-1000, 5000]` (38 bin CRN), evaluasi CRN `N=2.000` hari vs Threshold & Perfect Foresight | `power.rs`, `battery.rs`, `dp.rs`, `panel_e.rs` | `panel_e_day004e` | `day004f_electricity-prices-regimes-jumps.svg`, `day004f_battery-dispatch-and-strategies.svg`, `day004f_expanded-tail-grid-sensitivity.svg` | `power::tests`, `battery::tests`, `dp::tests` |
| **Panel F** | Valuasi proyek baterai 10 tahun (*resampling* DP-33 `1.000x365` hari/tahun), kurva dasar AFNS vs ketidakpastian fit `4b/4c/4d`, NPV pre-CAPEX, *break-even CAPEX*, skenario CAPEX/O&M ilustratif, durasi/konveksitas/KRD/PV01, *hedge* dua obligasi | `integration.rs` | `panel_f_integration` | `day004f_battery-afns-npv-krd.svg` | `integration::tests` |

## 2. Ringkasan Hasil Utama (`outputs/day004f/results.json`, `587,02s` pada `8` worker)

- **Seluruh 7 kriteria teknis (`c1..c7`) berstatus `tercapai`** (`all_technical_passed = true`).
- **Panel A (`4,86s`)**: galat integral NS/NSS `4,16e-17`, bootstrap par `2,78e-16`, Riccati RK4 `1,36e-12`, Van Loan `0`, jangkar Gaussian `T=5` `ll_diff = 2,42e-13`, `mean_diff = 6,02e-16`, `cov_diff = 0`; PCA `T=5000, R=100` lulus kriteria agregat (`cosine top-3 = [0,99995, 0,99927, 0,99745]`, perubahan tanda `[0, 1, 2]`, `noise_bulk_z = 0,46`).
- **Panel B (`46,33s`)**: bias `kp11` meluruh `0,3906 (10y) -> 0,2291 (20y) -> 0,1138 (40y) -> 0,0601 (80y)` (slope OLS `3,76 ± 0,58`, WLS `4,05 ± 0,48`); *infill* harian tidak menghilangkan bias `kp11` (`RMSE 0,2565`) sementara `lambda/sigma_eps` turun tajam (`0,0011` / `1,8e-6`).
- **Panel C (`158,94s`)**: `5.000/5.000` refit bootstrap *finite*; jackknife 20y menurunkan bias `kp11` (`+0,2164 -> -0,0072`) namun `22/50` hasil non-fisik (`<= 0`); bootstrap persentil `kp11` cakupan `8,0%`, `theta1` `28,0%`.
- **Panel D (`376,49s`)**: `10.000/10.000` refit *finite*; pada `kp11`, `Basic` mencapai cakupan `94,0%` tetapi `98/100` batas bawah non-fisik (`0/100` cakupan valid secara fisik), `Log-Basic` `0/100` non-fisik dengan cakupan `44,0%`; pada `theta1`, `Studentized` menaikkan cakupan ke `73,0%` (median lebar `0,0406`, `IQR 0,0695`, mean `23,68` akibat outlier regularisasi Hessian).
- **Panel E (`0,41s`)**: DP-33 mean profit `250,45 ± 5,05 EUR/hari` vs DP-17 `244,92 ± 5,10` vs Threshold `127,47` vs Perfect Foresight `502,19`. Audit *expanded-tail* `[-1000, 5000]` (38 bin CRN) menurunkan massa tail atas dari `1,40%` menjadi `0,00%` dan bawah dari `0,43%` menjadi `0,0057%`, dengan kesepakatan aksi evaluasi `100%` terhadap DP-33 (`diff = 0,00 EUR/hari`).
- **Panel F (`0,003s`)**: NPV pre-CAPEX proyek 10 tahun pada kurva AFNS dasar = **`734.424,80 EUR`** (break-even CAPEX sama), *resampling* arus kas `733.934,43 ± 10.942,57 EUR` (`Q05–Q95: [715.899, 751.671]`). Ketidakpastian parameter AFNS (fit `4b/4c/4d`) menghasilkan SD NPV `56,6k–79,4k EUR` (`~5–7x` lebih besar daripada variabilitas arus kas harian). Durasi efektif `5,1399` tahun (`rel_err = 2,44e-11`), konveksitas `34,5878` (`rel_err = 5,25e-8`), `PV01 = 377,49 EUR/bp`, dan *hedge* 5y+10y menetralkan KRD 5y & 10y hingga `< 4e-12`.

## 3. Cara Build, Test, dan Run

```bash
cd projects/rust/day004f_energy_transition
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --release
cargo run --release --bin run          # Runner penuh (~587s pada 8 worker)
cargo run --release --bin run -- --quick  # Mode cepat ke temp dir (~7s)
```

## 4. Artefak Terkait

- JSON keluaran: [`outputs/day004f/results.json`](../../../outputs/day004f/results.json)
- 7 grafik SVG: [`outputs/figures/day004f_*.svg`](../../../outputs/figures/)
- Penjelasan proyek 11 bagian: [`day004f_rust-afns-energy-storage.md`](../../../docs/project-explanations/day004f_rust-afns-energy-storage.md)
