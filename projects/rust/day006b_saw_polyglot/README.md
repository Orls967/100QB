# Hari 6b (`day006b`) — Self-Avoiding Walk (SAW) Poliglot: Rust + Go + Python

Simulasi dan enumerasi eksak **jalan menghindari-diri (*self-avoiding walk*, SAW)** pada kisi persegi dua dimensi $\mathbb{Z}^2$ menggunakan arsitektur poliglot tiga bahasa dengan pembagian peran fungsional:

- **Rust (`rust/`)**: inti komputasi tanpa *crate* eksternal (`#![forbid(unsafe_code)]`) untuk enumerasi *depth-first search* (DFS) bereduksi simetri dihedral $D_4$ ($\times 8$), pekerja awalan (*prefix worker*) berbasis pipa JSON, serta pencuplikan Markov Chain Monte Carlo (MCMC) algoritma **pivot Madras–Sokal** dengan pembangkit bilangan acak `xoshiro256++` buatan sendiri.
- **Go (`go/`)**: orkestrator konkuren pustaka standar (`os/exec`, pipa `stdin`/`stdout`, `context.WithTimeout`) yang memecah enumerasi pada awalan panjang $K = 7$ ($543$ tugas berbobot $4$), mengelola kolam proses pekerja Rust persisten, menjumlahkan bilangan bulat `uint64` secara deterministik dalam urutan tetap, serta menguji invarian jumlah pekerja (`1` vs `N`).
- **Python (`py/`)**: enumerasi *brute-force* murni tanpa simetri ($n \le 12$) sebagai jangkar independen, validator bebas-tabrakan, regresi metode rasio ($\hat{\mu}, \hat{\gamma}$) dan regresi skala log-log ($\hat{\nu}$), pembuatan berkas hasil [`results.json`](../../../outputs/day006b/results.json), serta empat visualisasi (`PNG` + `SVG`).

## Cara Menjalankan

Dari akar proyek maupun akar repositori, seluruh kompilasi dan alur kerja dapat dijalankan dengan satu skrip:

```bash
./projects/rust/day006b_saw_polyglot/run_all.sh
```

Menjalankan rangkaian pengujian unit per bahasa:

```bash
cargo test --release --manifest-path projects/rust/day006b_saw_polyglot/rust/Cargo.toml
(cd projects/rust/day006b_saw_polyglot/go && go test -v ./...)
(cd projects/rust/day006b_saw_polyglot/py && ../../../.venv/bin/python -m unittest -v)
```

## Ringkasan Hasil Utama (`seed = 20261124`)

| Besaran / Kriteria | Nilai Terukur | Acuan / Ambang | Status |
|---|---|---|---|
| $c_{12}$ (*brute-force* Python = Rust = Go) | `324.932` | Persis sama ($n \le 12$) | `tercapai` |
| $c_{18}$ (Go *sharded* $K=7$ vs Rust 1-proses) | `124.658.732` | Identik $n = 1..18$ (`0` selisih vs A001411) | `tercapai` |
| Konstanta konektivitas $\hat{\mu}$ ($n \in [11, 18]$) | `2,62900 ± 0,01760` | $\|\hat{\mu} - 2,638\| < 0,03$ | `tercapai` |
| Eksponen entropik $\hat{\gamma}$ ($n \in [11, 18]$) | `1,38338 ± 0,09593` | $\|\hat{\gamma} - 1,34375\| < 0,08$ | `tercapai` |
| Eksponen metrik Flory $\hat{\nu}$ ($N \in \{50,100,200,400\}$) | `0,74468 ± 0,00102` | $\|\hat{\nu} - 0,75\| < 0,03$ | `tercapai` |
| Validitas sampel pivot $N = 300$ | `100 / 100` valid | `0` tabrakan diri | `tercapai` |
| Waktu dinding sesi penuh (`c10_ten_minute_wall_clock`) | `666 detik` (`420 detik` pipeline komputasi) | `<= 600 detik` dari `T0` (`15:40:51` → push `15:51:57` WITA) | `belum tercapai` |

## Artefak Keluaran & Visualisasi

- Data hasil lengkap: [`outputs/day006b/results.json`](../../../outputs/day006b/results.json)
- Grafik sampel bertumpuk: [`day006b_saw_samples.svg`](../../../outputs/figures/day006b_saw_samples.svg)
- Grafik penskalaan log-log: [`day006b_saw_scaling.svg`](../../../outputs/figures/day006b_saw_scaling.svg)
- Grafik metode rasio: [`day006b_saw_ratio.svg`](../../../outputs/figures/day006b_saw_ratio.svg)
- Peta panas titik ujung $n = 14$: [`day006b_saw_endpoint_heatmap.svg`](../../../outputs/figures/day006b_saw_endpoint_heatmap.svg)
- Penjelasan teknis lengkap: [`day006b_sim_saw-polyglot.md`](../../../docs/project-explanations/day006b_sim_saw-polyglot.md)
