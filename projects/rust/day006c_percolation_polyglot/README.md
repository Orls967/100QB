# Hari 6c (`day006c`) — Perkolasi Bond Kisi Persegi Poliglot: Rust + Go + Python

Implementasi arsitektur poliglot tiga bahasa untuk **perkolasi bond di kisi persegi** pada persegi panjang swa-dual $L \times (L + 1)$ ($N = 2L^2 - 1$ ikatan):

- **Rust (`rust/`)**: biner `percolate` (`#![forbid(unsafe_code)]`, nol *crate* eksternal) dengan algoritma mikrokanonik **Newman–Ziff**, `Xoshiro256Pp` + pengacakan Fisher–Yates bebas-bias Lemire, serta dua struktur Union-Find (*union by size* + *path compression*) untuk penyeberangan simpul maya kiri–kanan ($n^*$) dan ukuran klaster fisik terbesar pada $n = \lfloor N/2 \rfloor$.
- **Go (`go/`)**: orkestrator `orchestrate` pustaka standar yang membagi pekerjaan ke potongan tetap $(L, \text{chunk})$ dengan *seed* `SplitMix64` independen dari jumlah pekerja, mengelola kolam proses `percolate job` persisten lewat pipa JSON di bawah `context.WithTimeout`, dan mereduksi histogram `uint64` secara deterministik.
- **Python (`py/`)**: enumerasi rasional eksak `fractions.Fraction` pada kisi kecil ([`percolation.py`](py/percolation.py)), konvolusi bobot binomial berbasis `math.lgamma`, regresi penskalaan kritis ($\hat{\nu}, \hat{D}$), pembuatan grafik, dan pencatat waktu [`finalize.py`](py/finalize.py).

## Cara Build dan Menjalankan

Dari akar repositori (tanpa memerlukan utilitas `timeout` eksternal yang tak portabel):

```bash
./projects/rust/day006c_percolation_polyglot/run_all.sh
```

Pengujian per bahasa:

```bash
cargo test --release --manifest-path projects/rust/day006c_percolation_polyglot/rust/Cargo.toml
(cd projects/rust/day006c_percolation_polyglot/go && go test -v ./...)
(cd projects/rust/day006c_percolation_polyglot/py && ../../../../.venv/bin/python -m unittest -v)
```

## Ringkasan Status Sesi (`seed = 20261124`)

| Kriteria | Status | Bukti Terukur |
|---|---|---|
| `c1_exact_self_dual_half` | `sebagian` | Enumerasi eksak $2 \times 3$ ($64/128 = 1/2$) lulus di `cargo test --release`; skrip Python belum dieksekusi saat `LOCK` 9:00 |
| `c2_nz_vs_exact_small_lattice` | `sebagian` | NZ $2 \times 3$ ($40.000$ sampel) cocok $P(1/2)=0{,}5$ di `cargo test`; grid 9 titik Python belum dieksekusi |
| `c3_pcross_half_all_L` – `c6_figures_valid` | `belum diketahui` | Pipeline penuh `run_all.sh` dihentikan oleh disiplin kunci waktu `dl.sh` pada menit ke-9 |
| `c7_worker_invariance` | `tercapai` | `TestDeterministicReductionOrder` (`1` vs `4` *worker* pada $L \in \{8,16\}$) identik bit-demi-bit (`PASS`) |
| `c_wall_clock` | `belum tercapai` | Melewati ambang `570` detik (lihat [`outputs/day006c/results.json`](../../../outputs/day006c/results.json)) |

Penjelasan lengkap: [`day006c_sim_percolation-polyglot.md`](../../../docs/project-explanations/day006c_sim_percolation-polyglot.md).
