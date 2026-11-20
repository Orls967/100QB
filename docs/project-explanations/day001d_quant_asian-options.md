# Penjelasan Proyek — day001d: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (Rust)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001d` (Hari 1, Sesi d)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` (off-roadmap) · **Jenis**: `quant` · **Bahasa**: Rust (`rust`)
- **Folder Proyek**: [projects/rust/day001d_asian_options](../../rust/day001d_asian_options) · **Metadata**: [projects/rust/day001d_asian_options/meta.json](../../rust/day001d_asian_options/meta.json)
- **Rujukan Silang**: Implementasi kembar dari [day001c (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) dengan spesifikasi parameter identik ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42$).

## 2. Masalah dan Tujuan
Bagaimana mereplikasi mesin penilaian opsi Asia aritmetik (Naive MC, Antithetic Variates, Control Variate Kemna–Vorst, serta Delta Pathwise & Bump-and-Revalue CRN) di dalam ekosistem **Rust** secara deterministik dan tervalidasi silang secara statistik terhadap Python `day001c` (anchor Black–Scholes `10.450583572185565`, Kemna–Vorst `5.5655088313493559`)?

## 3. Apa yang Dikerjakan Proyek
Program biner/runner di `projects/rust/day001d_asian_options` menjalankan simulasi $N=200.000$ jalur GBM ($M=252$ langkah), mengevaluasi sanity check Eropa, ketiga estimator harga opsi Asia, kedua estimator Delta, dan kurva konvergensi $N \in [1.000, 200.000]$, lalu menulis hasil ke [outputs/day001d/results.json](../../outputs/day001d/results.json) serta grafik SVG ke `outputs/figures/day001d_*.svg`.

## 4. Uraian Kode
- [projects/rust/day001d_asian_options/Cargo.toml](../../rust/day001d_asian_options/Cargo.toml): crate Rust edisi 2021 dengan `rayon`, `rand_chacha` (`ChaCha8Rng`), `statrs`, `serde_json`, dan `plotters`.
- [projects/rust/day001d_asian_options/src/lib.rs](../../rust/day001d_asian_options/src/lib.rs): `OptionParams`, `bs_call_price`, `geometric_asian_call_price`, simulator GBM paralel `rayon`, estimator Naive/Antithetic/Control Variate, `asian_delta_pathwise`, dan `asian_delta_bump`.
- [projects/rust/day001d_asian_options/src/main.rs](../../rust/day001d_asian_options/src/main.rs): runner CLI penghasil `outputs/day001d/results.json` dan grafik SVG.

## 5. Metode dan Alur Kerja
Arsitektur spesifik **Rust**: menggunakan ChaCha8Rng (`rand_chacha`) dengan seeding per jalur deterministik `seed + path_id` di atas work-stealing thread pool `rayon`. Setiap worker menghitung statistik cukup lokal pada jalur yang ditugaskan, lalu direduksi secara deterministik.

## 6. Output yang Ditargetkan
Cocok dengan anchor Black–Scholes (`10.450583572185565`) dan Kemna–Vorst (`5.5655088313493559`), MC Eropa dalam 3 SE dari BS, VRF Antithetic $\approx 2\times$ dan Control Variate $> 1200\times$, harga CV konsisten lintas bahasa di sekitar `~5.782`, Delta pathwise vs bump dalam 3 SE, dan slope log-log SE $\approx -0.50$.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Anchor tertutup BS & Kemna–Vorst vs Python `1c` | `tercapai` | Di [outputs/day001d/results.json](../../outputs/day001d/results.json): `bs_call = 10.450583572185565`, `y_exact = 5.565508831349343` |
| Sanity check MC Eropa (`z < 3`) | `tercapai` | MC Eropa = `10.470029 ± 0.032991`, `z_score = 0.5894 < 3.0` |
| Estimator Naive, Antithetic & Control Variate | `tercapai` | Naive = `5.796308 ± 0.017942`, Antithetic = `5.771260 ± 0.012407` (`VRF = 2.09x`), Control Variate = `5.782043 ± 0.000496` (`VRF = 1311.2x`, `beta = 1.035248`) |
| Konsistensi silang harga CV terhadap Python `1c` (`5.782468`) | `tercapai` | Selisih `-0.000425` konsisten dengan SE `0.000496` |
| Delta Pathwise vs Bump-and-Revalue CRN (`z < 3`) | `tercapai` | Pathwise = `0.591213 ± 0.001179`, Bump = `0.591020 ± 0.001157`, `z_score = 0.1166 < 3.0` |
| Slope konvergensi log-log SE ($\approx -0.50$) | `tercapai` | Slopes: Naive `-0.5010`, Antithetic `-0.4960`, Control Variate `-0.4690` |
| Verifikasi test suite Rust | `tercapai` | `cargo test --release -q` dijalankan pada sesi `day003c`: 5/5 unit test lulus dalam 0.07s; waktu inti 200k jalur = `379.8 ms` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi. **FAKTA**: Test suite `projects/rust/day001d_asian_options` telah dijalankan ulang pada sesi `day003c` dan lulus bersih.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Arsitektur Konkurensi & Memori Rust**: Trait iterator paralel `rayon::prelude::ParallelIterator` mereduksi akumulator `(sum, sum_sq, sum_xy)` secara *lock-free* dan deterministik karena setiap indeks jalur menurunkan stream `ChaCha8Rng` sendiri.
2. **Determinisme PRNG Lintas Worker**: Pembagian stream bilangan acak berbasis indeks jalur/chunk memastikan hasil simulasi tidak bergantung pada penjadwalan thread/proses OS.
3. **Transformasi Box–Muller Normal Standar**: Konversi pasangan seragam $(U_1, U_2) \in (0,1)^2$ menjadi normal standar independen $Z = \sqrt{-2\ln U_1}\cos(2\pi U_2)$.
4. **Komputasi CDF Normal $\Phi(x)$ via `erfc`**: Identitas $\Phi(x) = \frac{1}{2}\operatorname{erfc}(-x/\sqrt{2})$ menjaga akurasi di ekor distribusi.
5. **Simulasi Eksak Log-GBM Diskrit**: Akumulasi log-harga $\ln S_{t_k} = \ln S_{t_{k-1}} + (r - \frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ bebas bias diskretisasi waktu.
6. **Rata-Rata Geometrik via Penjumlahan Logaritma**: $\bar{S}_G = \exp(\frac{1}{M}\sum_{k=1}^M \ln S_{t_k})$ mencegah *overflow/underflow* perkalian $M=252$ harga saham.
7. **Koreksi Control Variate Kemna–Vorst**: Memanfaatkan korelasi $\rho > 0.9996$ antara payoff aritmetik dan geometrik pada jalur yang sama untuk mencapai VRF $> 1300\times$.
8. **Antithetic Variates pada Ruang Jalur**: Pembalikan tanda vektor shock $-Z$ pada langkah yang sama memangkas pembangkitan normal menjadi setengahnya.
9. **Estimator Delta Pathwise (IPA)**: $\Delta_{\text{PW}} = e^{-rT} \mathbf{1}\{\bar{S}_A > K\} (\bar{S}_A / S_0)$ dihitung dalam satu sapuan bersama harga opsi.
10. **Common Random Numbers (CRN) pada Beda Hingga**: Penggunaan ulang stream PRNG yang identik untuk $S_0(1+h)$ dan $S_0(1-h)$ menjaga `z_score` selisih terhadap Pathwise di bawah `0.2` SE.

## 10. Kesimpulan dan Batasan
Proyek `day001d` menunjukkan bahwa spesifikasi kuantitatif opsi Asia di `day001c` dapat direplikasi di **Rust** dengan kesesuaian rumus tertutup 15 digit dan konsistensi statistik penuh.

## 11. Referensi dan Jejak Bukti
- Rujukan induk: [day001c_quant_opsi-asia-monte-carlo-variance.md](day001c_quant_opsi-asia-monte-carlo-variance.md).
- Artefak: [projects/rust/day001d_asian_options/README.md](../../rust/day001d_asian_options/README.md), [outputs/day001d/results.json](../../outputs/day001d/results.json), [day001d_asian-se-convergence.svg](../../outputs/figures/day001d_asian-se-convergence.svg), [day001d_asian-ci-comparison.svg](../../outputs/figures/day001d_asian-ci-comparison.svg).
