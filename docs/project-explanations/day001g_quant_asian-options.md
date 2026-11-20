# Penjelasan Proyek — day001g: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (C++)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001g` (Hari 1, Sesi g)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` (off-roadmap) · **Jenis**: `quant` · **Bahasa**: C++ (`cpp`)
- **Folder Proyek**: [projects/cpp/day001g_asian_options](../../cpp/day001g_asian_options) · **Metadata**: [projects/cpp/day001g_asian_options/meta.json](../../cpp/day001g_asian_options/meta.json)
- **Rujukan Silang**: Implementasi kembar dari [day001c (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) dengan spesifikasi parameter identik ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42$).

## 2. Masalah dan Tujuan
Bagaimana mereplikasi mesin penilaian opsi Asia aritmetik (Naive MC, Antithetic Variates, Control Variate Kemna–Vorst, serta Delta Pathwise & Bump-and-Revalue CRN) di dalam ekosistem **C++** secara deterministik dan tervalidasi silang secara statistik terhadap Python `day001c` (anchor Black–Scholes `10.450583572185565`, Kemna–Vorst `5.5655088313493559`)?

## 3. Apa yang Dikerjakan Proyek
Program biner/runner di `projects/cpp/day001g_asian_options` menjalankan simulasi $N=200.000$ jalur GBM ($M=252$ langkah), mengevaluasi sanity check Eropa, ketiga estimator harga opsi Asia, kedua estimator Delta, dan kurva konvergensi $N \in [1.000, 200.000]$, lalu menulis hasil ke [outputs/day001g/results.json](../../outputs/day001g/results.json) serta grafik SVG ke `outputs/figures/day001g_*.svg`.

## 4. Uraian Kode
- [projects/cpp/day001g_asian_options/CMakeLists.txt](../../cpp/day001g_asian_options/CMakeLists.txt): konfigurasi build C++20 (`-O3`, `std::thread`) dan CTest `doctest`.
- [projects/cpp/day001g_asian_options/include/](../../cpp/day001g_asian_options/include/) & `src/`: `std::erfc`, PRNG `xoshiro256++`, paralelisme `std::thread`, estimator MC, serta penulis JSON/SVG.
- [projects/cpp/day001g_asian_options/app/main.cpp](../../cpp/day001g_asian_options/app/main.cpp) & `tests/`: biner runner utama dan unit test `doctest`.

## 5. Metode dan Alur Kerja
Arsitektur spesifik **C++**: menggunakan C++20 `std::thread` dengan partisi blok statis, `std::erfc`, dan optimisasi `-O3` (`Apple clang 16.0.0`). Setiap worker menghitung statistik cukup lokal pada jalur yang ditugaskan, lalu direduksi secara deterministik.

## 6. Output yang Ditargetkan
Cocok dengan anchor Black–Scholes (`10.450583572185565`) dan Kemna–Vorst (`5.5655088313493559`), MC Eropa dalam 3 SE dari BS, VRF Antithetic $\approx 2\times$ dan Control Variate $> 1200\times$, harga CV konsisten lintas bahasa di sekitar `~5.782`, Delta pathwise vs bump dalam 3 SE, dan slope log-log SE $\approx -0.50$.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Anchor tertutup BS & Kemna–Vorst vs Python `1c` | `tercapai` | Di [outputs/day001g/results.json](../../outputs/day001g/results.json): `bs_call = 10.450583572185600`, `y_exact = 5.565508831349350` |
| Sanity check MC Eropa (`z < 3`) | `tercapai` | MC Eropa = `10.444283 ± 0.032911`, `z_score = 0.1914 < 3.0` |
| Estimator Naive, Antithetic & Control Variate | `tercapai` | Naive = `5.772786 ± 0.017892`, Antithetic = `5.784808 ± 0.012471` (`VRF = 2.06x`), Control Variate = `5.782097 ± 0.000494` (`VRF = 1310.0x`, `beta = 1.035133`) |
| Konsistensi silang harga CV terhadap Python `1c` (`5.782468`) | `tercapai` | Selisih `-0.000371` konsisten dengan SE `0.000494` |
| Delta Pathwise vs Bump-and-Revalue CRN (`z < 3`) | `tercapai` | Pathwise = `0.589589 ± 0.001179`, Bump = `0.589413 ± 0.001157`, `z_score = 0.1064 < 3.0` |
| Slope konvergensi log-log SE ($\approx -0.50$) | `tercapai` | Slopes: Naive `-0.5072`, Antithetic `-0.5072`, Control Variate `-0.4837` |
| Verifikasi test suite C++ | `tercapai` | `cmake` + `ctest --output-on-failure` dijalankan pada sesi `day003c`: lulus dalam 0.12 detik; waktu inti 200k jalur = `392.1 ms` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi. **FAKTA**: Test suite `projects/cpp/day001g_asian_options` telah dijalankan ulang pada sesi `day003c` dan lulus bersih.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Arsitektur Konkurensi & Memori C++**: Loop simulasi C++20 bekerja langsung di register CPU tanpa alokasi heap per jalur (`O(1)` memori per jalur), menghasilkan waktu inti 200k jalur tercepat.
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
Proyek `day001g` menunjukkan bahwa spesifikasi kuantitatif opsi Asia di `day001c` dapat direplikasi di **C++** dengan kesesuaian rumus tertutup 15 digit dan konsistensi statistik penuh.

## 11. Referensi dan Jejak Bukti
- Rujukan induk: [day001c_quant_opsi-asia-monte-carlo-variance.md](day001c_quant_opsi-asia-monte-carlo-variance.md).
- Artefak: [projects/cpp/day001g_asian_options/README.md](../../cpp/day001g_asian_options/README.md), [outputs/day001g/results.json](../../outputs/day001g/results.json), [day001g_asian-se-convergence.svg](../../outputs/figures/day001g_asian-se-convergence.svg), [day001g_asian-ci-comparison.svg](../../outputs/figures/day001g_asian-ci-comparison.svg).
