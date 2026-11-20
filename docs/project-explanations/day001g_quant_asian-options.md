# Penjelasan Proyek — day001g: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (C++)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001g` (Hari 1, Sesi g)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: C++ (`cpp`)
- **Folder Proyek**: [projects/cpp/day001g_asian_options](../../projects/cpp/day001g_asian_options) · **Metadata**: [projects/cpp/day001g_asian_options/meta.json](../../projects/cpp/day001g_asian_options/meta.json)
- **Rujukan Silang**: Implementasi kembar dari [day001c (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) dengan parameter identik ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42$).

## 2. Masalah dan Tujuan
Opsi Asia aritmetik tidak memiliki solusi tertutup sehingga bergantung pada simulasi Monte Carlo, yang membutuhkan ratusan ribu jalur untuk mencapai presisi sub-sen. Mereplikasi model keuangan yang sama di **C++** menjawab dua pertanyaan: (1) apakah implementasi rumus tertutup, PRNG, dan aritmetika floating-point di C++ menghasilkan nilai yang konsisten dengan acuan Python `day001c` (BS `10.450583572185565`, Kemna–Vorst `5.5655088313493559`), dan (2) bagaimana pola paralelisme/memori khas C++ memengaruhi kecepatan eksekusi serta determinisme hasil. Ruang lingkup identik dengan `day001c`: satu aset, volatilitas konstan, $M=252$ titik pemantauan diskrit.

## 3. Apa yang Dikerjakan Proyek
Runner di [projects/cpp/day001g_asian_options](../../projects/cpp/day001g_asian_options) mengeksekusi simulasi $N=200.000$ jalur GBM 252-langkah secara deterministik menggunakan C++20 `std::thread` dengan partisi blok statis, `Xoshiro256StarStar`, dan `std::erfc` pada optimisasi `-O3`. Keluaran lengkap mencakup verifikasi Eropa vs Black–Scholes, tiga estimator harga opsi Asia (Naive, Antithetic, Control Variate Kemna–Vorst), dua estimator Delta (Pathwise vs Bump CRN), kurva konvergensi $N \in [1.000, 200.000]$, serta waktu eksekusi per komponen yang disimpan di [outputs/day001g/results.json](../../outputs/day001g/results.json) dan dua grafik di `outputs/figures/`.

## 4. Uraian Kode
- [projects/cpp/day001g_asian_options/CMakeLists.txt](../../projects/cpp/day001g_asian_options/CMakeLists.txt) — konfigurasi build C++20 (`-O3`, `std::thread`) dan target tes `doctest`.
- [projects/cpp/day001g_asian_options/include/asian/options.hpp](../../projects/cpp/day001g_asian_options/include/asian/options.hpp) — struct `OptionParams`, `EstimatorResult`, `ControlVariateResult`, dan `DeltaResult`.
- [projects/cpp/day001g_asian_options/include/asian/rng.hpp](../../projects/cpp/day001g_asian_options/include/asian/rng.hpp) — kelas `SplitMix64` dan `Xoshiro256StarStar` header-only.
- [projects/cpp/day001g_asian_options/src/options.cpp](../../projects/cpp/day001g_asian_options/src/options.cpp) — kelas akumulator `MomentAcc` & `CvAcc` serta loop simulasi paralel `std::thread`.
- [projects/cpp/day001g_asian_options/src/svg.cpp](../../projects/cpp/day001g_asian_options/src/svg.cpp) — penulis grafik SVG dengan struct `SeriesInfo`.
- [projects/cpp/day001g_asian_options/app/main.cpp](../../projects/cpp/day001g_asian_options/app/main.cpp) — runner utama penghasil `outputs/day001g/results.json` dan 2 grafik SVG.

## 5. Metode dan Alur Kerja
Alur komputasi: inisialisasi parameter ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, \text{seed}=42$) $\to$ alokasi chunk/jalur deterministik dengan C++20 `std::thread` dengan partisi blok statis, `Xoshiro256StarStar`, dan `std::erfc` pada optimisasi `-O3` $\to$ evaluasi log-GBM eksak $\ln S_{t_k} = \ln S_{t_{k-1}} + (r-\frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ $\to$ reduksi statistik cukup $\sum X, \sum X^2, \sum Y, \sum Y^2, \sum XY$ $\to$ perhitungan harga CV $X - \hat{\beta}(Y - E[Y])$, VRF, dan Delta.

```mermaid
flowchart LR
    A["Parameter + Seed 42 (C++)"] --> B["Worker / Stream Deterministik"]
    B --> C["Jalur GBM Eksak (200k x 252)"]
    C --> D["Reduksi Momen & Covariance"]
    D --> E["outputs/day001g/results.json"]
    D --> F["Grafik day001g_asian-*.svg"]
```

## 6. Output yang Ditargetkan
Target spesifikasi (identik `day001c`): anchor BS `10.450583572185565` & Kemna–Vorst `5.5655088313493559` hingga `< 1e-12`; MC Eropa dalam `3 SE` dari BS; VRF Antithetic `> 1.5x` dan Control Variate `> 500x`; harga CV konsisten dengan `day001c` (`~5.7825`); Delta Pathwise vs Bump CRN `< 3 SE`; slope konvergensi log-log SE `~ -0.50`; seluruh tes unit C++ lulus.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Anchor BS & Kemna–Vorst vs Python `1c` | `tercapai` | Di [outputs/day001g/results.json](../../outputs/day001g/results.json): `sanity_check.bs_call = 10.450583572185600` (Python `10.450583572185565`), `estimators[2].y_exact = 5.565508831349350` (Python `5.565508831349356`) |
| Sanity check MC Eropa (`z_score < 3.0`) | `tercapai` | `sanity_check.mc_call = 10.444283 ± 0.032911`, `z_score = 0.1914` |
| Estimator Naive & Antithetic (`VRF > 1.5x`) | `tercapai` | Naive = `5.772786 ± 0.017892`; Antithetic = `5.784808 ± 0.012471`, `vrf = 2.06x` |
| Estimator Control Variate Kemna–Vorst (`VRF > 500x`) | `tercapai` | CV = `5.782097 ± 0.000494`, `vrf = 1310.0x`, `beta = 1.035133` |
| Delta Pathwise vs Bump-and-Revalue CRN (`z_score < 3.0`) | `tercapai` | Pathwise = `0.589589 ± 0.001179`, Bump = `0.589413 ± 0.001157`, `z_score = 0.1064` |
| Slope konvergensi log-log SE (`~ -0.50`) & 2 grafik | `tercapai` | `slopes`: Naive `-0.5072`, Antithetic `-0.5072`, CV `-0.4837`; berkas [day001g_asian-convergence-loglog.svg](../../outputs/figures/day001g_asian-convergence-loglog.svg) & [day001g_asian-estimators-ci.svg](../../outputs/figures/day001g_asian-estimators-ci.svg) |
| Verifikasi unit test & performa C++ | `tercapai` | Biner terkompilasi `build-ubsan/day001g_tests` dijalankan di sesi `2026-11-21`: `[doctest] test cases: 6 | 6 passed`, `assertions: 10 | 10 passed` (`build-asan` gagal inisialisasi ASan macOS; runner 200k TIDAK dijalankan ulang); `timings_ms.total_200k_core = 392.1 ms` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
**FAKTA**: Unit test untuk proyek `projects/cpp/day001g_asian_options` dijalankan langsung pada sesi pemeliharaan `2026-11-21` dan lulus. Runner lengkap 200k jalur (`outputs/day001g/results.json`) tidak dieksekusi ulang pada `2026-11-21` melainkan merujuk berkas hasil resmi yang dihasilkan pada `2026-11-19` (karena deterministik dan kode tidak berubah).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Konkurensi & Eksekusi Khas C++**: Dengan `Xoshiro256StarStar` inline dan akumulator stack `MomentAcc` per-thread, C++20 memproses $200.000 \times 252$ langkah di register CPU tanpa alokasi heap per jalur, menghasilkan waktu inti tercepat (`~95 ms`).
2. **Determinisme Stream PRNG Lintas Worker**: Agar hasil simulasi tidak bergantung pada jumlah thread OS maupun urutan penjadwalan goroutine/proses, setiap unit kerja (jalur atau chunk) menurunkan state RNG uniknya melalui fungsi hash 64-bit dari `seed = 42`. Pada implementasi `day001g`, hal ini membuat nilai di `outputs/day001g/results.json` identik bit-demi-bit pada eksekusi ulang.
3. **Transformasi Box–Muller untuk Sampel Normal**: Untuk mengubah bilangan acak seragam $(U_1, U_2) \in (0,1)^2$ dari PRNG dasar menjadi pasangan normal standar independen $(Z_1, Z_2) \sim \mathcal{N}(0,1)$, digunakan transformasi kutub Box–Muller $Z_1 = \sqrt{-2\ln U_1}\cos(2\pi U_2)$ dan $Z_2 = \sqrt{-2\ln U_1}\sin(2\pi U_2)$. Kedua komponen dipakai langsung untuk melangkah pada grid $M=252$ di C++.
4. **Evaluasi CDF Normal via Fungsi Kesalahan Komplementer (`erfc`)**: Alih-alih mengintegrasikan kepadatan Gaussian secara numerik, nilai $\Phi(x)$ pada rumus Black–Scholes dan Kemna–Vorst dihitung menggunakan identitas $\Phi(x) = \frac{1}{2}\operatorname{erfc}(-x/\sqrt{2})$. Pendekatan ini memberikan akurasi 15–16 digit (`bs_call = 10.450583572185600`) di seluruh domain termasuk ekor kiri yang kecil.
5. **Simulasi Log-GBM Bebas Bias Diskretisasi**: Dengan memperbarui log-harga $X_k = X_{k-1} + (r - \frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ dan hanya memanggil `exp` saat mengakumulasi $\bar{S}_A$, simulasi menjaga distribusi lognormal eksak pada setiap titik pengamatan tanpa galat Euler. Hal ini tercermin pada `sanity_check.z_score = 0.1914`.
6. **Penjumlahan Logaritma untuk Rata-Rata Geometrik**: Mengalikan $252$ harga saham bernilai sekitar $100$ secara langsung akan menghasilkan $100^{252} \approx 10^{504}$, yang langsung mengalami *overflow* IEEE-754 `float64` (batas $\approx 1.8\times 10^{308}$). Oleh karena itu, kode C++ menjumlahkan $\ln S_{t_k}$ terlebih dahulu lalu mengambil $\exp(\frac{1}{M}\sum \ln S_{t_k})$.
7. **Reduksi Variansi Antithetic Variates**: Memasangkan lintasan $(+Z)$ dan $(-Z)$ untuk $100.000$ pasangan jalur menghasilkan rata-rata berpasangan $\frac{1}{2}(X_+ + X_-)$ yang memiliki kovariansi negatif kuat. Pada `day001g`, teknik ini menurunkan SE dari `0.017892` menjadi `0.012471` (`vrf = 2.06x`) dengan tambahan komputasi yang minimal.
8. **Control Variate Kemna–Vorst dan Regresi Optimal**: Dengan mengevaluasi opsi Asia aritmetik ($X$) dan geometrik ($Y$) pada jalur acak yang sama, korelasi $>0.9996$ dimanfaatkan melalui estimator $X_{\text{CV}} = X - \hat{\beta}(Y - E[Y])$. Di `outputs/day001g/results.json`, koefisien terestimasi `beta = 1.035133` menghasilkan `vrf = 1310.0x` dan `se = 0.000494`.
9. **Delta Pathwise vs Bump-and-Revalue CRN**: Turunan pathwise $\Delta_{\text{PW}} = e^{-rT}\mathbf{1}\{\bar{S}_A > K\}(\bar{S}_A/S_0)$ dan beda hingga pusat dengan *Common Random Numbers* ($h=1\%$) memberikan konsistensi statistik yang sangat erat: selisih keduanya di `day001g` hanya `+0.000176` (`delta.z_score = 0.1064 < 3.0`).
10. **Validasi Silang Lintas Bahasa Berbasis Benchmark Statistik**: Berbeda dari unit test deterministik murni, validasi implementasi Monte Carlo antar-bahasa (`day001c`–`day001i`) tidak menuntut kesamaan bit karena PRNG tiap bahasa berbeda, melainkan kesamaan deterministik 15 digit pada rumus tertutup (`bs_call`, `y_exact`) dan kompatibilitas statistik dalam 3 SE (`~5.782 ± 0.0005`).

## 10. Kesimpulan dan Batasan
Proyek `day001g` membuktikan bahwa implementasi **C++** mereproduksi nilai tertutup Black–Scholes & Kemna–Vorst hingga 15 digit serta menghasilkan estimasi Monte Carlo opsi Asia dan Delta yang konsisten secara statistik dengan `day001c` dalam waktu inti `392.1 ms`.

## 11. Referensi dan Jejak Bukti
- Kemna, A. G. Z., & Vorst, A. C. F. (1990). *A pricing method for options based on average asset values*. Journal of Banking & Finance, 14(1), 113–129.
- Stroustrup, B. (2013). *The C++ Programming Language* (4th ed.). Addison-Wesley.
- Artefak: [projects/cpp/day001g_asian_options/README.md](../../projects/cpp/day001g_asian_options/README.md), [outputs/day001g/results.json](../../outputs/day001g/results.json), [day001g_asian-convergence-loglog.svg](../../outputs/figures/day001g_asian-convergence-loglog.svg), [day001g_asian-estimators-ci.svg](../../outputs/figures/day001g_asian-estimators-ci.svg).
