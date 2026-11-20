# Penjelasan Proyek — day001i: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (Ruby)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001i` (Hari 1, Sesi i)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Ruby (`ruby`)
- **Folder Proyek**: [projects/ruby/day001i_asian_options](../../projects/ruby/day001i_asian_options) · **Metadata**: [projects/ruby/day001i_asian_options/meta.json](../../projects/ruby/day001i_asian_options/meta.json)
- **Rujukan Silang**: Implementasi kembar dari [day001c (Python)](day001c_quant_opsi-asia-monte-carlo-variance.md) dengan parameter identik ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42$).

## 2. Masalah dan Tujuan
Opsi Asia aritmetik tidak memiliki solusi tertutup sehingga bergantung pada simulasi Monte Carlo, yang membutuhkan ratusan ribu jalur untuk mencapai presisi sub-sen. Mereplikasi model keuangan yang sama di **Ruby** menjawab dua pertanyaan: (1) apakah implementasi rumus tertutup, PRNG, dan aritmetika floating-point di Ruby menghasilkan nilai yang konsisten dengan acuan Python `day001c` (BS `10.450583572185565`, Kemna–Vorst `5.5655088313493559`), dan (2) bagaimana pola paralelisme/memori khas Ruby memengaruhi kecepatan eksekusi serta determinisme hasil. Ruang lingkup identik dengan `day001c`: satu aset, volatilitas konstan, $M=252$ titik pemantauan diskrit.

## 3. Apa yang Dikerjakan Proyek
Runner di [projects/ruby/day001i_asian_options](../../projects/ruby/day001i_asian_options) mengeksekusi simulasi $N=200.000$ jalur GBM 252-langkah secara deterministik menggunakan Ruby 4.0.9 dengan `RubyVM::YJIT.enable` dan paralelisme multi-proses `Process.fork` melewati *Global VM Lock* (GVL). Keluaran lengkap mencakup verifikasi Eropa vs Black–Scholes, tiga estimator harga opsi Asia (Naive, Antithetic, Control Variate Kemna–Vorst), dua estimator Delta (Pathwise vs Bump CRN), kurva konvergensi $N \in [1.000, 200.000]$, serta waktu eksekusi per komponen yang disimpan di [outputs/day001i/results.json](../../outputs/day001i/results.json) dan dua grafik di `outputs/figures/`.

## 4. Uraian Kode
- [projects/ruby/day001i_asian_options/Gemfile](../../projects/ruby/day001i_asian_options/Gemfile) & [Rakefile](../../projects/ruby/day001i_asian_options/Rakefile) — proyek Ruby murni tanpa gem runtime eksternal, menggunakan `minitest` dan `rubocop`.
- [projects/ruby/day001i_asian_options/lib/asian_options/black_scholes.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/black_scholes.rb) — fungsi `bs_call`, `bs_put`, `geometric_asian_call`, dan `geometric_log_moments` via `Math.erfc`.
- [projects/ruby/day001i_asian_options/lib/asian_options/chunk_simulator.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/chunk_simulator.rb) — fungsi `run_european`, `run_naive`, `run_antithetic`, dan `run_control_variate`.
- [projects/ruby/day001i_asian_options/lib/asian_options/estimators.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/estimators.rb) & [greeks.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/greeks.rb) — orkestrasi estimator harga dan `delta_pathwise` vs `delta_bump`.
- [projects/ruby/day001i_asian_options/lib/asian_options/parallel.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/parallel.rb) — `map_chunks`, `run_forked`, dan `spawn_worker` via `Process.fork` + `IO.pipe` untuk menembus GVL CRuby.
- [projects/ruby/day001i_asian_options/lib/asian_options/convergence_svg.rb](../../projects/ruby/day001i_asian_options/lib/asian_options/convergence_svg.rb) — kelas `ConvergenceChart` penulis SVG mandiri.

## 5. Metode dan Alur Kerja
Alur komputasi: inisialisasi parameter ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, \text{seed}=42$) $\to$ alokasi chunk/jalur deterministik dengan Ruby 4.0.9 dengan `RubyVM::YJIT.enable` dan paralelisme multi-proses `Process.fork` melewati *Global VM Lock* (GVL) $\to$ evaluasi log-GBM eksak $\ln S_{t_k} = \ln S_{t_{k-1}} + (r-\frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ $\to$ reduksi statistik cukup $\sum X, \sum X^2, \sum Y, \sum Y^2, \sum XY$ $\to$ perhitungan harga CV $X - \hat{\beta}(Y - E[Y])$, VRF, dan Delta.

```mermaid
flowchart LR
    A["Parameter + Seed 42 (Ruby)"] --> B["Worker / Stream Deterministik"]
    B --> C["Jalur GBM Eksak (200k x 252)"]
    C --> D["Reduksi Momen & Covariance"]
    D --> E["outputs/day001i/results.json"]
    D --> F["Grafik day001i_asian-*.svg"]
```

## 6. Output yang Ditargetkan
Target spesifikasi (identik `day001c`): anchor BS `10.450583572185565` & Kemna–Vorst `5.5655088313493559` hingga `< 1e-12`; MC Eropa dalam `3 SE` dari BS; VRF Antithetic `> 1.5x` dan Control Variate `> 500x`; harga CV konsisten dengan `day001c` (`~5.7825`); Delta Pathwise vs Bump CRN `< 3 SE`; slope konvergensi log-log SE `~ -0.50`; seluruh tes unit Ruby lulus.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Anchor BS & Kemna–Vorst vs Python `1c` | `tercapai` | Di [outputs/day001i/results.json](../../outputs/day001i/results.json): `sanity_check.bs_call = 10.450583572185565` (Python `10.450583572185565`), `estimators[2].y_exact = 5.565508831349343` (Python `5.565508831349356`) |
| Sanity check MC Eropa (`z_score < 3.0`) | `tercapai` | `sanity_check.mc_call = 10.446569 ± 0.032887`, `z_score = 0.1221` |
| Estimator Naive & Antithetic (`VRF > 1.5x`) | `tercapai` | Naive = `5.776979 ± 0.017855`; Antithetic = `5.782943 ± 0.012394`, `vrf = 2.08x` |
| Estimator Control Variate Kemna–Vorst (`VRF > 500x`) | `tercapai` | CV = `5.782149 ± 0.000490`, `vrf = 1329.3x`, `beta = 1.035177` |
| Delta Pathwise vs Bump-and-Revalue CRN (`z_score < 3.0`) | `tercapai` | Pathwise = `0.589374 ± 0.001179`, Bump = `0.589333 ± 0.001157`, `z_score = 0.0248` |
| Slope konvergensi log-log SE (`~ -0.50`) & 2 grafik | `tercapai` | `slopes`: Naive `-0.4981`, Antithetic `-0.4903`, CV `-0.4905`; berkas [day001i_asian-convergence-loglog.svg](../../outputs/figures/day001i_asian-convergence-loglog.svg) & [day001i_asian-estimators-ci.svg](../../outputs/figures/day001i_asian-estimators-ci.svg) |
| Verifikasi unit test & performa Ruby | `tercapai` | `bundle exec rake test` dijalankan di sesi `2026-11-21`: `7 runs, 40 assertions, 0 failures, 0 errors` (`5.76s`); runner `bin/run` TIDAK dijalankan ulang di sesi ini; `timings_ms.total_200k_core = 4328.4 ms` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
**FAKTA**: Unit test untuk proyek `projects/ruby/day001i_asian_options` dijalankan langsung pada sesi pemeliharaan `2026-11-21` dan lulus. Runner lengkap 200k jalur (`outputs/day001i/results.json`) tidak dieksekusi ulang pada `2026-11-21` melainkan merujuk berkas hasil resmi yang dihasilkan pada `2026-11-19` (karena deterministik dan kode tidak berubah).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Konkurensi & Eksekusi Khas Ruby**: Karena CRuby mengunci eksekusi bytecode paralel lewat GVL, `parallel.rb` mem-fork proses anak per core (`run_forked`) dan mengirimkan hasil terurut via `IO.pipe`, sementara YJIT mengompilasi loop simulasi panas ke kode mesin natif.
2. **Determinisme Stream PRNG Lintas Worker**: Agar hasil simulasi tidak bergantung pada jumlah thread OS maupun urutan penjadwalan goroutine/proses, setiap unit kerja (jalur atau chunk) menurunkan state RNG uniknya melalui fungsi hash 64-bit dari `seed = 42`. Pada implementasi `day001i`, hal ini membuat nilai di `outputs/day001i/results.json` identik bit-demi-bit pada eksekusi ulang.
3. **Transformasi Box–Muller untuk Sampel Normal**: Untuk mengubah bilangan acak seragam $(U_1, U_2) \in (0,1)^2$ dari PRNG dasar menjadi pasangan normal standar independen $(Z_1, Z_2) \sim \mathcal{N}(0,1)$, digunakan transformasi kutub Box–Muller $Z_1 = \sqrt{-2\ln U_1}\cos(2\pi U_2)$ dan $Z_2 = \sqrt{-2\ln U_1}\sin(2\pi U_2)$. Kedua komponen dipakai langsung untuk melangkah pada grid $M=252$ di Ruby.
4. **Evaluasi CDF Normal via Fungsi Kesalahan Komplementer (`erfc`)**: Alih-alih mengintegrasikan kepadatan Gaussian secara numerik, nilai $\Phi(x)$ pada rumus Black–Scholes dan Kemna–Vorst dihitung menggunakan identitas $\Phi(x) = \frac{1}{2}\operatorname{erfc}(-x/\sqrt{2})$. Pendekatan ini memberikan akurasi 15–16 digit (`bs_call = 10.450583572185565`) di seluruh domain termasuk ekor kiri yang kecil.
5. **Simulasi Log-GBM Bebas Bias Diskretisasi**: Dengan memperbarui log-harga $X_k = X_{k-1} + (r - \frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k$ dan hanya memanggil `exp` saat mengakumulasi $\bar{S}_A$, simulasi menjaga distribusi lognormal eksak pada setiap titik pengamatan tanpa galat Euler. Hal ini tercermin pada `sanity_check.z_score = 0.1221`.
6. **Penjumlahan Logaritma untuk Rata-Rata Geometrik**: Mengalikan $252$ harga saham bernilai sekitar $100$ secara langsung akan menghasilkan $100^{252} \approx 10^{504}$, yang langsung mengalami *overflow* IEEE-754 `float64` (batas $\approx 1.8\times 10^{308}$). Oleh karena itu, kode Ruby menjumlahkan $\ln S_{t_k}$ terlebih dahulu lalu mengambil $\exp(\frac{1}{M}\sum \ln S_{t_k})$.
7. **Reduksi Variansi Antithetic Variates**: Memasangkan lintasan $(+Z)$ dan $(-Z)$ untuk $100.000$ pasangan jalur menghasilkan rata-rata berpasangan $\frac{1}{2}(X_+ + X_-)$ yang memiliki kovariansi negatif kuat. Pada `day001i`, teknik ini menurunkan SE dari `0.017855` menjadi `0.012394` (`vrf = 2.08x`) dengan tambahan komputasi yang minimal.
8. **Control Variate Kemna–Vorst dan Regresi Optimal**: Dengan mengevaluasi opsi Asia aritmetik ($X$) dan geometrik ($Y$) pada jalur acak yang sama, korelasi $>0.9996$ dimanfaatkan melalui estimator $X_{\text{CV}} = X - \hat{\beta}(Y - E[Y])$. Di `outputs/day001i/results.json`, koefisien terestimasi `beta = 1.035177` menghasilkan `vrf = 1329.3x` dan `se = 0.000490`.
9. **Delta Pathwise vs Bump-and-Revalue CRN**: Turunan pathwise $\Delta_{\text{PW}} = e^{-rT}\mathbf{1}\{\bar{S}_A > K\}(\bar{S}_A/S_0)$ dan beda hingga pusat dengan *Common Random Numbers* ($h=1\%$) memberikan konsistensi statistik yang sangat erat: selisih keduanya di `day001i` hanya `+0.000041` (`delta.z_score = 0.0248 < 3.0`).
10. **Validasi Silang Lintas Bahasa Berbasis Benchmark Statistik**: Berbeda dari unit test deterministik murni, validasi implementasi Monte Carlo antar-bahasa (`day001c`–`day001i`) tidak menuntut kesamaan bit karena PRNG tiap bahasa berbeda, melainkan kesamaan deterministik 15 digit pada rumus tertutup (`bs_call`, `y_exact`) dan kompatibilitas statistik dalam 3 SE (`~5.782 ± 0.0005`).

## 10. Kesimpulan dan Batasan
Proyek `day001i` membuktikan bahwa implementasi **Ruby** mereproduksi nilai tertutup Black–Scholes & Kemna–Vorst hingga 15 digit serta menghasilkan estimasi Monte Carlo opsi Asia dan Delta yang konsisten secara statistik dengan `day001c` dalam waktu inti `4328.4 ms`.

## 11. Referensi dan Jejak Bukti
- Kemna, A. G. Z., & Vorst, A. C. F. (1990). *A pricing method for options based on average asset values*. Journal of Banking & Finance, 14(1), 113–129.
- Thomas, D., Fowler, C., & Hunt, A. (2013). *Programming Ruby* & dokumentasi resmi Ruby YJIT.
- Artefak: [projects/ruby/day001i_asian_options/README.md](../../projects/ruby/day001i_asian_options/README.md), [outputs/day001i/results.json](../../outputs/day001i/results.json), [day001i_asian-convergence-loglog.svg](../../outputs/figures/day001i_asian-convergence-loglog.svg), [day001i_asian-estimators-ci.svg](../../outputs/figures/day001i_asian-estimators-ci.svg).
