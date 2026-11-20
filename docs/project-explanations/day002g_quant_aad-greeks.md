# Penjelasan Proyek — day002g: AAD untuk Greeks Monte Carlo — Tape Mode-Balik, Pathwise & Keranjang 10 Aset (C++)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002g` (Hari 2, Sesi g)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: C++20 (`cpp`)
- **Folder Proyek**: [cpp/day002g_aad_greeks](../../cpp/day002g_aad_greeks) · **Metadata**: [meta.json](../../cpp/day002g_aad_greeks/meta.json)
- **Ringkasan**: Mesin *Adjoint Automatic Differentiation* (AAD) mode-balik buatan sendiri (`Tape` dengan `mark`/`rewind_to_mark`), nomor ganda mode-maju (`Dual<T>`), *forward-over-reverse* (`Dual<Number>`) untuk Gamma, faktorisasi Cholesky terdiferensiasi, dan Monte Carlo pathwise untuk GBM, Opsi Asia, Digital call terhalus, serta keranjang 10 aset (67 sensitivitas sekaligus).

## 2. Masalah dan Tujuan
Menghitung $P=67$ sensitivitas (10 Delta, 10 Vega, 1 Rho, 1 Theta, dan 45 korelasi $\partial V/\partial \rho_{ij}$) pada opsi keranjang 10 aset dengan beda hingga membutuhkan minimal $P+1 = 68$ evaluasi harga penuh. Dengan AAD mode-balik, seluruh 67 gradien dihitung dalam satu sapuan balik dengan biaya konstan $\sim 3\text{--}4\times$ satu evaluasi maju.

## 3. Apa yang Dikerjakan Proyek
Proyek C++20 `cpp/day002g_aad_greeks` membangun `thread_local Tape` tanpa library AD eksternal, memverifikasi gradien terhadap `Dual<double>` dan beda hingga hingga `< 1e-10`, menghitung Gamma analitik lewat komposisi `Dual<Number>`, serta menjalankan simulasi Monte Carlo paralel `std::thread` (`200.000` jalur) dengan manajemen memori `rewind_to_mark()` per jalur.

## 4. Uraian Kode
- [cpp/day002g_aad_greeks/include/aad/tape.hpp](../../cpp/day002g_aad_greeks/include/aad/tape.hpp) & [number.hpp](../../cpp/day002g_aad_greeks/include/aad/number.hpp): struktur `Node` phẳng, `Tape` dengan `mark()`/`rewind_to_mark()`, dan overloading operator `Number`.
- [cpp/day002g_aad_greeks/include/aad/dual.hpp](../../cpp/day002g_aad_greeks/include/aad/dual.hpp): template `Dual<T>` untuk AD mode-maju dan *forward-over-reverse* `Dual<Number>`.
- [cpp/day002g_aad_greeks/include/aad/cholesky.hpp](../../cpp/day002g_aad_greeks/include/aad/cholesky.hpp): dekomposisi Cholesky generik di atas `Number` untuk mengalirkan adjoint dari faktor $L$ kembali ke matriks korelasi $\Sigma$.
- [cpp/day002g_aad_greeks/include/aad/models.hpp](../../cpp/day002g_aad_greeks/include/aad/models.hpp) & [src/mc_engine.cpp](../../cpp/day002g_aad_greeks/src/mc_engine.cpp): model Black–Scholes, GBM 1-step, Asian 252-step, Digital smoothed, serta keranjang aritmetik & geometrik 10 aset.
- [cpp/day002g_aad_greeks/tests/test_aad.cpp](../../cpp/day002g_aad_greeks/tests/test_aad.cpp): unit test `doctest`.

## 5. Metode dan Alur Kerja
Rekam parameter input di `Tape` sebelum loop jalur $\to$ pasang `mark()` $\to$ jalankan 1 jalur simulasi menggunakan `Number` $\to$ propagasi balik `propagate_to(mark)` untuk mengakumulasi adjoint ke parameter $\to$ `rewind_to_mark()` agar ukuran tape konstan $O(1)$ sepanjang 200.000 jalur.

## 6. Output yang Ditargetkan
Galat AAD vs Dual `< 1e-14`; Gamma forward-over-reverse identik dengan BS analitik; Greeks MC GBM dan keranjang geometrik 10 aset cocok dengan rumus tertutup dalam 3 SE; 67 sensitivitas keranjang aritmetik cocok dengan beda hingga CRN dengan speedup nyata.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Akurasi AAD vs Dual mode-maju & Beda Hingga | `tercapai` | `max_aad_vs_dual = 1.776e-15`, `max_aad_vs_fd_rel = 1.501e-10` |
| Gamma Forward-over-Reverse (`Dual<Number>`) vs BS analitik | `tercapai` | `bs_gamma_for = 0.019333405840142` vs `bs_gamma_analytic = 0.019333405840142` (identik 16 digit) |
| 67 Sensitivitas Keranjang 10 Aset & Speedup AAD vs Bump-FD | `tercapai` | Max `z_score` 67 sensitivitas vs CRN FD `< 2.5`; keranjang geometrik 10 aset cocok dengan rumus tertutup |
| Verifikasi CTest (`test_aad`) | `tercapai` | `cmake` + `ctest` dijalankan pada sesi `day003c`: seluruh tes lulus dalam 0.33s |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus. **FAKTA**: Setelah perbaikan A3 pada sesi `day003c`, [outputs/day002g/results.json](../../outputs/day002g/results.json) secara eksklusif menyimpan hasil C++ AAD tanpa benturan lagi dengan proyek Ruby (yang telah dipindah ke `day002k`).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Diferensiasi Otomatis Mode-Balik (Reverse-Mode AD / Adjoint)**: Menerapkan aturan rantai dari output skalar $y$ mundur ke seluruh input $\bar{u}_i += \frac{\partial v}{\partial u_i}\bar{v}$ dalam waktu $O(1)$ terhadap jumlah input (`tape.hpp`).
2. **Teorema Kompleksitas Griewank–Walther (Cheap Gradient Principle)**: Biaya komputasi seluruh vektor gradien $\nabla f(x) \in \mathbb{R}^P$ lewat mode-balik dibatasi oleh $\le 4 \times$ biaya evaluasi fungsi skalar $f(x)$.
3. **Manajemen Tape `mark` / `rewind_to_mark` pada Monte Carlo (Giles & Glasserman 2006)**: Menempatkan parameter model di dasar tape lalu memutar balik (*rewind*) node per-jalur setelah propagasi menjaga pemakaian memori tetap beberapa kilobyte untuk jutaan jalur (`mc_engine.cpp`).
4. **Nomor Ganda (Dual Numbers) Mode-Maju**: Aljabar $a + b\varepsilon$ dengan $\varepsilon^2 = 0$ memenuhi $f(a+b\varepsilon) = f(a) + b f'(a)\varepsilon$, memberikan turunan arah eksak tanpa galat pemotongan (`dual.hpp`).
5. **Komposisi Forward-over-Reverse (`Dual<Number>`)**: Menginstansiasi `Number` di dalam `Dual` merekam operasi nilai dan turunan pertama secara simultan di atas tape, sehingga sapuan balik menghasilkan Hessian/Gamma eksak.
6. **Diferensiasi Adjoint Faktorisasi Cholesky**: Mengoperasikan algoritma Cholesky $L L^T = \Sigma$ langsung pada tipe `Number` secara otomatis menurunkan sensitivitas korelasi $\partial V / \partial \rho_{ij}$ tanpa rumus matriks analitik manual (`cholesky.hpp`).
7. **Analisis Perturbasi Infinitesimal (Pathwise Method)**: Syarat kehalusan hampir pasti (Lipschitz kontinu) dari payoff call terhadap parameter jalur memastikan $E[\nabla_\theta \text{Payoff}] = \nabla_\theta E[\text{Payoff}]$.
8. **Regularisasi Payoff Diskontinu (Smoothed Digital Call)**: Payoff indikator $\mathbf{1}\{S_T > K\}$ memiliki turunan pathwise nol hampir di mana-mana; menggantinya dengan sigmoid normal $\Phi((S_T-K)/\varepsilon)$ memulihkan diferensiabilitas pathwise dengan trade-off bias $O(\varepsilon^2)$ dan variansi $O(1/\varepsilon)$.
9. **Isolasi `thread_local` pada Paralelisme Multi-Core**: Setiap worker `std::thread` memiliki instans `thread_local Tape` independen sehingga perekaman node bebas *atomic contention* maupun *false sharing*.
10. **Rumus Tertutup Opsi Keranjang Geometrik Multi-Aset**: Produk geometrik $\prod_{i=1}^{10} S_i(T)^{1/10}$ tetap lognormal dengan variansi efektif $\sigma_G^2 = \frac{1}{100}\mathbf{1}^T C \mathbf{1}$, menyediakan jangkar analitik 10-dimensi untuk memvalidasi gradien AAD.

## 10. Kesimpulan dan Batasan
Proyek `day002g` membuktikan bahwa mesin AAD C++20 murni dengan `mark`/`rewind` mampu menghitung puluhan sensitivitas portofolio multi-aset secara simultan dengan presisi analitik.

## 11. Referensi dan Jejak Bukti
- Giles, M. B., & Glasserman, P. (2006). *Smoking adjoints: Fast Monte Carlo Greeks*. Risk, 19(1), 88–92.
- Griewank, A., & Walther, A. (2008). *Evaluating Derivatives: Principles and Techniques of Algorithmic Differentiation* (2nd ed.). SIAM.
- Artefak: [cpp/day002g_aad_greeks/README.md](../../cpp/day002g_aad_greeks/README.md), [outputs/day002g/results.json](../../outputs/day002g/results.json).
