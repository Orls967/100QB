# Penjelasan Proyek — day002h: Proses Titik Hawkes & Mesin Limit Order Book Berkinerja Tinggi (C++)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002h` (Hari 2, Sesi h)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: C++20 (`cpp`)
- **Folder Proyek**: [projects/cpp/day002h_hawkes_orderbook](../../cpp/day002h_hawkes_orderbook) · **Metadata**: [meta.json](../../cpp/day002h_hawkes_orderbook/meta.json)
- **Ringkasan**: Proses titik self-exciting Hawkes (simulasi Ogata thinning vs eksak Dassios–Zhao, teorema penskalaan waktu Meyer–Papangelou dengan uji KS & Ljung–Box, MLE rekursif $O(n)$ + BFGS $R=200$), mesin *limit order book* (LOB) nol-alokasi heap (>23 juta operasi/detik) vs mesin referensi, serta simulasi mikrostruktur pasar gabungan.

## 2. Masalah dan Tujuan
Kedatangan order di bursa elektronik bersifat mengelompok (*self-exciting*) dan buku order harus memproses jutaan pesan *add/cancel/execute* per detik dengan latensi sub-mikrodetik tanpa alokasi memori dinamis.

## 3. Apa yang Dikerjakan Proyek
Proyek C++20 `projects/cpp/day002h_hawkes_orderbook` mengimplementasikan simulasi dan inferensi MLE $O(n)$ untuk proses Hawkes eksponensial, membangun mesin pencocokan LOB price-time priority menggunakan *pool allocator*, antrean intrusif doubly-linked list, dan bitmap hierarkis 64-bit, memverifikasi kesamaan bit-demi-bit terhadap mesin referensi pada 1.000.000 order, serta mensimulasikan aliran order Hawkes di atas LOB untuk mereproduksi fakta stilistik pasar.

## 4. Uraian Kode
- [projects/cpp/day002h_hawkes_orderbook/include/](../../cpp/day002h_hawkes_orderbook/include/) & `src/`: modul proses Hawkes (Ogata, Dassios–Zhao, kompensator $\Lambda(t_i)$, log-likelihood $O(n)$ + gradien analitik + BFGS), mesin LOB cepat (pool + bitmap 2-level) vs LOB referensi `std::map`, dan simulator pasar gabungan.
- [projects/cpp/day002h_hawkes_orderbook/app/main.cpp](../../cpp/day002h_hawkes_orderbook/app/main.cpp) & `tests/test_hawkes_ob.cpp`: runner benchmark dan suite tes `doctest`.

## 5. Metode dan Alur Kerja
Simulasi Hawkes Ogata vs Dassios–Zhao $\to$ transformasi waktu $\tau_i = \Lambda(t_i)$ & uji KS/Ljung–Box $\to$ estimasi MLE $O(n)$ dengan BFGS pada $R=200$ replikasi paralel $\to$ benchmark 1.000.000 pesan LOB (0 alokasi heap setelah warmup) vs mesin referensi $\to$ analisis *fat tails*, tanda autokorelasi order, dan kecekungan *market impact*.

## 6. Output yang Ditargetkan
Laju stasioner dan faktor Fano cocok dengan teori $\mu/(1-\alpha/\beta)$ dan $1/(1-\alpha/\beta)^2$; uji KS & Ljung–Box menerima parameter benar dan menolak parameter salah; MLE $R=200$ tak-bias dengan cakupan CI ~95%; LOB cepat identik 100% dengan referensi pada 1 juta pesan dengan 0 alokasi heap.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Simulasi Hawkes (Ogata vs Dassios–Zhao) & Teorema Meyer–Papangelou | `tercapai` | `rate_dz = 1.0105` (teori `1.0`), `ks_og_dz_p = 0.106`; pada parameter benar `ks_p = 0.255`, `lb_p = 0.786`; pada parameter salah `ks_p = 9.17e-54` |
| MLE Rekursif $O(n)$ + BFGS ($R=200$ replikasi) | `tercapai` | `mean_mu = 0.50004` (sejati `0.5`, cakupan `95.5%`) |
| Mesin Limit Order Book 1.000.000 Pesan vs Referensi | `tercapai` | Kecocokan state & trade 100%, `0` alokasi heap setelah warmup, throughput tinggi di [outputs/day002h/results.json](../../outputs/day002h/results.json) |
| Verifikasi CTest (`test_hawkes_ob`) | `tercapai` | `cmake` + `ctest` dijalankan pada sesi `day003c`: seluruh tes lulus dalam 8.21s |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `ctest` serta `outputs/day002h/results.json`.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Proses Titik Self-Exciting Hawkes (1971)**: Intensitas bersyarat $\lambda^*(t) = \mu + \sum_{t_i < t} \alpha e^{-\beta(t - t_i)}$ dengan rasio percabangan $n^* = \alpha/\beta < 1$ dan laju stasioner $\bar{\lambda} = \mu/(1 - n^*)$.
2. **Faktor Fano & Pengelompokan Kejadian (Clustering)**: Dalam limit jendela panjang, $\mathrm{Var}(N_T)/E[N_T] \to 1/(1 - \alpha/\beta)^2 > 1$ (overdispersi relatif terhadap Poisson homogen).
3. **Algoritma Ogata Modified Thinning (1981)**: Memanfaatkan batas atas lokal monoton $\lambda^*(t) \le \lambda^*(t_i^+)$ untuk membangkitkan kandidat eksponensial dan menerima dengan peluang $\lambda^*(t^*)/\bar{\lambda}$.
4. **Simulasi Eksak Dekomposisi Lompatan Dassios–Zhao (2013)**: Mengambil sampel langsung waktu tunggu minimum antara kedatangan imigran Poisson homogen $\mu$ dan waktu peluruhan eksponensial non-homogen tanpa penolakan (*rejection-free*).
5. **Rekursi $O(n)$ Log-Likelihood & Kompensator Ozaki (1979)**: Mendefinisikan $A_i = e^{-\beta(t_i - t_{i-1})}(1 + A_{i-1})$ mereduksi biaya evaluasi log-likelihood $\ell(\mu,\alpha,\beta)$ beserta gradien analitiknya dari $O(n^2)$ menjadi $O(n)$.
6. **Teorema Penskalaan Waktu Random (Meyer 1971; Papangelou 1972)**: Interval kompensator $\Delta\Lambda_i = \int_{t_{i-1}}^{t_i}\lambda^*(s)ds$ adalah variabel acak $\mathrm{Exp}(1)$ iid jika dan hanya jika model intensitas benar; diuji dengan Kolmogorov–Smirnov dan Ljung–Box.
7. **Quasi-Newton BFGS dengan Reparameterisasi Tak-Berkendala**: Pemetaan logaritmik/logit menjamin $\mu > 0, \beta > \alpha > 0$ selama pembaruan aproksimasi Hessian invers BFGS.
8. **Arsitektur Limit Order Book Price-Time Priority (FIFO)**: Setiap level harga menyimpan antrean intrusif doubly-linked list dari node yang dialokasikan di *pre-allocated object pool*, menjadikan operasi *cancel by ID* dan *match at best* berwaktu $O(1)$.
9. **Indeks Harga Bitmap Hierarkis 2-Level 64-Bit**: Menggunakan instruksi perangkat keras `std::countr_zero` / `std::countl_zero` pada dua tingkat kata 64-bit untuk menemukan *best bid* dan *best ask* dalam $O(1)$ siklus CPU.
10. **Fakta Stilistik Mikrostruktur Pasar**: Interaksi kedatangan order Hawkes dengan pencocokan LOB secara endogen menghasilkan distribusi return berekor tebal (*excess kurtosis* $> 0$), autokorelasi tanda order positif, dan kurva dampak harga cekung.

## 10. Kesimpulan dan Batasan
Proyek `day002h` menggabungkan ekonometrika proses titik waktu-kontinu dengan rekayasa sistem pencocokan buku order berlatensi rendah.

## 11. Referensi dan Jejak Bukti
- Hawkes, A. G. (1971). *Spectra of some self-exciting and mutually exciting point processes*. Biometrika, 58(1), 83–90.
- Dassios, A., & Zhao, H. (2013). *Exact simulation of Hawkes process with exponentially decaying intensity*. Electronic Communications in Probability, 18, 1–13.
- Artefak: [projects/cpp/day002h_hawkes_orderbook/README.md](../../cpp/day002h_hawkes_orderbook/README.md), [outputs/day002h/results.json](../../outputs/day002h/results.json).
