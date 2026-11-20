# Penjelasan Proyek — day002i: Model Hull–White 1-Faktor — Pohon Trinomial, Jamshidian, Kalibrasi & Bermudan LSM (C++)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002i` (Hari 2, Sesi i)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: C++20 (`cpp`)
- **Folder Proyek**: [cpp/day002i_hull_white](../../cpp/day002i_hull_white) · **Metadata**: [meta.json](../../cpp/day002i_hull_white/meta.json)
- **Ringkasan**: Model suku bunga pendek Hull–White 1-faktor di atas kurva awal Nelson–Siegel: rumus analitik obligasi $P(t,T)$ dan opsi ZBP/ZBC/caplet, pohon trinomial dua-tahap Hull–White (1994) dengan induksi maju Arrow–Debreu, dekomposisi Jamshidian untuk swaption Eropa vs pohon vs Monte Carlo OU, kalibrasi $(a,\sigma)$ ($R=100$), swaption Bermudan 10y (pohon vs LSM dua-pass), dan DV01.

## 2. Masalah dan Tujuan
Bagaimana membangun model suku bunga pendek bebas-arbitrase yang cocok persis dengan kurva diskonto awal $P^{M}(0,T)$ hingga presisi mesin, menilai opsi kupon dan swaption Eropa secara analitik lewat trik Jamshidian (1989), serta menilai swaption Bermudan baik dengan pohon trinomial maupun *Least-Squares Monte Carlo* (LSM)?

## 3. Apa yang Dikerjakan Proyek
Proyek C++20 `cpp/day002i_hull_white` mengimplementasikan seluruh pipeline derivatif suku bunga Hull–White dalam C++20 murni tanpa pustaka eksternal, menjalankan verifikasi konvergensi pohon $m \in \{12..192\}$ + ekstrapolasi Richardson, kalibrasi paralel 100 replikasi, dan menulis [outputs/day002i/results.json](../../outputs/day002i/results.json) serta 4 grafik SVG.

## 4. Uraian Kode
- [cpp/day002i_hull_white/include/](../../cpp/day002i_hull_white/include/) & `src/`: kurva Nelson–Siegel, rumus tertutup $A(t,T), B(t,T)$, ZBP/ZBC, solver Brent untuk $r^*$ Jamshidian, pohon trinomial Hull–White dengan harga keadaan Arrow–Debreu $Q_{i,j}$, simulator eksak OU + control variate, kalibrasi Nelder–Mead, dan Bermudan LSM.
- [cpp/day002i_hull_white/app/main.cpp](../../cpp/day002i_hull_white/app/main.cpp) & `tests/test_hw.cpp`: runner utama dan suite tes `doctest`.

## 5. Metode dan Alur Kerja
Konstruksi kurva Nelson–Siegel $\to$ kalibrasi pergeseran deterministik $\alpha_i$ pohon trinomial via induksi maju Arrow–Debreu $\to$ pricing ZBP/ZBC/Cap 5y analitik vs pohon $\to$ dekomposisi Jamshidian swaption Eropa vs pohon & MC $\to$ kalibrasi $(a,\sigma)$ $\to$ induksi mundur Bermudan di pohon vs regresi polinomial LSM dua-pass & sensitivitas DV01.

## 6. Output yang Ditargetkan
Galat kurva diskonto di pohon `< 1.5e-15`; konvergensi pohon + Richardson cocok dengan ZBP/ZBC/Cap dan swaption Jamshidian; MC swaption dalam 3 SE dari Jamshidian; harga Bermudan $\ge \max(\text{Eropa})$ dan LSM konsisten dengan pohon; identitas DV01 bucket vs paralel terpenuhi.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Kalibrasi kurva pohon Arrow–Debreu & Konvergensi ZBP/ZBC/Cap | `tercapai` | Galat kurva `< 1.5e-15` di semua $m$; Richardson(192,96) ZBP = `0.01278125` vs analitik `0.01277295` |
| Swaption Eropa: Dekomposisi Jamshidian vs Pohon vs MC | `tercapai` | Harga pohon dan MC antitetik + control variate cocok dengan rumus tertutup Jamshidian dalam batas toleransi/SE |
| Kalibrasi $(a, \sigma)$ bebas-noise & $R=100$ replikasi bernoise | `tercapai` | Pemulihan parameter eksak pada data bersih dan distribusi tak-bias pada 100 replikasi di `bagian3_calibration` |
| Swaption Bermudan 10y (Pohon vs Two-Pass LSM) & DV01 | `tercapai` | Bermudan melebihi batas bawah maksimum swaption Eropa individual; LSM dua-pass berada dalam SE dari harga pohon; jumlah DV01 bucket cocok dengan pergeseran paralel |
| Verifikasi CTest (`test_hw`) | `tercapai` | `cmake` + `ctest` dijalankan pada sesi `day003c`: seluruh tes lulus dalam 0.73s |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `ctest` serta `outputs/day002i/results.json`.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Model Suku Bunga Pendek Hull–White (1990) / Extended Vasicek**: $dr_t = (\theta(t) - a r_t)dt + \sigma dW_t$ dengan fungsi deterministik $\theta(t)$ yang dipilih agar cocok dengan struktur jangka waktu awal $P^M(0,T)$.
2. **Struktur Term Afin $P(t,T) = A(t,T)\exp(-B(t,T)r_t)$**: Dengan $B(t,T) = (1 - e^{-a(T-t)})/a$, harga obligasi zero-coupon pada waktu mendatang merupakan fungsi eksponensial monoton menurun dari $r_t$.
3. **Parametrisasi Kurva Hasil Nelson–Siegel (1987)**: $y(T) = \beta_0 + \beta_1\frac{1-e^{-T/\tau}}{T/\tau} + \beta_2\left(\frac{1-e^{-T/\tau}}{T/\tau} - e^{-T/\tau}\right)$ memberikan kurva hasil dan *instantaneous forward rate* $f^M(0,T)$ yang mulus.
4. **Rumus Tertutup Opsi Obligasi Zero-Coupon (ZBP/ZBC)**: Di bawah ukuran *T-forward* $\mathbb{Q}^T$, harga obligasi $P(T,S)/P(t,T)$ berdistribusi lognormal dengan volatilitas terintegrasi $\sigma_P = \sigma\sqrt{\frac{1-e^{-2a(T-t)}}{2a}} B(T,S)$.
5. **Dekomposisi Jamshidian (1989) untuk Opsi Obligasi Kupon & Swaption**: Karena setiap harga obligasi $P(T_0, T_k; r)$ monoton terhadap suku bunga pendek tunggal $r_{T_0}$, payoff $(\sum c_k P(T_0,T_k;r) - 1)_+$ terurai secara **eksak** menjadi portofolio opsi obligasi zero-coupon dengan strike $K_k = P(T_0,T_k; r^*)$.
6. **Pohon Trinomial Dua-Tahap Hull–White (1994)**: Tahap 1 membangun kisi trinomial simetris untuk proses OU terpusat $x_t$ (dengan percabangan normal/atas/bawah pada $|j| = j_{\max}$ agar probabilitas tetap di $[0,1]$); Tahap 2 menggeser setiap langkah waktu dengan $\alpha_i$.
7. **Induksi Maju Harga Keadaan Arrow–Debreu $Q_{i,j}$**: Mengakumulasi nilai kini kontingen $Q_{i,j}$ dari akar pohon memungkinkan penentuan $\alpha_i$ secara analitik dalam satu sapuan maju sehingga $P_{\text{tree}}(0, t_i) = P^M(0, t_i)$ hingga presisi mesin (`~1e-15`).
8. **Ukuran T-Forward & Simulasi Eksak Ornstein–Uhlenbeck**: Pergantian numeraire ke $P(t,T_0)$ menyederhanakan ekspektasi diskonto stokastik menjadi ekspektasi di bawah drift forward analitik.
9. **Two-Pass Least-Squares Monte Carlo (Longstaff–Schwartz 2001) untuk Bermudan**: Pass 1 mengestimasi koefisien regresi nilai kelanjutan pada jalur *in-sample*; Pass 2 menerapkan aturan latihan tersebut pada jalur independen *out-of-sample* untuk menghasilkan batas bawah tak-bias.
10. **Sensitivitas Kurva DV01 (Dollar Value of 01) Bucket vs Paralel**: Mengukur perubahan nilai instrumen terhadap perturbasi $+1\text{ bp}$ pada tiap bucket tenor kurva awal dan memverifikasi aditivitas linearnya terhadap pergeseran paralel seluruh kurva.

## 10. Kesimpulan dan Batasan
Proyek `day002i` menyajikan kerangka penilaian derivatif suku bunga Hull–White yang konsisten antara rumus analitik, pohon trinomial Arrow–Debreu, dan simulasi Monte Carlo.

## 11. Referensi dan Jejak Bukti
- Hull, J., & White, A. (1994). *Numerical procedures for implementing term structure models I: Single-factor models*. Journal of Derivatives, 2(1), 7–16.
- Jamshidian, F. (1989). *An exact bond option formula*. Journal of Finance, 44(1), 205–209.
- Artefak: [cpp/day002i_hull_white/README.md](../../cpp/day002i_hull_white/README.md), [outputs/day002i/results.json](../../outputs/day002i/results.json).
