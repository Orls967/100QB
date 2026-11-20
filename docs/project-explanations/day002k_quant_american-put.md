# Penjelasan Proyek — day002k: Opsi Put Amerika — CRR Binomial, LSM Monte Carlo & PSOR PDE (Ruby)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002k` (Hari 2, Sesi k; direlabel dari `2g` pada sesi `day003c` A3 untuk menghilangkan benturan ID dengan C++ AAD)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Ruby (`ruby`)
- **Folder Proyek**: [projects/ruby/day002k_american_put](../../ruby/day002k_american_put) · **Metadata**: [meta.json](../../ruby/day002k_american_put/meta.json)
- **Ringkasan**: Tiga metode independen untuk menilai *American Put Option* ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1$) dalam Ruby murni tanpa gem runtime: pohon binomial Cox–Ross–Rubinstein ($N=4.000$), *Least-Squares Monte Carlo* Longstaff–Schwartz ($100.000$ jalur, basis Laguerre/polinomial derajat 3), dan PDE Crank–Nicolson dengan *Projected Successive Over-Relaxation* (PSOR $400\times 400$).

## 2. Masalah dan Tujuan
Opsi Put Amerika memberikan hak latihan lebih awal kapan saja sebelum maturitas ($V(S,t) \ge (K-S)_+$), menciptakan masalah batas bebas (*free-boundary / optimal stopping problem*) yang tidak memiliki solusi analitik tertutup.

## 3. Apa yang Dikerjakan Proyek
Proyek Ruby `projects/ruby/day002k_american_put` menghitung harga, *early exercise premium* (EEP), dan Greeks opsi put Amerika dengan tiga paradigma numerik berbeda (pohon diskrit CRR, simulasi mundur LSM paralel `Process.fork`, dan masalah komplementaritas linear PDE dengan solver iteratif PSOR), lalu memvalidasi silang ketiganya di [outputs/day002k/results.json](../../outputs/day002k/results.json).

## 4. Uraian Kode
- [projects/ruby/day002k_american_put/lib/american_put/crr.rb](../../ruby/day002k_american_put/lib/american_put/crr.rb): pohon binomial CRR $N=4.000$ beserta ekstraksi Delta, Gamma, dan Theta dari simpul awal.
- [projects/ruby/day002k_american_put/lib/american_put/lsm.rb](../../ruby/day002k_american_put/lib/american_put/lsm.rb) & [lsm_simulator.rb](../../ruby/day002k_american_put/lib/american_put/lsm_simulator.rb): simulasi jalur antitetik dan regresi nilai kelanjutan Longstaff–Schwartz pada jalur *in-the-money*.
- [projects/ruby/day002k_american_put/lib/american_put/psor.rb](../../ruby/day002k_american_put/lib/american_put/psor.rb) & [psor_grid.rb](../../ruby/day002k_american_put/lib/american_put/psor_grid.rb): diskretisasi PDE Crank–Nicolson dan iterasi PSOR dengan proyeksi $\max(V_i^{(k)}, g_i)$.
- [projects/ruby/day002k_american_put/lib/american_put/runner.rb](../../ruby/day002k_american_put/lib/american_put/runner.rb) & [test/test_american_put.rb](../../ruby/day002k_american_put/test/test_american_put.rb): runner benchmark dan suite tes `minitest`.

## 5. Metode dan Alur Kerja
Evaluasi acuan Black–Scholes Put Eropa (`5.57352602`) $\to$ induksi mundur CRR ($N=4.000$) $\to$ penyelesaian LCP Crank–Nicolson dengan PSOR ($400\times 400$, $\omega=1.3$) $\to$ simulasi LSM ($100.000$ jalur, 50 langkah latihan) $\to$ uji konsistensi silang ketiga metode dan pemetaan profil harga terhadap spot $S_0 \in [70, 130]$.

## 6. Output yang Ditargetkan
Selisih harga CRR($4000$) vs PSOR($400\times 400$) `< 5e-3`; harga LSM berada dalam 3 SE dari CRR dan PSOR; harga Amerika selalu $\ge$ harga Eropa dan $\ge (K-S_0)_+$; Delta dan Gamma CRR vs PSOR cocok hingga `< 1e-4`.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Konsistensi CRR ($N=4000$) vs PSOR PDE ($400\times 400$) | `tercapai` | CRR = `6.090181` (Delta `-0.411073`, Gamma `0.022992`), PSOR = `6.087551` (Delta `-0.411052`, Gamma `0.022996`), selisih harga = `2.629637e-03` |
| Konsistensi LSM Monte Carlo vs CRR & PSOR (`z < 3`) | `tercapai` | LSM = `6.051915 ± 0.013566`, `lsm_vs_crr_z = 2.821 < 3.0`, `lsm_vs_psor_z = 2.627 < 3.0` |
| Early Exercise Premium (EEP) positif & batas bebas | `tercapai` | BS Eropa = `5.573526`, EEP CRR = `0.517155 > 0` |
| Verifikasi Minitest, RuboCop & Reproduksi Riwayat Git (A3) | `tercapai` | `bundle exec rake test` (6 runs, 44 assertions, 0 failures), `rubocop` (19 files, no offenses), dan `bin/run` dijalankan pada sesi `day003c`: 195/203 field numerik non-timing **identik bit-demi-bit** dengan commit `31af7ca` |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus. **FAKTA**: Pada sesi `day003c` (Bagian A3), proyek ini dipindahkan dari `day002g` ke `day002k` dan dijalankan ulang penuh; seluruh harga, SE, Greeks, dan kedua grafik SVG identik dengan riwayat git sebelum tertimpa.

## 9. Landasan Teori (Minimal 10 Poin)
1. **Masalah Penghentian Optimal (Optimal Stopping) Snell Envelope**: $V_0 = \sup_{\tau \in \mathcal{T}_{[0,T]}} E^{\mathbb{Q}}[e^{-r\tau}(K - S_\tau)_+]$ di mana batas latihan optimal $S^*(t)$ memisahkan daerah kelanjutan dari daerah latihan seketika.
2. **Ketergantungan Early Exercise Put pada Suku Bunga $r > 0$**: Untuk put tanpa dividen, nilai waktu bunga $rK$ yang diterima saat latihan dini dapat melebihi nilai asuransi volatilitas ketika $S_t$ cukup rendah ($S_t \le S^*(t)$).
3. **Pohon Binomial Cox–Ross–Rubinstein (1979) dengan Proyeksi Dinamis**: $u = e^{\sigma\sqrt{\Delta t}}, d = 1/u, p = (e^{r\Delta t}-d)/(u-d)$ dengan pembaruan mundur $V_{i,j} = \max((K - S_{i,j})_+, e^{-r\Delta t}(p V_{i+1,j+1} + (1-p)V_{i+1,j}))$ (`crr.rb`).
4. **Ekstraksi Greeks langsung dari Lapisan Awal Pohon Binomial**: Delta dihitung dari selisih simpul $t_1$, Gamma dari tiga simpul $t_2$, dan Theta dari simpul tengah $V_{2,1} - V_{0,0}$ tanpa re-pricing pohon.
5. **Algoritma Longstaff–Schwartz (2001) Least-Squares Monte Carlo**: Mengaproksimasi ekspektasi nilai kelanjutan $E[e^{-r\Delta t}V_{t+1} | S_t]$ melalui regresi kuadrat terkecil terhadap basis polinomial hanya pada jalur yang sedang *in-the-money* (`lsm.rb`).
6. **Sifat Batas Bawah (Low-Bias) Aturan Latihan Sub-Optimal LSM**: Karena aturan penghentian hasil regresi merupakan strategi latihan yang valid (namun sub-optimal terhadap batas sejati) dan dievaluasi dengan langkah diskrit $M=50$, harga LSM (`6.0519`) berada sedikit di bawah nilai waktu-kontinu (`6.0902`).
7. **Formulasi Masalah Komplementaritas Linear (LCP) PDE Amerika**: Pada setiap langkah waktu Crank–Nicolson, vektor solusi $v$ harus memenuhi $v \ge g$, $(A v - b) \ge 0$, dan $(v - g)^T (A v - b) = 0$ dengan $g_i = (K - S_i)_+$ (`psor.rb`).
8. **Algoritma Projected Successive Over-Relaxation (PSOR — Cryer 1971)**: Memodifikasi pembaruan Gauss–Seidel tridiagonal dengan relaksasi $\omega \in (1,2)$ dan proyeksi langsung $v_i^{(k+1)} = \max\big(g_i, (1-\omega)v_i^{(k)} + \frac{\omega}{A_{ii}}(b_i - \dots)\big)$, konvergen dalam rata-rata `~15.5` sapuan per langkah waktu.
9. **Kondisi Smooth Pasting di Batas Latihan $S^*(t)$**: Solusi nilai put Amerika menempel secara tangensial pada garis payoff intrinsik: $V(S^*(t), t) = K - S^*(t)$ dan $\frac{\partial V}{\partial S}(S^*(t), t) = -1$.
10. **Dekomposisi Early Exercise Premium (EEP)**: $V_{\text{Am}}(S_0) = V_{\text{Eu}}(S_0) + \text{EEP}(S_0)$ memungkinkan verifikasi batas bawah $V_{\text{Am}} \ge V_{\text{Eu}}$ di seluruh profil spot $S_0 \in [70, 130]$.

## 10. Kesimpulan dan Batasan
Proyek `day002k` memvalidasi kesesuaian tiga metode numerik utama penilaian opsi Amerika dalam Ruby murni.

## 11. Referensi dan Jejak Bukti
- Longstaff, F. A., & Schwartz, E. S. (2001). *Valuing American options by simulation: A simple least-squares approach*. Review of Financial Studies, 14(1), 113–147.
- Wilmott, P., Dewynne, J., & Howison, S. (1993). *Option Pricing: Mathematical Models and Computation*. Oxford Financial Press.
- Artefak: [projects/ruby/day002k_american_put/README.md](../../ruby/day002k_american_put/README.md), [outputs/day002k/results.json](../../outputs/day002k/results.json).
