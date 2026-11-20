# Penjelasan Proyek — day001c: Opsi Asia — Monte Carlo, Variance Reduction & Greeks (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day001c` (Hari 1, Sesi c)
- **Tanggal**: 2026-11-19 (WITA) · **Topik No**: `bonus` (off-roadmap) · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb](../../days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb)
- **Ringkasan**: Penilaian *arithmetic Asian call option* ($S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000$) dengan Naive MC, Antithetic Variates, Control Variate Kemna–Vorst, dan estimasi Delta (Pathwise vs Bump-and-Revalue CRN).

## 2. Masalah dan Tujuan
Rata-rata aritmetik dari variabel lognormal tidak berdistribusi lognormal sehingga opsi Asia aritmetik tidak memiliki rumus tertutup. Selain itu, Monte Carlo naif membutuhkan jutaan jalur untuk mencapai presisi di bawah 1 sen. Proyek ini menerapkan *antithetic variates* dan *control variate* berbasis opsi Asia geometrik tertutup (Kemna–Vorst 1990) untuk menekan variansi lebih dari $1300\times$.

## 3. Apa yang Dikerjakan Proyek
Modul `src/qb100/options.py` dan notebook `day001c` mensimulasikan jalur GBM eksak 252 langkah, memvalidasi mesin MC terhadap Black–Scholes Eropa (`10.450583572185565`), membandingkan tiga estimator harga opsi Asia aritmetik dan dua estimator Delta, serta memetakan konvergensi SE pada skala log-log.

## 4. Uraian Kode
- [src/qb100/options.py](../../src/qb100/options.py): fungsi `bs_call_price`, `geometric_asian_call_price`, `simulate_gbm_paths`, `price_asian_naive`, `price_asian_antithetic`, `price_asian_control_variate`, `asian_delta_pathwise`, dan `asian_delta_bump`.
- [days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb](../../days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb): 14 sel notebook eksperimen dan visualisasi.
- [tests/test_options.py](../../tests/test_options.py): 5 unit test pytest untuk rumus tertutup BS/Kemna–Vorst, reduksi variansi, dan konsistensi Delta.

## 5. Metode dan Alur Kerja
Parameter tetap $(S_0=100, K=100, r=0.05, \sigma=0.2, T=1, M=252, N=200.000, \text{seed}=42)$ $\to$ simulasi eksak $S_{t_k} = S_0 \exp(\sum_{j=1}^k [(r-\frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t} Z_j])$ $\to$ hitung payoff aritmetik $X = e^{-rT}(\bar{S}_A - K)_+$ dan geometrik $Y = e^{-rT}(\bar{S}_G - K)_+$ $\to$ koreksi control variate $X_{\text{CV}} = X - \hat{\beta}(Y - E[Y])$ $\to$ hitung Delta pathwise dan finite-difference CRN.

## 6. Output yang Ditargetkan
Anchor BS Eropa = `10.450584`, Kemna–Vorst geometrik = `5.565509`; MC Eropa dalam 3 SE dari BS; VRF Antithetic > 1.5x dan Control Variate > 500x; Delta pathwise dan bump-and-revalue cocok dalam 3 SE; slope log-log SE $\approx -0.50$.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Sanity check MC Eropa vs Black–Scholes | `tercapai` | Output sel 5: BS = `10.450584`, MC = `10.504703 ± 0.033110` (selisih `+0.054119`, `1.63 SE < 3 SE`) |
| Harga Asian Call & Variance Reduction Factor (VRF) | `tercapai` | Output sel 6: Naive = `5.804757 ± 0.017935`, Antithetic = `5.791835 ± 0.012451` (VRF `2.08x`), Control Variate = `5.782468 ± 0.000495` (VRF `1313x`, $\hat{\beta}=1.035673$, $E[Y]=5.565509$) |
| Delta Pathwise vs Bump-and-Revalue (CRN) | `tercapai` | Output sel 7: Pathwise = `0.589628 ± 0.001179`, Bump = `0.589588 ± 0.001158` (selisih `+0.000040`, `0.02 SE`) |
| Konvergensi log-log SE vs $N$ | `tercapai` | Output sel 8: slope Naive `-0.5087`, Antithetic `-0.5022`, Control Variate `-0.4991`; grafik [day001c_asian-se-convergence.png](../../outputs/figures/day001c_asian-se-convergence.png) & [day001c_asian-ci-comparison.png](../../outputs/figures/day001c_asian-ci-comparison.png) |
| Unit test `tests/test_options.py` | `tercapai` | `pytest -q` dijalankan di sesi ini: 5/5 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus. Implementasi Python `day001c` ini juga menjadi rujukan validasi silang untuk enam implementasi bahasa lain (`day001d`–`day001i`).

## 9. Landasan Teori (Minimal 10 Poin)
1. **Gerak Brown Geometrik (GBM) Diskretisasi Eksak**: Solusi SDE $dS_t = r S_t dt + \sigma S_t dW_t$ di bawah ukuran risiko-netral $\mathbb{Q}$ adalah $S_{t_k} = S_{t_{k-1}}\exp((r-\frac{1}{2}\sigma^2)\Delta t + \sigma\sqrt{\Delta t}Z_k)$, bebas bias diskretisasi waktu pada `simulate_gbm_paths`.
2. **Penilaian Risiko-Netral (Risk-Neutral Valuation)**: Harga opsi adalah ekspektasi terdiskonto $V_0 = e^{-rT} E^{\mathbb{Q}}[\text{Payoff}]$ di bawah ukuran martingale ekuivalen.
3. **Rumus Black–Scholes Eropa (Anchor)**: $C_{\text{BS}} = S_0 N(d_1) - K e^{-rT} N(d_2) = 10.450583572185565$ pada `bs_call_price`, menjadi jangkar validasi generator jalur GBM.
4. **Opsi Asia Aritmetik vs Geometrik**: Rata-rata aritmetik $\bar{S}_A = \frac{1}{M}\sum_{k=1}^M S_{t_k} \ge \bar{S}_G = (\prod_{k=1}^M S_{t_k})^{1/M}$ (ketaksamaan AM–GM), sehingga harga call aritmetik (`~5.7825`) selalu lebih tinggi dari geometrik (`5.5655`).
5. **Rumus Tertutup Kemna–Vorst (1990) untuk Opsi Asia Geometrik**: Karena $\ln \bar{S}_G$ adalah kombinasi linear variabel normal, $\bar{S}_G$ tetap lognormal dengan $\sigma_G = \sigma \sqrt{\frac{(M+1)(2M+1)}{6M^2}}$ dan drift efektif $\mu_G = (r - \frac{1}{2}\sigma^2)\frac{M+1}{2M} + \frac{1}{2}\sigma_G^2$, memberikan $E[Y] = 5.565508831349356$ pada `geometric_asian_call_price`.
6. **Antithetic Variates**: Memasangkan jalur $(Z_1,\dots,Z_M)$ dan $(-Z_1,\dots,-Z_M)$ menghasilkan kovariansi negatif karena payoff monoton terhadap shock, menggandakan efisiensi (`VRF ~ 2.08x`) pada `price_asian_antithetic`.
7. **Control Variates & Koefisien Regresi Optimal**: Estimator $X_{\text{CV}} = X - \beta^*(Y - E[Y])$ dengan $\beta^* = \mathrm{Cov}(X,Y)/\mathrm{Var}(Y)$ mengurangi variansi sebesar faktor $1/(1 - \rho_{XY}^2)$; karena $\rho_{XY} > 0.9996$, VRF mencapai `~1313x` pada `price_asian_control_variate`.
8. **Pathwise Derivative (IPA — Infinitesimal Perturbation Analysis)**: Untuk payoff Lipschitz $(\bar{S}_A - K)_+$, turunan dapat dimasukkan ke dalam ekspektasi: $\Delta_{\text{PW}} = e^{-rT} \mathbf{1}\{\bar{S}_A > K\} \frac{\bar{S}_A}{S_0}$ karena $\partial S_{t_k}/\partial S_0 = S_{t_k}/S_0$ (`asian_delta_pathwise`).
9. **Finite-Difference Delta dengan Common Random Numbers (CRN)**: Menggunakan matriks bilangan acak $Z$ yang sama pada $S_0(1+h)$ dan $S_0(1-h)$ membuat variansi selisih menjadi $O(h^2)$, menjaga SE tetap kecil (`0.001158`) pada `asian_delta_bump`.
10. **Variance Reduction Factor (VRF)**: Rasio $\mathrm{Var}(\hat{\theta}_{\text{naive}})/\mathrm{Var}(\hat{\theta}_{\text{VR}})$ mengukur penghematan jumlah jalur simulasi yang dibutuhkan untuk mencapai lebar CI 95% yang sama.

## 10. Kesimpulan dan Batasan
Proyek `day001c` membuktikan efektivitas control variate analitik Kemna–Vorst (VRF > 1300x) dan kesetaraan estimator Delta pathwise dengan beda hingga CRN.

## 11. Referensi dan Jejak Bukti
- Kemna, A. G. Z., & Vorst, A. C. F. (1990). *A pricing method for options based on average asset values*. Journal of Banking & Finance, 14(1), 113–129.
- Broadie, M., & Glasserman, P. (1996). *Estimating security price derivatives using simulation*. Management Science, 42(2), 269–285.
- Artefak: [day001c_quant_opsi-asia-monte-carlo-variance.ipynb](../../days/day001c_quant_opsi-asia-monte-carlo-variance.ipynb), [src/qb100/options.py](../../src/qb100/options.py).
