# Penjelasan Proyek — day002e: Model Heston — Pricing Fourier, Kalibrasi & Monte Carlo QE (Python)

## 1. Identitas dan Ringkasan Proyek
- **ID Progress**: `day002e` (Hari 2, Sesi e)
- **Tanggal**: 2026-11-20 (WITA) · **Topik No**: `bonus` · **Jenis**: `quant` · **Bahasa**: Python
- **Notebook**: [days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb](../../days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb)
- **Ringkasan**: Model volatilitas stokastik Heston (1993) dengan pelanggaran kondisi Feller ($2\kappa\theta = 0.16 < \xi^2 = 0.25$): fungsi karakteristik formasi Albrecher (*little Heston trap*), pricing Fourier Gil–Pelaez vs Fang–Oosterlee COS, kalibrasi & analisis identifiabilitas SVD Jacobian, serta simulasi skema *Quadratic-Exponential* (QE) Andersen (2008) vs Euler *full truncation*.

## 2. Masalah dan Tujuan
Ketika kondisi Feller $2\kappa\theta \ge \xi^2$ dilanggar (kasus umum di pasar ekuitas), proses variansi CIR menyentuh nol sehingga skema Euler naif menjadi bias/negatif, sementara cabang logaritma kompleks pada rumus fungsi karakteristik Heston asli mengalami diskontinuitas cabang (*branch-cut discontinuity*).

## 3. Apa yang Dikerjakan Proyek
Paket `src/qb100/heston/` dan notebook `day002e` mengimplementasikan fungsi karakteristik Heston yang bebas *branch-cut*, membandingkan kuadratur Gil–Pelaez dan deret Fourier-Cosine (COS), mengkalibrasi 5 parameter Heston pada permukaan *implied volatility* beserta analisis nilai singular Jacobian dan studi MC 100 replikasi, serta membandingkan skema simulasi QE Andersen dengan control variate martingale terhadap Euler full-truncation.

## 4. Uraian Kode
- [src/qb100/heston/cf.py](../../src/qb100/heston/cf.py): `HestonParams`, `DEFAULT_PARAMS`, `heston_cf` (formulasi Albrecher et al. 2007), dan `heston_cumulants`.
- [src/qb100/heston/pricing.py](../../src/qb100/heston/pricing.py): `heston_gil_pelaez`, `heston_cos`, `bs_price`, `implied_vol`, dan `market_grid`.
- [src/qb100/heston/calibration.py](../../src/qb100/heston/calibration.py): `calibrate_heston`, `identifiability_from_jacobian`, dan `run_noisy_calibration_mc`.
- [src/qb100/heston/simulation.py](../../src/qb100/heston/simulation.py): `simulate_qe_paths`, `simulate_euler_ft_paths`, `price_from_paths`, `fourier_greeks`, dan `qe_greeks_crn`.
- [tests/test_heston_pricing.py](../../tests/test_heston_pricing.py) & [tests/test_heston_calib_mc.py](../../tests/test_heston_calib_mc.py): 6 suite unit test pytest.

## 5. Metode dan Alur Kerja
Evaluasi $\phi_T(u)$ formasi Albrecher $\to$ inversi Gil–Pelaez & ekspansi COS ($N=256$ suku, rentang pemotongan dari kumulan) $\to$ inversi IV & kalibrasi non-linear least squares (multi-start + SVD Jacobian) $\to$ simulasi variansi CIR dengan pencocokan dua momen QE Andersen (kuadratik Gaussian untuk $\psi \le 1.5$, eksponensial + atom di 0 untuk $\psi > 1.5$) + skema log-spot bebas kebocoran martingale.

## 6. Output yang Ditargetkan
Gil–Pelaez dan COS(256) cocok hingga `< 1e-8`; paritas put-call `< 1e-8`; syarat martingale $\phi(-i) = S_0 e^{rT}$ hingga `< 1e-12`; limit $\xi \to 0$ menuju Black–Scholes; kalibrasi noise-free RMSE `< 1e-10`; harga QE-64 dengan control variate cocok dengan COS dalam 2 SE dan menjaga $v_t \ge 0$ secara eksak.

## 7. Status Pencapaian

| Tujuan / Kriteria | Status | Bukti Spesifik (Artefak & Verifikasi) |
|---|---|---|
| Kesetaraan Gil–Pelaez vs COS(256) & Martingale check | `tercapai` | Output sel 3: `max |GP - COS(256)| = 2.501e-9`, `max |Put-Call Parity| = 6.594e-9`, `max |phi(-i) - S0*e^(rT)| = 9.948e-14` |
| Kalibrasi noise-free, multi-start & 100 replikasi noisy | `tercapai` | Output sel 5: RMSE noise-free `1.626e-11`; 100 replikasi noisy selesai dalam `46.29 s` dengan analisis korelasi $\kappa$–$\xi$ dan kondisi Jacobian |
| Simulasi QE Andersen vs Fourier COS & Euler-FT | `tercapai` | Output sel 7: COS ATM = `8.929410`, QE-64 CV = `8.928157 ± 0.012536` (`VRF = 3.44x`, `min(v) = 0.0`), fraksi cabang $\psi > 1.5 = 3.21\%$ |
| Unit test `test_heston_pricing.py` & `test_heston_calib_mc.py` | `tercapai` | `pytest -q` dijalankan pada sesi `day003c`: 6/6 tes lulus |

## 8. Jika Belum Tercapai atau Belum Diketahui (Fakta vs Hipotesis)
Seluruh kriteria lulus dan terverifikasi lewat `pytest -q` serta output notebook.

## 9. Landasan Teori (Minimal 10 Poin)
1. **SDE Volatilitas Stokastik Heston (1993)**: $dS_t = r S_t dt + \sqrt{v_t} S_t dW_t^S$ dan $dv_t = \kappa(\theta - v_t)dt + \xi\sqrt{v_t}dW_t^v$ dengan korelasi $d\langle W^S, W^v\rangle_t = \rho dt$.
2. **Kondisi Feller untuk Proses CIR**: Titik $v=0$ tidak pernah tersentuh jika dan hanya jika $2\kappa\theta \ge \xi^2$; dengan parameter proyek ($2\kappa\theta = 0.16 < \xi^2 = 0.25$), proses menyentuh nol secara berulang.
3. **Formulasi Little Heston Trap (Albrecher et al. 2007)**: Menulis faktor eksponensial fungsi karakteristik dengan $g = (b - \rho\xi i u - d)/(b - \rho\xi i u + d)$ memastikan $|g e^{-dT}| < 1$ sehingga $\ln((1 - g e^{-dT})/(1-g))$ tidak pernah melintasi *branch cut* negatif (`heston_cf`).
4. **Rumus Inversi Fourier Gil–Pelaez (1951)**: Peluang latihan risiko-netral $P_j = \frac{1}{2} + \frac{1}{\pi}\int_0^\infty \operatorname{Re}\left[\frac{e^{-iu\ln K}\psi_j(u)}{iu}\right]du$ dihitung lewat kuadratur adaptif (`heston_gil_pelaez`).
5. **Metode Fourier-Cosine (COS) Fang–Oosterlee (2008)**: Mengekspansi kepadatan $\ln(S_T/K)$ dalam deret kosinus pada interval $[a,b]$ yang ditentukan oleh kumulan analitik $(c_1, c_2, c_4)$, memberikan konvergensi eksponensial dalam $0.40$ ms (`heston_cos`).
6. **Efek Leverage ($\rho < 0$) & Kemiringan Implied Volatility**: Korelasi negatif $\rho = -0.7$ menghasilkan ekor kiri distribusi return yang tebal (*negative skew*), sedangkan $\xi$ mengendalikan kelengkungan (*vol-of-vol smile*).
7. **Skema Quadratic-Exponential (QE) Andersen (2008)**: Mencocokkan mean dan variansi bersyarat eksak $v_{t+\Delta}$ dari distribusi non-central chi-square: memakai transformasi kuadratik normal $a(b+Z)^2$ saat $\psi = s^2/m^2 \le 1.5$ dan campuran atom di nol plus eksponensial saat $\psi > 1.5$ (`simulate_qe_paths`).
8. **Koreksi Drift Martingale pada Diskretisasi Log-Spot**: Mengganti aproksimasi integral $\int v_s ds$ dengan konstanta $K_0^*$ yang mengkompensasi fungsi pembangkit momen skema QE menjaga $E[S_T] = S_0 e^{rT}$ secara eksak.
9. **Skema Euler Full-Truncation (Lord et al. 2010)**: Memotong variansi menjadi $v_t^+ = \max(v_t, 0)$ di dalam drift dan difusi mengurangi bias dibanding refleksi, namun tetap memiliki bias langkah waktu lebih besar daripada QE (`simulate_euler_ft_paths`).
10. **Identifiabilitas Lokal & Dekomposisi SVD Jacobian**: Nilai singular terkecil dari matriks Jacobian $\partial \sigma_{\text{imp}} / \partial \theta$ mengungkap kolinearitas struktural antara kecepatan *mean-reversion* $\kappa$ dan *vol-of-vol* $\xi$ pada maturitas terbatas (`identifiability_from_jacobian`).

## 10. Kesimpulan dan Batasan
Proyek `day002e` memberikan implementasi referensi Heston yang lengkap dari pricing Fourier sub-milidetik hingga simulasi QE Andersen pada rezim Feller terlanggar.

## 11. Referensi dan Jejak Bukti
- Heston, S. L. (1993). *A closed-form solution for options with stochastic volatility*. Review of Financial Studies, 6(2), 327–343.
- Fang, F., & Oosterlee, C. W. (2008). *A novel pricing method for European options based on Fourier-cosine series expansions*. SIAM Journal on Scientific Computing, 31(2), 826–848.
- Andersen, L. (2008). *Simple and efficient simulation of the Heston stochastic volatility model*. Journal of Computational Finance, 11(3), 1–42.
- Artefak: [day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb](../../days/day002e_quant_model-heston-pricing-fourier-kalibrasi.ipynb), [src/qb100/heston/](../../src/qb100/heston/__init__.py).
