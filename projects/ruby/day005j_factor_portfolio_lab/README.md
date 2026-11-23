# Day 5j — Laboratorium Model Faktor & Optimisasi Portofolio (Ruby murni)

Proyek Ruby murni (`>= 3.3`, nol gem numerik) yang mengimplementasikan rantai penuh *quantitative equity & multi-asset portfolio management* dalam 12 bagian terstruktur (`36` modul di `lib/factor_portfolio_lab/`, `83` unit/integrasi tes dengan `712` asersi, `42` kriteria verifikasi kuantitatif `c1`–`c42`, dan `26` grafik SVG di `outputs/figures/day005j_*.svg`).

## Cakupan 12 Bagian

1. **Factor World & Spiked PCA (`factor_world.rb`, `pca.rb`, `spiked_study.rb`)**: transisi fase Baik–Ben Arous–Péché (BBP) dan tepi Marchenko–Pastur ($\lambda_\pm = \sigma^2(1 \pm \sqrt{\gamma})^2$, $R = 500$).
2. **Regresi Deret Waktu & Uji GRS (`ols.rb`, `time_series_regression.rb`, `grs.rb`, `grs_study.rb`)**: uji Gibbons–Ross–Shanken $F(N, T-N-K)$ eksak vs Wald asimtotik $\chi^2(N)$ ($R = 5000$).
3. **Fama–MacBeth & Koreksi Shanken (`fama_macbeth.rb`, `shanken.rb`, `fm_study.rb`)**: estimasi dua-tahap premi faktor dan koreksi *errors-in-variables* Shanken ($R = 1000$).
4. **Risk Parity (`risk_parity.rb`, `rp_study.rb`)**: *Equal Risk Contribution* melalui Newton Spinu (konvergensi 6 iterasi, deviasi `< 1e-16`) vs *Cyclic Coordinate Descent*.
5. **Black–Litterman (`black_litterman.rb`, `bl_study.rb`)**: *implied equilibrium* $\pi = \delta \Sigma w_{\mathrm{mkt}}$, kalibrasi He–Litterman & Idzorek, dan verifikasi paritas Bayes *importance sampling* $R = 500{.}000$.
6. **Markowitz QP Terkendala (`constraints.rb`, `markowitz.rb`, `mv_study.rb`)**: proyeksi *capped simplex*, *Projected Barzilai–Borwein + active-set KKT polish* (residual `< 1e-11`, `45` titik frontier).
7. **Michaud Resampling & Shrinkage (`shrinkage.rb`, `resampling.rb`)**: *Ledoit–Wolf* analitik, *James–Stein* mean, dan frontier *bootstrap* Michaud ($R = 500$).
8. **Backtest Rolling Dua-Rezim (`strategy.rb`, `metrics.rb`, `backtest.rb`, `backtest_study.rb`)**: `8` strategi seragam pada $R = 300$ dunia dua-rezim tanpa *look-ahead* dan biaya transaksi $c \in \{0, 5, 20\}$ bps.
9. **Atribusi Risiko Euler & Salah Spesifikasi $K$ (`attribution.rb`)**: dekomposisi Euler eksak (`err = 2.08e-17`) dan kebocoran variansi sistematis ke idio saat $K_{\mathrm{fit}} < K_{\mathrm{true}}$.
10. **Stress Test & Reverse Stress Test (`stress.rb`)**: skenario faktor linear + gamma kuadratik dan *closed-form Mahalanobis reverse stress test* di metrik $\Sigma_f^{-1}$.
11. **Ketahanan Ekor Berat Student-$t(4)$ (`robustness.rb`)**: dampak *tail risk* terhadap kalibrasi GRS dan CVaR 95% portofolio.
12. **Skalabilitas Komputasi (`complexity.rb`)**: profil waktu $N \in \{10, 25, 50, 100, 200\}$ dan eksponen log-log empiris.

## Menjalankan Proyek

```bash
cd projects/ruby/day005j_factor_portfolio_lab
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Ringkasan Hasil Utama (`MAIN_SEED = 20261123`)

- **Waktu eksekusi penuh (`bin/run`, 10 pekerja YJIT ON)**: `364.21s` (`39/42` kriteria `tercapai`, `26` grafik SVG).
- **BBP Spiked PCA (`R=500`, $\gamma=0{,}25$)**: pada $\ell=0{,}3 < \sqrt{\gamma}=0{,}5$, $\hat{\lambda}_1 = 2{,}168$ (teori `2,250`) dan alignment `0,096`; pada $\ell=2{,}0$, $\hat{\lambda}_1 = 3{,}346$ (teori `3,375`) dan alignment `0,832` (teori `0,833`).
- **GRS vs Wald (`R=5000`, $N=50, T=120$)**: penolakan GRS $F$ `0,0522 +- 0,0031` vs Wald $\chi^2$ `0,8482` (*over-rejection* masif).
- **Black–Litterman (`R=500.000`)**: `max_mu_diff = 5.41e-5`, `max_cov_diff = 9.57e-7`, `ESS = 288,497`.
- **Markowitz QP**: `max_kkt_long_only = 4.44e-16`, `max_kkt_capped = 9.98e-12`.
- **Kriteria `belum tercapai` (jujur, 3/42)**: `c30` (*James–Stein* L1 `1,7023` menyamai *plug-in* `1,7023` pada $\kappa=0{,}5, T=120$), `c32` (`MinVarShrink` ann_vol `0,26368` sedikit di atas `MinVarSample` `0,26358` pada dunia dua-rezim), dan `c33` (`RiskParitySample` Sharpe `-0,0375` vs `EqualWeight` `-0,0266` saat rezim *bear* ber-drift negatif).

Penjelasan teknis lengkap tersedia di [day005j_quant_factor-portfolio-lab.md](../../../docs/project-explanations/day005j_quant_factor-portfolio-lab.md).
