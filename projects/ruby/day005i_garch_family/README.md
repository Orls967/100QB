# Day 005i — Keluarga GARCH, Inferensi QMLE Sandwich, VaR/ES, DCC(1,1), HAR-RV & Salah Spesifikasi (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** yang mengimplementasikan secara end-to-end pemodelan volatilitas bersyarat univariat hingga multivariat dan portofolio 5-aset: `GARCH(1,1)`, `GJR-GARCH(1,1)`, `EGARCH(1,1)`, `EWMA`, `Normal`/`Student-t`/`Hansen Skew-t`, optimasi BFGS ber-skala, kovariansi *Sandwich* Bollerslev–Wooldridge, uji batas $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_1$, backtesting VaR/ES (Kupiec, Christoffersen, Acerbi–Székely $Z_2$, FHS), seleksi QLIKE/Diebold–Mariano/NIC, DCC(1,1) dua-tahap Engle, Corsi HAR-RV intraday, efek semu IGARCH pada *structural break*, dan stres korelasi portofolio 5-aset. Seluruh data bersifat sintetis dan proyek ini bukan nasihat investasi.

## Intuisi & Rumus Kunci

1. **Keluarga Volatilitas Bersyarat**:
   - **GARCH(1,1)**: $h_t = \omega + \alpha \varepsilon_{t-1}^2 + \beta h_{t-1}$, $\mathbb{E}[h_t] = \frac{\omega}{1 - \alpha - \beta}$, kurtosis $\kappa_r = \frac{3(1 - (\alpha+\beta)^2)}{1 - (\alpha+\beta)^2 - 2\alpha^2}$ saat $3\alpha^2 + 2\alpha\beta + \beta^2 < 1$.
   - **GJR-GARCH(1,1)**: $h_t = \omega + (\alpha + \gamma \mathbf{1}_{\{\varepsilon_{t-1} < 0\}})\varepsilon_{t-1}^2 + \beta h_{t-1}$, persistensi $\alpha + \frac{1}{2}\gamma + \beta$.
   - **EGARCH(1,1)**: $\ln h_t = \omega + \alpha(|z_{t-1}| - \sqrt{2/\pi}) + \gamma z_{t-1} + \beta \ln h_{t-1}$, $h_t > 0$ otomatis.
2. **QMLE & Sandwich SE (Bollerslev–Wooldridge, 1992)**:
   Saat ekor return *Student-$t(5)$* ($\kappa_z = 9 \neq 3$) diestimasi memakai *quasi-likelihood* Gaussian, matriks informasi Hessian $A^{-1}$ meremehkan ketidakpastian dan gagal menjaga cakupan 95% (`0,720`), sedangkan kovariansi *Sandwich* $A^{-1} B A^{-1}$ menaikkan rasio SE `1,676x` dan memulihkan cakupan (`0,908`).
3. **Skor Eksak & Backtesting VaR/ES**:
   Model GARCH-$t$ meraih efisiensi MSE `2,47x` atas QMLE pada $\alpha$; uji batas $1/\nu = 0$ maupun $\gamma = 0$ mengikuti campuran $\frac{1}{2}\chi^2_0 + \frac{1}{2}\chi^2_1$ (nilai kritis `2,7055`). Backtesting Acerbi–Székely $Z_2$ memisahkan kegagalan ekor ES GARCH-Normal (`Z2 = -0,949`) dari kalibrasi bersih GARCH-$t$ (`+0,143`) dan FHS (`-0,031`).
4. **Multivariat DCC(1,1), HAR-RV Intraday & Salah Spesifikasi**:
   DCC dua-tahap memulihkan $(a=0,0512,\ b=0,8907)$ (`200/200` konvergen) dan menjaga portofolio terkalibrasi saat korelasi tinggi; HAR-RV memanfaatkan $RV_t$ (*unbiased* terhadap $IV_t$ dengan galat $O(1/M)$) sehingga mengungguli GARCH harian (DM $t=5,51,\ p=3,6\times 10^{-8}$); *break* struktural pada $\omega$ menggelembungkan persistensi GARCH ke `0,9977` yang pulih ke `0,8435` saat dipisah dan terdeteksi `100%` oleh CUSUM Inclán–Tiao.

## Instalasi, Lint, Test & Eksekusi

```bash
cd projects/ruby/day005i_garch_family
bundle install --local
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Struktur Kode (`26` Modul `lib/garch_family/` & `44` Unit/Integrasi Tes)

- [lib/garch_family/rng.rb](lib/garch_family/rng.rb) & [parallel.rb](lib/garch_family/parallel.rb) — `SplitMix64`, Gamma/Student-$t$/Skew-$t$ inversi, serta `Parallel.map_chunks` & `Parallel.map_items` berbasis `Process.fork` + `IO.pipe` + `Marshal`.
- [lib/garch_family/special.rb](lib/garch_family/special.rb) & [linalg.rb](lib/garch_family/linalg.rb) — Beta tak-lengkap, $\chi^2$ & $\frac{1}{2}\chi^2_0+\frac{1}{2}\chi^2_1$, kuadratur Gauss–Legendre, Gauss–Jordan, serta OLS + Newey–West HAC.
- [lib/garch_family/innovations.rb](lib/garch_family/innovations.rb), [vol_model.rb](lib/garch_family/vol_model.rb), [garch11.rb](lib/garch_family/garch11.rb), [gjr_garch11.rb](lib/garch_family/gjr_garch11.rb), [egarch11.rb](lib/garch_family/egarch11.rb), [ewma.rb](lib/garch_family/ewma.rb) — Antarmuka seragam model volatilitas beserta filter, simulasi, transformasi, dan turunan skor analitik rekursif.
- [lib/garch_family/optim.rb](lib/garch_family/optim.rb), [inference.rb](lib/garch_family/inference.rb), [moments_study.rb](lib/garch_family/moments_study.rb), [mle.rb](lib/garch_family/mle.rb), [innovation_study.rb](lib/garch_family/innovation_study.rb) — BFGS berkondisi awal, Hessian/OPG/*Sandwich*, momen orde ke-4, studi $\sqrt{T}$ QMLE ($T=500..2500$), efisiensi $t$-MLE, dan uji batas.
- [lib/garch_family/risk_backtest.rb](lib/garch_family/risk_backtest.rb), [model_comparison.rb](lib/garch_family/model_comparison.rb), [dcc.rb](lib/garch_family/dcc.rb), [dcc_studies.rb](lib/garch_family/dcc_studies.rb), [har_rv.rb](lib/garch_family/har_rv.rb), [misspec_portfolio.rb](lib/garch_family/misspec_portfolio.rb) — Backtest VaR/ES *rolling*, Patton QLIKE & DM, DCC(1,1), Corsi HAR-RV, *break* struktural IGARCH semu, dan stres portofolio 5-aset.
- [lib/garch_family/json_out.rb](lib/garch_family/json_out.rb), [svg.rb](lib/garch_family/svg.rb), [figures.rb](lib/garch_family/figures.rb), [criteria.rb](lib/garch_family/criteria.rb), [runner.rb](lib/garch_family/runner.rb) — Penulis JSON, generator 13 grafik SVG, evaluator kriteria `c1`–`c22`, dan orkestrator `Runner`.
- `test/test_*.rb` (`6` berkas, `44` tes, `258` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `107,92s`)

| Bagian | Metrik Utama | Hasil Terukur (`20/22` `tercapai`, `c13 = sebagian`, `c6 = belum tercapai`) |
|---|---|---|
| 1. Momen & Skor Analitik | Variansi `15` grid, kurtosis, ACF, forecast | `max_var_z = 1,95`, `max_kurt_z = 1,95`, `max_acf_z = 2,35`, `fc_err = 0,57%`, skor FD `< 1e-5` |
| 2. QMLE & Sandwich SE | `RMSE` Normal ($R=1000$), `t(5)` Naif vs Sandwich | `RMSE(alpha) 0,0385 -> 0,0256 -> 0,0154`; `t(5)`: `cov_sw` `0,834 -> 0,881 -> 0,913` ($T=2500$, `c6 = belum tercapai` vs `0,95 +- 3 SE`, namun `extras` $T=10\mathrm{k}$ mencapai `0,9380 +- 0,0076`, `1,58 SE`) |
| 3. Student-$t$ & Uji Batas | $\hat{\nu}$, efisiensi MSE, $\frac{1}{2}\chi^2_0+\frac{1}{2}\chi^2_1$ | `nu_med = 5,048`, `eff(alpha) = 2,47x`; Null `zero = 0,577`, `mix = 0,0433`, kuasa `1,000` |
| 4. Backtest VaR & ES | Rolling 1% ($T=3000$), $Z_2$, MC `R=500` | `garch_n: 1,80%, Z2=-0,949`; `garch_t: 0,87%, Z2=+0,143`; `fhs: 1,03%, Z2=-0,031` |
| 5. Seleksi, DM, NIC & Leverage | AIC/BIC, QLIKE, DM, MC $\gamma=0$ | `BIC_diff = +7,28`, `NIC = 1,668`, `DM p = 0,076` (`sebagian`; `extras T=4k p = 0,0075`); LR `size=0,035`, `pwr=0,98` |
| 6. DCC(1,1) Multivariat | Reduksi CCC, recovery `R=200`, VaR 5% | `ccc_err = 0,0`, `corr = 0,963`; `a = 0,0512`, `b = 0,8907`; port `DCC p_cc = 0,812` vs `Uncorr 2,9e-23` |
| 7. Intraday HAR-RV vs Benchmark | *Infill* $M=13..390$, OOS QLIKE setara info | `max_z = 1,66`, `MSE(13)/MSE(390) = 34,59x`; `HAR = 0,0117` vs `GARCH = 0,0319`, `RV-Naif = 0,0260`, `EWMA-RV = 0,0153`, `GARCH-X = 0,0174` (`DM p <= 9,3e-7`) |
| 8. Break Struktural | IGARCH semu vs split & CUSUM `R=300` | `full = 0,9972 +- 0,0001` (`[0,9917..0,9999]`) vs `split = 0,8277 +- 0,0038` (`med = 0,8383`), CUSUM `100%` (`tau = 1514,9`) |
| 9. Portofolio 5-Aset & Stres | Euler identity, OOS VaR, Stres `2.5x` | `euler_err = 8,67e-19`, `OOS 5% p_uc = 0,0695`, `div_drop = 0,313`, `var1_amp = 3,23x` |

Artefak keluaran: [results.json](../../../outputs/day005i/results.json), `13` grafik SVG `outputs/figures/day005i_*.svg`, dan dokumen penjelasan [day005i_quant_garch-family.md](../../../docs/project-explanations/day005i_quant_garch-family.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan pembelajaran kuantitatif dan tidak merupakan nasihat investasi.
