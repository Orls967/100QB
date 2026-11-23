# Day 005h — Mesin Backtest, Statistik Sharpe, Bias Seleksi, Pengujian Ganda & PBO (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** yang membangun mesin backtest kausal beserta studi statistik dan kuantifikasi *overfitting* pada 300 konfigurasi strategi perdagangan (*SMA Crossover*, *Momentum*, dan *Mean-Reversion Z-Score*). Seluruh data bersifat sintetis dan proyek ini tidak merupakan nasihat investasi.

## Intuisi & Rumus Kunci

1. **Anchor AR(1) Tanpa Look-Ahead**: Untuk pasar AR(1) stasioner $r_t = \phi r_{t-1} + \varepsilon_t$ ($\operatorname{SD}(r_t) = \sigma$, $\phi = 0{,}1$) dan strategi kausal $s_t = \operatorname{sign}(r_{t-1})$, ekspektasi return tanpa biaya adalah:
   $$\mathbb{E}[r_t s_t] = \phi\,\mathbb{E}[|r_{t-1}|] = \phi\,\sigma\sqrt{\frac{2}{\pi}} \approx 0{,}0007979.$$
2. **Metode Delta SE Sharpe & Koreksi Non-Gaussian / Autokorelasi**:
   $$\mathrm{SE}(\widehat{\mathrm{SR}}) = \sqrt{\frac{1 - \hat{\gamma}_3 \widehat{\mathrm{SR}} + \frac{\hat{\gamma}_4 - 1}{4}\widehat{\mathrm{SR}}^2}{T}}, \qquad \eta(q) = \frac{q}{q + 2\sum_{k=1}^{q-1}(q - k)\phi^k}.$$
   Anualisasi naif $\sqrt{q}$ melebih-lebihkan Sharpe saat autokorelasi $\phi > 0$, dan *stationary bootstrap* (Politis–Romano, 1994) memulihkan cakupan interval kepercayaan 95%.
3. **Bias Seleksi, Deflated Sharpe Ratio (DSR) & Pengujian Ganda**:
   $$\mathbb{E}\!\left[\max_{i \le N} Z_i\right] = \int_{-\infty}^{\infty} x N\phi(x)\Phi(x)^{N-1}\,dx \approx (1 - \gamma)\Phi^{-1}\!\left(1 - \frac{1}{N}\right) + \gamma\Phi^{-1}\!\left(1 - \frac{1}{Ne}\right).$$
   Under pure null $H_0$, Benjamini–Hochberg menjamin $\mathrm{FDR} = \pi_0 \alpha$, sedangkan Bonferroni/Holm dan *Romano–Wolf step-down bootstrap* mengontrol $\mathrm{FWER} \le \alpha$.
4. **CSCV / PBO & Validasi Bersih**: Enumerasi eksak seluruh $\binom{16}{8} = 12.870$ pembagian *In-Sample / Out-of-Sample* menghasilkan $\mathrm{PBO} = \mathbb{P}(\lambda \le 0) = 0{,}50$ untuk *noise* simetris yang dapat ditukar, dan *Purged + Embargoed K-Fold* menghilangkan kebocoran korelasi pada label maju tumpang-tindih $h = 10$.

## Instalasi, Lint, Test & Eksekusi

```bash
cd projects/ruby/day005h_backtest_overfitting
bundle install --local
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Struktur Kode

- [lib/backtest_overfitting/rng.rb](lib/backtest_overfitting/rng.rb) & [parallel.rb](lib/backtest_overfitting/parallel.rb) — `SplitMix64`, sampler Normal / Student-$t$ / Geometrik, dan paralelisme deterministik `Process.fork` + `IO.pipe` + `Marshal`.
- [lib/backtest_overfitting/special.rb](lib/backtest_overfitting/special.rb) — `norm_cdf`, `norm_inv` (Halley `< 1e-13`), `binom` exact `Integer`, `student_t_pvalue`, `gauss_legendre`, `pearson_corr`, `ols_fit`, dan `ranks`.
- [lib/backtest_overfitting/data_generator.rb](lib/backtest_overfitting/data_generator.rb), [strategy.rb](lib/backtest_overfitting/strategy.rb), [backtest_engine.rb](lib/backtest_overfitting/backtest_engine.rb), [metrics.rb](lib/backtest_overfitting/metrics.rb) — Generator sintetis (`iid`, `AR(1)`, rezim 2-keadaan, `t_5`, `lognormal`, `GBM`), pustaka 300 strategi kausal (`SMA`, `MOM`, `MRZ`), mesin backtest beserta uji invarian/tanpa *look-ahead*, dan ringkasan metrik risiko.
- [lib/backtest_overfitting/sharpe_stats.rb](lib/backtest_overfitting/sharpe_stats.rb), [selection_bias.rb](lib/backtest_overfitting/selection_bias.rb), [multiple_testing.rb](lib/backtest_overfitting/multiple_testing.rb) — Studi MC $R=20\mathrm{k}$ SE & cakupan Sharpe, faktor Lo $\eta(q)$ & *stationary bootstrap*, $\mathbb{E}[\max]$ & kalibrasi DSR ($R=2\mathrm{k}$), serta Bonferroni / Holm / BH / Romano–Wolf ($R=2\mathrm{k}$).
- [lib/backtest_overfitting/cscv_pbo.rb](lib/backtest_overfitting/cscv_pbo.rb) & [purged_cv.rb](lib/backtest_overfitting/purged_cv.rb) — CSCV eksak $12.870$ pembagian, studi biaya transaksi $c \in \{0, 5, 20\}$ bps, dan *Purged + Embargoed K-Fold* ($h = 10$).
- [lib/backtest_overfitting/runner.rb](lib/backtest_overfitting/runner.rb), [criteria.rb](lib/backtest_overfitting/criteria.rb), [figures.rb](lib/backtest_overfitting/figures.rb), [json_out.rb](lib/backtest_overfitting/json_out.rb), [svg.rb](lib/backtest_overfitting/svg.rb) — Orkestrator, 11 kriteria `c1`–`c11`, penulis JSON dan 6 grafik SVG.
- [test/test_core.rb](test/test_core.rb), [test/test_engine.rb](test/test_engine.rb), [test/test_sharpe_mt.rb](test/test_sharpe_mt.rb), [test/test_cscv_purged.rb](test/test_cscv_purged.rb) — `31` tes unit + integrasi (`13.163` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `20,62s`)

| Bagian | Metrik Utama | Hasil Terukur (`11/11` `tercapai`) |
|---|---|---|
| 1. Anchor & Mesin | Ekuitas, no-lookahead, AR(1), GBM | `eq_err = 0`, `zero_cost_err = 0`, `ar1_z = 0,65`, `gbm_z = 0,55` |
| 2. Statistik Sharpe | Cakupan 95%, Lo $\eta(q)$, Bootstrap | Gauss cov `0,948–0,951`; `AR(1)`: Naif `0,846` $\to$ Stationary Boot `0,914` |
| 3. Bias Seleksi & DSR | $\mathbb{E}[\max]$, $N_{\mathrm{eff}}$, Kalibrasi | Exact vs approx `< 0,018` ($N=300$); Pustaka $\bar{\rho}=0{,}229 \implies N_{\mathrm{eff}}=231{,}5$; DSR pass `0,001` vs Naif `1,000` |
| 4. Pengujian Ganda | FWER, BH FDR $= \pi_0\alpha$, Romano–Wolf | Bonf FWER `0,041`, BH FDR `0,0405` ($\pi_0\alpha = 0{,}045$), Power `0,573` vs `0,259`; RW ($\rho=0{,}4$, `extras`): FWER `0,046`, Power `0,601` > Bonf `0,560` |
| 5. CSCV / PBO & Biaya | $12.870$ pembagian, $c \in \{0, 5, 20\}$ bps | Noise PBO `0,464 +- 0,036` (`0,462` di $N=300$ `extras`), Sinyal PBO `0,000`; `MOM(1)` terbaik `0,112` di `0 bps` runtuh diganti `SMA(10,40)` di `20 bps` |
| 6. Validasi Bersih | Korelasi OOS ($h=10$, fitur noise) | Naive Shuffled CV `0,3741` (**bocor**) vs Purged+Embargo `-0,0037 +- 0,0095` (`z = 0,39`) |

Artefak keluaran: [results.json](../../../outputs/day005h/results.json), 6 grafik SVG `outputs/figures/day005h_*.svg`, dan dokumen penjelasan [day005h_quant_backtest-overfitting.md](../../../docs/project-explanations/day005h_quant_backtest-overfitting.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
