# Day 005f — Valuasi Tranche CDO: LHP Vasicek, Rekursi ASB, Base Correlation & Exponential Twisting (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** menggunakan **topik cadangan** (*valuasi tranche CDO: LHP Vasicek vs rekursi Andersen–Sidenius–Basu vs MC dengan exponential twisting*) setelah topik utama (*Dupire local vol / SVI*) terdeteksi telah dikerjakan penuh pada [`day003a`](../../../docs/project-explanations/day003a_quant_volatilitas-lokal-dupire-ssvi.md) dan model Heston pada [`day002e`](../../../docs/project-explanations/day002e_quant_model-heston-pricing-fourier-kalibrasi.md). Seluruh data bersifat sintetis (125 nama, tenor 5Y, kupon kuartalan, $r = 0{,}02$).

## Intuisi & Rumus Kunci

1. **Rekursi Andersen–Sidenius–Basu (ASB, 2003)**: Bersyarat pada faktor umum $M \sim \mathcal{N}(0, 1)$, kegagalan tiap obligor independen dengan peluang $p_i(t \mid M) = \Phi\!\left(\frac{\Phi^{-1}(p_i(t)) - \sqrt{\rho}\,M}{\sqrt{1-\rho}}\right)$. Distribusi kerugian portofolio bersyarat dihitung eksak via konvolusi DP pada grid satuan kerugian lalu diintegralkan atas $M$ dengan kuadratur Gauss–Legendre 64-titik.
2. **Large Homogeneous Pool (LHP Vasicek)**: Pada batas $N \to \infty$ homogen, CDF kerugian memiliki bentuk tertutup $F_{\mathrm{LHP}}(x) = \Phi\!\left(\frac{\sqrt{1-\rho}\,\Phi^{-1}(x/(1-R)) - \Phi^{-1}(p)}{\sqrt{\rho}}\right)$, dan $\mathbb{E}[\min(L, D)]$ dapat dievaluasi via ortan bivariat Gaussian.
3. **Fair Spread & Konservasi Tranche**: $\mathrm{PV}_{\mathrm{prot}} = \sum_m DF(t_m)(\mathrm{EL}(t_m) - \mathrm{EL}(t_{m-1}))$, $\mathrm{PV}_{\mathrm{prem}} = \sum_m \Delta t_m DF(t_m)((D-A) - \frac{1}{2}(\mathrm{EL}(t_m)+\mathrm{EL}(t_{m-1})))$, dengan $\sum_k \mathrm{PV}_{\mathrm{prot},k} = \mathrm{PV}_{\mathrm{prot,port}}$.
4. **Exponential Twisting (Glasserman–Li, 2005)**: Untuk *tranche senior* yang jarang tembus, pergeseran faktor $M \sim \mathcal{N}(\mu^*, 1)$ digabung dengan *exponential twisting* Bernoulli $q_i(\theta, M) = \frac{p_i e^{\theta \ell}}{1 - p_i + p_i e^{\theta \ell}}$ menargetkan ambang *attachment* $A$.

## Instalasi, Lint, Test & Eksekusi

```bash
cd projects/ruby/day005f_cdo_tranche
bundle install --local
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Struktur Kode

- [lib/cdo_tranche/rng.rb](lib/cdo_tranche/rng.rb) & [parallel.rb](lib/cdo_tranche/parallel.rb) — `SplitMix64`, sampler normal/Gamma/$\chi^2$, dan `Process.fork` + `IO.pipe` + `Marshal` deterministik.
- [lib/cdo_tranche/special.rb](lib/cdo_tranche/special.rb) & [models.rb](lib/cdo_tranche/models.rb) — CDF/kuantil normal (AS241) & $t$, Gauss–Legendre, ortan bivariat Plackett, portofolio 125-nama homogen/heterogen.
- [lib/cdo_tranche/asb.rb](lib/cdo_tranche/asb.rb) & [lhp.rb](lib/cdo_tranche/lhp.rb) — Rekursi ASB, binomial eksak, dan bentuk tertutup LHP Vasicek.
- [lib/cdo_tranche/pricing.rb](lib/cdo_tranche/pricing.rb) & [base_corr.rb](lib/cdo_tranche/base_corr.rb) — *Fair spread* 6 *tranche*, studi granularitas $N$, *base correlation* vs *compound correlation*.
- [lib/cdo_tranche/importance.rb](lib/cdo_tranche/importance.rb) & [model_risk.rb](lib/cdo_tranche/model_risk.rb) — *Exponential twisting* 2-tahap ($R=100\mathrm{k}$) dan *model risk* Gaussian vs $t(\nu=4)$ vs *random recovery*.
- [lib/cdo_tranche/runner.rb](lib/cdo_tranche/runner.rb), [criteria.rb](lib/cdo_tranche/criteria.rb), [figures.rb](lib/cdo_tranche/figures.rb), [json_out.rb](lib/cdo_tranche/json_out.rb), [svg.rb](lib/cdo_tranche/svg.rb) — Runner, kriteria `c1`–`c9`, dan 5 grafik SVG.
- [test/test_core.rb](test/test_core.rb) & [test/test_studies.rb](test/test_studies.rb) — `10` tes unit + integrasi (`31` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `24,59s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. Anchor ASB & LHP | ASB vs binomial & LHP vs bivariat | `max_asb_err = 7,18e-16`, `max_lhp_err = 2,52e-10`, galat LHP turun `6,99e-4` ($N=25$) $\to$ `8,48e-5` ($N=500$) |
| 2. Valuasi 6 Tranche | Konservasi PV proteksi | `homog_pv_sum_err = 6,94e-18`, `hetero_pv_sum_err = 6,94e-18` |
| 3. Base & Compound Corr | Kalibrasi ke pasar $t(\nu=4)$ | `base_corr` `[0,455; 0,438; 0,432; 0,424; 0,410]`; *mezzanine* `0,3969`, *senior* `[]` (tidak ada akar) |
| 4. Importance Sampling | *Super-senior* `30-100%` ($R=100\mathrm{k}$) | `max |is_z| = 2,64`, penghematan variansi **`309,1x`** (`VR = 353,5x` pada `extras` $R=400\mathrm{k}$) |
| 5. Model Risk | *Spread* senior `15-30%` | Gaussian `27,44 bps` vs $t(\nu=4)$ `74,31 bps` vs *Random Recovery* `63,46 bps` |

Artefak keluaran: [results.json](../../../outputs/day005f/results.json), 5 grafik SVG `outputs/figures/day005f_*.svg`, dan dokumen penjelasan [day005f_quant_cdo-tranche.md](../../../docs/project-explanations/day005f_quant_cdo-tranche.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
