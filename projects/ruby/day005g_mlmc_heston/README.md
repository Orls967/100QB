# Day 005g — Multilevel Monte Carlo (Giles) untuk CIR & Heston (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** menggunakan **topik cadangan** (*Multilevel Monte Carlo (Giles) untuk CIR/Heston dengan skema Milstein/full-truncation: laju kuat, hukum kompleksitas, dan estimasi biaya $\varepsilon^{-2}$ vs $\varepsilon^{-3}$*) setelah topik utama (*Hull–White & Bermudan swaption*) terdeteksi telah dikerjakan penuh pada [`day002i`](../../../docs/project-explanations/day002i_quant_hull-white.md). Seluruh data bersifat sintetis.

## Intuisi & Rumus Kunci

1. **Identitas Teleskopik MLMC (Giles, 2008)**: Untuk estimasi $\mathbb{E}[P_L]$ pada grid paling halus $L$:
   $$\mathbb{E}[P_L] = \mathbb{E}[P_0] + \sum_{\ell=1}^L \mathbb{E}[P_\ell - P_{\ell-1}],$$
   di mana $(P_\ell, P_{\ell-1})$ disimulasikan dengan inkremen Brown yang sama ($\Delta W_j^c = \Delta W_{2j}^f + \Delta W_{2j+1}^f$).
2. **Alokasi Sampel Optimal & Teorema Kompleksitas**: Dengan variansi level $V_\ell = \operatorname{Var}(P_\ell - P_{\ell-1}) \propto 2^{-\beta \ell}$, bias lemah $2^{-\alpha\ell}$, dan biaya per jalur $C_\ell \propto 2^{\gamma \ell}$ ($\gamma = 1$), alokasi optimal $N_\ell \propto \varepsilon^{-2}\sqrt{V_\ell / C_\ell}$ memprediksi kompleksitas asimtotik:
   - $\beta > 1$ (GBM Milstein, $\beta = 1{,}945$): $O(\varepsilon^{-2})$.
   - $\beta = 1$ (GBM Euler, $\beta = 0{,}982$): $O(\varepsilon^{-2}(\ln \varepsilon)^2)$.
   - $\beta < 1$ (Heston *Full-Truncation* saat kondisi Feller dilanggar, $\beta = 0{,}821$, $\alpha = 1{,}039$): $O(\varepsilon^{-2 - (1-\beta)/\alpha}) = O(\varepsilon^{-2{,}172})$ asimtotik, atau kemiringan pra-asimtotik diskret `2,475` pada $\varepsilon \in [0{,}01, 0{,}08]$.
3. **Koppel Antitetik Giles–Szpruch**: Menukar pasangan inkremen halus $(\Delta W_{2k}, \Delta W_{2k+1})$ meniadakan kesalahan *Lévy area* antar-koordinat; untuk skema Euler tanpa suku kuadratik Milstein, batas atas reduksi variansi adalah `1,0` (1D simetris GBM), meluruh menuju `1,0` (1D CIR dengan bias truncation kecil di $v=0$), dan $R_{\max} = \frac{2}{1+\rho^2} = 1{,}3423$ pada Heston 2D ($\rho = -0{,}7$). Terukur pada $\ell = 1..6$: GBM `1,0000`, Heston `1,2744 -> 1,1869`, CIR `1,0900 -> 1,0611` (`c6` tetap `belum tercapai` terhadap batas a priori `<= 1,04` untuk CIR pada $\ell=6$ akibat penolakan/pemotongan batas $v_n^+=0$).

## Instalasi, Lint, Test & Eksekusi

```bash
cd projects/ruby/day005g_mlmc_heston
bundle install --local
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Struktur Kode

- [lib/mlmc_heston/rng.rb](lib/mlmc_heston/rng.rb) & [parallel.rb](lib/mlmc_heston/parallel.rb) — `SplitMix64`, Box–Muller normal, dan `Process.fork` + `IO.pipe` + `Marshal` deterministik.
- [lib/mlmc_heston/special.rb](lib/mlmc_heston/special.rb) & [heston_cf.rb](lib/mlmc_heston/heston_cf.rb) — `norm_cdf`, `bs_call`, Gauss–Legendre, aritmetika kompleks (`cexp`, `clog`, `csqrt`), dan Heston *little trap* CF + Gil–Pelaez.
- [lib/mlmc_heston/simulators.rb](lib/mlmc_heston/simulators.rb) — Momen CIR eksak dan simulator terkopel GBM (Euler/Milstein) & Heston (*Full-Truncation* European/Asian).
- [lib/mlmc_heston/mlmc.rb](lib/mlmc_heston/mlmc.rb), [complexity.rb](lib/mlmc_heston/complexity.rb), [antithetic.rb](lib/mlmc_heston/antithetic.rb) — Regresi laju $(\alpha, \beta)$, MLMC adaptif, prediksi teori Giles + titik silang, dan *antithetic coupling*.
- [lib/mlmc_heston/runner.rb](lib/mlmc_heston/runner.rb), [criteria.rb](lib/mlmc_heston/criteria.rb), [figures.rb](lib/mlmc_heston/figures.rb), [json_out.rb](lib/mlmc_heston/json_out.rb), [svg.rb](lib/mlmc_heston/svg.rb) — Runner, kriteria `c1`–`c8`, dan 5 grafik SVG.
- [test/test_core.rb](test/test_core.rb) & [test/test_studies.rb](test/test_studies.rb) — `22` tes unit + integrasi (`38` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `90,26s`)

| Bagian | Metrik Utama | Hasil Terukur (`7/8` `tercapai`, `c6 = belum tercapai`) |
|---|---|---|
| 1. Anchor Tertutup | BS ($\xi=0$), paritas, konvergensi CF | `xi_zero_err = 4,26e-14`, `parity_err = 7,11e-15`, `cf_quad_err = 8,23e-10` |
| 2. Laju Kuat ($\beta$) | GBM & Heston (Feller OK vs dilanggar) | `GBM Euler = 0,982`, `GBM Milstein = 1,945`, `Heston Feller OK = 1,276`, `Feller dilanggar = 0,821` |
| 3. Kompleksitas & Speedup | Terukur vs Teori Giles & titik silang | `GBM slope = 1,876` (teori `2,000`, speedup **`14,90x`**); `Heston slope = 2,471` (pra-asimtotik `2,475`, asimtotik `2,172`) vs `Std MC = 3,000`, namun pada $\varepsilon=0{,}01$ **MLMC Heston justru `5,44x` lebih lambat** (`speedup = 0,1837x`) dan di $\varepsilon=0{,}005$ masih **`0,4978x`** (*crossover* ekstrapolasi $\varepsilon^* \approx 5{,}05\times 10^{-4}$..$1{,}59\times 10^{-3}$) |
| 4. Antithetic Giles–Szpruch | Multi-model Euler $\ell = 1..6$ & batas teori | `GBM = 1,0000`, `CIR: 1,0900 -> 1,0611`, `Heston: 1,2744 -> 1,1869` ($\le 1,3423$); `c6 = belum tercapai` karena `CIR_l6 = 1,0611 > 1,04` (batas Feller $v_n^+=0$) |

Artefak keluaran: [results.json](../../../outputs/day005g/results.json), 5 grafik SVG `outputs/figures/day005g_*.svg`, dan dokumen penjelasan [day005g_quant_mlmc-heston.md](../../../docs/project-explanations/day005g_quant_mlmc-heston.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
