# Day 005g — Multilevel Monte Carlo (Giles) untuk CIR & Heston (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** menggunakan **topik cadangan** (*Multilevel Monte Carlo (Giles) untuk CIR/Heston dengan skema Milstein/full-truncation: laju kuat, hukum kompleksitas, dan estimasi biaya $\varepsilon^{-2}$ vs $\varepsilon^{-3}$*) setelah topik utama (*Hull–White & Bermudan swaption*) terdeteksi telah dikerjakan penuh pada [`day002i`](../../../docs/project-explanations/day002i_quant_hull-white.md). Seluruh data bersifat sintetis.

## Intuisi & Rumus Kunci

1. **Identitas Teleskopik MLMC (Giles, 2008)**: Untuk estimasi $\mathbb{E}[P_L]$ pada grid paling halus $L$:
   $$\mathbb{E}[P_L] = \mathbb{E}[P_0] + \sum_{\ell=1}^L \mathbb{E}[P_\ell - P_{\ell-1}],$$
   di mana $(P_\ell, P_{\ell-1})$ disimulasikan dengan inkremen Brown yang sama ($\Delta W_j^c = \Delta W_{2j}^f + \Delta W_{2j+1}^f$).
2. **Alokasi Sampel Optimal & Teorema Kompleksitas**: Dengan variansi level $V_\ell = \operatorname{Var}(P_\ell - P_{\ell-1}) \propto 2^{-\beta \ell}$ dan biaya per jalur $C_\ell \propto 2^{\gamma \ell}$ ($\gamma = 1$), alokasi optimal $N_\ell \propto \varepsilon^{-2}\sqrt{V_\ell / C_\ell}$ menghasilkan total kompleksitas:
   - $\beta > 1$ (GBM Milstein, $\beta \approx 1{,}95$): $O(\varepsilon^{-2})$.
   - $\beta = 1$ (GBM Euler, $\beta \approx 0{,}98$): $O(\varepsilon^{-2}(\ln \varepsilon)^2)$.
   - $\beta < 1$ (Heston *Full-Truncation* saat kondisi Feller $2\kappa\theta \ge \xi^2$ dilanggar, $\beta \approx 0{,}82$): $O(\varepsilon^{-2 - (1-\beta)/\alpha})$.
3. **Koppel Antitetik Giles–Szpruch**: Menukar pasangan inkremen halus $(\Delta W_{2k}, \Delta W_{2k+1})$ menghasilkan dua jalur halus dengan rata-rata yang lebih dekat ke jalur kasar, menekan variansi level Heston tanpa koreksi *Lévy area*.

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
- [lib/mlmc_heston/mlmc.rb](lib/mlmc_heston/mlmc.rb), [complexity.rb](lib/mlmc_heston/complexity.rb), [antithetic.rb](lib/mlmc_heston/antithetic.rb) — Regresi laju $(\alpha, \beta)$, MLMC adaptif, studi kompleksitas vs MC standar, dan *antithetic coupling*.
- [lib/mlmc_heston/runner.rb](lib/mlmc_heston/runner.rb), [criteria.rb](lib/mlmc_heston/criteria.rb), [figures.rb](lib/mlmc_heston/figures.rb), [json_out.rb](lib/mlmc_heston/json_out.rb), [svg.rb](lib/mlmc_heston/svg.rb) — Runner, kriteria `c1`–`c8`, dan 5 grafik SVG.
- [test/test_core.rb](test/test_core.rb) & [test/test_studies.rb](test/test_studies.rb) — `10` tes unit + integrasi (`15` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `72,14s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. Anchor Tertutup | BS ($\xi=0$), paritas, konvergensi CF | `xi_zero_err = 4,26e-14`, `parity_err = 7,11e-15`, `cf_quad_err = 8,23e-10` |
| 2. Laju Kuat ($\beta$) | GBM & Heston (Feller OK vs dilanggar) | `GBM Euler = 0,982`, `GBM Milstein = 1,945`, `Heston Feller OK = 1,276`, `Feller dilanggar = 0,821` |
| 3. Kompleksitas & Akurasi | $\ln(\text{Cost})$ vs $\ln(1/\varepsilon)$ | `GBM Milstein slope = 1,876` (`14,9x` speedup), `Heston MLMC = 2,471` vs `Std MC = 3,000`; `z <= 0,65` |
| 4. Antithetic Giles–Szpruch | Pengurangan variansi $\ell = 3$ | `13,51 -> 10,78` (`1,253x`) |

Artefak keluaran: [results.json](../../../outputs/day005g/results.json), 5 grafik SVG `outputs/figures/day005g_*.svg`, dan dokumen penjelasan [day005g_quant_mlmc-heston.md](../../../docs/project-explanations/day005g_quant_mlmc-heston.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
