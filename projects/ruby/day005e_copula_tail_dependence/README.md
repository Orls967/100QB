# Day 005e — Kopula, Ketergantungan Ekor, & Alokasi Modal Euler (Ruby)

Proyek **Ruby murni ($\ge 3.3$, nol gem numerik)** yang mengimplementasikan pemodelan ketergantungan non-linear via kopula bivariat/multivariat, koefisien ketergantungan ekor, alokasi modal Euler (*Expected Shortfall*), estimasi *Canonical Maximum Likelihood* (CML) beserta studi salah spesifikasi, serta *importance sampling* kejadian gagal-bersama langka pada data sintetis.

## Intuisi & Rumus Kunci

1. **Teorema Sklar & Kopula Bivariat**: Setiap distribusi bersama kontinu $F(x_1, x_2) = C(F_1(x_1), F_2(x_2))$ dapat dipisahkan menjadi marginal dan struktur ketergantungan $C$:
   - **Gaussian($\rho$) & Student-$t(\nu, \rho)$**: $\tau = \frac{2}{\pi}\arcsin\rho$; Gaussian memiliki $\rho_S = \frac{6}{\pi}\arcsin(\rho/2)$ dan $\lambda_L = \lambda_U = 0$ ($|\rho|<1$), sedangkan $t(\nu,\rho)$ memiliki ketergantungan ekor simetris $\lambda_L = \lambda_U = 2\,t_{\nu+1}\!\left(-\sqrt{(\nu+1)\frac{1-\rho}{1+\rho}}\right)$.
   - **Clayton($\theta$)** (*frailty* Gamma Marshall–Olkin): $\tau = \frac{\theta}{\theta+2}$, $\lambda_L = 2^{-1/\theta}$, $\lambda_U = 0$.
   - **Gumbel($\theta$)** (*frailty* *positive-stable* CMS): $\tau = 1 - \frac{1}{\theta}$, $\lambda_U = 2 - 2^{1/\theta}$, $\lambda_L = 0$.
2. **Alokasi Modal Euler**: Untuk portofolio $L = -w^\top X$ dengan $X \sim \mathcal{N}(0, \Sigma)$, $\mathrm{ES}_\alpha = \sigma_p \frac{\phi(z_\alpha)}{1-\alpha}$ homogen derajat $1$ terhadap $w$, sehingga kontribusi Euler $c_i = w_i \frac{\partial \mathrm{ES}_\alpha}{\partial w_i} = \frac{w_i(\Sigma w)_i}{\sigma_p}\frac{\phi(z_\alpha)}{1-\alpha} = \mathbb{E}[L_i \mid L \ge \mathrm{VaR}_\alpha]$ memenuhi $\sum_{i=1}^d c_i = \mathrm{ES}_\alpha$.
3. **Importance Sampling Faktor Tunggal**: Untuk $Z_i = \sqrt{\rho}\,F + \sqrt{1-\rho}\,\varepsilon_i \le a$, pergeseran rerata analitik $\mu^* \approx \frac{d\sqrt{\rho}\,a}{1 - \rho + d\rho}$ (dipoles Newton) pada $F \sim \mathcal{N}(\mu^*, 1)$ digabungkan dengan *conditional MC* menghasilkan penghematan variansi ribuan hingga jutaan kali lipat.

## Instalasi, Lint, Test & Eksekusi

```bash
cd projects/ruby/day005e_copula_tail_dependence
bundle install --local
bundle exec rubocop
bundle exec rake test
bundle exec bin/run
bundle exec ruby bin/run_extras
```

## Struktur Kode

- [lib/copula_tail_dependence/rng.rb](lib/copula_tail_dependence/rng.rb) — `SplitMix64`, sampler polar normal, Marsaglia–Tsang Gamma, $\chi^2$, Student-$t$, dan *positive-stable* CMS.
- [lib/copula_tail_dependence/parallel.rb](lib/copula_tail_dependence/parallel.rb) — Paralelisme deterministik `Process.fork` + `IO.pipe` + `Marshal` (chunk tetap).
- [lib/copula_tail_dependence/special.rb](lib/copula_tail_dependence/special.rb) — PDF/CDF normal, kuantil Wichura AS241, beta tak-lengkap, CDF/kuantil $t$, Gauss–Legendre, dan ortan bivariat Plackett.
- [lib/copula_tail_dependence/linalg.rb](lib/copula_tail_dependence/linalg.rb) & [lib/copula_tail_dependence/stats.rb](lib/copula_tail_dependence/stats.rb) — Cholesky, produk matriks, statistik sampel, dan SE *jackknife*.
- [lib/copula_tail_dependence/copulas.rb](lib/copula_tail_dependence/copulas.rb) & [lib/copula_tail_dependence/dependence.rb](lib/copula_tail_dependence/dependence.rb) — Sampler & log-densitas kopula, Kendall $\tau$ $O(n\log n)$ (Fenwick), Spearman $\rho_S$, dan $\hat{\chi}(u)$.
- [lib/copula_tail_dependence/euler.rb](lib/copula_tail_dependence/euler.rb) — Alokasi Euler analitik, beda hingga, dan MC ($R=200\mathrm{k}$ / $1\mathrm{M}$).
- [lib/copula_tail_dependence/fit.rb](lib/copula_tail_dependence/fit.rb) — *Pseudo-observations*, CML Gaussian/Clayton/profil-$t$, studi $R=500$, dan salah spesifikasi ES99 ($R=200$).
- [lib/copula_tail_dependence/importance.rb](lib/copula_tail_dependence/importance.rb) — Kuadratur 1D & *importance sampling* gagal-bersama langka ($R=100\mathrm{k}$).
- [lib/copula_tail_dependence/runner.rb](lib/copula_tail_dependence/runner.rb), [json_out.rb](lib/copula_tail_dependence/json_out.rb), [svg.rb](lib/copula_tail_dependence/svg.rb) — Runner utama, penulis JSON, dan 7 grafik SVG.
- [test/test_core.rb](test/test_core.rb) & [test/test_studies.rb](test/test_studies.rb) — `12` tes unit + integrasi (`62` asersi).

## Ringkasan Hasil (`MAIN_SEED = 20261123`, Ruby `4.0.9`, YJIT `ON`, `10` worker, `9,20s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. Ketergantungan ($R=200\mathrm{k}$) | Kendall $\tau$, Spearman $\rho_S$, ortan | `max |tau_z| = 1,66`, `spearman_err = 8,8e-5`, `max_orthant_err = 5,55e-17` |
| 1b. Ekor $\chi(u)$ | Konvergensi ke $\lambda$ & $\chi(u)$ eksak | `Gauss max |z| = 0,23` vs $\chi(u)$ eksak; `t4/Clayton/Gumbel` $u=0{,}999$ `max |z| = 0,88` |
| 2. Alokasi Euler ($d=5$) | Eksak vs FD vs MC | `sum_err = 0`, `fd_err = 5,50e-12`, `max_mc_z = 0,53` (`ES99 = 0,372` Gauss vs `0,414` $t_4$ vs `0,551` Clayton) |
| 3. Estimasi & Salah Spesifikasi | CML vs $\tau$, LRT, ES99 | `n=5000 |bias| < 0,0008`; `LRT = 155,35` (Gaussian ditolak); `ES99 = 2,062` (Gauss) vs `2,228` ($t_4$, `-7,44%`) |
| 4. Importance Sampling ($R=100\mathrm{k}$) | $d=2, 5$; $\rho=0,5, 0,8$ | `max |is_z| = 2,45`, penghematan variansi `2.604x` hingga `5.165.715x` |

Artefak keluaran: [results.json](../../../outputs/day005e/results.json), 7 grafik SVG `outputs/figures/day005e_*.svg`, dan dokumen penjelasan [day005e_quant_copula-tail-dependence.md](../../../docs/project-explanations/day005e_quant_copula-tail-dependence.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
