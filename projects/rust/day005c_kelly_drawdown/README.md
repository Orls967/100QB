# Day 005c — Kriteria Kelly dan Risiko Drawdown (Rust)

Proyek Rust mandiri (*Rust stable, edition 2021*, tanpa crate eksternal) yang membedah **kriteria Kelly** dan **risiko drawdown** pada taruhan biner diskret maupun pasar kontinu (GBM), lengkap dengan rumus peluang sentuh eksak, estimator jembatan Brownian tanpa bias, analisis Kelly fraksional, kerugian pertumbuhan akibat galat estimasi, penyusutan (*shrinkage*), serta simulasi ekor berat Student-$t$ ($\nu = 4$).

## Intuisi & Rumus Kunci

1. **Taruhan Biner Diskret**: Untuk peluang menang $p = 0{,}55$ dan *odds* bersih $b = 1$, laju pertumbuhan log adalah $g(f) = p\ln(1+bf) + (1-p)\ln(1-f)$ dengan fraksi Kelly $f^* = p - (1-p)/b = 0{,}10$ dan titik impas $f_{\mathrm{be}} \approx 0{,}198668$.
2. **Pasar Kontinu (GBM) & Peluang Sentuh Eksak**: Dengan $\mu = 0{,}08$, $\sigma = 0{,}20$, $r = 0{,}03$, log-kekayaan mengikuti $\mathrm{BM}(g(f), s(f))$ dengan $g(f) = r + f(\mu - r) - \frac{1}{2}f^2\sigma^2$, $s(f) = f\sigma$, dan $f^* = (\mu - r)/\sigma^2 = 1{,}25$.
   - Horizon tak hingga: $P(\inf_{t\ge 0} W_t \le x W_0) = x^{2g/s^2}$ (untuk $f = f^*$, $2g^*/s^{*2} = 1 + 2r/(f^{*2}\sigma^2)$, sehingga saat $r = 0$ diperoleh $P = x$).
   - Horizon hingga $T$: $P(\min_{t \le T} W_t \le x W_0) = \Phi\!\left(\frac{a - gT}{s\sqrt{T}}\right) + e^{2ga/s^2}\Phi\!\left(\frac{a + gT}{s\sqrt{T}}\right)$, $a = \ln x$.
3. **Jembatan Brownian Tanpa Bias**: Pemantauan diskret memiliki bias negatif $O(\sqrt{\Delta})$; probabilitas melintas antar langkah $\exp\!\left(-\frac{2(y_i - a)(y_{i+1} - a)}{s^2 \Delta}\right)$ menghilangkan bias diskretisasi.
4. **Kelly Fraksional**: Untuk $f = c f^*$, rasio pertumbuhan berlebih adalah $(g(cf^*) - r)/(g^* - r) = c(2 - c)$. Setengah-Kelly ($c = 0{,}5$) mempertahankan $75\%$ pertumbuhan berlebih sambil menekan peluang sentuh $0{,}5$ dari $25{,}70\%$ menjadi $0{,}87\%$.
5. **Galat Estimasi & Penyusutan**: Jika $\hat{\mu} \sim \mathcal{N}(\mu, \sigma^2 / T_{\mathrm{est}})$, kerugian pertumbuhan rata-rata $\mathbb{E}[g(f^*) - g(\hat{f})] = \frac{1}{2 T_{\mathrm{est}}}$ (invarian terhadap $\mu$ dan $\sigma$). Penyusutan optimal $c^* = \frac{f^{*2}}{f^{*2} + 1/(\sigma^2 T_{\mathrm{est}})}$ memulihkan pertumbuhan dan menekan *drawdown*.

## Build, Run & Test

```bash
cargo fmt --manifest-path projects/rust/day005c_kelly_drawdown/Cargo.toml --all -- --check
cargo clippy --manifest-path projects/rust/day005c_kelly_drawdown/Cargo.toml --all-targets --all-features -- -D warnings
cargo test --release --manifest-path projects/rust/day005c_kelly_drawdown/Cargo.toml
cargo run --release --manifest-path projects/rust/day005c_kelly_drawdown/Cargo.toml --bin run
```

## Struktur Kode

- [src/rng.rs](src/rng.rs) — `SplitMix64`, `xoshiro256++`, sampler uniform, Bernoulli, normal, Gamma, dan Student-$t$
- [src/special.rs](src/special.rs) — PDF/CDF normal (Cody) dan invers CDF normal
- [src/stats.rs](src/stats.rs) — `mean`, `variance_sample`, `sd`, `se`, `quantile`
- [src/parallel.rs](src/parallel.rs) — Eksekusi paralel deterministik (`std::thread::scope`, reduksi blok tetap `64`)
- [src/kelly_binary.rs](src/kelly_binary.rs) — Bagian 1: Kelly biner diskret, optimisasi numerik, simulasi $N=1000, R=100\mathrm{k}$
- [src/gbm.rs](src/gbm.rs) — Parameter GBM, `drift_log`, `vol_log`, `kelly_fraction`
- [src/hitting.rs](src/hitting.rs) — Bagian 2: rumus sentuh eksak & estimator jembatan Brownian ($R = 10^6$)
- [src/fractional.rs](src/fractional.rs) — Bagian 3: kurva $c(2-c)$, peta $(c, T)$, dan akar numerik vs tertutup
- [src/estimation.rs](src/estimation.rs) — Bagian 4.1–4.2: kerugian $1/(2T_{\mathrm{est}})$ dan penyusutan optimal
- [src/heavy_tail.rs](src/heavy_tail.rs) — Bagian 4.3–4.4: return harian Student-$t(4)$ dan tabel ringkasan praktis
- [src/jsonout.rs](src/jsonout.rs) & [src/svg.rs](src/svg.rs) — Penulis JSON 17-digit dan grafik SVG
- [src/bin/run.rs](src/bin/run.rs) — Runner penuh (`MAIN_SEED = 20261123`)
- [tests/integration_tests.rs](tests/integration_tests.rs) — Tes integrasi invariansi worker, determinisme, dan anchor 3 SE

## Ringkasan Hasil (`MAIN_SEED = 20261123`, `10` worker, `46,80s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. Kelly Biner | $f^*$ & $f_{\mathrm{be}}$; $E[\ln W_{1000}]$ | $f^* = 0{,}10$ (`diff = 0`), $f_{\mathrm{be}} = 0{,}198668$; `max_z_log = 2,07` |
| 2. Sentuh GBM | Sel informatif ($33$ sel, $R=10^6$) | `33/33` bridge dalam `3 SE` (`max |z| = 2,80`), naif negatif `33/33`, kemiringan $\sqrt{\Delta}$ = `0,4946 ± 0,0168` |
| 3. Kelly Fraksional | $c(2-c)$ & akar | `max_err = 4,44e-16`; $c=0{,}5$: `75%` *excess growth*, $P_\infty(0{,}5) = 0{,}008729$ |
| 4. Galat Estimasi | $E[g(f^*) - g(\hat{f})] = 1/(2T_{\mathrm{est}})$ | $T_{\mathrm{est}} \in \{1,5,10,20\}$: `max |z| = 1,41` (invarian terhadap $\mu,\sigma$) |
| 4. Penyusutan | $c^* = 0{,}2381$ ($T_{\mathrm{est}}=5$, 20 thn) | $g_{\mathrm{shrink}} = 0{,}03751$ vs $g_{\mathrm{full}} = -0{,}03810$ (`+0,0756 ± 0,0004`); $P_{20}(0{,}5)$: `5,82%` vs `53,71%` |
| 4. Ekor Berat $t_4$ | $t_4$ harian vs Gaussian harian ($20$ thn) | $2f^*$: $\Delta g = -0{,}00631 \pm 0{,}00057$ (`z = -11,04`), `ruin = 0,387%` vs `0%` |

Tes: `8` unit test + `3` integration test (`cargo test --release` lolos semua).

Artefak keluaran: [results.json](../../../outputs/day005c/results.json), grafik `outputs/figures/day005c_*.svg`, serta dokumen penjelasan [day005c_quant_kelly-drawdown.md](../../../docs/project-explanations/day005c_quant_kelly-drawdown.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
