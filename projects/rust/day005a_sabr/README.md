# Day 005a — Model SABR: CEV Eksak, Ekspansi Hagan, Kepadatan Negatif, Kalibrasi & Delta Bartlett (Rust)

Proyek Rust mandiri (*zero external crates*) yang membedah model **SABR** (*Stochastic Alpha, Beta, Rho*) dari fondasi analitik hingga simulasi Monte Carlo, analisis arbitrase kepadatan, kalibrasi smile, dan lindung nilai delta Bartlett pada data sintetis:

$$dF_t = \alpha_t F_t^\beta\,dW_{1,t}, \qquad d\alpha_t = \nu \alpha_t\,dW_{2,t}, \qquad d\langle W_1, W_2 \rangle_t = \rho\,dt,$$

dengan batas *absorbing* di $F = 0$, $r = 0$, $F_0 = 0{,}04$, $\beta = 0{,}5$, $\alpha_0 = 0{,}0325$, $\nu = 0{,}4$, $\rho = -0{,}3$, $T \in \{1, 5, 10\}$ tahun, dan *strike* $K \in \{0{,}5\%; 1\%; 2\%; 3\%; 4\%; 5\%; 6\%; 8\%; 10\%\}$.

## Intuisi & Rumus Kunci

1. **Kasus $\nu = 0$ (CEV Eksak via BESQ)**: Transformasi $X_t = \frac{F_t^{2(1-\beta)}}{\alpha^2(1-\beta)^2}$ menghasilkan proses *squared Bessel* berdimensi $\delta = \frac{1-2\beta}{1-\beta}$. Harga *call* dengan penyerapan di $0$ adalah:
   $$C_{\mathrm{CEV}}(F_0, K) = F_0\,Q_{\chi'^2}\!\left(y;\; 2 + \frac{1}{1-\beta},\; x\right) - K\,F_{\chi'^2}\!\left(x;\; \frac{1}{1-\beta},\; y\right),$$
   dengan $x = \frac{F_0^{2(1-\beta)}}{\alpha^2(1-\beta)^2 T}$, $y = \frac{K^{2(1-\beta)}}{\alpha^2(1-\beta)^2 T}$, serta massa penyerapan $P(F_T = 0) = \Gamma_{\mathrm{inc}}\!\left(\frac{1}{2(1-\beta)}, \frac{x}{2}\right)$ (untuk $\beta = 1/2$, $e^{-2F_0/(\alpha^2 T)}$).
2. **Ekspansi Hagan (2002)**: Aproksimasi asimtotik volatilitas implisit Black $\sigma_B(K, F)$ ke orde $O(T)$, menangani $z \to 0$ secara stabil via deret Taylor $z/x(z) = 1 - \frac{\rho}{2}z + \frac{2-3\rho^2}{12}z^2 + O(z^3)$.
3. **Referensi $\rho = 0$ Kondisional**: Bersyarat pada jalur $\alpha$, $F_T$ berdistribusi CEV dengan variansi terintegrasi $I = \int_0^T \alpha_s^2\,ds$, sehingga $C = \mathbb{E}[C_{\mathrm{CEV}}(F_0, K;\, \alpha^2 T \to I)]$.
4. **Kepadatan Implisit Breeden–Litzenberger**: $q(K) = \partial^2 C / \partial K^2$. Karena Hagan mengabaikan massa penyerapan di $0$ dan kehilangan akurasi pada $T$ besar / $K$ kecil, $q_{\mathrm{Hagan}}(K)$ menjadi negatif di sayap bawah.
5. **Delta Bartlett (2006)**: Karena $d\alpha_t = \frac{\rho\nu}{F_t^\beta}\,dF_t + \nu\alpha_t\sqrt{1-\rho^2}\,dW_2^\perp$, delta *minimum-variance* adalah:
   $$\Delta_{\mathrm{Bartlett}} = \frac{\partial C}{\partial F} + \frac{\partial C}{\partial \alpha}\cdot\frac{\rho\nu}{F^\beta}.$$

## Build, Run & Test

```bash
cargo fmt --manifest-path projects/rust/day005a_sabr/Cargo.toml --all -- --check
cargo clippy --manifest-path projects/rust/day005a_sabr/Cargo.toml --all-targets --all-features -- -D warnings
cargo test --release --manifest-path projects/rust/day005a_sabr/Cargo.toml
cargo run --release --manifest-path projects/rust/day005a_sabr/Cargo.toml --bin run
cargo run --release --manifest-path projects/rust/day005a_sabr/Cargo.toml --bin run_extras
```

## Struktur Kode

- [src/rng.rs](src/rng.rs) — SplitMix64, `xoshiro256++`, normal, Poisson, Gamma, dan `ncx2`
- [src/special.rs](src/special.rs) — PDF/CDF normal, `lgamma`, gamma tak lengkap, CDF/SF $\chi'^2$ dengan rekurensi, kuadratur Gauss–Legendre
- [src/black.rs](src/black.rs) — Harga Black-76, delta, vega, dan inversi volatilitas implisit (Newton + bisection)
- [src/cev.rs](src/cev.rs) — Rumus tertutup CEV, solver PDE Crank–Nicolson 1D (grid `sinh` + Rannacher), MC transisi eksak, dan MC Euler
- [src/hagan.rs](src/hagan.rs) — Volatilitas implisit dan harga *call* Hagan (2002)
- [src/sabr_mc.rs](src/sabr_mc.rs) — MC kondisional $\rho = 0$ (Simpson) dan MC Euler penuh (*Absorbing* & *FullTruncation*)
- [src/density.rs](src/density.rs) — Kepadatan implisit Breeden–Litzenberger, deteksi massa negatif, dan harga digital/*butterfly*
- [src/calib.rs](src/calib.rs) — Optimisasi Nelder–Mead, kalibrasi *smile*, dan matriks korelasi Hessian numerik
- [src/hedge.rs](src/hedge.rs) — Delta Bartlett vs naif vs Black ATM vs CEV pada simulasi lindung nilai harian CRN
- [src/bin/run.rs](src/bin/run.rs) & [src/bin/run_extras.rs](src/bin/run_extras.rs) — Runner utama dan eksperimen tambahan

## Ringkasan Hasil (`MAIN_SEED = 20261123`, `8` worker, `222,91s`)

| Bagian | Metrik Utama | Hasil Terukur |
|---|---|---|
| 1. CEV ($\nu = 0$) | Eksak vs PDE CN vs MC Eksak/Euler | `max_pde_rel_err = 1,67e-8`, `max_mc_exact_z = 1,78`, `max_mc_euler_z = 1,63`; $P(F_{10}=0) = 0{,}000514$ vs `0,000475` |
| 2. Hagan $\nu \to 0$ | Galat IV vs CEV eksak | `0,001..0,20 bps` (ATM), `0,75..1,90 bps` ($K=1\%$) |
| 3. Referensi $\nu > 0$ | $\rho = 0$ kondisional vs Euler ($R=2\mathrm{M}$) | `cond = 26,2515 ± 0,0062 bps`, `z_abs(n=800) = 0,53`, `z_ft(n=800) = 0,25` |
| 3. Smile ($\rho=-0{,}3$) | Galat Hagan vs MC ($K \in [2\%, 8\%]$) | $T=1$: `max_abs = 7,00 bps`; $T=5$: `31,1 bps` ATM; $T=10$: `134,2 bps` ATM |
| 4. Kepadatan & Arbitrase | Massa negatif $\int q^- dK$ | $T=1$: `0`; $T=5$: `0,00156`; $T=10$: `0,10195` (`min_q = -88,59`) |
| 5. Kalibrasi & Bartlett | SD P&L harian ($T=1, R=20\mathrm{k}$ CRN) | Bartlett `6,318 bps` < Naif `6,598 bps` < Black `7,055 bps` < CEV `7,167 bps` (`var_diff = -3,62 ± 0,18 bps²`) |

Artefak keluaran: [results.json](../../../outputs/day005a/results.json), [extras.json](../../../outputs/day005a/extras.json), grafik `outputs/figures/day005a_*.svg`, serta dokumen penjelasan lengkap [day005a_quant_sabr.md](../../../docs/project-explanations/day005a_quant_sabr.md).

> Catatan: Seluruh data bersifat sintetis untuk tujuan riset/pembelajaran kuantitatif dan bukan nasihat investasi.
