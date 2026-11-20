# Hari 003b — Eksekusi Optimal (Rust Murni)

Implementasi Rust murni (`#![forbid(unsafe_code)]`, nol crate eksternal runtime) untuk tiga model eksekusi optimal:
1. **Almgren–Chriss (2000)**: *implementation shortfall* dengan dampak permanen ($\gamma$) dan sementara ($\eta$), QP kuadratik via Cholesky KKT, solusi tertutup `sinh` diskrit, limit kontinu, frontier efisien, dan simulasi Monte Carlo antitetik (`400.000` jalur).
2. **Dampak Transien Obizhaeva–Wang (2013) & Kernel Umum**: kernel eksponensial $G_{ij} = \lambda_{\mathrm{ow}} e^{-\rho |t_i - t_j|}$, solver analitik $O(N)$ vs Cholesky $O(N^3)$, kernel *power-law*, kernel berosilasi positif-definit dengan manipulasi harga (*round-trip*, Alfonsi–Schied–Slynko 2012), dan kernel non-PD.
3. **Likuiditas Stokastik (Dynamic Programming)**: rantai Markov 2-rezim (*likuid* / *ilikuid*), rekursi mundur kuadratik $V_k(y, s) = a_k(s) y^2$, verifikasi *brute-force* $N=5$, simulasi `400.000` jalur paralel deterministik, nilai adaptivitas vs strategi statis, dan sensitivitas persistensi.

## Bukti Aljabar Invers Tridiagonal Obizhaeva–Wang
Misalkan $a = e^{-\rho \Delta} \in (0, 1)$ dan $M_{ij} = a^{|i-j|}$ untuk $i, j = 1 \dots N$, sehingga $G = \lambda_{\mathrm{ow}} M$. Perhatikan matriks tridiagonal simetris $T$ dengan entri:
- $T_{11} = T_{NN} = \frac{1}{1-a^2}$,
- $T_{ii} = \frac{1+a^2}{1-a^2}$ untuk $2 \le i \le N-1$,
- $T_{i,i\pm 1} = \frac{-a}{1-a^2}$.

Untuk baris $i = 1$: $(T M)_{1j} = \frac{a^{j-1} - a \cdot a^{|2-j|}}{1-a^2}$. Jika $j = 1$, nilainya $\frac{1 - a^2}{1-a^2} = 1$; jika $j \ge 2$, $|2-j| = j-2$, sehingga $a^{j-1} - a^{j-1} = 0$. Untuk baris interior $2 \le i \le N-1$:
$$(T M)_{ij} = \frac{-a \cdot a^{|i-1-j|} + (1+a^2) a^{|i-j|} - a \cdot a^{|i+1-j|}}{1-a^2}.$$
Jika $j = i$, hasilnya $\frac{-a^2 + 1 + a^2 - a^2}{1-a^2} = 1$. Jika $j > i$, $|i-1-j| = j-i+1$ dan $|i+1-j| = j-i-1$, sehingga pembilangnya menjadi $(-a^2 + (1+a^2) - 1) a^{j-i} = 0$ (identik untuk $j < i$). Maka $T = M^{-1}$.

Mengalikan $M^{-1}$ dengan vektor $\mathbf{1}$ menghasilkan:
- $(M^{-1}\mathbf{1})_1 = (M^{-1}\mathbf{1})_N = \frac{1-a}{1-a^2} = \frac{1}{1+a}$,
- $(M^{-1}\mathbf{1})_i = \frac{1+a^2-2a}{1-a^2} = \frac{1-a}{1+a}$ untuk $2 \le i \le N-1$.

Menormalkan agar $\sum_{i=1}^N x_i^* = X$ langsung memberikan:
$$x_1^* = x_N^* = \frac{X}{2 + (N-2)(1-a)}, \qquad x_i^* = \frac{X(1-a)}{2 + (N-2)(1-a)} \quad (2 \le i \le N-1).$$

## Cara Build, Run, dan Test
```bash
cd rust/day003b_optimal_execution
cargo fmt --check
cargo clippy --all-targets -- -D warnings
cargo test --release
cargo run --release
```

## Struktur Kode
- `src/linalg.rs` — Cholesky, Thomas tridiagonal $O(N)$, rotasi Jacobi eigenvalue simetris.
- `src/rng.rs` — `Xoshiro256PlusPlus` (seeding SplitMix64) + Box–Muller normal.
- `src/parallel.rs` — paralelisme `std::thread` dengan chunk tetap (invarian terhadap jumlah thread).
- `src/stats.rs` — rata-rata, variansi, SE, dan *batch means* SE variansi.
- `src/json.rs` & `src/svg.rs` — penulis JSON (17 digit) dan grafik SVG buatan sendiri.
- `src/ac.rs` — model Almgren–Chriss (QP KKT, `sinh`, limit kontinu, frontier, MC).
- `src/transient.rs` — Obizhaeva–Wang, benchmark $N=2000$, power-law, osilasi, non-PD, gabungan.
- `src/dp_liquidity.rs` — DP kuadratik 2-rezim, verifikasi brute-force, simulasi, adaptivitas.
- `src/main.rs` — runner utama penghasil `outputs/day003b/results.json` dan 7 grafik SVG.

## Ringkasan Hasil
- **Almgren–Chriss**: QP vs `sinh` diskrit selisih relatif `4.15e-16`; MC `400.000` jalur: variansi `z_var = 0.641`; mean tanpa antitetik `z = 1.032` (optimal), `0.945` (TWAP); mean dengan antitetik adalah uji vakum (z ~ 1e-11).
- **Obizhaeva–Wang**: Anchor eksak selisih relatif `1.37e-15`; pada $N=2000$, Cholesky `0.625 s` vs analitik `4 us`.
- **DP Likuiditas**: Brute-force $N=5$ selisih `1.20e-16`; simulasi `400.000` jalur `z = 1.641`; penghematan adaptivitas `14.09%` terhadap strategi statis stasioner.
