# Hari 8j (`day008j`) — Diagonalisasi Eksak Rantai Spin-½ Heisenberg/XXZ (Rust Murni, ≤ 499 Baris)

Implementasi Rust murni (edisi 2021, tanpa crate luar, total `Cargo.toml + src/*.rs ≤ 499` baris setelah `cargo fmt`) untuk rantai periodik $N$ spin-$\tfrac{1}{2}$ model XXZ ($J=1, \Delta \in \{0, 0.5, 1, 2\}$):

- `selftest`: verifikasi nilai tangan $N=4, \Delta=1$ ($E_0 = -2$) dan matriks ikatan tunggal 2-spin ($E_0 = -3/4$).
- `ed`: energi dasar $E_0(N, \Delta)$ untuk $N \in \{4,6,\dots,24\}$ (Lanczos dengan reortogonalisasi penuh $\le 120$), pembanding Jacobi dense ($N \le 12$), anchor eksak fermion bebas Jordan–Wigner ($\Delta=0, N \le 24$), dan identitas internal ($\langle S_{\text{tot}}^2 \rangle$, $\langle S_i\cdot S_{i+1}\rangle = E_0/N$, $\sum_r \langle S_0\cdot S_r\rangle = 0$).
- `scale`: regresi CFT $E_0/N$ vs $1/N^2$ ($\Delta=1$) serta celah singlet–triplet $S^z=1$ vs $S^z=0$ ($\Delta=1, 2$).
- `entangle`: dekomposisi Schmidt per blok $S^z_A$ dan regresi Calabrese–Cardy $S(\ell)$ vs $\ln[(N/\pi)\sin(\pi\ell/N)]$.
- `marshall`: uji aturan tanda Marshall pada sublattice ganjil ($N \le 16$, semua $\Delta$).
- `pics` & `all`: 6 grafik SVG ke `../../../outputs/figures/day008j_*.svg`, CSV ringkas ke `out/`, serta `../../../outputs/day008j/results.json`.

```sh
make
make test
make lines
make run
```
