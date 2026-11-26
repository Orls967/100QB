# Hari 8j (`day008j`) — Diagonalisasi Eksak Rantai Spin-½ Heisenberg/XXZ (Rust Murni, ≤ 499 Baris)

Implementasi Rust murni (edisi 2021, **tanpa crate eksternal**, total `Cargo.toml + src/*.rs = 494` baris pasca-`cargo fmt` $\le 499$) untuk rantai periodik $N \in \{4,6,\dots,24\}$ spin-$\tfrac{1}{2}$ model XXZ ($J=1, \Delta \in \{0, 0.5, 1, 2\}$):

- `selftest`: verifikasi nilai tangan $N=4, \Delta=1$ ($E_0 = -2$) dan matriks ikatan tunggal 2-spin ($E_0 = -3/4$).
- `all`: menjalankan seluruh alur secara deterministik (`SEED = 20261126`):
  - Pembangkitan basis sektor $S^z$ dengan *Gosper's hack* + *combinadic ranking* $O(N)$, tanpa matriks jarang eksplisit (*on-the-fly* $H|v\rangle$).
  - Lanczos dengan reortogonalisasi MGS penuh ($\le 120$ vektor) dan pembanding Jacobi dense eksak ($N \le 12$), ditambah anchor fermion bebas Jordan–Wigner ($\Delta=0, N \le 24$).
  - Identitas internal ($\langle S_{\text{tot}}^2 \rangle = 0$, $\langle S_0\cdot S_1\rangle = E_0/N$, dan kemonotonan $E_0(\Delta)$).
  - Ekstrapolasi CFT $E_0/N$ vs $1/N^2$ ($\Delta=1$) serta celah singlet–triplet $S^z=1$ vs $S^z=0$ ($N \le 20, \Delta=1, 2$).
  - Entropi keterikatan von Neumann lewat dekomposisi Schmidt per blok $S^z_A$ ($N \in \{12,16,20\}$) dan regresi Calabrese–Cardy.
  - 6 grafik SVG ke `../../../outputs/figures/day008j_*.svg`, 5 CSV ringkas ke `out/`, tolok ukur kecepatan 5 ulangan, serta `../../../outputs/day008j/results.json`.

> Fitur `SHOULD` tambahan (aturan tanda Marshall detail/peta dan vektor korelasi penuh $\sum_r \langle S_0\cdot S_r\rangle$) dipangkas sesuai aturan batas baris $\le 499$ (`belum dijalankan` di `results.json`).

```sh
make test
make lines
make run
```
