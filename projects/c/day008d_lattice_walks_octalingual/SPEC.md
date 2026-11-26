# Spesifikasi Teknis: Jalan Acak Kisi Pólya Lintas 8 Bahasa (`day008d`)

Dokumen ini adalah kontrak tunggal yang diikuti secara identik oleh **8 bahasa** (**C11**, **C++20**, **Rust**, **Go**, **Ruby**, **Python**, **TypeScript**, dan **Jupyter Notebook `.ipynb`**). Setiap bahasa mengimplementasikan keseluruhan **12 tahap (`s01`–`s12` dan `all`)** secara mandiri dan menulis berkas keluaran bernama sama ke `out/<bahasa>/`.

---

## 1. Model Jalan Acak Kisi ($\mathbb{Z}^d$, $d \in \{1, 2, 3\}$)

- **Keadaan awal**: Titik asal $(0, \dots, 0) \in \mathbb{Z}^d$.
- **Aturan langkah**: Pada tiap langkah, tarik tepat satu bilangan acak $u_{01} \in [0, 1)$, hitung $k = \lfloor u_{01} \cdot 2d \rfloor \in \{0, \dots, 2d-1\}$ (dijaga $k = 2d-1$ bila pembulatan mencapai $2d$). Sumbu yang bergerak adalah $a = \lfloor k/2 \rfloor$, dengan perpindahan $+1$ jika $k \bmod 2 = 0$ dan $-1$ jika $k \bmod 2 = 1$.

---

## 2. RNG Bit-Identik Lintas 8 Bahasa & Hash FNV-1a

- **Benih utama**: `MASTER_SEED = 20261126` (`uint64`).
- **SplitMix64**:
  ```text
  splitmix64_mix(z):
      z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) mod 2^64
      z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) mod 2^64
      return (z ^ (z >> 31)) mod 2^64
  ```
- **xoshiro256++**: state `s[0..3]` diisi 4 panggilan berurutan `st = (st + 0x9E3779B97F4A7C15) mod 2^64; s[i] = splitmix64_mix(st)`.
  ```text
  rotl(x, k) = ((x << k) | (x >> (64 - k))) mod 2^64
  next_u64():
      res = (rotl((s[0] + s[3]) mod 2^64, 23) + s[0]) mod 2^64
      t = (s[1] << 17) mod 2^64
      s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]
      s[2] ^= t; s[3] = rotl(s[3], 45)
      return res
  u01() = (next_u64() >> 11) * 2^-53
  ```
- **Vektor acuan S01 (`Xoshiro256pp(20261126)`)**:
  - `x[0] = 0xA9920E0E7DAB3EF0`
  - `x[1] = 0x033B9CEB4B3BA7A0`
  - `x[2] = 0xE7473A6A39AFCEA6`
- **Benih per tugas**:
  ```text
  task_seed(d, stage, idx):
      h = 0x517CC1B727220A95
      h = splitmix64_mix(h ^ uint64(d))
      h = splitmix64_mix(h ^ uint64(stage))
      h = splitmix64_mix(h ^ uint64(idx))
      return splitmix64_mix(MASTER_SEED ^ h)
  ```
- **Hash FNV-1a 64-bit** (S02): Untuk tiap $d$, kumpulkan $8 \times d$ koordinat akhir `int32` little-endian (4 byte per koordinat, urutan jalan $w = 0 \dots 7$, sumbu $a = 0 \dots d-1$), lalu hash byte-demi-byte dengan basis `0xCBF29CE484222325` dan prima `0x100000001B3`.

---

## 3. Anchor Eksak (Dihitung Mandiri oleh Tiap Bahasa)

1. **MSD**: $\mathbb{E}[\|R_t\|^2] = t$ persis untuk semua $d \in \{1, 2, 3\}$.
2. **Peluang di asal $u_{2n}$ & fungsi pembaruan $F_d(T)$**:
   - $d = 1$: $u_0^{(1)} = 1$, $u_{2n}^{(1)} = u_{2(n-1)}^{(1)} \cdot \frac{2n - 1}{2n}$.
   - $d = 2$: $u_{2n}^{(2)} = \bigl(u_{2n}^{(1)}\bigr)^2$.
   - $d = 3$: $u_{2n}^{(3)} = u_{2n}^{(1)} \cdot w_n$ dengan $w_0 = 1$, $w_1 = \frac{1}{3}$, dan rekursi 3-suku eksak Domb–Joyce:
     $$w_{n+1} = \frac{(10n^2 + 10n + 3)\,w_n - n^2\,w_{n-1}}{9(n+1)^2}$$
     yang diverifikasi terhadap jumlah multinomial langsung $\sum_{i+j+k=n} \bigl(\frac{n!}{3^n i!j!k!}\bigr)^2$ untuk $n \le 4$ pada `S01` (`selftest`).
   - **Renewal**: $f_{2n} = u_{2n} - \sum_{k=1}^{n-1} f_{2k} u_{2(n-k)}$ (suku ganjil $0$), $F_d(T) = \sum_{2n \le T} f_{2n}$.
   - Nilai literatur pembanding (dari ingatan): $p_1 = p_2 = 1$ (rekuren), $p_3 \approx 0{,}340537$ (konstanta Watson, transien).
3. **Gambler's Ruin di $[0, N]$ mulai dari $x$**: $P(\text{capai } N) = \frac{x}{N}$, $\mathbb{E}[\tau] = x(N - x)$.
4. **Distribusi koordinat akhir $X_t$ (S08)**: $\mathrm{Var}(X_t) = t / d$, kurtosis eksak $\mathrm{Kurt}(X_t) = 3 + \frac{d - 3}{t}$.

---

## 4. Kontrak 12 Tahap (`out/<bahasa>/`)

- `s01` → `s01_selftest.txt`: verifikasi 3 vektor acuan RNG + multinomial $n \le 4$.
- `s02` → `s02_parity.csv` (`d,walk,x0,x1,x2,fnv1a_hex`): 8 jalan × 200 langkah per $d$, wajib bit-identik 8 bahasa.
- `s03` → `s03_msd.csv` (`d,t,msd,se,exact`): $t \in \{1, 2, 4, \dots, 4096\}$, $W$ jalan.
- `s04` → `s04_return_sim.csv` (`d,T,p_hat,se`): $T \in \{10, 100, 1000, 10000\}$, $W$ jalan.
- `s05` → `s05_return_exact.csv` (`d,T,F_exact`): renewal eksak hingga $T = 10^4$.
- `s06` → `s06_return_compare.csv` (`d,T,p_hat,se,F_exact,z_score,p_inf_est`): perbandingan S04 vs S05 + ekstrapolasi $T^{-1/2}$ untuk $d=3$.
- `s07` → `s07_ruin.csv` (`N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z`): $N=10$ ($x \in \{2,5,8\}$), $N=50$ ($x \in \{10,25,40\}$).
- `s08` → `s08_endpoint.csv` (`d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count`): $t=1024$, ringkasan $d=1,2,3$ + histogram 40 kotak $d=1$.
- `s09` → `s09_ci.csv` (`metric,d,estimate,boot_se,ci_lo,ci_hi`): bootstrap 200 ulang (RNG sama) pada `msd_1024` dan `return_10000`.
- `s10` → `s10_fit.csv` (`d,slope,slope_se,intercept`): regresi log–log MSD vs $t$.
- `s11` → `s11_bench.csv` (`d,steps_per_rep,reps,ns_med,ns_mad`): $\ge 7$ ulang, median $\pm$ MAD.
- `s12` → `figures/{msd,return_prob,ruin}.svg` (3 SVG tulis tangan) + `s12_summary.json`.

### Beban Kerja Per Bahasa
- **C, C++, Rust, Go**: $W = 200{.}000$ (`s03`, `s04`), `ruin_W = 20.000`, `ep_W = 20.000`, `bench_steps = 10.000.000`.
- **TypeScript (`node --experimental-strip-types`)**: $W = 20{.}000$, `ruin_W = 5.000`, `ep_W = 5.000`, `bench_steps = 2.000.000`.
- **Ruby & Python**: $W = 4{.}000$, `ruin_W = 2.000`, `ep_W = 4.000`, `bench_steps = 200.000`.
- **Notebook (`ipynb`)**: $W = 2{.}000$, `ruin_W = 1.000`, `ep_W = 2.000`, `bench_steps = 100.000`.

### Kriteria A Priori (`c1`–`c10`)
- `c1`: `s01`–`s02` bit-identik di seluruh 8 bahasa.
- `c2`: MSD $|\hat{m} - t| \le 4\text{ SE}$ di $\ge 95\%$ titik tiap bahasa.
- `c3`: Peluang kembali simulasi vs renewal eksak dalam $\le 4\text{ SE}$ di $\ge 95\%$ titik tiap bahasa.
- `c4`: Gambler's ruin dalam $\le 4\text{ SE}$ di $\ge 95\%$ titik tiap bahasa.
- `c5`: Kemiringan log–log MSD $\in [0{,}97, 1{,}03]$ di semua bahasa.
- `c6`: $F_3(10^4)$ eksak dalam $0{,}005$ dari $0{,}3405$; $F_1, F_2$ monoton naik menuju $1$ dengan $F_2(10^4) < 0{,}9$.
- `c7`: `all` selesai bersih di seluruh 8 bahasa.
- `c8`: Seluruh `selftest` lulus.
- `c9`: 3 SVG tiap bahasa + 3 grafik gabungan ada, tidak kosong, dan valid XML.
- `c10`: Kecepatan (`ns/langkah`, median $\pm$ MAD) dan jumlah baris kode per bahasa dilaporkan terukur ($\le 450$ baris per program).
