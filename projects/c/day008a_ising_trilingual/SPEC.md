# Spesifikasi Teknis: Model Ising 2D Lintas Tiga Bahasa (`day008a`)

Dokumen ini adalah kontrak tunggal (*single source of truth*) yang diikuti secara literal oleh implementasi **C11**, **C++20**, dan **Python** untuk menjamin **paritas bit-identik** pada beban uji paritas serta konsistensi analisis fisika statistik pada beban penuh.

---

## 1. Model Fisika & Konvensi Kisi

- **Geometri**: Kisi persegi 2D berukuran $L \times L$ dengan syarat batas periodik (*periodic boundary conditions*, PBC), jumlah situs $N = L^2$.
- **Variabel Spin**: $s_i \in \{+1, -1\}$ bertipe `int8` (`int8_t` di C/C++, `int` dalam daftar/array di Python).
- **Hamiltonian**: $J = 1$, medan luar $H = 0$, suhu invers $\beta = 1 / T$:
  $$E = -\sum_{y=0}^{L-1}\sum_{x=0}^{L-1} s_{y,x}\bigl(s_{y,(x+1)\bmod L} + s_{(y+1)\bmod L, x}\bigr) \in [-2N, 2N] \cap \mathbb{Z}$$
  $$M = \sum_{i=0}^{N-1} s_i \in [-N, N] \cap \mathbb{Z}$$
- **Pengindeksan Situs**: Urutan baris-utama (*row-major*):
  $$i = y \cdot L + x, \qquad x \in \{0, \dots, L-1\},\; y \in \{0, \dots, L-1\}$$
- **Tetangga Periodik** (urutan kanonik: **kanan, atas, kiri, bawah**):
  - `right(x, y)` $= y \cdot L + ((x + 1) \bmod L)$
  - `up(x, y)` $= ((y + 1) \bmod L) \cdot L + x$
  - `left(x, y)` $= y \cdot L + ((x + L - 1) \bmod L)$
  - `down(x, y)` $= ((y + L - 1) \bmod L) \cdot L + x$
- **Inisialisasi**:
  - **Beban Paritas**: Mulai dingin (*cold start*), seluruh $s_i = +1$.
  - **Beban Penuh**: Mulai panas (*hot start*), tiap $s_i$ ditarik acak seragam dari $\{+1, -1\}$ menggunakan `u01 < 0.5 ? 1 : -1` pada urutan indeks $i = 0, \dots, N-1$ memakai RNG tugas yang bersangkutan sebelum loop termalisasi.

---

## 2. Generator Bilangan Acak (RNG) Identik Lintas Bahasa

Seluruh bahasa menggunakan aritmetika unsigned 64-bit eksak (`uint64_t` di C/C++, `& 0xFFFFFFFFFFFFFFFF` di Python) — **dilarang** memakai `rand()`, `std::mt19937`, atau `numpy.random` pada simulasi.

### 2.1 SplitMix64 (Inisialisasi Benih)
```text
splitmix64_next(state):
    state = (state + 0x9E3779B97F4A7C15) mod 2^64
    z = state
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) mod 2^64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) mod 2^64
    return z ^ (z >> 31)
```
Fungsi pembaur (*hash*) deterministik `splitmix64_mix(x)`: menjalankan transformasi `z` di atas tanpa mengubah state eksternal.

### 2.2 xoshiro256++
State 4 kata 64-bit `s[0..3]` diisi dengan 4 panggilan berurutan `splitmix64_next(&seed)`.
```text
rotl(x, k) = ((x << k) | (x >> (64 - k))) mod 2^64

xoshiro256pp_next(s):
    result = (rotl((s[0] + s[3]) mod 2^64, 23) + s[0]) mod 2^64
    t = (s[1] << 17) mod 2^64
    s[2] ^= s[0]
    s[3] ^= s[1]
    s[1] ^= s[2]
    s[0] ^= s[3]
    s[2] ^= t
    s[3] = rotl(s[3], 45)
    return result
```

### 2.3 Konversi ke `u01` dan Bilangan Bulat `[0, n)`
- `u01`:
  $$u_{01} = (x \gg 11) \times 2^{-53} \in [0, 1)$$
  Di C/C++: `(double)(x >> 11) * 0x1.0p-53`. Di Python: `(x >> 11) * (2.0 ** -53)`.
- Bilangan bulat dalam $[0, n)$ (tanpa operator modulo `%`):
  $$k = \lfloor u_{01} \cdot n \rfloor$$
  Dengan penjagaan batas aman IEEE-754 bila pembulatan menghasilkan $k \ge n \implies k = n - 1$.

### 2.4 Skema Benih Per-Tugas
- `MASTER_SEED = 20261126` (`uint64_t`).
- ID Algoritma: `0 = Metropolis`, `1 = Wolff`, `2 = Swendsen-Wang`.
- Fungsi hash tugas:
  ```text
  task_hash(L, T_index, algo_id, replica):
      h = 0x517CC1B727220A95
      h = splitmix64_mix(h ^ uint64(L))
      h = splitmix64_mix(h ^ uint64(T_index))
      h = splitmix64_mix(h ^ uint64(algo_id))
      h = splitmix64_mix(h ^ uint64(replica))
      return h

  task_seed = splitmix64_mix(MASTER_SEED ^ task_hash(L, T_index, algo_id, replica))
  ```
Karena setiap tugas `(L, T_index, algo_id, replica)` memiliki benih independen yang hanya bergantung pada koordinat tugasnya, hasil eksekusi paralel (`pthread` / `std::thread`) untuk `1, 2, 4, N` *thread* dijamin **bit-identik**.

---

## 3. Urutan Operasi Algoritma

### 3.1 Metropolis (Sapuan Sekuensial Leksikografis)
- Tabel penerimaan dihitung sekali sebelum loop:
  $$p_4 = \exp(-4\beta), \qquad p_8 = \exp(-8\beta)$$
  Menggunakan `exp()` pustaka standar (`libm`). Tersedia pula mode opsional `--prob-file <path>` yang membaca `(p_4, p_8)` atau `p_add` dalam format teks heksadesimal IEEE-754 (`%a`) untuk mengisolasi perbedaan 1-ulp `exp()` lintas bahasa bila terjadi.
- Satu sapuan (*sweep*) mengunjungi $i = 0, 1, \dots, N-1$ berurutan:
  1. Hitung jumlah tetangga `sum_nn = s[right(i)] + s[up(i)] + s[left(i)] + s[down(i)]`.
  2. $\Delta E = 2 \cdot s_i \cdot \text{sum\_nn} \in \{-8, -4, 0, 4, 8\}$.
  3. **Selalu tarik tepat satu `u01`** per situs (`u = rng_u01()`), bahkan ketika $\Delta E \le 0$.
  4. Terima pembalikan $s_i \leftarrow -s_i$ jika $\Delta E \le 0$, atau $\Delta E = 4 \land u < p_4$, atau $\Delta E = 8 \land u < p_8$.

### 3.2 Algoritma Klaster Wolff
- Probabilitas aktivasi ikatan dihitung sekali: $p_{\text{add}} = 1 - \exp(-2\beta)$.
- Satu langkah klaster Wolff:
  1. Tarik satu `u01` untuk memilih situs benih: $i_0 = \lfloor u_{01} \cdot N \rfloor$.
  2. Catat `old_spin = s[i0]`, balik `s[i0] = -old_spin`, dorong `i0` ke tumpukan LIFO eksplisit.
  3. Selama tumpukan tidak kosong:
     - Pop situs `curr` dari puncak tumpukan.
     - Periksa 4 tetangga dalam urutan **(kanan, atas, kiri, bawah)**.
     - Untuk tiap tetangga `nb`: jika `s[nb] == old_spin` (searah dan belum masuk klaster), tarik satu `u01`; jika $u < p_{\text{add}}$, balik `s[nb] = -old_spin` dan dorong `nb` ke tumpukan.

### 3.3 Algoritma Klaster Swendsen–Wang
- Probabilitas aktivasi ikatan dihitung sekali: $p_{\text{add}} = 1 - \exp(-2\beta)$.
- Struktur **Union-Find** berukuran $N$: array `parent[0..N-1]` diinisialisasi `-1` (`parent[r] < 0` berarti akar dengan ukuran `-parent[r]`).
  - `find(x)`: *path halving* (`while (parent[x] >= 0) { if (parent[parent[x]] >= 0) parent[x] = parent[parent[x]]; x = parent[x]; }`).
  - `unite(a, b)`: *union by size* (gabungkan pohon lebih kecil ke pohon lebih besar; jika sama, jadikan akar `a` sebagai induk).
- Satu langkah Swendsen–Wang:
  1. Reset `parent[0..N-1] = -1`.
  2. Untuk $i = 0, 1, \dots, N-1$, periksa dua ikatan ke arah **(kanan, lalu atas)**:
     - Jika `s[i] == s[right(i)]`, tarik satu `u01`; jika $u < p_{\text{add}}$, `unite(i, right(i))`.
     - Jika `s[i] == s[up(i)]`, tarik satu `u01`; jika $u < p_{\text{add}}$, `unite(i, up(i))`.
  3. Untuk $i = 0, 1, \dots, N-1$: jika `parent[i] < 0` (indeks $i$ adalah akar klaster), tarik satu `u01` dan tetapkan `flip_root[i] = (u < 0.5) ? 1 : 0`.
  4. Untuk $i = 0, 1, \dots, N-1$: cari `r = find(i)`; jika `flip_root[r]`, balik `s[i] = -s[i]`.

---

## 4. Observabel, Hash Konfigurasi & Format CSV

- Observabel per langkah dihitung dalam aritmetika bilangan bulat eksak `(E_sum_int, M_sum_int)`, lalu dibagi $N$:
  $$e = \frac{E_{\text{sum}}}{N}, \qquad m = \frac{M_{\text{sum}}}{N}, \qquad |m|, \qquad m^2, \qquad m^4.$$
- **Hash FNV-1a 64-bit** atas konfigurasi spin akhir $s_0, \dots, s_{N-1}$ (diinterpretasikan sebagai byte tak bertanda `(uint8_t)s[i]`):
  - Basis offset: `0xCBF29CE484222325`, prima: `0x100000001B3`.
  - Untuk tiap $i = 0 \dots N-1$: `h ^= (uint8_t)s[i]; h = (h * 0x100000001B3) mod 2^64`.
- **Format CSV**: Pemisah koma `,`, header tetap, akhiran baris `\n`, bilangan pecahan diformat `%.17g`.

---

## 5. Beban Uji, Kriteria A Priori & Anchor Eksak

### 5.1 Beban Paritas (Wajib, Identik 100% di C, C++, Python)
- Parameter: $L = 8$, $\beta = 0{,}44$, mulai dingin ($s_i = +1$), `3000` sapuan (Metropolis) atau `3000` langkah klaster (Wolff, SW).
- Benih: `task_seed(8, 0, algo_id, 0)` (satu benih per algoritma).
- Keluaran: deret `(step, E_sum_int, M_sum_int)` dan `final_fnv1a_hex`.

### 5.2 Beban Penuh
- Ukuran kisi:
  - C & C++: $L \in \{8, 16, 32, 64, 128\}$.
  - Notebook Python: $L \in \{8, 16, 32\}$ penuh, $L = 64$ dikurangi (faktor pengurangan dicatat di `results.json`), $L = 128$ tidak dijalankan di notebook (dicatat eksplisit).
- Grid suhu (37 titik): 25 titik seragam $T \in [1{,}8, 3{,}4]$ digabung dan diurutkan dengan 12 titik rapat $T \in [2{,}15, 2{,}40]$ (total 37 nilai $T$ unik).
- Jumlah langkah per titik $(L, T)$:
  - Uji pilot memperkirakan $\tau_{\text{int}}$ pada $T_c$: termalisasi $\ge 20\,\tau_{\text{int}}$ (minimal `2000` langkah Wolff/SW, `8192` sapuan Metropolis).
  - Pengukuran: $\ge \max(50\,\tau_{\text{int}}, 2^{17})$ untuk Wolff/SW (disesuaikan untuk $L=128$ agar total runtime C/C++ $\le 15$ menit, hitungan nyata dicatat di `results.json`), dan $\ge 2^{18}$ sapuan untuk Metropolis pada $L \le 64$ (disesuaikan untuk $L=128$ jika diperlukan, dicatat di `results.json`).

### 5.3 Anchor Eksak
1. **Enumerasi penuh** $L = 4$ ($2^{16} = 65{.}536$ keadaan, 3 bahasa) dan $L = 5$ ($2^{25} = 33{.}554{.}432$ keadaan, C & C++) pada $\beta \in \{0{,}2, 0{,}35, 0{,}4406867935097715, 0{,}5, 0{,}7\}$: fungsi partisi $Z$, $\langle E \rangle / N$, $\langle m^2 \rangle$.
2. **Suhu kritis eksak**: $T_c = \frac{2}{\ln(1 + \sqrt{2})} \approx 2{,}269185314213022$.
3. **Energi Onsager $L \to \infty$**:
   $$e_{\text{Onsager}}(T) = -\coth(2\beta)\left[1 + \frac{2}{\pi}(2\tanh^2(2\beta) - 1) K(k_1)\right], \qquad k_1 = \frac{2\sinh(2\beta)}{\cosh^2(2\beta)}$$
   dengan $K(k_1)$ dihitung lewat **Arithmetic-Geometric Mean (AGM)** $\frac{\pi}{2\,\text{AGM}(1, \sqrt{1 - k_1^2})}$ pada $T \in \{2{,}0, 2{,}6, 3{,}2\}$, serta $e(T_c) = -\sqrt{2}$.
4. **Magnetisasi spontan Yang ($T < T_c$)**:
   $$m_{\text{Yang}}(T) = \left(1 - \sinh^{-4}(2\beta)\right)^{1/8}$$
   pada $T = 2{,}0$ (vs $L = 128$) dan $T = 2{,}1$ (vs $L = 64$).
5. **Konstanta pembanding dari ingatan literatur** (hanya pembanding longgar): $U_4^\star \approx 0{,}6107$, $\gamma/\nu = 7/4 = 1{,}75$, $\nu = 1$, $\beta/\nu = 1/8 = 0{,}125$, $z_{\text{Metropolis}} \approx 2{,}17$, $z_{\text{Wolff}} \approx 0{,}25$, $z_{\text{SW}} \approx 0{,}2$.

### 5.4 Kriteria A Priori (`c1`–`c14`)
- `c1`: Paritas bit-identik 3 bahasa pada 3 algoritma (0 selisih pada 3000 langkah + hash FNV-1a akhir identik).
- `c2`: Enumerasi $L=4$ (3 bahasa) & $L=5$ (C/C++) selisih relatif $\le 10^{-12}$, dan simulasi $L=4$ cocok eksak dalam $\le 4\,\text{SE}$.
- `c3`: Hasil simulasi C dan C++ identik bit-demi-bit untuk `1, 2, 4, N` *thread*.
- `c4`: Energi $L=128$ di $T \in \{2{,}0, 2{,}6, 3{,}2\}$ cocok Onsager dalam $\le 4\,\text{SE}$ (dengan catatan bias ukuran-hingga).
- `c5`: $|m|$ pada $T = 2{,}0$ ($L=128$) dan $T = 2{,}1$ ($L=64$) cocok Yang dalam $\le 4\,\text{SE}$ (dengan catatan bias ukuran-hingga).
- `c6`: Titik silang Binder $U_4(L, 2L)$ berada dalam $|T - T_c| \le 0{,}01$ untuk pasangan $L \ge 16$.
- `c7`: $|\hat{\gamma}/\hat{\nu} - 1{,}75| \le 0{,}10$.
- `c8`: $|\hat{\beta}/\hat{\nu} - 0{,}125| \le 0{,}03$.
- `c9`: $z_{\text{Wolff}} < z_{\text{Metropolis}} - 1{,}0$ dan $\tilde{z}_{\text{Metropolis}} \in [1{,}8, 2{,}5]$.
- `c10`: Notebook Python berjalan bersih dua kali dari nol dengan angka identik.
- `c11`: Seluruh tes (`make test` C, `make test` C++, dan sel `assert` notebook) lulus.
- `c12`: Seluruh 8 grafik `outputs/figures/day008a_*` ada, tidak kosong, dan SVG valid XML.
- `c13`: Ukuran kode dilaporkan dan setiap berkas sumber $\le 600$ baris setelah format.
- `c14`: Perbandingan kecepatan 3 bahasa (ns/pembaruan-spin, median $\pm$ MAD $\ge 7$ ulang) dilaporkan tanpa klaim tak terukur.
