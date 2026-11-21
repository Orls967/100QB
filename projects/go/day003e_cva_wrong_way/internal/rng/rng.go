package rng

import "math"

// MainSeed adalah benih utama deterministik sesi day003e (tanggal sesi 2026-11-21).
const MainSeed uint64 = 20261121

// SplitMix64 menghasilkan bilangan acak 64-bit dan memajukan status.
func SplitMix64(state *uint64) uint64 {
	*state += 0x9e3779b97f4a7c15
	z := *state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

// DeriveChunkSeed menurunkan benih unik per chunk secara deterministik
// dari benih utama, salt eksperimen, dan indeks chunk (independen terhadap jumlah worker).
func DeriveChunkSeed(mainSeed, salt, chunkIdx uint64) uint64 {
	x := mainSeed ^ (salt * 0x9e3779b97f4a7c15) ^ (chunkIdx*0xbf58476d1ce4e5b9 + 0x94d049bb133111eb)
	_ = SplitMix64(&x)
	return SplitMix64(&x)
}

func rotl(x uint64, k int) uint64 {
	return (x << k) | (x >> (64 - k))
}

// RNG mengimplementasikan generator xoshiro256++ dengan cache Box-Muller.
type RNG struct {
	s        [4]uint64
	hasSpare bool
	spare    float64
}

// New menginisialisasi RNG dari sebuah benih 64-bit melalui SplitMix64.
func New(seed uint64) *RNG {
	st := seed
	if st == 0 {
		st = 0xdeadbeefcafebabe
	}
	var s [4]uint64
	for i := 0; i < 4; i++ {
		s[i] = SplitMix64(&st)
	}
	return &RNG{s: s}
}

// Uint64 mengembalikan integer 64-bit acak seragam.
func (r *RNG) Uint64() uint64 {
	res := rotl(r.s[0]+r.s[3], 23) + r.s[0]
	t := r.s[1] << 17

	r.s[2] ^= r.s[0]
	r.s[3] ^= r.s[1]
	r.s[1] ^= r.s[2]
	r.s[0] ^= r.s[3]

	r.s[2] ^= t
	r.s[3] = rotl(r.s[3], 45)
	return res
}

// Float64 mengembalikan bilangan titik ambang seragam di [0, 1).
func (r *RNG) Float64() float64 {
	return float64(r.Uint64()>>11) * (1.0 / float64(uint64(1)<<53))
}

// Float64Open mengembalikan bilangan seragam di (0, 1) secara ketat (menghindari 0 dan 1).
func (r *RNG) Float64Open() float64 {
	return (float64(r.Uint64()>>11) + 0.5) * (1.0 / float64(uint64(1)<<53))
}

// ExpFloat64 mengembalikan sampel eksponensial standar Exp(1).
func (r *RNG) ExpFloat64() float64 {
	return -math.Log(r.Float64Open())
}

// NormFloat64 mengembalikan sampel normal standar N(0, 1) lewat metode Box-Muller.
func (r *RNG) NormFloat64() float64 {
	if r.hasSpare {
		r.hasSpare = false
		return r.spare
	}
	u1 := r.Float64Open()
	u2 := r.Float64()
	mag := math.Sqrt(-2.0 * math.Log(u1))
	ang := 2.0 * math.Pi * u2
	r.spare = mag * math.Sin(ang)
	r.hasSpare = true
	return mag * math.Cos(ang)
}

// Gamma mengembalikan sampel distribusi Gamma(shape, scale=1) via Marsaglia-Tsang (2000)
// dan penyesuaian Ahrens-Dieter untuk shape < 1.
func (r *RNG) Gamma(shape float64) float64 {
	if shape <= 0 {
		return 0
	}
	if shape < 1.0 {
		u := r.Float64Open()
		return r.Gamma(shape+1.0) * math.Pow(u, 1.0/shape)
	}
	d := shape - 1.0/3.0
	c := 1.0 / math.Sqrt(9.0*d)
	for {
		x := r.NormFloat64()
		v := 1.0 + c*x
		if v <= 0 {
			continue
		}
		v3 := v * v * v
		u := r.Float64Open()
		x2 := x * x
		if u < 1.0-0.0331*(x2*x2) {
			return d * v3
		}
		if math.Log(u) < 0.5*x2+d*(1.0-v3+math.Log(v3)) {
			return d * v3
		}
	}
}

// Beta mengembalikan sampel distribusi Beta(a, b) via rasio dua Gamma independen.
func (r *RNG) Beta(a, b float64) float64 {
	ga := r.Gamma(a)
	gb := r.Gamma(b)
	sum := ga + gb
	if sum == 0 {
		return 0.5
	}
	return ga / sum
}

func (r *RNG) binomial(n int, p float64) int {
	k := 0
	for n > 0 {
		if p <= 0 {
			return k
		}
		if p >= 1 {
			return k + n
		}
		if n <= 15 {
			for j := 0; j < n; j++ {
				if r.Float64() < p {
					k++
				}
			}
			return k
		}
		// Pemecahan Beta-Binomial eksak Ahrens-Dieter: orde ke-i dari n sampel U(0,1)
		// berdistribusi Beta(i, n + 1 - i).
		i := (n + 1) / 2
		y := r.Beta(float64(i), float64(n+1-i))
		if p >= y {
			k += i
			n -= i
			p = (p - y) / (1.0 - y)
		} else {
			n = i - 1
			p = p / y
		}
	}
	return k
}

// Poisson mengembalikan sampel distribusi Poisson(lam) secara eksak:
// algoritma Knuth untuk lam < 15, dan dekomposisi Gamma + Binomial Ahrens-Dieter untuk lam >= 15.
func (r *RNG) Poisson(lam float64) int {
	k := 0
	for lam > 0 {
		if lam < 15.0 {
			limit := math.Exp(-lam)
			prod := 1.0
			cnt := 0
			for {
				cnt++
				prod *= r.Float64Open()
				if prod <= limit {
					break
				}
			}
			return k + cnt - 1
		}
		m := int(0.875 * lam)
		if m < 1 {
			m = 1
		}
		x := r.Gamma(float64(m))
		if x <= lam {
			k += m
			lam -= x
		} else {
			return k + r.binomial(m-1, lam/x)
		}
	}
	return k
}
