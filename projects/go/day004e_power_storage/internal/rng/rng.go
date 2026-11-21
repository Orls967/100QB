package rng

import "math"

const (
	// MainSeed adalah benih utama deterministik tetap untuk day004e (2026-11-22).
	MainSeed uint64 = 2026112205
	// KernelSeed memisahkan generator kernel transisi DP dari lintasan evaluasi.
	KernelSeed uint64 = 2026112215
	// EvalSeed memisahkan generator lintasan evaluasi N=2.000 hari dari kernel DP.
	EvalSeed uint64 = 2026112225
	// VerifySeed digunakan untuk pengujian momen statistik model.
	VerifySeed uint64 = 2026112235
)

// SplitMix64 menghasilkan bilangan acak 64-bit dan memajukan status.
func SplitMix64(state *uint64) uint64 {
	*state += 0x9e3779b97f4a7c15
	z := *state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

// DeriveSeed menurunkan benih deterministik per (salt, index) yang bebas urutan worker goroutine.
func DeriveSeed(baseSeed, salt, index uint64) uint64 {
	x := baseSeed ^ (salt * 0x9e3779b97f4a7c15) ^ (index*0xbf58476d1ce4e5b9 + 0x94d049bb133111eb)
	_ = SplitMix64(&x)
	return SplitMix64(&x)
}

// DeriveStateSeed menurunkan benih deterministik unik untuk state DP (hour, regime, priceBin).
func DeriveStateSeed(baseSeed, salt, hour, regime, priceBin uint64) uint64 {
	x := baseSeed ^ (salt * 0x9e3779b97f4a7c15)
	x ^= (hour + 1) * 0xbf58476d1ce4e5b9
	x ^= (regime + 1) * 0x94d049bb133111eb
	x ^= (priceBin + 1) * 0xd6e8feb86659fd93
	_ = SplitMix64(&x)
	return SplitMix64(&x)
}

func rotl(x uint64, k int) uint64 {
	return (x << k) | (x >> (64 - k))
}

// RNG mengimplementasikan pembangkit xoshiro256++ buatan sendiri dengan cache Box-Muller.
type RNG struct {
	s        [4]uint64
	hasSpare bool
	spare    float64
}

// New menginisialisasi RNG dari benih 64-bit lewat SplitMix64.
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

// Uint64 mengembalikan integer 64-bit seragam.
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

// Float64 mengembalikan bilangan seragam di [0, 1).
func (r *RNG) Float64() float64 {
	return float64(r.Uint64()>>11) * (1.0 / float64(uint64(1)<<53))
}

// Float64Open mengembalikan bilangan seragam ketat di (0, 1) agar log(u) aman.
func (r *RNG) Float64Open() float64 {
	return (float64(r.Uint64()>>11) + 0.5) * (1.0 / float64(uint64(1)<<53))
}

// NormFloat64 menghasilkan sampel normal standar N(0, 1) via transformasi Box-Muller.
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

// ExpFloat64 menghasilkan sampel eksponensial dengan rata-rata mean (>= 0).
func (r *RNG) ExpFloat64(mean float64) float64 {
	if mean <= 0 {
		return 0
	}
	return -mean * math.Log(r.Float64Open())
}

// Poisson menghasilkan sampel cacah Poisson(lambda) secara eksak (algoritma Knuth untuk lambda kecil).
func (r *RNG) Poisson(lambda float64) int {
	if lambda <= 0 {
		return 0
	}
	if lambda < 30.0 {
		limit := math.Exp(-lambda)
		k := 0
		prod := 1.0
		for {
			k++
			prod *= r.Float64Open()
			if prod <= limit {
				return k - 1
			}
		}
	}
	// Pendekatan transformasi penolakan Atkinson untuk lambda besar bila diperlukan.
	c := 0.767 - 3.36/lambda
	beta := math.Pi / math.Sqrt(3.0*lambda)
	alpha := beta * lambda
	kConst := math.Log(c) - lambda - math.Log(beta)
	for {
		u := r.Float64Open()
		x := (alpha - math.Log((1.0-u)/u)) / beta
		if x < -0.5 {
			continue
		}
		n := int(math.Floor(x + 0.5))
		v := r.Float64Open()
		lhs := alpha - beta*x + math.Log(v/math.Pow(1.0+math.Exp(alpha-beta*x), 2))
		rhs := kConst + float64(n)*math.Log(lambda) - lgammaInt(n+1)
		if lhs <= rhs {
			return n
		}
	}
}

func lgammaInt(n int) float64 {
	v, _ := math.Lgamma(float64(n))
	return v
}
