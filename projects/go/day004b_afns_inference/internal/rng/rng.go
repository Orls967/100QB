package rng

import "math"

// MainSeed adalah benih utama deterministik khusus sesi day004b (2026-11-22).
const MainSeed uint64 = 2026112202

// SplitMix64 menghasilkan bilangan acak 64-bit dan memajukan status.
func SplitMix64(state *uint64) uint64 {
	*state += 0x9e3779b97f4a7c15
	z := *state
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb
	return z ^ (z >> 31)
}

// DeriveSeed menurunkan benih unik secara deterministik dari benih utama, salt eksperimen,
// dan indeks replikasi atau chunk (independen terhadap jumlah goroutine).
func DeriveSeed(mainSeed, salt, index uint64) uint64 {
	x := mainSeed ^ (salt * 0x9e3779b97f4a7c15) ^ (index*0xbf58476d1ce4e5b9 + 0x94d049bb133111eb)
	_ = SplitMix64(&x)
	return SplitMix64(&x)
}

func rotl(x uint64, k int) uint64 {
	return (x << k) | (x >> (64 - k))
}

// RNG mengimplementasikan generator xoshiro256++ buatan sendiri dengan cache Box-Muller.
type RNG struct {
	s        [4]uint64
	hasSpare bool
	spare    float64
}

// New menginisialisasi RNG dari benih 64-bit melalui SplitMix64.
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

// Float64 mengembalikan bilangan acak seragam di [0, 1).
func (r *RNG) Float64() float64 {
	return float64(r.Uint64()>>11) * (1.0 / float64(uint64(1)<<53))
}

// Float64Open mengembalikan bilangan seragam ketat di (0, 1) agar Box-Muller bebas log(0).
func (r *RNG) Float64Open() float64 {
	return (float64(r.Uint64()>>11) + 0.5) * (1.0 / float64(uint64(1)<<53))
}

// UniformRange mengembalikan sampel seragam di [lo, hi).
func (r *RNG) UniformRange(lo, hi float64) float64 {
	return lo + (hi-lo)*r.Float64()
}

// NormFloat64 menghasilkan sampel normal standar N(0, 1) lewat metode Box-Muller.
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
