package rng

import (
	"math"
	"math/bits"
)

// SplitMix64 menurunkan seed 64-bit secara deterministik.
func SplitMix64(x uint64) uint64 {
	x += 0x9e3779b97f4a7c15
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9
	x = (x ^ (x >> 27)) * 0x94d049bb133111eb
	return x ^ (x >> 31)
}

// RNG mengimplementasikan xoshiro256++ dengan cache Box-Muller.
type RNG struct {
	s        [4]uint64
	hasSpare bool
	spare    float64
}

// New membuat generator xoshiro256++ dari seed 64-bit.
func New(seed uint64) *RNG {
	var r RNG
	x := seed
	for i := 0; i < 4; i++ {
		x = SplitMix64(x)
		r.s[i] = x
	}
	return &r
}

// Uint64 menghasilkan bilangan acak uint64.
func (r *RNG) Uint64() uint64 {
	res := bits.RotateLeft64(r.s[0]+r.s[3], 23) + r.s[0]
	t := r.s[1] << 17
	r.s[2] ^= r.s[0]
	r.s[3] ^= r.s[1]
	r.s[1] ^= r.s[2]
	r.s[0] ^= r.s[3]
	r.s[2] ^= t
	r.s[3] = bits.RotateLeft64(r.s[3], 45)
	return res
}

// Float64 menghasilkan float64 seragam pada [0, 1).
func (r *RNG) Float64() float64 {
	return float64(r.Uint64()>>11) * (1.0 / float64(uint64(1)<<53))
}

// Intn menghasilkan bilangan bulat seragam pada [0, n).
func (r *RNG) Intn(n int) int {
	if n <= 1 {
		return 0
	}
	return int(r.Uint64() % uint64(n))
}

// NormFloat64 menghasilkan sampel N(0, 1) via Box-Muller.
func (r *RNG) NormFloat64() float64 {
	if r.hasSpare {
		r.hasSpare = false
		return r.spare
	}
	u1 := r.Float64()
	for u1 <= 1e-15 {
		u1 = r.Float64()
	}
	u2 := r.Float64()
	rad := math.Sqrt(-2.0 * math.Log(u1))
	theta := 2.0 * math.Pi * u2
	r.spare = rad * math.Sin(theta)
	r.hasSpare = true
	return rad * math.Cos(theta)
}

// Poisson menghasilkan sampel Poisson(lam) dengan algoritma Knuth.
func (r *RNG) Poisson(lam float64) int {
	if lam <= 0 {
		return 0
	}
	l := math.Exp(-lam)
	k := 0
	p := 1.0
	for {
		k++
		p *= r.Float64()
		if p <= l {
			return k - 1
		}
	}
}

// Perm menghasilkan permutasi acak Fisher-Yates dari 0..n-1.
func (r *RNG) Perm(n int) []int {
	p := make([]int, n)
	for i := 0; i < n; i++ {
		p[i] = i
	}
	for i := n - 1; i > 0; i-- {
		j := r.Intn(i + 1)
		p[i], p[j] = p[j], p[i]
	}
	return p
}
