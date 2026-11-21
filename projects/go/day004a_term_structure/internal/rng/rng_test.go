package rng

import (
	"math"
	"testing"
)

func TestXoshiroAndNormalMoments(t *testing.T) {
	rg := New(DeriveChunkSeed(MainSeed, 0x101, 0))
	const n = 200_000
	sumU, sumU2 := 0.0, 0.0
	sumZ, sumZ2, sumZ4 := 0.0, 0.0, 0.0
	for i := 0; i < n; i++ {
		u := rg.Float64()
		if u < 0 || u >= 1 {
			t.Fatalf("Float64 out of [0,1): %v", u)
		}
		sumU += u
		sumU2 += u * u

		z := rg.NormFloat64()
		sumZ += z
		z2 := z * z
		sumZ2 += z2
		sumZ4 += z2 * z2
	}
	fn := float64(n)
	meanU := sumU / fn
	varU := (sumU2 - sumU*sumU/fn) / (fn - 1.0)
	if math.Abs(meanU-0.5) > 3.0*math.Sqrt(1.0/(12.0*fn)) {
		t.Fatalf("Uniform mean out of 3 SE: got %v", meanU)
	}
	if math.Abs(varU-1.0/12.0) > 5e-3 {
		t.Fatalf("Uniform var out of tol: got %v", varU)
	}

	meanZ := sumZ / fn
	varZ := (sumZ2 - sumZ*sumZ/fn) / (fn - 1.0)
	kurtZ := (sumZ4 / fn) / (varZ * varZ)
	if math.Abs(meanZ) > 3.0/math.Sqrt(fn) {
		t.Fatalf("Normal mean out of 3 SE: got %v", meanZ)
	}
	if math.Abs(varZ-1.0) > 3.0*math.Sqrt(2.0/fn) {
		t.Fatalf("Normal var out of 3 SE: got %v", varZ)
	}
	if math.Abs(kurtZ-3.0) > 0.05 {
		t.Fatalf("Normal kurtosis far from 3: got %v", kurtZ)
	}

	// Determinisme bit-identik
	rgA := New(DeriveChunkSeed(MainSeed, 0x999, 7))
	rgB := New(DeriveChunkSeed(MainSeed, 0x999, 7))
	for i := 0; i < 100; i++ {
		if math.Float64bits(rgA.NormFloat64()) != math.Float64bits(rgB.NormFloat64()) {
			t.Fatalf("RNG not deterministic at step %d", i)
		}
	}
}
