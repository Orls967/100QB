package rng

import (
	"math"
	"testing"
)

func TestSeedDeterminismAndMoments(t *testing.T) {
	s1 := DeriveSeed(MainSeed, 0x101, 7)
	s2 := DeriveSeed(MainSeed, 0x101, 7)
	s3 := DeriveSeed(MainSeed, 0x101, 8)
	s4 := DeriveSeed(MainSeed, 0x102, 7)
	if s1 != s2 {
		t.Fatalf("expected identical derived seed, got %d vs %d", s1, s2)
	}
	if s1 == s3 || s1 == s4 || s3 == s4 {
		t.Fatalf("expected distinct derived seeds for different index/salt")
	}

	rA := New(s1)
	rB := New(s1)
	for i := 0; i < 64; i++ {
		if math.Float64bits(rA.NormFloat64()) != math.Float64bits(rB.NormFloat64()) {
			t.Fatalf("RNG sequence diverged at step %d", i)
		}
	}

	// Uji momen empiris N(0, 1) sederhana (20.000 sampel)
	rg := New(MainSeed)
	const n = 20000
	sum, sumSq := 0.0, 0.0
	for i := 0; i < n; i++ {
		z := rg.NormFloat64()
		sum += z
		sumSq += z * z
	}
	mean := sum / n
	variance := (sumSq - sum*sum/n) / float64(n-1)
	if math.Abs(mean) > 0.04 {
		t.Fatalf("empirical N(0,1) mean too large: %v", mean)
	}
	if math.Abs(variance-1.0) > 0.05 {
		t.Fatalf("empirical N(0,1) variance too far from 1: %v", variance)
	}
}
