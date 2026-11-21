package rng

import (
	"math"
	"testing"
)

func TestSeedDerivationBootstrapUniquenessAndMoments(t *testing.T) {
	s1 := DeriveSeed(MainSeed, 0x101, 7)
	s2 := DeriveSeed(MainSeed, 0x101, 7)
	s3 := DeriveSeed(MainSeed, 0x101, 8)
	if s1 != s2 || s1 == s3 {
		t.Fatalf("DeriveSeed invariance/uniqueness failure: %d %d %d", s1, s2, s3)
	}

	b12 := DeriveBootstrapSeed(MainSeed, 0x202, 1, 2)
	b21 := DeriveBootstrapSeed(MainSeed, 0x202, 2, 1)
	b00 := DeriveBootstrapSeed(MainSeed, 0x202, 0, 0)
	if b12 == b21 || b12 == b00 || b21 == b00 {
		t.Fatalf("expected distinct bootstrap seeds across (outer, boot): %d %d %d", b12, b21, b00)
	}
	if b12 != DeriveBootstrapSeed(MainSeed, 0x202, 1, 2) {
		t.Fatalf("expected deterministic bootstrap seed repeatability")
	}

	rA := New(b12)
	rB := New(b12)
	for i := 0; i < 64; i++ {
		if math.Float64bits(rA.NormFloat64()) != math.Float64bits(rB.NormFloat64()) {
			t.Fatalf("RNG sequence diverged at draw %d", i)
		}
	}

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
	if math.Abs(mean) > 0.04 || math.Abs(variance-1.0) > 0.05 {
		t.Fatalf("unexpected empirical N(0,1) mean=%v var=%v", mean, variance)
	}
}
