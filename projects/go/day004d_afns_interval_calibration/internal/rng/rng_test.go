package rng

import (
	"math"
	"testing"
)

func TestSeedDerivationBootstrapUniquenessAndMoments(t *testing.T) {
	s1 := DeriveSeed(MainSeed, 0x4D01, 5)
	s2 := DeriveSeed(MainSeed, 0x4D01, 5)
	s3 := DeriveSeed(MainSeed, 0x4D01, 6)
	if s1 != s2 || s1 == s3 {
		t.Fatalf("unexpected DeriveSeed behavior: s1=%d s2=%d s3=%d", s1, s2, s3)
	}
	b00 := DeriveBootstrapSeed(MainSeed, 0x4D02, 0, 0)
	b01 := DeriveBootstrapSeed(MainSeed, 0x4D02, 0, 1)
	b10 := DeriveBootstrapSeed(MainSeed, 0x4D02, 1, 0)
	if b00 == b01 || b00 == b10 || b01 == b10 {
		t.Fatalf("bootstrap seeds collided")
	}

	rg := New(s1)
	const n = 20000
	sum, sumSq := 0.0, 0.0
	for i := 0; i < n; i++ {
		z := rg.NormFloat64()
		sum += z
		sumSq += z * z
	}
	mean := sum / n
	variance := sumSq/n - mean*mean
	if math.Abs(mean) > 0.03 || math.Abs(variance-1.0) > 0.04 {
		t.Fatalf("unexpected normal moments: mean=%g var=%g", mean, variance)
	}
}
