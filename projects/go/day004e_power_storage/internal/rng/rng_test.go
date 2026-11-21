package rng

import (
	"math"
	"testing"
)

func TestDeterminismAndDerivedSeeds(t *testing.T) {
	s1 := DeriveSeed(MainSeed, 0x101, 42)
	s2 := DeriveSeed(MainSeed, 0x101, 42)
	s3 := DeriveSeed(MainSeed, 0x101, 43)
	if s1 != s2 || s1 == s3 {
		t.Fatalf("unexpected DeriveSeed results: s1=%d s2=%d s3=%d", s1, s2, s3)
	}
	st1 := DeriveStateSeed(KernelSeed, 0x201, 12, 1, 5)
	st2 := DeriveStateSeed(KernelSeed, 0x201, 12, 2, 5)
	if st1 == st2 {
		t.Fatalf("expected distinct state seeds for different regimes")
	}

	rA := New(s1)
	rB := New(s2)
	for i := 0; i < 64; i++ {
		if rA.Uint64() != rB.Uint64() {
			t.Fatalf("RNG sequence diverged at step %d", i)
		}
	}
}

func TestDistributionMoments(t *testing.T) {
	const n = 200000
	rg := New(DeriveSeed(VerifySeed, 0x301, 0))

	// Normal(0, 1)
	var sumN, sumN2 float64
	for i := 0; i < n; i++ {
		z := rg.NormFloat64()
		sumN += z
		sumN2 += z * z
	}
	meanN := sumN / n
	varN := (sumN2 - n*meanN*meanN) / (n - 1)
	if math.Abs(meanN) > 4.0/math.Sqrt(n) {
		t.Fatalf("normal mean out of 4 SE: %.5f", meanN)
	}
	if math.Abs(varN-1.0) > 4.0*math.Sqrt(2.0/n) {
		t.Fatalf("normal variance out of 4 SE: %.5f", varN)
	}

	// Exp(mean = 120)
	const expMean = 120.0
	var sumE, sumE2 float64
	for i := 0; i < n; i++ {
		e := rg.ExpFloat64(expMean)
		if e < 0 {
			t.Fatalf("negative exponential sample: %f", e)
		}
		sumE += e
		sumE2 += e * e
	}
	meanE := sumE / n
	varE := (sumE2 - n*meanE*meanE) / (n - 1)
	if math.Abs(meanE-expMean) > 4.0*expMean/math.Sqrt(n) {
		t.Fatalf("exp mean out of 4 SE: %.5f", meanE)
	}
	theoVarE := expMean * expMean
	seVarE := theoVarE * math.Sqrt(8.0/n)
	if math.Abs(varE-theoVarE) > 4.0*seVarE {
		t.Fatalf("exp variance out of 4 SE: %.2f vs %.2f", varE, theoVarE)
	}

	// Poisson(lambda = 0.08)
	const lam = 0.08
	var sumP, sumP2 float64
	for i := 0; i < n; i++ {
		k := float64(rg.Poisson(lam))
		sumP += k
		sumP2 += k * k
	}
	meanP := sumP / n
	varP := (sumP2 - n*meanP*meanP) / (n - 1)
	if math.Abs(meanP-lam) > 4.0*math.Sqrt(lam/n) {
		t.Fatalf("poisson mean out of 4 SE: %.6f vs %.6f", meanP, lam)
	}
	if math.Abs(varP-lam) > 5.0*math.Sqrt((lam+2.0*lam*lam)/n) {
		t.Fatalf("poisson variance out of tolerance: %.6f vs %.6f", varP, lam)
	}
}
