package rng

import (
	"math"
	"testing"
)

func TestNormalAndGammaMoments(t *testing.T) {
	rg := New(MainSeed)
	const n = 100000
	sumZ, sqZ := 0.0, 0.0
	for i := 0; i < n; i++ {
		z := rg.NormFloat64()
		sumZ += z
		sqZ += z * z
	}
	meanZ := sumZ / n
	varZ := sqZ/n - meanZ*meanZ
	if math.Abs(meanZ) > 0.01 {
		t.Fatalf("normal mean out of tolerance: %v", meanZ)
	}
	if math.Abs(varZ-1.0) > 0.015 {
		t.Fatalf("normal variance out of tolerance: %v", varZ)
	}

	// Uji Gamma(shape=0.3) < 1 dan Gamma(shape=4.5) >= 1
	for _, shape := range []float64{0.3, 4.5} {
		s, sq := 0.0, 0.0
		for i := 0; i < n; i++ {
			g := rg.Gamma(shape)
			if g < 0 || math.IsNaN(g) {
				t.Fatalf("invalid gamma sample: %v", g)
			}
			s += g
			sq += g * g
		}
		m := s / n
		v := sq/n - m*m
		seM := math.Sqrt(shape / n)
		if math.Abs(m-shape) > 4.0*seM {
			t.Fatalf("gamma(%.1f) mean=%v, expected=%.1f", shape, m, shape)
		}
		if math.Abs(v-shape)/shape > 0.04 {
			t.Fatalf("gamma(%.1f) var=%v, expected=%.1f", shape, v, shape)
		}
	}
}

func TestPoissonMoments(t *testing.T) {
	rg := New(MainSeed ^ 12345)
	const n = 60000
	for _, lam := range []float64{0.5, 5.0, 25.0, 200.0} {
		s, sq := 0.0, 0.0
		for i := 0; i < n; i++ {
			k := rg.Poisson(lam)
			if k < 0 {
				t.Fatalf("negative poisson: %d", k)
			}
			fk := float64(k)
			s += fk
			sq += fk * fk
		}
		m := s / n
		v := sq/n - m*m
		se := math.Sqrt(lam / n)
		if math.Abs(m-lam) > 4.0*se {
			t.Fatalf("poisson(%.1f) mean=%v (SE=%v)", lam, m, se)
		}
		if math.Abs(v-lam)/lam > 0.04 {
			t.Fatalf("poisson(%.1f) var=%v", lam, v)
		}
	}
}
