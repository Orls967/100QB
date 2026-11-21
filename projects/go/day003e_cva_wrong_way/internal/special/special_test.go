package special

import (
	"math"
	"testing"
)

func TestGaussHermiteStandardNormalMoments(t *testing.T) {
	nodes := GaussHermiteStandardNormal(32)
	if len(nodes) != 32 {
		t.Fatalf("expected 32 nodes, got %d", len(nodes))
	}
	wSum := 0.0
	for _, nd := range nodes {
		wSum += nd.W
	}
	if math.Abs(wSum-1.0) > 1e-14 {
		t.Fatalf("weights do not sum to 1: %.17g", wSum)
	}

	// Eksak untuk polinomial hingga derajat 2n-1 = 63; cek sampai derajat 20.
	// Besaran skala ekspektasi mutlak E[|Z|^k] digunakan agar uji momen ganjil bebas dari ilusi pembatalan nol.
	doubleFact := 1.0
	absScale := 1.0
	for k := 0; k <= 20; k++ {
		moment := 0.0
		for _, nd := range nodes {
			moment += nd.W * math.Pow(nd.X, float64(k))
		}
		expected := 0.0
		if k == 0 {
			expected = 1.0
			absScale = 1.0
		} else if k%2 == 0 {
			doubleFact *= float64(k - 1)
			expected = doubleFact
			absScale = doubleFact
		} else {
			absScale *= float64(k) * 0.8
		}
		tol := 1e-12 * (1.0 + absScale)
		if math.Abs(moment-expected) > tol {
			t.Fatalf("moment k=%d got %.17g want %.17g (tol=%.3g)", k, moment, expected, tol)
		}
	}
}

func TestNormQuantileRoundTrip(t *testing.T) {
	for _, p := range []float64{1e-10, 1e-5, 0.01, 0.25, 0.5, 0.95, 0.99, 1 - 1e-10} {
		x := NormQuantile(p)
		got := NormCDF(x)
		if math.Abs(got-p) > 1e-14 {
			t.Fatalf("p=%v x=%v roundtrip=%v", p, x, got)
		}
	}
}
