package rng

import (
	"math"
	"testing"
)

func TestSamplersMoments(t *testing.T) {
	rg := New(MainSeed)
	const n = 100000

	// Eksponensial(rate=2), Erlang(2, rate=2), Weibull(0.5, 0.5), Lognormal(0, 0.5)
	cases := []struct {
		name     string
		theoMean float64
		theoVar  float64
		draw     func() float64
	}{
		{"Exp(2)", 0.5, 0.25, func() float64 { return rg.ExpRate(2.0) }},
		{"Erlang(2,2)", 1.0, 0.5, func() float64 { return rg.Erlang(2, 2.0) }},
		{"Pareto(2.5,1.5)", 1.0, 5.0, func() float64 { return rg.ParetoLomax(2.5, 1.5) }},
		{"Weibull(0.5,0.5)", 1.0, 5.0, func() float64 { return rg.Weibull(0.5, 0.5) }},
		{"Lognormal(0,0.5)", math.Exp(0.125), (math.Exp(0.25) - 1.0) * math.Exp(0.25), func() float64 { return rg.Lognormal(0, 0.5) }},
	}

	for _, tc := range cases {
		s, sq := 0.0, 0.0
		for i := 0; i < n; i++ {
			x := tc.draw()
			s += x
			sq += x * x
		}
		m := s / n
		se := math.Sqrt(tc.theoVar / n)
		if math.Abs(m-tc.theoMean) > 4.0*se {
			t.Fatalf("%s mean=%v want=%v se=%v", tc.name, m, tc.theoMean, se)
		}
	}

	// Poisson(10) dan Poisson(120)
	for _, lam := range []float64{10.0, 120.0} {
		s, sq := 0.0, 0.0
		for i := 0; i < 40000; i++ {
			k := float64(rg.Poisson(lam))
			s += k
			sq += k * k
		}
		m := s / 40000
		se := math.Sqrt(lam / 40000)
		if math.Abs(m-lam) > 4.0*se {
			t.Fatalf("Poisson(%.0f) mean=%v se=%v", lam, m, se)
		}
	}
}
