package options

import (
	"math"
	"testing"
)

func defaultParams() Params {
	return Params{
		S0:     100.0,
		K:      100.0,
		R:      0.05,
		Sigma:  0.20,
		T:      1.0,
		NSteps: 252,
	}
}

func TestBSPutCallParity(t *testing.T) {
	cases := []struct {
		name string
		k    float64
	}{
		{"ATM", 100.0},
		{"OTM_call", 105.0},
		{"ITM_call", 90.0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			p := defaultParams()
			p.K = tc.k
			c := BSCall(p)
			put := BSPut(p)
			rhs := p.S0 - p.K*math.Exp(-p.R*p.T)
			if math.Abs((c-put)-rhs) > 1e-10 {
				t.Fatalf("put-call parity mismatch: got %.12f, want %.12f", c-put, rhs)
			}
		})
	}
}

func TestGeometricAsianMatchesLargeMC(t *testing.T) {
	p := defaultParams()
	exact, err := GeometricAsianCall(p)
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	const nPaths = 50_000
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPaths, 4, func(c, length int) momentAcc {
		rng := chunkRNG(42, c)
		var acc momentAcc
		for range length {
			_, _, geomAvg := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
			pay := math.Max(geomAvg-p.K, 0) * disc
			acc.sum += pay
			acc.sumSq += pay * pay
		}
		return acc
	})
	mc := summarizeMoments(nPaths, sum, sumSq)
	if math.Abs(mc.Price-exact) >= 3.0*mc.SE {
		t.Fatalf("MC %.6f not within 3 SE (%.6f) of exact %.6f", mc.Price, mc.SE, exact)
	}
}

func TestControlVariateBeatsNaiveSE(t *testing.T) {
	p := defaultParams()
	naive, err := PriceAsianNaive(p, 20_000, 42)
	if err != nil {
		t.Fatalf("naive error: %v", err)
	}
	cv, err := PriceAsianControlVariate(p, 20_000, 42)
	if err != nil {
		t.Fatalf("cv error: %v", err)
	}
	if cv.SE >= naive.SE/10.0 {
		t.Fatalf("expected CV SE (%.6f) << Naive SE (%.6f)", cv.SE, naive.SE)
	}
}

func TestDeterministicSameSeedAndWorkers(t *testing.T) {
	p := defaultParams()
	p.NSteps = 50
	a, err := PriceAsianControlVariate(p, 10_000, 7, 1)
	if err != nil {
		t.Fatalf("1-worker error: %v", err)
	}
	b, err := PriceAsianControlVariate(p, 10_000, 7, 8)
	if err != nil {
		t.Fatalf("8-worker error: %v", err)
	}
	if a.Price != b.Price || a.SE != b.SE {
		t.Fatalf("non-deterministic across worker counts: 1w=%+v vs 8w=%+v", a, b)
	}
}

func TestDeltaPathwiseMatchesBump(t *testing.T) {
	p := defaultParams()
	pw, err := DeltaAsianPathwise(p, 50_000, 42)
	if err != nil {
		t.Fatalf("pathwise error: %v", err)
	}
	bmp, err := DeltaAsianBump(p, 50_000, 0.01, 42)
	if err != nil {
		t.Fatalf("bump error: %v", err)
	}
	seComb := math.Hypot(pw.SE, bmp.SE)
	if math.Abs(pw.Delta-bmp.Delta) >= 3.0*seComb {
		t.Fatalf("delta mismatch: pw=%.6f vs bump=%.6f (seComb=%.6f)", pw.Delta, bmp.Delta, seComb)
	}
}
