package cir

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
)

func TestCIRFellerAffineRK4AndMoments(t *testing.T) {
	p := DefaultParams()
	lhs, rhs, ok := p.FellerCheck()
	if !ok || lhs <= rhs {
		t.Fatalf("Feller condition expected satisfied: lhs=%v rhs=%v", lhs, rhs)
	}

	for _, T := range []float64{1, 2, 5, 10} {
		qAnal := p.SurvivalProb(T)
		qRK4 := p.SurvivalRiccatiRK4(T, 2000)
		diff := math.Abs(qAnal - qRK4)
		if diff >= 1e-10 {
			t.Fatalf("T=%.1f analytic=%v rk4=%v diff=%v >= 1e-10", T, qAnal, qRK4, diff)
		}
	}

	// Uji momen bersyarat transisi eksak CIR (Poisson-Gamma).
	const dt = 0.25
	const n = 100000
	lam0 := 0.035
	expMean, expVar := p.TransitionMoments(lam0, dt)
	es := p.NewExactStepPrecomp(dt)
	rg := rng.New(rng.MainSeed)
	s, sq := 0.0, 0.0
	for i := 0; i < n; i++ {
		val := es.SampleNext(lam0, rg)
		s += val
		sq += val * val
	}
	m := s / n
	v := sq/n - m*m
	se := math.Sqrt(expVar / n)
	if math.Abs(m-expMean) > 4.0*se {
		t.Fatalf("exact CIR transition mean=%v want=%v se=%v", m, expMean, se)
	}
	if math.Abs(v-expVar)/expVar > 0.03 {
		t.Fatalf("exact CIR transition var=%v want=%v", v, expVar)
	}
}

func TestCIRSurvivalMCAndCDS(t *testing.T) {
	p := DefaultParams()
	rows := p.VerifySurvivalMC([]float64{1, 5, 10}, []int{20, 40}, 40000, 8, 4, rng.MainSeed)
	for _, r := range rows {
		if r.ExactCIRSE > 0 {
			z := math.Abs(r.ExactCIRMC-r.AnalyticQ) / r.ExactCIRSE
			if z > 3.5 {
				t.Fatalf("Survival Exact CIR mismatch at T=%.1f m=%d: z=%v", r.Horizon, r.SubStepsPerYear, z)
			}
		}
		if r.EulerFTSE > 0 {
			z := math.Abs(r.EulerFTMC-r.AnalyticQ) / r.EulerFTSE
			if z > 3.5 {
				t.Fatalf("Survival Euler FT mismatch at T=%.1f m=%d: z=%v", r.Horizon, r.SubStepsPerYear, z)
			}
		}
	}

	cds := p.VerifyParCDS([]float64{1, 3, 5, 10}, 0.03, 40000, 8, 4, rng.MainSeed)
	for _, c := range cds {
		if c.ZScoreProt > 3.5 {
			t.Fatalf("CDS T=%.1f PVProt MC vs anal z=%v > 3.5", c.Maturity, c.ZScoreProt)
		}
	}
}
