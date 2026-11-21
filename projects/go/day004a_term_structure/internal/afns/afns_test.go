package afns

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
)

func TestRiccatiClosedFormVsRK4AndNSLoadings(t *testing.T) {
	p := DefaultParams()
	taus := []float64{0.01, 0.05, 0.25, 1.0, 5.0, 10.0, 20.0, 30.0}
	prevAdj := 0.0
	maxBErr, maxAErr, maxAdjErr := 0.0, 0.0, 0.0

	for _, tau := range taus {
		h := p.YieldLoadings(tau)
		if math.Abs(h[0]-1.0) > 1e-15 || math.Abs(h[1]-ns.G1(tau, p.Lambda)) > 1e-15 || math.Abs(h[2]-ns.G2(tau, p.Lambda)) > 1e-15 {
			t.Fatalf("YieldLoadings mismatch vs NS at tau=%g", tau)
		}

		bCl := p.BVec(tau)
		adjCl, aCl := p.YieldAdjustment(tau)
		if adjCl >= 0 {
			t.Fatalf("expected negative yield adjustment at tau=%g, got %g", tau, adjCl)
		}
		if adjCl >= prevAdj {
			t.Fatalf("expected |adjustment| strictly increasing with tau=%g: got %g vs prev %g", tau, adjCl, prevAdj)
		}
		prevAdj = adjCl

		steps := int(math.Ceil(tau * 400))
		if steps < 400 {
			steps = 400
		}
		bRK, aRK, bRich, aRich := p.SolveRiccatiRK4(tau, steps)
		for i := 0; i < 3; i++ {
			if d := math.Abs(bCl[i] - bRich[i]); d > maxBErr {
				maxBErr = d
			}
			if d := math.Abs(bCl[i] - bRK[i]); d > maxBErr {
				maxBErr = d
			}
		}
		if d := math.Abs(aCl - aRK); d > maxAErr {
			maxAErr = d
		}
		if d := math.Abs(aCl - aRich); d > maxAErr {
			maxAErr = d
		}
		if d := math.Abs(adjCl - (-aRich / tau)); d > maxAdjErr {
			maxAdjErr = d
		}
	}
	t.Logf("Riccati vs RK4/Richardson maxBErr=%.3e maxAErr=%.3e maxYieldAdjErr=%.3e", maxBErr, maxAErr, maxAdjErr)
	if maxBErr >= 1e-10 || maxAErr >= 1e-10 || maxAdjErr >= 1e-10 {
		t.Fatalf("Riccati closed-form vs RK4 error >= 1e-10: B=%g A=%g adj=%g", maxBErr, maxAErr, maxAdjErr)
	}
	_, _, _, adj30 := p.YieldAdjustmentComponents(30.0)
	t.Logf("30y yield adjustment = %.2f bps", adj30*10_000.0)
}

func TestTransitionPMatrixExpVanLoanGLAndMC(t *testing.T) {
	// 1. Uji segitiga bawah penuh KP dengan elemen luar-diagonal tak-nol terhadap MatExpSeries, GL, dan Van Loan
	fullP := DefaultParams()
	fullP.KP = [3][3]float64{
		{0.05, 0, 0},
		{0.12, 0.5, 0},
		{-0.08, 0.25, 0.4},
	}
	dt := 0.75
	phiCl := fullP.TransitionMatrixClosed(dt)
	negKdt := lag.Zeros(3, 3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			negKdt[i][j] = -fullP.KP[i][j] * dt
		}
	}
	phiSer := lag.MatExpSeries(negKdt, 35)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if d := math.Abs(phiCl[i][j] - phiSer[i][j]); d > 1e-14 {
				t.Fatalf("Lower-tri Phi(%g) mismatch at (%d,%d): %v vs %v", dt, i, j, phiCl[i][j], phiSer[i][j])
			}
		}
	}
	qGL := fullP.TransitionCovGL(dt, 32)
	qVL := fullP.TransitionCovVanLoan(dt)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			if d := math.Abs(qGL[i][j] - qVL[i][j]); d > 1e-14 {
				t.Fatalf("Q_dt GL vs VanLoan mismatch at (%d,%d): %v vs %v", i, j, qGL[i][j], qVL[i][j])
			}
		}
	}

	// 2. Uji diagonal default + MC transisi + MC premi jangka waktu + determinisme pekerja
	p := DefaultParams()
	x0 := [3]float64{0.052, -0.015, 0.008}
	mc1 := p.VerifyTransitionMC(x0, 1.0/52.0, 40_000, 8, 1, rng.MainSeed)
	mc4 := p.VerifyTransitionMC(x0, 1.0/52.0, 40_000, 8, 4, rng.MainSeed)
	if math.Float64bits(mc1.MaxZScore) != math.Float64bits(mc4.MaxZScore) {
		t.Fatalf("TransitionMC worker invariance violated")
	}
	if mc1.MaxZScore > 3.5 {
		t.Fatalf("TransitionMC z-score too large: %v", mc1.MaxZScore)
	}

	tp := p.VerifyTermPremiumMC(10.0, 40, 20_000, 8, 4, rng.MainSeed)
	t.Logf("10y Term Premium closed=%.2f bps vs MC=%.2f+-%.2f bps (z=%.2f)",
		tp.ClosedFormTPBps, tp.MCTermPremium*10_000.0, tp.MCTermPremiumSE*10_000.0, tp.ZScore)
	if tp.ZScore > 3.0 {
		t.Fatalf("TermPremium MC z-score too large: %v", tp.ZScore)
	}
}
