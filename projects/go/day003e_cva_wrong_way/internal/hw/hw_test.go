package hw

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
)

func TestHullWhiteMartingaleAndJamshidian(t *testing.T) {
	p := DefaultParams()
	for _, tc := range [][2]float64{{2.0, 10.0}, {5.0, 10.0}, {8.0, 10.0}} {
		lh1, lh2, p0T, err1, err2 := p.MartingaleCheckGH(tc[0], tc[1], 32)
		if err1 > 1e-14 || err2 > 1e-14 {
			t.Fatalf("martingale GH failed at (t=%.1f, T=%.1f): lh1=%v lh2=%v p0T=%v err1=%v err2=%v",
				tc[0], tc[1], lh1, lh2, p0T, err1, err2)
		}
		mcMean, mcSE, _, z := p.MartingaleCheckMC(tc[0], tc[1], 100000, 16, 4, rng.MainSeed)
		if math.Abs(mcMean-p0T) >= 5e-4 || z > 3.5 {
			t.Fatalf("martingale MC failed at (t=%.1f, T=%.1f): mc=%v se=%v p0T=%v z=%v",
				tc[0], tc[1], mcMean, mcSE, p0T, z)
		}
	}

	payDates := []float64{3, 4, 5, 6, 7, 8, 9, 10}
	kFwd := p.ForwardSwapRate(2.0, payDates)
	strikes := []float64{kFwd - 0.01, kFwd, kFwd + 0.01}
	labels := []string{"ITM", "ATM", "OTM"}
	res := p.VerifySwaption(2.0, payDates, strikes, labels, 32, 100000, 16, 4, rng.MainSeed)
	for _, r := range res {
		if math.Abs(r.CheckSumCiXi-1.0) > 1e-12 {
			t.Fatalf("Jamshidian sum c_i X_i = %v != 1", r.CheckSumCiXi)
		}
		if r.AbsErrJamVsGH > 1e-8 {
			t.Fatalf("Swaption %s Jamshidian=%v GH2D=%v diff=%v > 1e-8", r.Label, r.Jamshidian, r.QuadGH2D, r.AbsErrJamVsGH)
		}
		if r.ZScoreMC > 3.0 {
			t.Fatalf("Swaption %s MC z-score=%.3f > 3.0", r.Label, r.ZScoreMC)
		}
	}
}

func TestExposureGHvsMCAndWorkerDeterminism(t *testing.T) {
	p := DefaultParams()
	swap := p.NewAnnualSwap(10, 0, true)
	tGrid := make([]float64, 41)
	for j := range tGrid {
		tGrid[j] = float64(j) * 0.25
	}
	pts1, epeGH1, epeMC1 := p.ComputeSwapExposureProfile(swap, tGrid, 32, 100000, 8, 1, rng.MainSeed)
	pts4, epeGH4, epeMC4 := p.ComputeSwapExposureProfile(swap, tGrid, 32, 100000, 8, 4, rng.MainSeed)

	if math.Float64bits(epeGH1) != math.Float64bits(epeGH4) || math.Float64bits(epeMC1) != math.Float64bits(epeMC4) {
		t.Fatalf("worker count changed bits: EPE (%v vs %v)", epeMC1, epeMC4)
	}
	maxZ := 0.0
	for j := range pts1 {
		if math.Float64bits(pts1[j].EEMC) != math.Float64bits(pts4[j].EEMC) ||
			math.Float64bits(pts1[j].PFE95MC) != math.Float64bits(pts4[j].PFE95MC) {
			t.Fatalf("worker determinism failed at t=%.2f", pts1[j].T)
		}
		if math.Abs(pts1[j].EEDiscGH-pts1[j].EEDiscGH2D) > 1e-13 {
			t.Fatalf("1D vs 2D quad mismatch at t=%.2f: %v vs %v", pts1[j].T, pts1[j].EEDiscGH, pts1[j].EEDiscGH2D)
		}
		if pts1[j].EESEMC > 0 {
			z := math.Abs(pts1[j].EEMC-pts1[j].EEGH) / pts1[j].EESEMC
			if z > maxZ {
				maxZ = z
			}
			if z > 4.0 {
				t.Fatalf("EE MC vs GH mismatch at t=%.2f: z=%v (MC=%v Anal=%v SE=%v)", pts1[j].T, z, pts1[j].EEMC, pts1[j].EEGH, pts1[j].EESEMC)
			}
		}
	}
	t.Logf("Max z-score across 41 quarterly dates = %.3f", maxZ)
}
