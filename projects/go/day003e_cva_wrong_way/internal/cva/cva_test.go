package cva

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/cir"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/hw"
	"github.com/Orls967/100QB/go/day003e_cva_wrong_way/internal/rng"
)

func TestAnalyticCVASensitivityAndMC(t *testing.T) {
	hwP := hw.DefaultParams()
	cirP := cir.DefaultParams()
	swap := hwP.NewAnnualSwap(10, 0, true)
	tGrid := QuarterlyGrid(10)

	anal := ComputeAnalyticCVAAndSensitivities(hwP, cirP, swap, tGrid, 32)
	if anal.DLambda0AbsErr >= 1e-6 {
		t.Fatalf("dCVA/dLambda0 analytic vs central FD err=%v >= 1e-6", anal.DLambda0AbsErr)
	}
	if anal.CVAAnalytic <= 0 || anal.CS01Lambda0_1bp <= 0 || anal.CS01Theta_1bp <= 0 {
		t.Fatalf("unexpected non-positive CVA or CS01: %+v", anal)
	}

	rows := VerifyIndependentCVAMC(hwP, cirP, swap, tGrid, anal.CVAAnalytic, []int{20, 40, 80}, 40000, 8, 4, rng.MainSeed)
	for i, r := range rows {
		t.Logf("m=%d anal=%.8f tau=%.8f(SE=%.3g,z=%.2f) rb=%.8f(SE=%.3g,z=%.2f) pairedBias320=%.3g intRMSE=%.3g",
			r.SubStepsPerYear, anal.CVAAnalytic, r.CVATauBucket, r.SETauBucket, r.ZScoreTauBucket,
			r.CVARaoBlackwell, r.SERaoBlackwell, r.ZScoreRB, r.PairedBiasVs320, r.IntRMSEVs320)
		if r.ZScoreTauBucket > 3.0 {
			t.Fatalf("m=%d tau bucket CVA z-score=%v > 3.0", r.SubStepsPerYear, r.ZScoreTauBucket)
		}
		if i > 0 && r.IntRMSEVs320 >= rows[i-1].IntRMSEVs320 {
			t.Fatalf("intensity integral RMSE vs m=320 did not decrease: m=%d (%v) >= m=%d (%v)",
				r.SubStepsPerYear, r.IntRMSEVs320, rows[i-1].SubStepsPerYear, rows[i-1].IntRMSEVs320)
		}
	}
	if rows[len(rows)-1].ZScoreRB > 3.5 {
		t.Fatalf("m=%d RB CVA z-score=%v > 3.5", rows[len(rows)-1].SubStepsPerYear, rows[len(rows)-1].ZScoreRB)
	}
}

func TestWWRNettingAndCollateralMonotonicity(t *testing.T) {
	hwP := hw.DefaultParams()
	cirP := cir.DefaultParams()
	swap10Pay := hwP.NewAnnualSwap(10, 0, true)
	swap10Rec := hwP.NewAnnualSwap(10, 0, false)
	swap5Pay := hwP.NewAnnualSwap(5, 0, true)
	tGrid := QuarterlyGrid(10)

	rhos := []float64{-0.8, -0.4, 0, 0.4, 0.8}
	wwr1 := ComputeWWRCurve(hwP, cirP, swap10Pay, tGrid, rhos, 40, 20000, 8, 1, rng.MainSeed)
	wwr4 := ComputeWWRCurve(hwP, cirP, swap10Pay, tGrid, rhos, 40, 20000, 8, 4, rng.MainSeed)
	for i := range wwr1 {
		if math.Float64bits(wwr1[i].CVARaoBlackwell) != math.Float64bits(wwr4[i].CVARaoBlackwell) {
			t.Fatalf("WWR worker determinism failed at rho=%.1f", wwr1[i].Rho)
		}
		if i > 0 && wwr1[i].CVARaoBlackwell <= wwr1[i-1].CVARaoBlackwell {
			t.Fatalf("WWR not strictly increasing for payer swap: rho=%.1f (%v) vs rho=%.1f (%v)",
				wwr1[i-1].Rho, wwr1[i-1].CVARaoBlackwell, wwr1[i].Rho, wwr1[i].CVARaoBlackwell)
		}
	}

	net := ComputeNetting(hwP, cirP, swap5Pay, swap10Rec, swap10Rec, tGrid, []float64{-0.8, 0, 0.8}, 32, 40, 20000, 8, 4, rng.MainSeed)
	if net.MaxPointwiseViol > 1e-15 {
		t.Fatalf("subadditivity violated pointwise: %v", net.MaxPointwiseViol)
	}
	if net.OppositeNetMaxEE > 1e-15 || math.Abs(net.OppositeNetCVA) > 1e-15 {
		t.Fatalf("opposite netting identity failed: ee=%v cva=%v", net.OppositeNetMaxEE, net.OppositeNetCVA)
	}
	for _, r := range net.ByRho {
		if r.NettedCVA > r.StandaloneSum {
			t.Fatalf("netted CVA exceeded standalone sum at rho=%.1f: %v > %v", r.Rho, r.NettedCVA, r.StandaloneSum)
		}
	}

	// Uji kolateral: delta = 0.04 across H = [0, 0.01, 0.05, Inf] dan MPOR -> 0 pada H = 0.
	colH := ComputeCollateralCVA(hwP, cirP, swap10Pay, tGrid, 0.04, []float64{0, 0.01, 0.05, math.Inf(1)}, []float64{0, 0.8}, 40, 20000, 8, 4, rng.MainSeed)
	for _, rho := range []float64{0, 0.8} {
		var prev float64
		first := true
		for _, r := range colH {
			if r.Rho != rho {
				continue
			}
			if !first && r.CVA <= prev {
				t.Fatalf("collateral CVA not increasing with threshold at rho=%.1f: prev=%v curr=%v", rho, prev, r.CVA)
			}
			prev = r.CVA
			first = false
		}
	}

	var prevMPOR float64
	for idx, d := range []float64{0.04, 0.01, 0.0025} {
		rows := ComputeCollateralCVA(hwP, cirP, swap10Pay, tGrid, d, []float64{0}, []float64{0}, 40, 20000, 8, 4, rng.MainSeed)
		if idx > 0 && rows[0].CVA >= prevMPOR {
			t.Fatalf("MPOR convergence failed: delta=%v cva=%v >= prev=%v", d, rows[0].CVA, prevMPOR)
		}
		prevMPOR = rows[0].CVA
	}
}
