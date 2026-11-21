package panjer

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/stats"
)

func TestPanjerBracketAndOrder(t *testing.T) {
	hList := []float64{0.1, 0.05, 0.02, 0.01}
	uList := []float64{0, 1, 5, 10, 20, 50}

	for _, d := range claims.LightTailedDists() {
		rows := EvaluateGrid(d, hList, uList)
		widthsU10 := make([]float64, 0, len(hList))
		logH := make([]float64, 0, len(hList))
		for _, r := range rows {
			if !r.Bracketed {
				t.Fatalf("%s not bracketed at h=%g u=%g: [%v, %v] vs exact %v", d.Name, r.H, r.U, r.PsiLow, r.PsiHigh, r.PsiExact)
			}
			if r.U == 10.0 {
				widthsU10 = append(widthsU10, math.Log(r.Width))
				logH = append(logH, math.Log(r.H))
			}
		}
		slope := stats.LinearSlope(logH, widthsU10)
		if math.Abs(slope-1.0) > 0.08 {
			t.Fatalf("%s discretization order slope=%v not near 1.0", d.Name, slope)
		}
		t.Logf("%s Panjer width order at u=10: slope=%.4f", d.Name, slope)
	}

	// Uji ekor berat stabil hingga u = 200 dengan h = 0.05 (cepat di unit test).
	for _, d := range claims.HeavyTailedDists() {
		ser := ComputeBounds(claims.Rho, d.LadderSurvival, 200.0, 0.05)
		lo, _, hi, _, loCdf, _ := ser.At(200.0)
		if lo <= 0 || lo > hi {
			t.Fatalf("%s positive-tail Panjer failed at u=200: lo=%v hi=%v", d.Name, lo, hi)
		}
		t.Logf("%s at u=200 (h=0.05): tail=[%.6g, %.6g], cdfSumLow=%.6g", d.Name, lo, hi, loCdf)
	}
}
