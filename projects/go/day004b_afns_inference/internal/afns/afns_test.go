package afns

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004b_afns_inference/internal/rng"
)

func TestG1G2AndYieldAdjustmentSmallU(t *testing.T) {
	p := DefaultParams()
	for _, tau := range []float64{1e-6, 5e-5, 0.05, 0.25, 1.0, 10.0, 30.0} {
		if g1 := G1(tau, p.Lambda); g1 <= 0 || g1 > 1.0000001 {
			t.Fatalf("invalid G1(%g)=%g", tau, g1)
		}
		if g2 := G2(tau, p.Lambda); g2 < 0 || g2 > 0.35 {
			t.Fatalf("invalid G2(%g)=%g", tau, g2)
		}
	}
	// Kontinuitas batas Gauss-Legendre u < 0.05 vs rumus tertutup di u = 0.05
	pLeft := p
	pRight := p
	tauLo := (0.05 - 1e-8) / p.Lambda
	tauHi := (0.05 + 1e-8) / p.Lambda
	_, _, _, adjLo := pLeft.YieldAdjustmentComponents(tauLo)
	_, _, _, adjHi := pRight.YieldAdjustmentComponents(tauHi)
	if math.Abs(adjLo-adjHi) > 1e-11 {
		t.Fatalf("YieldAdjustment discontinuity across u=0.05: %v vs %v", adjLo, adjHi)
	}
}

func TestDiscreteLyapunovStationarity(t *testing.T) {
	p := DefaultParams()
	qInf := p.StationaryCov()
	for _, dt := range []float64{1.0 / 12.0, 1.0 / 52.0, 1.0 / 252.0, 1.0} {
		phi := p.TransitionMatrixClosed(dt)
		qDt := p.TransitionCovClosed(dt)
		for i := 0; i < 3; i++ {
			for j := 0; j < 3; j++ {
				lhs := qDt[i][j]
				for a := 0; a < 3; a++ {
					for b := 0; b < 3; b++ {
						lhs += phi[i][a] * qInf[a][b] * phi[j][b]
					}
				}
				if diff := math.Abs(lhs - qInf[i][j]); diff > 1e-14 {
					t.Fatalf("discrete Lyapunov mismatch for dt=%g at (%d,%d): diff=%g", dt, i, j, diff)
				}
			}
		}
	}
}

func TestSimulateStationaryDeterminism(t *testing.T) {
	p := DefaultParams()
	taus := []float64{0.25, 2.0, 10.0}
	r1 := rng.New(rng.DeriveSeed(rng.MainSeed, 0x11, 0))
	r2 := rng.New(rng.DeriveSeed(rng.MainSeed, 0x11, 0))
	_, y1 := SimulateStationaryPath(p, taus, 25, 1.0/52.0, 0.0005, r1)
	_, y2 := SimulateStationaryPath(p, taus, 25, 1.0/52.0, 0.0005, r2)
	for i := range y1 {
		for j := range y1[i] {
			if math.Float64bits(y1[i][j]) != math.Float64bits(y2[i][j]) {
				t.Fatalf("path mismatch at (%d,%d)", i, j)
			}
		}
	}
}
