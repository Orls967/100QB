package ns

import (
	"math"
	"testing"
)

func TestYieldForwardIdentityLimitsAndSmallTauStability(t *testing.T) {
	curves := []struct {
		name string
		p    Params
	}{
		{"NS", DefaultNS()},
		{"NSS", DefaultNSS()},
	}
	taus := []float64{0.1, 1.0, 5.0, 30.0}

	for _, c := range curves {
		for _, tau := range taus {
			y := c.p.Yield(tau)
			avgF := c.p.ForwardIntegralAvg(tau, 32)
			err := math.Abs(y - avgF)
			if err >= 1e-12 {
				t.Fatalf("%s tau=%g: |y - (1/tau)int f| = %g >= 1e-12", c.name, tau, err)
			}
		}

		// Limit tau -> 0: y -> L + S
		yTiny := c.p.Yield(1e-12)
		if math.Abs(yTiny-(c.p.L+c.p.S)) > 1e-11 {
			t.Fatalf("%s limit tau->0 mismatch: got %v, want %v", c.name, yTiny, c.p.L+c.p.S)
		}
		// Limit tau -> inf: y -> L
		yHuge := c.p.Yield(1e9)
		if math.Abs(yHuge-c.p.L) > 1e-8 {
			t.Fatalf("%s limit tau->inf mismatch: got %v, want %v", c.name, yHuge, c.p.L)
		}
	}

	// Demonstrasi kehilangan digit pada evaluasi naif vs stabilitas deret Taylor
	p := DefaultNS()
	rows := p.CompareSmallTauStability([]float64{1e-4, 1e-8, 1e-12, 1e-14, 1e-16})
	for _, r := range rows {
		t.Logf("tau=%.0e | g1RelErr=%.3e | g2RelErr=%.3e | yAbsErr=%.3e", r.Tau, r.G1RelErr, r.G2RelErr, r.YieldAbsErr)
	}
	// Pada tau=1e-12, versi naif kehilangan ~4 digit pada g1 dan >100% pada g2, sedangkan versi stabil mulus ke L+S
	if rows[2].G1RelErr < 1e-6 || rows[2].G2RelErr < 0.1 {
		t.Fatalf("expected naive evaluation to lose digits at tau=1e-12, got g1RelErr=%g g2RelErr=%g", rows[2].G1RelErr, rows[2].G2RelErr)
	}
	if math.Abs(rows[4].YieldStable-(p.L+p.S)) > 1e-15 {
		t.Fatalf("stable yield at tau=1e-16 failed: got %v, want %v", rows[4].YieldStable, p.L+p.S)
	}
}
