package closedform

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/lundberg"
)

func TestLundbergAndClosedForm(t *testing.T) {
	for _, d := range claims.LightTailedDists() {
		mean, _, cv := d.Moments()
		if math.Abs(mean-1.0) > 1e-14 {
			t.Fatalf("%s mean=%v != 1", d.Name, mean)
		}
		if d.Kind == claims.KindHyperexp2 && cv <= 1.0 {
			t.Fatalf("Hyperexp2 CV=%v <= 1", cv)
		}

		resL, err := lundberg.Solve(d)
		if err != nil {
			t.Fatalf("lundberg.Solve(%s): %v", d.Name, err)
		}
		cf, err := Derive(d)
		if err != nil {
			t.Fatalf("closedform.Derive(%s): %v", d.Name, err)
		}

		if math.Abs(cf.Eval(0)-claims.Rho) > 1e-14 {
			t.Fatalf("%s psi(0)=%v != rho=%v", d.Name, cf.Eval(0), claims.Rho)
		}
		if math.Abs(resL.R-cf.R1) >= 1e-12 {
			t.Fatalf("%s Brent R=%v vs Laplace R1=%v", d.Name, resL.R, cf.R1)
		}
		if math.Abs(resL.C-cf.C1) >= 1e-12 {
			t.Fatalf("%s Cramér-Lundberg C=%v vs residue C1=%v", d.Name, resL.C, cf.C1)
		}
		if d.Kind == claims.KindExponential {
			rExact := claims.Theta / ((1.0 + claims.Theta) * claims.Mu)
			if math.Abs(resL.R-rExact) >= 1e-12 {
				t.Fatalf("Exp Brent R=%v != exact %v", resL.R, rExact)
			}
		}

		// Monoton menurun & batas Lundberg psi(u) <= e^{-R u}
		prev := 1.0
		for _, u := range []float64{0, 1, 5, 10, 20, 50, 100} {
			v := cf.Eval(u)
			if v >= prev || v > lundberg.Bound(resL.R, u)+1e-15 {
				t.Fatalf("%s failed monotonicity/Lundberg bound at u=%v: psi=%v bound=%v", d.Name, u, v, lundberg.Bound(resL.R, u))
			}
			prev = v
		}
		ratio100 := cf.Eval(100.0) / (resL.C * math.Exp(-resL.R*100.0))
		if math.Abs(ratio100-1.0) > 1e-8 {
			t.Fatalf("%s Cramér-Lundberg asymptotic ratio(100)=%v != 1", d.Name, ratio100)
		}
	}

	for _, d := range claims.HeavyTailedDists() {
		if _, err := lundberg.Solve(d); err == nil {
			t.Fatalf("expected error from lundberg.Solve for heavy tail %s", d.Name)
		}
	}
}
