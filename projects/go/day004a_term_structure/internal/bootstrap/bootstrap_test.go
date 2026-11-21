package bootstrap

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
)

func TestBootstrapExactRecovery(t *testing.T) {
	mats := SemiAnnualGrid(30.0)
	for _, tc := range []struct {
		name  string
		curve ns.Params
	}{
		{"NS", ns.DefaultNS()},
		{"NSS", ns.DefaultNSS()},
	} {
		par, trueDF, trueY := ParRatesFromCurve(tc.curve, mats, 0.5)
		recDF, recY, err := BootstrapFromParRates(mats, par, 0.5)
		if err != nil {
			t.Fatalf("%s bootstrap failed: %v", tc.name, err)
		}
		maxDFErr := 0.0
		maxYErr := 0.0
		for i := range mats {
			if e := math.Abs(recDF[i] - trueDF[i]); e > maxDFErr {
				maxDFErr = e
			}
			if e := math.Abs(recY[i] - trueY[i]); e > maxYErr {
				maxYErr = e
			}
		}
		t.Logf("%s bootstrap maxDFErr=%.3e maxYErr=%.3e", tc.name, maxDFErr, maxYErr)
		if maxDFErr >= 1e-12 || maxYErr >= 1e-12 {
			t.Fatalf("%s bootstrap recovery exceeded 1e-12: dfErr=%g yErr=%g", tc.name, maxDFErr, maxYErr)
		}
	}
}
