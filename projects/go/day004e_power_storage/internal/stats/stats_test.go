package stats

import (
	"math"
	"testing"
)

func TestQuantileVaRESAndPairedDiff(t *testing.T) {
	xs := make([]float64, 100)
	for i := 0; i < 100; i++ {
		xs[i] = float64(i + 1) // 1..100
	}
	if med := Quantile(xs, 0.5); math.Abs(med-50.5) > 1e-12 {
		t.Fatalf("median mismatch: %v", med)
	}
	// VaR 5% pada 1..100: pos = 0.05*99 = 4.95 -> 0.05*5 + 0.95*6 = 5.95.
	// ES 5%: k = ceil(0.05*100) = 5 -> mean(1,2,3,4,5) = 3.0.
	v, es := VaRAndES(xs, 0.05)
	if math.Abs(v-5.95) > 1e-12 || math.Abs(es-3.0) > 1e-12 {
		t.Fatalf("VaR/ES mismatch: var=%v es=%v", v, es)
	}

	a := []float64{10, 12, 14, 16}
	b := []float64{8, 10, 11, 13}
	// diffs = [2, 2, 3, 3] -> mean=2.5, sd=sqrt(1/3), se=sqrt(1/12)
	pd := PairedDiffCI(a, b)
	expSD := math.Sqrt(1.0 / 3.0)
	expSE := expSD / 2.0
	if pd.N != 4 || math.Abs(pd.MeanDiff-2.5) > 1e-12 || math.Abs(pd.SDDiff-expSD) > 1e-12 || math.Abs(pd.MCSE-expSE) > 1e-12 {
		t.Fatalf("PairedDiffCI mismatch: %+v", pd)
	}
	if pd.PositiveFrac != 1.0 {
		t.Fatalf("expected PositiveFrac=1.0, got %v", pd.PositiveFrac)
	}
}
