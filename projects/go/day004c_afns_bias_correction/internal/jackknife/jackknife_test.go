package jackknife

import (
	"math"
	"testing"
)

func TestJackknifeAlgebraicFirstOrderBiasCancellation(t *testing.T) {
	// Uji aljabar murni atas asumsi ekspansi bias c/H + d/H^2:
	// Ini HANYA membuktikan identitas aljabar apabila asumsi c/H berlaku, bukan bahwa
	// MLE AFNS nonlinier memenuhi ekspansi tersebut secara eksak.
	trueVal := 0.05
	c := 4.0
	for _, H := range []float64{20.0, 80.0} {
		psiFull := trueVal + c/H
		psiH1 := trueVal + c/(0.5*H)
		psiH2 := trueVal + c/(0.5*H)
		jk := CombineSplit(psiFull, psiH1, psiH2)
		if err := math.Abs(jk - trueVal); err > 1e-14 {
			t.Fatalf("O(1/H) exact cancellation failed for H=%g: jk=%v diff=%g", H, jk, err)
		}

		// Dengan suku orde dua d/H^2:
		// 2*(psi + c/H + d/H^2) - (psi + 2c/H + 4d/H^2) = psi - 2*d/H^2.
		d := 12.0
		psiFull2 := trueVal + c/H + d/(H*H)
		psiH12 := trueVal + c/(0.5*H) + d/((0.5*H)*(0.5*H))
		psiH22 := psiH12
		jk2 := CombineSplit(psiFull2, psiH12, psiH22)
		expected2 := trueVal - 2.0*d/(H*H)
		if err := math.Abs(jk2 - expected2); err > 1e-14 {
			t.Fatalf("second-order O(1/H^2) identity failed for H=%g: got=%v want=%v", H, jk2, expected2)
		}
	}

	// Uji tidak ada clamping: hasil negatif kp11 tetap dipertahankan dan ditandai NonPhysical & OutOfBounds
	full := [11]float64{0.6, 0.02, 0.5, 0.4, 0.045, -0.02, 0.01, 0.006, 0.01, 0.02, 0.0005}
	h1 := [11]float64{0.6, 0.10, 0.5, 0.4, 0.045, -0.02, 0.01, 0.006, 0.01, 0.02, 0.0005}
	h2 := [11]float64{0.6, 0.10, 0.5, 0.4, 0.045, -0.02, 0.01, 0.006, 0.01, 0.02, 0.0005}
	res := CombineVec(full, h1, h2) // kp11 = 2*0.02 - 0.10 = -0.06
	if math.Abs(res.EstParams[1]-(-0.06)) > 1e-12 || !res.OutOfBoundsByParam[1] || !res.NonPhysicalByParam[1] {
		t.Fatalf("expected unclamped negative kp11 flagged out-of-bounds & non-physical, got %+v", res)
	}
}
