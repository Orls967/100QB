package interval

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/stats"
)

func TestIntervalMethodsExactSyntheticAndInvalidHandling(t *testing.T) {
	// Uji data sintetis B=100 pada parameter kp11 (index 1, wajib positif, trueVal=0.05)
	// psiHat = 0.25, seHat = 0.10
	// bootEsts = 0.25 + 0.20 + 0.004*i (i=0..99) -> min=0.45, max=0.846
	psiHat := 0.25
	seHat := 0.10
	trueVal := 0.05
	bootEsts := make([]float64, 100)
	bootSEs := make([]float64, 100)
	for i := 0; i < 100; i++ {
		bootEsts[i] = 0.45 + 0.004*float64(i)
		bootSEs[i] = 0.12
	}
	qLo := stats.Quantile(bootEsts, 0.025)
	qHi := stats.Quantile(bootEsts, 0.975)

	w := EvalWald(1, psiHat, seHat, trueVal, true)
	if !w.Valid || math.Abs(w.Lower-(psiHat-stats.Z95*seHat)) > 1e-12 || w.Covered {
		t.Fatalf("unexpected Wald result: %+v", w)
	}

	p := EvalPercentile(1, trueVal, true, bootEsts, 100)
	if !p.Valid || math.Abs(p.Lower-qLo) > 1e-12 || math.Abs(p.Upper-qHi) > 1e-12 || p.Covered {
		t.Fatalf("unexpected Percentile result: %+v", p)
	}

	b := EvalBasic(1, psiHat, trueVal, true, bootEsts, 100)
	expBasicLo := 2.0*psiHat - qHi // ~ 0.50 - 0.8361 = -0.3361 (negatif, tidak di-clamp!)
	expBasicHi := 2.0*psiHat - qLo // ~ 0.50 - 0.4599 = 0.0401
	if !b.Valid || math.Abs(b.Lower-expBasicLo) > 1e-12 || math.Abs(b.Upper-expBasicHi) > 1e-12 {
		t.Fatalf("unexpected Basic result: %+v", b)
	}
	if !b.NonPhysical || !b.OutOfBounds || b.Lower >= 0 {
		t.Fatalf("expected Basic negative lower bound preserved & flagged NonPhysical/OutOfBounds: %+v", b)
	}

	// Studentized: t* = (bootEst - 0.25) / 0.12 -> q_0.025(t*) = (qLo - 0.25)/0.12
	st := EvalStudentized(1, psiHat, seHat, trueVal, true, bootEsts, bootSEs, 100)
	expStLo := psiHat - ((qHi-psiHat)/0.12)*seHat
	expStHi := psiHat - ((qLo-psiHat)/0.12)*seHat
	if !st.Valid || math.Abs(st.Lower-expStLo) > 1e-12 || math.Abs(st.Upper-expStHi) > 1e-12 {
		t.Fatalf("unexpected Studentized result: %+v", st)
	}

	// LogBasic: selalu positif (interpolasi kuantil di skala log)
	lb := EvalLogBasic(1, psiHat, trueVal, true, bootEsts, 100)
	logDraws := make([]float64, 100)
	for i := 0; i < 100; i++ {
		logDraws[i] = math.Log(bootEsts[i])
	}
	expLbLo := math.Exp(2.0*math.Log(psiHat) - stats.Quantile(logDraws, 0.975))
	expLbHi := math.Exp(2.0*math.Log(psiHat) - stats.Quantile(logDraws, 0.025))
	if !lb.Valid || lb.Lower <= 0 || lb.NonPhysical || math.Abs(lb.Lower-expLbLo) > 1e-12 || math.Abs(lb.Upper-expLbHi) > 1e-12 {
		t.Fatalf("unexpected LogBasic result: %+v", lb)
	}

	// Kasus SE* nonpositif (>5 draw) -> Studentized invalid sementara Percentile masih valid
	badSEs := make([]float64, 100)
	copy(badSEs, bootSEs)
	for i := 0; i < 6; i++ {
		badSEs[i] = 0.0
	}
	stBad := EvalStudentized(1, psiHat, seHat, trueVal, true, bootEsts, badSEs, 100)
	if stBad.Valid || stBad.Covered || stBad.InvalidReason != ReasonInsufficientValidSE {
		t.Fatalf("expected Studentized invalid due to SE*<=0, got %+v", stBad)
	}
}
