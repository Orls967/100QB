package risk

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
)

func TestBondRiskKRDPartitionAndTwoBondHedge(t *testing.T) {
	curve := ns.DefaultNS()
	const bumpDur = 1e-6
	const bumpConv = 1e-4
	m := AnalyzeBondRisk(curve, 10.0, 0.04, 2.0, bumpDur, bumpConv)
	if m.EffDurationZeroRelErr >= 1e-6 || m.EffConvexityZeroRelErr >= 1e-6 ||
		m.ModDurationYTMRelErr >= 1e-6 || m.ConvexityYTMRelErr >= 1e-6 {
		t.Fatalf("Duration/convexity closed vs FD error >= 1e-6: %+v", m)
	}
	t.Logf("10y 4%% Bond: P=%.6f YTM=%.6f ModDur=%.6f (relErr=%.2e) Conv=%.6f (relErr=%.2e)",
		m.Price, m.YTM, m.ModDurationYTMClosed, m.ModDurationYTMRelErr, m.ConvexityYTMClosed, m.ConvexityYTMRelErr)

	times, cfs := BondCashFlows(10.0, 0.04, 2.0)
	krd := ComputeKeyRateDurations(curve, times, cfs, StandardKeyRateNodes, bumpDur)
	if krd.SumVsParFDErr >= 1e-8 || krd.SumVsParClErr >= 1e-12 {
		t.Fatalf("KRD sum vs parallel duration error too high: fd=%g cl=%g", krd.SumVsParFDErr, krd.SumVsParClErr)
	}
	t.Logf("KRD sum=%.10f vs effPar=%.10f (fdErr=%.2e clErr=%.2e)",
		krd.SumKRDFD, krd.EffParallelFD, krd.SumVsParFDErr, krd.SumVsParClErr)

	sw := AnalyzeSwapPV01(curve, 10.0, 2.0)
	if math.Abs(sw.PayerCurveDV01Ex-sw.PayerCurveDV01FD) > 1e-8 {
		t.Fatalf("Swap DV01 exact vs FD mismatch: %g vs %g", sw.PayerCurveDV01Ex, sw.PayerCurveDV01FD)
	}

	hedge, err := HedgeTwoKRDs(curve, StandardKeyRateNodes, 3, 5, bumpDur)
	if err != nil {
		t.Fatalf("HedgeTwoKRDs failed: %v", err)
	}
	if math.Abs(hedge.ResidualNode1) > 1e-13 || math.Abs(hedge.ResidualNode2) > 1e-13 {
		t.Fatalf("Hedge residual not zero at neutralized nodes: %g, %g", hedge.ResidualNode1, hedge.ResidualNode2)
	}
	t.Logf("Hedge weights (5y4%%, 10y2%%) = (%.6f, %.6f) | residuals=%.2e, %.2e",
		hedge.WeightHedge1, hedge.WeightHedge2, hedge.ResidualNode1, hedge.ResidualNode2)
}
