package pca

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/kalman"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
)

func TestPCAEigenvaluesCosineSignPatternsAndWorkers(t *testing.T) {
	spec := kalman.DefaultSpec()
	res1, err := RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 2000, 20, 4, 1, rng.MainSeed)
	if err != nil {
		t.Fatalf("RunPCAStudy failed: %v", err)
	}
	res4, err := RunPCAStudy(spec, kalman.StandardTenors, kalman.WeeklyDt, 2000, 20, 4, 4, rng.MainSeed)
	if err != nil {
		t.Fatalf("RunPCAStudy 4 workers failed: %v", err)
	}
	for k := range res1.EmpMeanEigenvalues {
		if math.Float64bits(res1.EmpMeanEigenvalues[k]) != math.Float64bits(res4.EmpMeanEigenvalues[k]) {
			t.Fatalf("PCA worker invariance violated at k=%d", k)
		}
	}

	if res1.TheoSignChanges != [3]int{0, 1, 2} {
		t.Fatalf("unexpected theoretical sign pattern: %v (want [0 1 2])", res1.TheoSignChanges)
	}
	if res1.EmpSignChanges != [3]int{0, 1, 2} {
		t.Fatalf("unexpected empirical sign pattern: %v (want [0 1 2])", res1.EmpSignChanges)
	}
	for k := 0; k < 3; k++ {
		if res1.EmpMeanTop3Cos[k] <= 0.99 {
			t.Fatalf("factor %d cosine similarity %.6f <= 0.99", k+1, res1.EmpMeanTop3Cos[k])
		}
		if res1.ZScoreSampleMean[k] > 3.0 {
			t.Fatalf("factor %d eigenvalue mean z-score %.2f > 3.0", k+1, res1.ZScoreSampleMean[k])
		}
	}
	if res1.NoiseBulkZScore > 3.0 {
		t.Fatalf("noise bulk mean z-score %.2f > 3.0 (emp=%g vs theo=%g)", res1.NoiseBulkZScore, res1.NoiseBulkEmpMean, res1.NoiseBulkTheoMean)
	}
	acc := res1.EvaluateAcceptance()
	if !acc.AggregatePassed || acc.Status != "tercapai" {
		t.Fatalf("expected aggregate PCA acceptance to pass, got %+v", acc)
	}
	t.Logf("PCA top3 var=%.2f%% | signChanges=%v | cos=[%.5f, %.5f, %.5f] | bulkZ=%.2f | diagAll10=%v | agg=%v (%s)",
		res1.TheoTop3CumVar*100.0, res1.TheoSignChanges,
		res1.EmpMeanTop3Cos[0], res1.EmpMeanTop3Cos[1], res1.EmpMeanTop3Cos[2],
		res1.NoiseBulkZScore, acc.DiagnosticAll10SortedWithin3SE, acc.AggregatePassed, acc.Status)
}

func TestPCAAcceptanceDiagnosticSortedNoiseSplitting(t *testing.T) {
	// Simulasikan pembelahan statistik urutan PC4-PC10 (z > 3.0) ketika sinyal top-3, noise bulk,
	// kosinus > 0.99, dan pola tanda [0, 1, 2] seluruhnya memenuhi syarat.
	base := PCAReplicationSummary{
		ZScoreSampleMean: []float64{1.32, 0.72, 0.64, 38.90, 27.15, 16.85, 1.09, 16.79, 33.60, 45.82},
		NoiseBulkZScore:  0.28,
		EmpMeanTop3Cos:   [3]float64{0.99995, 0.99923, 0.99748},
		TheoSignChanges:  [3]int{0, 1, 2},
		EmpSignChanges:   [3]int{0, 1, 2},
	}
	acc := base.EvaluateAcceptance()
	if acc.DiagnosticAll10SortedWithin3SE {
		t.Fatalf("expected DiagnosticAll10SortedWithin3SE=false when sorted noise eigenvalues split")
	}
	if !acc.AggregatePassed || acc.Status != "tercapai" {
		t.Fatalf("expected AggregatePassed=true and Status=tercapai despite sorted noise split, got %+v", acc)
	}

	// Pastikan kriteria agregat tetap gagal apabila rata-rata noise bulk melampaui 3 SE.
	badBulk := base
	badBulk.NoiseBulkZScore = 3.5
	accBad := badBulk.EvaluateAcceptance()
	if accBad.AggregatePassed || accBad.Status != "sebagian" {
		t.Fatalf("expected AggregatePassed=false and Status=sebagian when bulkZ > 3, got %+v", accBad)
	}
}
