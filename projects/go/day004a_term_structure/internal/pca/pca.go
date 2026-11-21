package pca

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/afns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/kalman"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/parallel"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/stats"
)

// TheoreticalCov menghitung matriks kovarians teoretis perubahan yield mingguan:
// Sigma_dy = H * [(Phi - I) Q_inf (Phi - I)^T + Q_dt] * H^T + 2 * sigma_eps^2 * I_M,
// sekaligus versi pendekatan inovasi H * Q_dt * H^T + 2 * sigma_eps^2 * I_M.
func TheoreticalCov(spec kalman.ModelSpec, taus []float64, dt float64) (exactCov, approxCov [][]float64) {
	m := len(taus)
	hMat := lag.Zeros(m, 3)
	for i, tau := range taus {
		row := spec.AFNS.YieldLoadings(tau)
		hMat[i][0], hMat[i][1], hMat[i][2] = row[0], row[1], row[2]
	}
	phi := spec.AFNS.TransitionMatrixClosed(dt)
	qDt := spec.AFNS.TransitionCovClosed(dt)
	qInf := spec.AFNS.StationaryCov()

	var phiMinusI [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			phiMinusI[i][j] = phi[i][j]
			if i == j {
				phiMinusI[i][j] -= 1.0
			}
		}
	}

	var tmp [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			tmp[i][j] = phiMinusI[i][0]*qInf[0][j] + phiMinusI[i][1]*qInf[1][j] + phiMinusI[i][2]*qInf[2][j]
		}
	}
	sigmaDx := lag.Zeros(3, 3)
	sigmaDxApprox := lag.Zeros(3, 3)
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			sigmaDxApprox[i][j] = qDt[i][j]
			sigmaDx[i][j] = qDt[i][j] + tmp[i][0]*phiMinusI[j][0] + tmp[i][1]*phiMinusI[j][1] + tmp[i][2]*phiMinusI[j][2]
		}
	}

	hT := lag.Transpose(hMat)
	exactCov = lag.MatMul(hMat, lag.MatMul(sigmaDx, hT))
	approxCov = lag.MatMul(hMat, lag.MatMul(sigmaDxApprox, hT))
	noiseVar := 2.0 * spec.SigmaEps * spec.SigmaEps
	for i := 0; i < m; i++ {
		exactCov[i][i] += noiseVar
		approxCov[i][i] += noiseVar
	}
	return exactCov, approxCov
}

// CountSignChanges menghitung banyaknya perubahan tanda di sepanjang komponen vektor eigen.
func CountSignChanges(v []float64) int {
	cnt := 0
	prevSign := 0
	for _, x := range v {
		s := 0
		if x > 1e-12 {
			s = 1
		} else if x < -1e-12 {
			s = -1
		}
		if s == 0 {
			continue
		}
		if prevSign != 0 && s != prevSign {
			cnt++
		}
		prevSign = s
	}
	return cnt
}

func orientEigenvectors(vecs [][]float64) {
	if len(vecs) >= 1 {
		// Faktor 1 (Level): orientasi agar rata-rata positif
		if stats.Mean(vecs[0]) < 0 {
			for i := range vecs[0] {
				vecs[0][i] = -vecs[0][i]
			}
		}
	}
	if len(vecs) >= 2 && len(vecs[1]) >= 2 {
		// Faktor 2 (Slope): orientasi agar naik dari tenor pendek ke panjang
		if vecs[1][len(vecs[1])-1] < vecs[1][0] {
			for i := range vecs[1] {
				vecs[1][i] = -vecs[1][i]
			}
		}
	}
	if len(vecs) >= 3 && len(vecs[2]) >= 1 {
		// Faktor 3 (Kurvatur): orientasi agar ujung pendek positif
		if vecs[2][0] < 0 {
			for i := range vecs[2] {
				vecs[2][i] = -vecs[2][i]
			}
		}
	}
}

func alignVectorSign(v, ref []float64) []float64 {
	out := make([]float64, len(v))
	copy(out, v)
	if lag.Dot(out, ref) < 0 {
		for i := range out {
			out[i] = -out[i]
		}
	}
	return out
}

// SampleYieldDiffCov menghitung matriks kovarians sampel dari perubahan yield Delta y_t = y_t - y_{t-1}.
func SampleYieldDiffCov(yields [][]float64) [][]float64 {
	tSteps := len(yields)
	m := len(yields[0])
	nDiff := tSteps - 1
	means := make([]float64, m)
	for t := 1; t < tSteps; t++ {
		for j := 0; j < m; j++ {
			means[j] += yields[t][j] - yields[t-1][j]
		}
	}
	for j := 0; j < m; j++ {
		means[j] /= float64(nDiff)
	}
	cov := lag.Zeros(m, m)
	for t := 1; t < tSteps; t++ {
		for i := 0; i < m; i++ {
			di := (yields[t][i] - yields[t-1][i]) - means[i]
			for j := i; j < m; j++ {
				dj := (yields[t][j] - yields[t-1][j]) - means[j]
				cov[i][j] += di * dj
			}
		}
	}
	denom := float64(nDiff - 1)
	for i := 0; i < m; i++ {
		for j := i; j < m; j++ {
			val := cov[i][j] / denom
			cov[i][j] = val
			cov[j][i] = val
		}
	}
	return cov
}

// PCAReplicationSummary menyimpan perbandingan penuh eigen-dekomposisi teoretis vs empiris (T, R).
type PCAReplicationSummary struct {
	Tenors             []float64
	TSteps             int
	NumReps            int
	TheoEigenvalues    []float64
	ApproxEigenvalues  []float64
	TheoVarExplained   []float64
	TheoTop3CumVar     float64
	TheoEigenvectors   [][]float64
	TheoSignChanges    [3]int
	EmpMeanEigenvalues []float64
	EmpSDEigenvalues   []float64 // SE satu sampel T langkah (SD lintas R replikasi)
	EmpMeanSE          []float64 // SE rata-rata R replikasi = SD / sqrt(R)
	ZScoreSampleMean   []float64 // |EmpMean - Theo| / EmpMeanSE
	SingleRunEigenvals []float64 // replikasi pertama sebagai contoh 1 lintasan T
	SingleRunZScores   []float64 // |SingleRun - Theo| / EmpSDEigenvalues
	NoiseBulkTheoMean  float64
	NoiseBulkEmpMean   float64
	NoiseBulkMeanSE    float64
	NoiseBulkZScore    float64
	SingleRunTop3Cos   [3]float64
	EmpMeanTop3Cos     [3]float64
	EmpMinTop3Cos      [3]float64
	EmpMeanTop3Vecs    [][]float64
	EmpSignChanges     [3]int
	EmpVarExplained    []float64
	EmpTop3CumVar      float64
}

// PCAAcceptance merangkum evaluasi kriteria PCA yang membedakan statistik sinyal/bulk teridentifikasi
// dari diagnostik eigenvalue mode derau terdegenerasi yang sudah diurutkan.
type PCAAcceptance struct {
	Top3EigenvaluesWithin3SE       bool
	NoiseBulkMeanWithin3SE         bool
	DiagnosticAll10SortedWithin3SE bool
	Top3CosineGt0p99               bool
	SignPatterns012                bool
	AggregatePassed                bool
	Status                         string
}

// EvaluateAcceptance menghitung status agregat PCA berdasarkan bagian yang teridentifikasi
// (PC1-PC3 dalam 3 SE, rata-rata noise bulk PC4-PC10 dalam 3 SE, kosinus top-3 > 0.99, pola tanda [0, 1, 2]),
// sambil tetap mencatat DiagnosticAll10SortedWithin3SE sebagai diagnostik pemisahan Wishart/Marchenko-Pastur.
func (s PCAReplicationSummary) EvaluateAcceptance() PCAAcceptance {
	top3OK := len(s.ZScoreSampleMean) >= 3
	for k := 0; k < 3 && k < len(s.ZScoreSampleMean); k++ {
		if s.ZScoreSampleMean[k] > 3.0 {
			top3OK = false
		}
	}
	bulkOK := s.NoiseBulkZScore <= 3.0
	all10OK := len(s.ZScoreSampleMean) > 0
	for _, z := range s.ZScoreSampleMean {
		if z > 3.0 {
			all10OK = false
		}
	}
	cosOK := s.EmpMeanTop3Cos[0] > 0.99 && s.EmpMeanTop3Cos[1] > 0.99 && s.EmpMeanTop3Cos[2] > 0.99
	signOK := s.TheoSignChanges == [3]int{0, 1, 2} && s.EmpSignChanges == [3]int{0, 1, 2}

	agg := top3OK && bulkOK && cosOK && signOK
	status := "belum tercapai"
	if agg {
		status = "tercapai"
	} else if top3OK || bulkOK || cosOK || signOK {
		status = "sebagian"
	}

	return PCAAcceptance{
		Top3EigenvaluesWithin3SE:       top3OK,
		NoiseBulkMeanWithin3SE:         bulkOK,
		DiagnosticAll10SortedWithin3SE: all10OK,
		Top3CosineGt0p99:               cosOK,
		SignPatterns012:                signOK,
		AggregatePassed:                agg,
		Status:                         status,
	}
}

// RunPCAStudy menghitung eigen-dekomposisi teoretis Jacobi dan membandingkan dengan R simulasi
// sepanjang T minggu secara paralel deterministik.
func RunPCAStudy(spec kalman.ModelSpec, taus []float64, dt float64, tSteps, nRep, nChunks, nWorkers int, seed uint64) (PCAReplicationSummary, error) {
	m := len(taus)
	theoCov, approxCov := TheoreticalCov(spec, taus, dt)
	theoVals, theoVecs, _, err := lag.SymmetricEigenJacobi(theoCov, 100, 1e-14)
	if err != nil {
		return PCAReplicationSummary{}, fmt.Errorf("theoretical Jacobi failed: %w", err)
	}
	approxVals, _, _, err := lag.SymmetricEigenJacobi(approxCov, 100, 1e-14)
	if err != nil {
		return PCAReplicationSummary{}, fmt.Errorf("approx Jacobi failed: %w", err)
	}
	orientEigenvectors(theoVecs)

	theoTotalVar := 0.0
	for _, ev := range theoVals {
		theoTotalVar += ev
	}
	theoVarExp := make([]float64, m)
	for i, ev := range theoVals {
		theoVarExp[i] = ev / theoTotalVar
	}
	theoTop3Cum := theoVarExp[0] + theoVarExp[1] + theoVarExp[2]
	theoSigns := [3]int{
		CountSignChanges(theoVecs[0]),
		CountSignChanges(theoVecs[1]),
		CountSignChanges(theoVecs[2]),
	}

	type repPCA struct {
		vals []float64
		vecs [3][]float64
		cos  [3]float64
	}
	perChunk := nRep / nChunks
	rem := nRep % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) []repPCA {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		out := make([]repPCA, cnt)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x4101, uint64(c)))
		for r := 0; r < cnt; r++ {
			_, yields := afns.SimulateStationaryPath(spec.AFNS, taus, tSteps, dt, spec.SigmaEps, rg)
			sCov := SampleYieldDiffCov(yields)
			ev, eVecs, _, jErr := lag.SymmetricEigenJacobi(sCov, 100, 1e-13)
			if jErr != nil {
				panic(jErr)
			}
			var top3 [3][]float64
			var cos [3]float64
			for k := 0; k < 3; k++ {
				top3[k] = alignVectorSign(eVecs[k], theoVecs[k])
				cos[k] = stats.CosineSimilarity(top3[k], theoVecs[k])
			}
			out[r] = repPCA{vals: ev, vecs: top3, cos: cos}
		}
		return out
	})

	all := make([]repPCA, 0, nRep)
	for _, ch := range chunks {
		all = append(all, ch...)
	}

	empMeanVals := make([]float64, m)
	empSDVals := make([]float64, m)
	empMeanSE := make([]float64, m)
	zMean := make([]float64, m)
	zSingle := make([]float64, m)
	for k := 0; k < m; k++ {
		series := make([]float64, nRep)
		for r := 0; r < nRep; r++ {
			series[r] = all[r].vals[k]
		}
		empMeanVals[k] = stats.Mean(series)
		empSDVals[k] = stats.SD(series)
		empMeanSE[k] = empSDVals[k] / math.Sqrt(float64(nRep))
		zMean[k] = math.Abs(empMeanVals[k]-theoVals[k]) / math.Max(empMeanSE[k], 1e-30)
		zSingle[k] = math.Abs(all[0].vals[k]-theoVals[k]) / math.Max(empSDVals[k], 1e-30)
	}

	// Hitung rata-rata bulk noise PC4..m per replikasi (invarian rotasi tanpa bias pengurutan)
	noiseBulkTheo := 0.0
	for k := 3; k < m; k++ {
		noiseBulkTheo += theoVals[k]
	}
	noiseBulkTheo /= float64(m - 3)
	noiseBulkSeries := make([]float64, nRep)
	for r := 0; r < nRep; r++ {
		s := 0.0
		for k := 3; k < m; k++ {
			s += all[r].vals[k]
		}
		noiseBulkSeries[r] = s / float64(m-3)
	}
	noiseBulkEmp := stats.Mean(noiseBulkSeries)
	noiseBulkSE := stats.SE(noiseBulkSeries)
	noiseBulkZ := math.Abs(noiseBulkEmp-noiseBulkTheo) / math.Max(noiseBulkSE, 1e-30)

	empMeanTotal := 0.0
	for _, ev := range empMeanVals {
		empMeanTotal += ev
	}
	empVarExp := make([]float64, m)
	for i, ev := range empMeanVals {
		empVarExp[i] = ev / empMeanTotal
	}
	empTop3Cum := empVarExp[0] + empVarExp[1] + empVarExp[2]

	empMeanVecs := lag.Zeros(3, m)
	var empMeanCos, empMinCos [3]float64
	for k := 0; k < 3; k++ {
		empMinCos[k] = 1.0
		cosSum := 0.0
		for r := 0; r < nRep; r++ {
			cosSum += all[r].cos[k]
			if all[r].cos[k] < empMinCos[k] {
				empMinCos[k] = all[r].cos[k]
			}
			for j := 0; j < m; j++ {
				empMeanVecs[k][j] += all[r].vecs[k][j]
			}
		}
		empMeanCos[k] = cosSum / float64(nRep)
		norm := math.Sqrt(lag.Dot(empMeanVecs[k], empMeanVecs[k]))
		for j := 0; j < m; j++ {
			empMeanVecs[k][j] /= norm
		}
	}
	empSigns := [3]int{
		CountSignChanges(empMeanVecs[0]),
		CountSignChanges(empMeanVecs[1]),
		CountSignChanges(empMeanVecs[2]),
	}

	return PCAReplicationSummary{
		Tenors:             taus,
		TSteps:             tSteps,
		NumReps:            nRep,
		TheoEigenvalues:    theoVals,
		ApproxEigenvalues:  approxVals,
		TheoVarExplained:   theoVarExp,
		TheoTop3CumVar:     theoTop3Cum,
		TheoEigenvectors:   theoVecs,
		TheoSignChanges:    theoSigns,
		EmpMeanEigenvalues: empMeanVals,
		EmpSDEigenvalues:   empSDVals,
		EmpMeanSE:          empMeanSE,
		ZScoreSampleMean:   zMean,
		SingleRunEigenvals: all[0].vals,
		SingleRunZScores:   zSingle,
		NoiseBulkTheoMean:  noiseBulkTheo,
		NoiseBulkEmpMean:   noiseBulkEmp,
		NoiseBulkMeanSE:    noiseBulkSE,
		NoiseBulkZScore:    noiseBulkZ,
		SingleRunTop3Cos:   all[0].cos,
		EmpMeanTop3Cos:     empMeanCos,
		EmpMinTop3Cos:      empMinCos,
		EmpMeanTop3Vecs:    empMeanVecs,
		EmpSignChanges:     empSigns,
		EmpVarExplained:    empVarExp,
		EmpTop3CumVar:      empTop3Cum,
	}, nil
}
