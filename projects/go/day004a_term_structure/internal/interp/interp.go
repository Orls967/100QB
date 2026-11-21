package interp

import (
	"fmt"
	"math"
	"sort"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/bootstrap"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/parallel"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/stats"
)

// LogLinearP menginterpolasi faktor diskonto secara log-linear atas tau.
type LogLinearP struct {
	taus []float64
	lnP  []float64
}

// NewLogLinearP membangun interpolator log-linear pada P(tau) dengan jangkar P(0)=1.
func NewLogLinearP(taus, dfs []float64) (*LogLinearP, error) {
	if len(taus) == 0 || len(taus) != len(dfs) {
		return nil, fmt.Errorf("NewLogLinearP: invalid dimensions")
	}
	t := make([]float64, len(taus)+1)
	lp := make([]float64, len(dfs)+1)
	t[0] = 0
	lp[0] = 0
	for i := range taus {
		if taus[i] <= t[i] || dfs[i] <= 0 {
			return nil, fmt.Errorf("NewLogLinearP: non-increasing tau or non-positive df at %d", i)
		}
		t[i+1] = taus[i]
		lp[i+1] = math.Log(dfs[i])
	}
	return &LogLinearP{taus: t, lnP: lp}, nil
}

// Eval mengembalikan (yield kontinu y(tau), forward seketika f(tau)).
func (ll *LogLinearP) Eval(tau float64) (y, f float64) {
	n := len(ll.taus)
	if tau <= 0 {
		f0 := -(ll.lnP[1] - ll.lnP[0]) / (ll.taus[1] - ll.taus[0])
		return f0, f0
	}
	idx := sort.SearchFloat64s(ll.taus, tau) - 1
	if idx < 0 {
		idx = 0
	}
	if idx >= n-1 {
		idx = n - 2
	}
	dt := ll.taus[idx+1] - ll.taus[idx]
	slope := (ll.lnP[idx+1] - ll.lnP[idx]) / dt
	lnPTau := ll.lnP[idx] + slope*(tau-ll.taus[idx])
	f = -slope
	y = -lnPTau / tau
	return y, f
}

// NaturalCubicSplineY menginterpolasi zero yield y(tau) memakai spline kubik natural buatan sendiri.
type NaturalCubicSplineY struct {
	taus   []float64
	yields []float64
	m2     []float64 // turunan kedua S''(tau_i)
}

// NewNaturalCubicSplineY membangun spline kubik natural (S” = 0 di ujung) via solver tridiagonal Thomas.
func NewNaturalCubicSplineY(taus, yields []float64) (*NaturalCubicSplineY, error) {
	n := len(taus)
	if n < 3 || len(yields) != n {
		return nil, fmt.Errorf("NewNaturalCubicSplineY: need at least 3 knots")
	}
	t := make([]float64, n)
	y := make([]float64, n)
	copy(t, taus)
	copy(y, yields)
	h := make([]float64, n-1)
	for i := 0; i < n-1; i++ {
		h[i] = t[i+1] - t[i]
		if h[i] <= 0 {
			return nil, fmt.Errorf("NewNaturalCubicSplineY: knots must be strictly increasing")
		}
	}

	m := n - 2
	diag := make([]float64, m)
	rhs := make([]float64, m)
	sub := make([]float64, m-1)
	sup := make([]float64, m-1)
	for row := 0; row < m; row++ {
		i := row + 1
		diag[row] = 2.0 * (h[i-1] + h[i])
		rhs[row] = 6.0 * ((y[i+1]-y[i])/h[i] - (y[i]-y[i-1])/h[i-1])
		if row > 0 {
			sub[row-1] = h[i-1]
		}
		if row < m-1 {
			sup[row] = h[i]
		}
	}
	innerM2, err := lag.SolveTridiagonal(sub, diag, sup, rhs)
	if err != nil {
		return nil, err
	}
	m2 := make([]float64, n)
	copy(m2[1:n-1], innerM2)
	return &NaturalCubicSplineY{taus: t, yields: y, m2: m2}, nil
}

// Eval mengembalikan (yield kontinu y(tau), forward seketika f(tau) = y(tau) + tau * y'(tau)).
func (sp *NaturalCubicSplineY) Eval(tau float64) (y, f float64) {
	n := len(sp.taus)
	idx := sort.SearchFloat64s(sp.taus, tau) - 1
	if idx < 0 {
		idx = 0
	}
	if idx >= n-1 {
		idx = n - 2
	}
	hi := sp.taus[idx+1] - sp.taus[idx]
	a := (sp.taus[idx+1] - tau) / hi
	b := (tau - sp.taus[idx]) / hi
	y = a*sp.yields[idx] + b*sp.yields[idx+1] +
		(hi*hi/6.0)*((a*a*a-a)*sp.m2[idx]+(b*b*b-b)*sp.m2[idx+1])
	dy := (sp.yields[idx+1]-sp.yields[idx])/hi +
		(hi/6.0)*(-(3.0*a*a-1.0)*sp.m2[idx]+(3.0*b*b-1.0)*sp.m2[idx+1])
	f = y + tau*dy
	return y, f
}

// MethodMetrics menyimpan ringkasan satu metode interpolasi lintas R replikasi.
type MethodMetrics struct {
	NegForwardMeanCount float64
	NegForwardMaxCount  int
	NegForwardFracReps  float64
	ForwardTVMean       float64
	ForwardTVSD         float64
	ForwardJumpRMSMean  float64
	ForwardJumpRMSSD    float64
	ForwardD2RMSMean    float64
	ForwardRMSEMean     float64
	ForwardRMSESD       float64
	YieldRMSEMean       float64
	YieldRMSESD         float64
	SampleRepZeroGrid   []float64
	SampleRepLogLinF    []float64
	SampleRepSplineF    []float64
	SampleRepLogLinY    []float64
	SampleRepSplineY    []float64
	SampleRepEvalTaus   []float64
}

// NoiseExperimentSummary merangkum perbandingan stabilitas LogLinearP vs NaturalCubicSplineY.
type NoiseExperimentSummary struct {
	NumReps     int
	SigmaPar    float64
	NumKnots    int
	NumEvalGrid int
	LogLinearP  MethodMetrics
	CubicSpline MethodMetrics
}

type repRow struct {
	llNegCount int
	llTV       float64
	llJumpRMS  float64
	llD2RMS    float64
	llFRMSE    float64
	llYRMSE    float64
	spNegCount int
	spTV       float64
	spJumpRMS  float64
	spD2RMS    float64
	spFRMSE    float64
	spYRMSE    float64

	llF []float64
	llY []float64
	spF []float64
	spY []float64
}

func evalMetricsOnGrid(eval func(float64) (float64, float64), evalTaus, trueY, trueF []float64) (negCount int, tv, jumpRMS, d2RMS, fRMSE, yRMSE float64, ys, fs []float64) {
	m := len(evalTaus)
	ys = make([]float64, m)
	fs = make([]float64, m)
	for j, tau := range evalTaus {
		y, f := eval(tau)
		ys[j] = y
		fs[j] = f
		if f < 0 {
			negCount++
		}
	}
	for j := 1; j < m; j++ {
		df := fs[j] - fs[j-1]
		tv += math.Abs(df)
		jumpRMS += df * df
	}
	jumpRMS = math.Sqrt(jumpRMS / float64(m-1))
	for j := 1; j < m-1; j++ {
		d2 := fs[j+1] - 2.0*fs[j] + fs[j-1]
		d2RMS += d2 * d2
	}
	d2RMS = math.Sqrt(d2RMS / float64(m-2))
	fRMSE = stats.RMSE(fs, trueF)
	yRMSE = stats.RMSE(ys, trueY)
	return negCount, tv, jumpRMS, d2RMS, fRMSE, yRMSE, ys, fs
}

// RunNoisyBootstrapComparison mengevaluasi stabilitas kurva forward untuk LogLinearP dan NaturalCubicSplineY
// dari data par semi-tahunan 0.5..30 tahun dengan derau i.i.d. N(0, sigmaPar^2).
func RunNoisyBootstrapComparison(curve ns.Params, sigmaPar float64, nRep, nChunks, nWorkers int, seed uint64) NoiseExperimentSummary {
	mats := bootstrap.SemiAnnualGrid(30.0)
	parTrue, _, _ := bootstrap.ParRatesFromCurve(curve, mats, 0.5)

	// Grid evaluasi halus [0.5, 30.0] langkah 0.05 (591 titik)
	const nEval = 591
	evalTaus := make([]float64, nEval)
	trueY := make([]float64, nEval)
	trueF := make([]float64, nEval)
	for j := 0; j < nEval; j++ {
		tau := 0.5 + 0.05*float64(j)
		evalTaus[j] = tau
		trueY[j] = curve.Yield(tau)
		trueF[j] = curve.Forward(tau)
	}

	perChunk := nRep / nChunks
	rem := nRep % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) []repRow {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		out := make([]repRow, cnt)
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x1301, uint64(c)))
		noisyPar := make([]float64, len(mats))
		for r := 0; r < cnt; r++ {
			for k := range mats {
				noisyPar[k] = parTrue[k] + sigmaPar*rg.NormFloat64()
			}
			dfs, zy, err := bootstrap.BootstrapFromParRates(mats, noisyPar, 0.5)
			if err != nil {
				panic(err)
			}
			ll, err := NewLogLinearP(mats, dfs)
			if err != nil {
				panic(err)
			}
			sp, err := NewNaturalCubicSplineY(mats, zy)
			if err != nil {
				panic(err)
			}
			llNeg, llTV, llJump, llD2, llFR, llYR, llY, llF := evalMetricsOnGrid(ll.Eval, evalTaus, trueY, trueF)
			spNeg, spTV, spJump, spD2, spFR, spYR, spY, spF := evalMetricsOnGrid(sp.Eval, evalTaus, trueY, trueF)
			row := repRow{
				llNegCount: llNeg, llTV: llTV, llJumpRMS: llJump, llD2RMS: llD2, llFRMSE: llFR, llYRMSE: llYR,
				spNegCount: spNeg, spTV: spTV, spJumpRMS: spJump, spD2RMS: spD2, spFRMSE: spFR, spYRMSE: spYR,
			}
			if c == 0 && r == 0 {
				row.llF, row.llY, row.spF, row.spY = llF, llY, spF, spY
			}
			out[r] = row
		}
		return out
	})

	llNegs := make([]float64, 0, nRep)
	llTVs := make([]float64, 0, nRep)
	llJumps := make([]float64, 0, nRep)
	llD2s := make([]float64, 0, nRep)
	llFRs := make([]float64, 0, nRep)
	llYRs := make([]float64, 0, nRep)

	spNegs := make([]float64, 0, nRep)
	spTVs := make([]float64, 0, nRep)
	spJumps := make([]float64, 0, nRep)
	spD2s := make([]float64, 0, nRep)
	spFRs := make([]float64, 0, nRep)
	spYRs := make([]float64, 0, nRep)

	llMaxNeg, llHasNeg := 0, 0
	spMaxNeg, spHasNeg := 0, 0
	var sampleLLF, sampleLLY, sampleSPF, sampleSPY []float64
	for cIdx, ch := range chunks {
		for rIdx, row := range ch {
			if cIdx == 0 && rIdx == 0 {
				sampleLLF = row.llF
				sampleLLY = row.llY
				sampleSPF = row.spF
				sampleSPY = row.spY
			}
			llNegs = append(llNegs, float64(row.llNegCount))
			if row.llNegCount > llMaxNeg {
				llMaxNeg = row.llNegCount
			}
			if row.llNegCount > 0 {
				llHasNeg++
			}
			llTVs = append(llTVs, row.llTV)
			llJumps = append(llJumps, row.llJumpRMS)
			llD2s = append(llD2s, row.llD2RMS)
			llFRs = append(llFRs, row.llFRMSE)
			llYRs = append(llYRs, row.llYRMSE)

			spNegs = append(spNegs, float64(row.spNegCount))
			if row.spNegCount > spMaxNeg {
				spMaxNeg = row.spNegCount
			}
			if row.spNegCount > 0 {
				spHasNeg++
			}
			spTVs = append(spTVs, row.spTV)
			spJumps = append(spJumps, row.spJumpRMS)
			spD2s = append(spD2s, row.spD2RMS)
			spFRs = append(spFRs, row.spFRMSE)
			spYRs = append(spYRs, row.spYRMSE)
		}
	}

	return NoiseExperimentSummary{
		NumReps:     nRep,
		SigmaPar:    sigmaPar,
		NumKnots:    len(mats),
		NumEvalGrid: nEval,
		LogLinearP: MethodMetrics{
			NegForwardMeanCount: stats.Mean(llNegs),
			NegForwardMaxCount:  llMaxNeg,
			NegForwardFracReps:  float64(llHasNeg) / float64(nRep),
			ForwardTVMean:       stats.Mean(llTVs),
			ForwardTVSD:         stats.SD(llTVs),
			ForwardJumpRMSMean:  stats.Mean(llJumps),
			ForwardJumpRMSSD:    stats.SD(llJumps),
			ForwardD2RMSMean:    stats.Mean(llD2s),
			ForwardRMSEMean:     stats.Mean(llFRs),
			ForwardRMSESD:       stats.SD(llFRs),
			YieldRMSEMean:       stats.Mean(llYRs),
			YieldRMSESD:         stats.SD(llYRs),
			SampleRepEvalTaus:   evalTaus,
			SampleRepLogLinF:    sampleLLF,
			SampleRepLogLinY:    sampleLLY,
		},
		CubicSpline: MethodMetrics{
			NegForwardMeanCount: stats.Mean(spNegs),
			NegForwardMaxCount:  spMaxNeg,
			NegForwardFracReps:  float64(spHasNeg) / float64(nRep),
			ForwardTVMean:       stats.Mean(spTVs),
			ForwardTVSD:         stats.SD(spTVs),
			ForwardJumpRMSMean:  stats.Mean(spJumps),
			ForwardJumpRMSSD:    stats.SD(spJumps),
			ForwardD2RMSMean:    stats.Mean(spD2s),
			ForwardRMSEMean:     stats.Mean(spFRs),
			ForwardRMSESD:       stats.SD(spFRs),
			YieldRMSEMean:       stats.Mean(spYRs),
			YieldRMSESD:         stats.SD(spYRs),
			SampleRepEvalTaus:   evalTaus,
			SampleRepSplineF:    sampleSPF,
			SampleRepSplineY:    sampleSPY,
		},
	}
}
