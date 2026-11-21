package bootstrap

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/afns"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/interval"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/optim"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/parallel"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/rng"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/stats"
)

// TargetParamIndices adalah indeks 4 parameter fokus inferensi & kalibrasi interval:
// kp11 (1), kp22 (2), kp33 (3), theta1 (4).
var TargetParamIndices = [4]int{1, 2, 3, 4}

// PerturbBootstrapGuess memberi gangguan deterministik kecil di sekitar fit luar
// sambil menjaga titik awal strictly interior terhadap kotak batas [lo, hi].
func PerturbBootstrapGuess(outerFit kalman.ModelSpec, rg *rng.RNG) kalman.ModelSpec {
	v := kalman.SpecToVec(outerFit)
	var initV [11]float64
	for i := 0; i < 11; i++ {
		mult := rg.UniformRange(0.92, 1.08)
		if math.Abs(mult-1.0) < 0.02 {
			if mult < 1.0 {
				mult = 0.95
			} else {
				mult = 1.05
			}
		}
		lo := kalman.LowerBoundsOriginal[i]
		hi := kalman.UpperBoundsOriginal[i]
		var cand float64
		if kalman.IsPositiveParam[i] {
			cand = v[i] * mult
			minInterior := lo * 1.05
			maxInterior := hi * 0.95
			if cand < minInterior {
				cand = minInterior
			} else if cand > maxInterior {
				cand = maxInterior
			}
		} else {
			if math.Abs(v[i]) < 0.005 {
				shift := rg.UniformRange(-0.002, 0.002)
				cand = v[i] + shift
			} else {
				cand = v[i] * mult
			}
			minInterior := lo + 0.005
			maxInterior := hi - 0.005
			if cand < minInterior {
				cand = minInterior
			} else if cand > maxInterior {
				cand = maxInterior
			}
		}
		initV[i] = cand
	}
	return kalman.VecToSpec(initV, outerFit.UseAFNSAdj)
}

// RunOuterFits membangkitkan R dataset luar independen dari distribusi stasioner P
// lalu menjalankan MLE penuh 11 parameter beserta Hessian SE secara paralel deterministik.
func RunOuterFits(trueSpec kalman.ModelSpec, taus []float64, dt float64, tSteps, nRep, nWorkers int, mainSeed, salt uint64) []kalman.MLEResult {
	return parallel.RunChunks(nRep, nWorkers, func(r int) kalman.MLEResult {
		rg := rng.New(rng.DeriveSeed(mainSeed, salt, uint64(r)))
		_, y := afns.SimulateStationaryPath(trueSpec.AFNS, taus, tSteps, dt, trueSpec.SigmaEps, rg)
		initSpec := kalman.PerturbInitialGuess(trueSpec, rg)
		return kalman.EstimateMLE(initSpec, trueSpec, taus, dt, y)
	})
}

// InnerBootDraw mencatat hasil MLE dan Hessian SE* satu refit bootstrap dalam.
type InnerBootDraw struct {
	OuterRepIdx            int
	BootIdx                int
	Finite                 bool
	Converged              bool
	TerminationReason      string
	HitBound               bool
	HitBoundByParam        [11]bool
	BestStage              int
	IterationsBest         int
	IterationsCum          int
	GradNormOpt            float64
	ScaledGradNormOpt      float64
	MinAlphaTested         float64
	LastStepNorm           float64
	CandidateBitIdentical  bool
	LastObjChange          float64
	ObjBitIdentical        bool
	LineSearchFailureCause string
	EstLogLik              float64
	EstParams              [11]float64
	HessianSE              [11]float64
	HessianMinEig          float64
	HessianRegularized     bool
}

// OptAggregateDiag merangkum statistik penghentian optimizer, batas aktif, dan pencarian garis.
type OptAggregateDiag struct {
	TotalFits               int
	FiniteCount             int
	NonFiniteCount          int
	ConvergedCount          int
	ConvergedPct            float64
	ScaledGradLe1e6Count    int
	ScaledGradLe1e6Pct      float64
	BoundHitCount           int
	BoundHitPct             float64
	HessianRegularizedCount int
	GradTolCount            int
	BoundaryStationaryCount int
	LineSearchFailedCount   int
	MaxIterCount            int
	InvalidObjCount         int
	LSObjBitIdenticalCount  int
	LSArmijoNotMetCount     int
	LSStepUnderflowCount    int
	LSInvalidObjCount       int
}

// MethodSummary merangkum kalibrasi satu metode interval pada satu parameter sepanjang R replikasi luar.
type MethodSummary struct {
	Method               string
	Applicable           bool
	TotalOuterReps       int
	ValidCount           int
	InvalidCount         int
	CoveredCount         int
	Coverage95           float64
	WilsonLower          float64
	WilsonUpper          float64
	MeanWidth            float64
	MedianWidth          float64
	NonPhysicalCount     int
	NonPhysicalPct       float64
	OutOfBoundsCount     int
	OutOfBoundsPct       float64
	InvalidOuterNonFin   int
	InvalidOuterSENonPos int
	InvalidOuterNonPosLg int
	InvalidInsuffFinite  int
	InvalidInsuffValidSE int
	InvalidInsuffPos     int
}

// OuterRepCalibration menyimpan seluruh hasil interval dan ringkasan bootstrap pada satu dataset luar r.
type OuterRepCalibration struct {
	OuterRepIdx        int
	OuterMLE           kalman.MLEResult
	FiniteDrawsCount   int
	InvalidDrawsCount  int
	BootBias           [11]float64
	BootCorrected      [11]float64
	WaldIntervals      [11]interval.IntervalResult
	PctIntervals       [11]interval.IntervalResult
	BasicIntervals     [11]interval.IntervalResult
	StudIntervals      [11]interval.IntervalResult
	LogBasicIntervals  [11]interval.IntervalResult
	InnerDraws         []InnerBootDraw
}

// ParamCalibrationSummary merangkum bias dan perbandingan kelima metode interval untuk satu parameter.
type ParamCalibrationSummary struct {
	ParamIdx                 int
	Name                     string
	TrueValue                float64
	TotalOuterReps           int
	ValidCorrectedOuterCount int
	RawMean                  float64
	RawBias                  float64
	RawBiasMCSE              float64
	RawEmpiricalSD           float64
	RawRMSE                  float64
	RawBoundHitCount         int
	RawBoundHitPct           float64
	BootMeanConditionalBias  float64
	BootCorrectedMean        float64
	BootCorrectedBias        float64
	BootCorrectedBiasMCSE    float64
	BootCorrectedEmpiricalSD float64
	BootCorrectedRMSE        float64
	RMSERatio                float64
	Wald                     MethodSummary
	Percentile               MethodSummary
	Basic                    MethodSummary
	Studentized              MethodSummary
	LogBasic                 MethodSummary
}

// StratificationSummary mengukur sensitivitas lebar dan cakupan theta1 terhadap kp11 luar.
type StratificationSummary struct {
	MedianOuterKP11       float64
	LowKP11Count          int
	HighKP11Count         int
	BoundHitOuterCount    int
	LowKP11Theta1WaldCov  float64
	HighKP11Theta1WaldCov float64
	LowKP11Theta1StudCov  float64
	HighKP11Theta1StudCov float64
	LowKP11Theta1StudW    float64
	HighKP11Theta1StudW   float64
	LowKP11Theta1WaldW    float64
	HighKP11Theta1WaldW   float64
	CorrKP11WithWaldW     float64
	CorrKP11WithStudW     float64
}

// CalibrationStudyResult menyimpan seluruh hasil studi kalibrasi interval R x B.
type CalibrationStudyResult struct {
	HorizonYears              float64
	Dt                        float64
	TSteps                    int
	NumOuterReps              int
	NumInnerBoot              int
	TotalInnerRefits          int
	MinRequiredValidPerOuter  int
	MinValidFractionThreshold float64
	NumWorkers                int
	OuterElapsedSeconds       float64
	InnerElapsedSeconds       float64
	OuterOptDiag              OptAggregateDiag
	InnerOptDiag              OptAggregateDiag
	Params                    [11]ParamCalibrationSummary
	Stratification            StratificationSummary
	OuterReplications         []OuterRepCalibration
}

func makeNaNDraw(r, b int) InnerBootDraw {
	var nanV [11]float64
	var nanSE [11]float64
	for i := 0; i < 11; i++ {
		nanV[i] = math.NaN()
		nanSE[i] = math.NaN()
	}
	return InnerBootDraw{
		OuterRepIdx:            r,
		BootIdx:                b,
		Finite:                 false,
		Converged:              false,
		TerminationReason:      optim.TermInvalidObjective,
		LineSearchFailureCause: optim.LSCauseInvalidObjective,
		EstLogLik:              math.NaN(),
		EstParams:              nanV,
		HessianSE:              nanSE,
	}
}

func summarizeMethod(method string, applicable bool, results []interval.IntervalResult) MethodSummary {
	n := len(results)
	ms := MethodSummary{
		Method:         method,
		Applicable:     applicable,
		TotalOuterReps: n,
	}
	if !applicable || n == 0 {
		return ms
	}
	validFlags := make([]bool, n)
	covFlags := make([]bool, n)
	widths := make([]float64, 0, n)
	for r, ir := range results {
		validFlags[r] = ir.Valid
		covFlags[r] = ir.Covered
		if ir.Valid {
			widths = append(widths, ir.Width)
			if ir.NonPhysical {
				ms.NonPhysicalCount++
			}
			if ir.OutOfBounds {
				ms.OutOfBoundsCount++
			}
		} else {
			switch ir.InvalidReason {
			case interval.ReasonOuterNonFinite:
				ms.InvalidOuterNonFin++
			case interval.ReasonOuterSENonPositive:
				ms.InvalidOuterSENonPos++
			case interval.ReasonOuterNonPositiveLog:
				ms.InvalidOuterNonPosLg++
			case interval.ReasonInsufficientFinite:
				ms.InvalidInsuffFinite++
			case interval.ReasonInsufficientValidSE:
				ms.InvalidInsuffValidSE++
			case interval.ReasonInsufficientPositive:
				ms.InvalidInsuffPos++
			}
		}
	}
	vCnt, invCnt, covCnt, tot, rate, wLo, wHi := stats.CoverageWithValidity(validFlags, covFlags)
	ms.ValidCount = vCnt
	ms.InvalidCount = invCnt
	ms.CoveredCount = covCnt
	ms.Coverage95 = rate
	ms.WilsonLower = wLo
	ms.WilsonUpper = wHi
	ms.MeanWidth = stats.Mean(widths)
	if len(widths) > 0 {
		ms.MedianWidth = stats.Quantile(widths, 0.5)
	}
	ft := math.Max(float64(tot), 1.0)
	ms.NonPhysicalPct = 100.0 * float64(ms.NonPhysicalCount) / ft
	ms.OutOfBoundsPct = 100.0 * float64(ms.OutOfBoundsCount) / ft
	return ms
}

func pearsonCorr(xs, ys []float64) float64 {
	if len(xs) < 2 || len(xs) != len(ys) {
		return 0
	}
	mx, my := stats.Mean(xs), stats.Mean(ys)
	sxy, sxx, syy := 0.0, 0.0, 0.0
	for i := range xs {
		dx := xs[i] - mx
		dy := ys[i] - my
		sxy += dx * dy
		sxx += dx * dx
		syy += dy * dy
	}
	if sxx <= 0 || syy <= 0 {
		return 0
	}
	return sxy / math.Sqrt(sxx*syy)
}

func finalizeOptDiag(d *OptAggregateDiag) {
	fn := math.Max(float64(d.TotalFits), 1.0)
	d.ConvergedPct = 100.0 * float64(d.ConvergedCount) / fn
	d.ScaledGradLe1e6Pct = 100.0 * float64(d.ScaledGradLe1e6Count) / fn
	d.BoundHitPct = 100.0 * float64(d.BoundHitCount) / fn
}

func accumulateOptDiag(d *OptAggregateDiag, finite, converged, hitBound, hessReg bool, sgrad float64, termReason, lsCause string, objBitIdent bool) {
	d.TotalFits++
	if finite {
		d.FiniteCount++
	} else {
		d.NonFiniteCount++
	}
	if converged {
		d.ConvergedCount++
	}
	if finite && sgrad <= 1e-6 {
		d.ScaledGradLe1e6Count++
	}
	if hitBound {
		d.BoundHitCount++
	}
	if hessReg {
		d.HessianRegularizedCount++
	}
	switch termReason {
	case optim.TermGradientTol:
		d.GradTolCount++
	case optim.TermBoundaryStationary:
		d.BoundaryStationaryCount++
	case optim.TermLineSearchFailed:
		d.LineSearchFailedCount++
		if objBitIdent {
			d.LSObjBitIdenticalCount++
		}
		switch lsCause {
		case optim.LSCauseArmijoNotMet:
			d.LSArmijoNotMetCount++
		case optim.LSCauseStepUnderflow:
			d.LSStepUnderflowCount++
		case optim.LSCauseInvalidObjective:
			d.LSInvalidObjCount++
		}
	case optim.TermMaxIterations:
		d.MaxIterCount++
	case optim.TermInvalidObjective:
		d.InvalidObjCount++
	}
}

// RunCalibrationStudy mengeksekusi R*B simulasi bootstrap parametrik beserta Hessian SE*
// lalu mengevaluasi kelima metode interval secara berpasangan pada dataset luar dan draw bootstrap yang sama.
func RunCalibrationStudy(
	trueSpec kalman.ModelSpec,
	taus []float64,
	horizonYears, dt float64,
	tSteps int,
	outerFits []kalman.MLEResult,
	outerElapsed float64,
	numBoot, nWorkers int,
	mainSeed, bootSalt uint64,
) CalibrationStudyResult {
	t0 := time.Now()
	nRep := len(outerFits)
	totalJobs := nRep * numBoot

	draws := parallel.RunChunks(totalJobs, nWorkers, func(jobIdx int) InnerBootDraw {
		r := jobIdx / numBoot
		b := jobIdx % numBoot
		outer := outerFits[r]
		if !outer.Finite {
			return makeNaNDraw(r, b)
		}
		fitSpec := kalman.VecToSpec(outer.EstParams, trueSpec.UseAFNSAdj)
		rg := rng.New(rng.DeriveBootstrapSeed(mainSeed, bootSalt, uint64(r), uint64(b)))
		_, bootYields := afns.SimulateStationaryPath(fitSpec.AFNS, taus, tSteps, dt, fitSpec.SigmaEps, rg)
		bootInit := PerturbBootstrapGuess(fitSpec, rg)

		res := kalman.EstimateMLEWithOptions(bootInit, fitSpec, taus, dt, bootYields, true)
		return InnerBootDraw{
			OuterRepIdx:            r,
			BootIdx:                b,
			Finite:                 res.Finite,
			Converged:              res.Converged,
			TerminationReason:      res.TerminationReason,
			HitBound:               res.HitBound,
			HitBoundByParam:        res.HitBoundByParam,
			BestStage:              res.BestStage,
			IterationsBest:         res.BestStageIterations,
			IterationsCum:          res.CumulativeIterations,
			GradNormOpt:            res.GradNormOpt,
			ScaledGradNormOpt:      res.ScaledGradNormOpt,
			MinAlphaTested:         res.MinAlphaTested,
			LastStepNorm:           res.LastStepNorm,
			CandidateBitIdentical:  res.CandidateBitIdentical,
			LastObjChange:          res.LastObjChange,
			ObjBitIdentical:        res.ObjBitIdentical,
			LineSearchFailureCause: res.LineSearchFailureCause,
			EstLogLik:              res.EstLogLik,
			EstParams:              res.EstParams,
			HessianSE:              res.HessianSE,
			HessianMinEig:          res.HessianMinEig,
			HessianRegularized:     res.HessianRegularized,
		}
	})

	trueV := kalman.SpecToVec(trueSpec)
	minReq := interval.MinRequiredValid(numBoot)

	var outerDiag, innerDiag OptAggregateDiag
	outerReps := make([]OuterRepCalibration, nRep)

	for r := 0; r < nRep; r++ {
		outer := outerFits[r]
		accumulateOptDiag(&outerDiag, outer.Finite, outer.Converged, outer.HitBound, outer.HessianRegularized,
			outer.ScaledGradNormOpt, outer.TerminationReason, outer.LineSearchFailureCause, outer.ObjBitIdentical)

		slice := draws[r*numBoot : (r+1)*numBoot]
		rep := OuterRepCalibration{
			OuterRepIdx: r,
			OuterMLE:    outer,
			InnerDraws:  slice,
		}
		for _, d := range slice {
			if d.Finite {
				rep.FiniteDrawsCount++
			} else {
				rep.InvalidDrawsCount++
			}
			accumulateOptDiag(&innerDiag, d.Finite, d.Converged, d.HitBound, d.HessianRegularized,
				d.ScaledGradNormOpt, d.TerminationReason, d.LineSearchFailureCause, d.ObjBitIdentical)
		}

		bootVals := make([]float64, numBoot)
		bootSEs := make([]float64, numBoot)
		for i := 0; i < 11; i++ {
			finiteVals := make([]float64, 0, numBoot)
			for b := 0; b < numBoot; b++ {
				if slice[b].Finite {
					bootVals[b] = slice[b].EstParams[i]
					bootSEs[b] = slice[b].HessianSE[i]
					finiteVals = append(finiteVals, slice[b].EstParams[i])
				} else {
					bootVals[b] = math.NaN()
					bootSEs[b] = math.NaN()
				}
			}
			if outer.Finite && len(finiteVals) >= minReq {
				bBias := stats.Mean(finiteVals) - outer.EstParams[i]
				rep.BootBias[i] = bBias
				rep.BootCorrected[i] = outer.EstParams[i] - bBias
			} else {
				rep.BootBias[i] = math.NaN()
				rep.BootCorrected[i] = math.NaN()
			}
			rep.WaldIntervals[i] = interval.EvalWald(i, outer.EstParams[i], outer.HessianSE[i], trueV[i], outer.Finite)
			rep.PctIntervals[i] = interval.EvalPercentile(i, trueV[i], outer.Finite, bootVals, numBoot)
			rep.BasicIntervals[i] = interval.EvalBasic(i, outer.EstParams[i], trueV[i], outer.Finite, bootVals, numBoot)
			rep.StudIntervals[i] = interval.EvalStudentized(i, outer.EstParams[i], outer.HessianSE[i], trueV[i], outer.Finite, bootVals, bootSEs, numBoot)
			if kalman.IsPositiveParam[i] {
				rep.LogBasicIntervals[i] = interval.EvalLogBasic(i, outer.EstParams[i], trueV[i], outer.Finite, bootVals, numBoot)
			}
		}
		outerReps[r] = rep
	}

	finalizeOptDiag(&outerDiag)
	finalizeOptDiag(&innerDiag)

	var paramSummaries [11]ParamCalibrationSummary
	fRep := math.Max(float64(nRep), 1.0)
	for i := 0; i < 11; i++ {
		rawVals := make([]float64, 0, nRep)
		corrVals := make([]float64, 0, nRep)
		condBiases := make([]float64, 0, nRep)
		waldRes := make([]interval.IntervalResult, nRep)
		pctRes := make([]interval.IntervalResult, nRep)
		basicRes := make([]interval.IntervalResult, nRep)
		studRes := make([]interval.IntervalResult, nRep)
		logRes := make([]interval.IntervalResult, nRep)
		hitRawCnt := 0

		for r := 0; r < nRep; r++ {
			rep := outerReps[r]
			if rep.OuterMLE.Finite {
				rawVals = append(rawVals, rep.OuterMLE.EstParams[i])
			}
			if rep.OuterMLE.HitBoundByParam[i] {
				hitRawCnt++
			}
			if !math.IsNaN(rep.BootCorrected[i]) {
				corrVals = append(corrVals, rep.BootCorrected[i])
				condBiases = append(condBiases, rep.BootBias[i])
			}
			waldRes[r] = rep.WaldIntervals[i]
			pctRes[r] = rep.PctIntervals[i]
			basicRes[r] = rep.BasicIntervals[i]
			studRes[r] = rep.StudIntervals[i]
			logRes[r] = rep.LogBasicIntervals[i]
		}

		rawMean := stats.Mean(rawVals)
		rawRMSE := stats.RMSE(rawVals, trueV[i])
		corrMean := stats.Mean(corrVals)
		corrRMSE := stats.RMSE(corrVals, trueV[i])
		rmseRatio := 0.0
		if rawRMSE > 0 {
			rmseRatio = corrRMSE / rawRMSE
		}

		paramSummaries[i] = ParamCalibrationSummary{
			ParamIdx:                 i,
			Name:                     kalman.ParamNames[i],
			TrueValue:                trueV[i],
			TotalOuterReps:           nRep,
			ValidCorrectedOuterCount: len(corrVals),
			RawMean:                  rawMean,
			RawBias:                  rawMean - trueV[i],
			RawBiasMCSE:              stats.SE(rawVals),
			RawEmpiricalSD:           stats.SD(rawVals),
			RawRMSE:                  rawRMSE,
			RawBoundHitCount:         hitRawCnt,
			RawBoundHitPct:           100.0 * float64(hitRawCnt) / fRep,
			BootMeanConditionalBias:  stats.Mean(condBiases),
			BootCorrectedMean:        corrMean,
			BootCorrectedBias:        corrMean - trueV[i],
			BootCorrectedBiasMCSE:    stats.SE(corrVals),
			BootCorrectedEmpiricalSD: stats.SD(corrVals),
			BootCorrectedRMSE:        corrRMSE,
			RMSERatio:                rmseRatio,
			Wald:                     summarizeMethod(interval.MethodWald, true, waldRes),
			Percentile:               summarizeMethod(interval.MethodPercentile, true, pctRes),
			Basic:                    summarizeMethod(interval.MethodBasic, true, basicRes),
			Studentized:              summarizeMethod(interval.MethodStudentized, true, studRes),
			LogBasic:                 summarizeMethod(interval.MethodLogBasic, kalman.IsPositiveParam[i], logRes),
		}
	}

	// Stratifikasi sensitivitas theta1 menurut kp11_hat luar (<= median vs > median)
	kp11Vals := make([]float64, 0, nRep)
	for _, rep := range outerReps {
		if rep.OuterMLE.Finite {
			kp11Vals = append(kp11Vals, rep.OuterMLE.EstParams[1])
		}
	}
	medKP11 := stats.Quantile(kp11Vals, 0.5)
	var strat StratificationSummary
	strat.MedianOuterKP11 = medKP11
	var loWCov, hiWCov, loSCov, hiSCov int
	var loWW, hiWW, loSW, hiSW []float64
	var xKP, yWW, ySW []float64
	for _, rep := range outerReps {
		if rep.OuterMLE.HitBound {
			strat.BoundHitOuterCount++
		}
		if !rep.OuterMLE.Finite {
			continue
		}
		kp := rep.OuterMLE.EstParams[1]
		wi := rep.WaldIntervals[4]
		si := rep.StudIntervals[4]
		if wi.Valid && si.Valid {
			xKP = append(xKP, kp)
			yWW = append(yWW, wi.Width)
			ySW = append(ySW, si.Width)
		}
		if kp <= medKP11 {
			strat.LowKP11Count++
			if wi.Covered {
				loWCov++
			}
			if si.Covered {
				loSCov++
			}
			if wi.Valid {
				loWW = append(loWW, wi.Width)
			}
			if si.Valid {
				loSW = append(loSW, si.Width)
			}
		} else {
			strat.HighKP11Count++
			if wi.Covered {
				hiWCov++
			}
			if si.Covered {
				hiSCov++
			}
			if wi.Valid {
				hiWW = append(hiWW, wi.Width)
			}
			if si.Valid {
				hiSW = append(hiSW, si.Width)
			}
		}
	}
	if strat.LowKP11Count > 0 {
		strat.LowKP11Theta1WaldCov = float64(loWCov) / float64(strat.LowKP11Count)
		strat.LowKP11Theta1StudCov = float64(loSCov) / float64(strat.LowKP11Count)
	}
	if strat.HighKP11Count > 0 {
		strat.HighKP11Theta1WaldCov = float64(hiWCov) / float64(strat.HighKP11Count)
		strat.HighKP11Theta1StudCov = float64(hiSCov) / float64(strat.HighKP11Count)
	}
	strat.LowKP11Theta1WaldW = stats.Mean(loWW)
	strat.HighKP11Theta1WaldW = stats.Mean(hiWW)
	strat.LowKP11Theta1StudW = stats.Mean(loSW)
	strat.HighKP11Theta1StudW = stats.Mean(hiSW)
	strat.CorrKP11WithWaldW = pearsonCorr(xKP, yWW)
	strat.CorrKP11WithStudW = pearsonCorr(xKP, ySW)

	return CalibrationStudyResult{
		HorizonYears:              horizonYears,
		Dt:                        dt,
		TSteps:                    tSteps,
		NumOuterReps:              nRep,
		NumInnerBoot:              numBoot,
		TotalInnerRefits:          totalJobs,
		MinRequiredValidPerOuter:  minReq,
		MinValidFractionThreshold: interval.MinValidFraction,
		NumWorkers:                nWorkers,
		OuterElapsedSeconds:       outerElapsed,
		InnerElapsedSeconds:       time.Since(t0).Seconds(),
		OuterOptDiag:              outerDiag,
		InnerOptDiag:              innerDiag,
		Params:                    paramSummaries,
		Stratification:            strat,
		OuterReplications:         outerReps,
	}
}
