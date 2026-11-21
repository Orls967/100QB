package bootstrap

import (
	"math"
	"time"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/afns"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/optim"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/parallel"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/stats"
)

// InnerBootDraw mencatat hasil penuh beserta diagnostik penghentian dan pencarian garis satu refit bootstrap dalam.
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
	Stage1Iterations       int
	Stage2Iterations       int
	GradNormOpt            float64
	ScaledGradNormOpt      float64
	EstLogLik              float64
	EstParams              [11]float64
	MinAlphaTested         float64
	LastStepNorm           float64
	CandidateBitIdentical  bool
	LastObjChange          float64
	ObjBitIdentical        bool
	LineSearchFailureCause string
}

// OuterRepBootstrapResult merangkum distribusi bootstrap (B refit) pada satu dataset luar r.
type OuterRepBootstrapResult struct {
	OuterRepIdx             int
	OuterEst                [11]float64
	OuterHessianSE          [11]float64
	TotalBootDraws          int
	FiniteDrawsCount        int
	InvalidDrawsCount       int
	ConvergedCount          int
	ScaledGradLe1e6Count    int
	BoundHitCount           int
	GradTolCount            int
	BoundaryStationaryCount int
	LineSearchFailedCount   int
	MaxIterCount            int
	InvalidObjCount         int
	IntervalValid           bool
	Params                  [11]OuterParamBootstrap
	InnerDraws              []InnerBootDraw
}

// BootParamSummary membandingkan MLE mentah vs bootstrap terkoreksi serta cakupan Wald vs persentil 95%.
type BootParamSummary struct {
	ParamIdx                   int
	Name                       string
	TrueValue                  float64
	TotalOuterReps             int
	ValidCorrectedOuterCount   int
	RawMean                    float64
	RawBias                    float64
	RawBiasMCSE                float64
	RawEmpiricalSD             float64
	RawRMSE                    float64
	BootMeanOfConditionalBias  float64
	BootCorrectedMean          float64
	BootCorrectedBias          float64
	BootCorrectedBiasMCSE      float64
	BootCorrectedEmpiricalSD   float64
	BootCorrectedRMSE          float64
	PairedBiasDiff             float64
	PairedBiasDiffMCSE         float64
	AbsBiasReduction           float64
	RMSEDiff                   float64
	RMSERatio                  float64
	OutOfBoundsCount           int
	OutOfBoundsPct             float64
	NonPhysicalCount           int
	NonPhysicalPct             float64
	WaldCoverage95             float64
	WaldWilsonLower            float64
	WaldWilsonUpper            float64
	WaldCoveredCount           int
	BootPercentileCoverage95   float64
	BootPercentileWilsonLower  float64
	BootPercentileWilsonUpper  float64
	BootPercentileCoveredCount int
	BootPercentileValidOuter   int
	BootPercentileInvalidOuter int
	BootPercentileMeanWidth    float64
	BootPercentileMedianWidth  float64
}

// ParametricBootstrapStudyResult menyimpan seluruh keluaran studi bootstrap parametrik R x B.
type ParametricBootstrapStudyResult struct {
	HorizonYears                 float64
	Dt                           float64
	TSteps                       int
	NumOuterReps                 int
	NumInnerBoot                 int
	TotalRefits                  int
	MinRequiredFinitePerOuter    int
	MinValidFractionThreshold    float64
	NumWorkers                   int
	ElapsedSeconds               float64
	TotalFiniteRefits            int
	TotalNonFiniteRefits         int
	TotalConvergedRefits         int
	TotalScaledGradLe1e6Refits   int
	TotalBoundHitRefits          int
	TotalGradTolRefits           int
	TotalBoundaryStationaryRefit int
	TotalLineSearchFailedRefits  int
	TotalMaxIterRefits           int
	TotalInvalidObjRefits        int
	TotalLSObjBitIdentical       int
	TotalLSCandBitIdentical      int
	TotalLSArmijoNotMet          int
	TotalLSStepUnderflow         int
	TotalLSInvalidObj            int
	ValidOuterIntervals          int
	InvalidOuterIntervals        int
	Params                       [11]BootParamSummary
	OuterReplications            []OuterRepBootstrapResult
}

func makeNaNDraw(r, b int) InnerBootDraw {
	var nanV [11]float64
	for i := range nanV {
		nanV[i] = math.NaN()
	}
	return InnerBootDraw{
		OuterRepIdx:       r,
		BootIdx:           b,
		Finite:            false,
		Converged:         false,
		TerminationReason: optim.TermInvalidObjective,
		EstLogLik:         math.NaN(),
		EstParams:         nanV,
	}
}

// RunParametricBootstrapStudy menjalankan bootstrap parametrik pada setiap dari R dataset luar
// dengan B simulasi & refit penuh 11 parameter secara paralel deterministik (chunk index = r*B + b).
func RunParametricBootstrapStudy(
	trueSpec kalman.ModelSpec,
	taus []float64,
	horizonYears, dt float64,
	tSteps int,
	outerFits []kalman.MLEResult,
	numBoot, nWorkers int,
	mainSeed, salt uint64,
) ParametricBootstrapStudyResult {
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
		rg := rng.New(rng.DeriveBootstrapSeed(mainSeed, salt, uint64(r), uint64(b)))
		_, bootYields := afns.SimulateStationaryPath(fitSpec.AFNS, taus, tSteps, dt, fitSpec.SigmaEps, rg)
		bootInit := PerturbBootstrapGuess(fitSpec, rg)

		res := kalman.EstimateMLEWithOptions(bootInit, fitSpec, taus, dt, bootYields, false)
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
			Stage1Iterations:       res.Stage1.Iterations,
			Stage2Iterations:       res.Stage2.Iterations,
			GradNormOpt:            res.GradNormOpt,
			ScaledGradNormOpt:      res.ScaledGradNormOpt,
			EstLogLik:              res.EstLogLik,
			EstParams:              res.EstParams,
			MinAlphaTested:         res.MinAlphaTested,
			LastStepNorm:           res.LastStepNorm,
			CandidateBitIdentical:  res.CandidateBitIdentical,
			LastObjChange:          res.LastObjChange,
			ObjBitIdentical:        res.ObjBitIdentical,
			LineSearchFailureCause: res.LineSearchFailureCause,
		}
	})

	trueV := kalman.SpecToVec(trueSpec)
	minReq := MinRequiredFinite(numBoot)

	outerReps := make([]OuterRepBootstrapResult, nRep)
	resOut := ParametricBootstrapStudyResult{
		HorizonYears:              horizonYears,
		Dt:                        dt,
		TSteps:                    tSteps,
		NumOuterReps:              nRep,
		NumInnerBoot:              numBoot,
		TotalRefits:               totalJobs,
		MinRequiredFinitePerOuter: minReq,
		MinValidFractionThreshold: MinValidFraction,
		NumWorkers:                nWorkers,
	}

	for r := 0; r < nRep; r++ {
		outer := outerFits[r]
		slice := draws[r*numBoot : (r+1)*numBoot]
		rep := OuterRepBootstrapResult{
			OuterRepIdx:    r,
			OuterEst:       outer.EstParams,
			OuterHessianSE: outer.HessianSE,
			TotalBootDraws: numBoot,
			InnerDraws:     slice,
		}
		for _, d := range slice {
			if d.Finite {
				rep.FiniteDrawsCount++
				resOut.TotalFiniteRefits++
			} else {
				rep.InvalidDrawsCount++
				resOut.TotalNonFiniteRefits++
			}
			if d.Converged {
				rep.ConvergedCount++
				resOut.TotalConvergedRefits++
			}
			if d.Finite && d.ScaledGradNormOpt <= 1e-6 {
				rep.ScaledGradLe1e6Count++
				resOut.TotalScaledGradLe1e6Refits++
			}
			if d.HitBound {
				rep.BoundHitCount++
				resOut.TotalBoundHitRefits++
			}
			switch d.TerminationReason {
			case optim.TermGradientTol:
				rep.GradTolCount++
				resOut.TotalGradTolRefits++
			case optim.TermBoundaryStationary:
				rep.BoundaryStationaryCount++
				resOut.TotalBoundaryStationaryRefit++
			case optim.TermLineSearchFailed:
				rep.LineSearchFailedCount++
				resOut.TotalLineSearchFailedRefits++
				if d.ObjBitIdentical {
					resOut.TotalLSObjBitIdentical++
				}
				if d.CandidateBitIdentical {
					resOut.TotalLSCandBitIdentical++
				}
				switch d.LineSearchFailureCause {
				case optim.LSCauseArmijoNotMet:
					resOut.TotalLSArmijoNotMet++
				case optim.LSCauseStepUnderflow:
					resOut.TotalLSStepUnderflow++
				case optim.LSCauseInvalidObjective:
					resOut.TotalLSInvalidObj++
				}
			case optim.TermMaxIterations:
				rep.MaxIterCount++
				resOut.TotalMaxIterRefits++
			case optim.TermInvalidObjective:
				rep.InvalidObjCount++
				resOut.TotalInvalidObjRefits++
			}
		}
		rep.IntervalValid = outer.Finite && rep.FiniteDrawsCount >= minReq
		if rep.IntervalValid {
			resOut.ValidOuterIntervals++
		} else {
			resOut.InvalidOuterIntervals++
		}

		paramVals := make([]float64, numBoot)
		for i := 0; i < 11; i++ {
			for b := 0; b < numBoot; b++ {
				if slice[b].Finite {
					paramVals[b] = slice[b].EstParams[i]
				} else {
					paramVals[b] = math.NaN()
				}
			}
			rep.Params[i] = EvalParamBootstrap(i, outer.EstParams[i], trueV[i], paramVals, numBoot)
		}
		outerReps[r] = rep
	}

	fRep := math.Max(float64(nRep), 1.0)
	for i := 0; i < 11; i++ {
		rawVals := make([]float64, 0, nRep)
		corrVals := make([]float64, 0, nRep)
		diffVals := make([]float64, 0, nRep)
		condBiases := make([]float64, 0, nRep)
		widths := make([]float64, 0, nRep)

		waldValid := make([]bool, nRep)
		waldCov := make([]bool, nRep)
		bootValid := make([]bool, nRep)
		bootCov := make([]bool, nRep)
		oobCnt, nonPhysCnt := 0, 0

		for r := 0; r < nRep; r++ {
			outer := outerFits[r]
			bp := outerReps[r].Params[i]
			if outer.Finite {
				rawVals = append(rawVals, outer.EstParams[i])
			}
			se := outer.HessianSE[i]
			waldValid[r] = outer.Finite && se > 0
			waldCov[r] = outer.Finite && se > 0 && math.Abs(outer.EstParams[i]-trueV[i]) <= stats.Z95*se

			bootValid[r] = outerReps[r].IntervalValid && bp.IntervalValid
			bootCov[r] = bp.Percentile95Covered
			if outerReps[r].IntervalValid && bp.IntervalValid && !math.IsNaN(bp.BootCorrected) {
				corrVals = append(corrVals, bp.BootCorrected)
				diffVals = append(diffVals, bp.BootCorrected-outer.EstParams[i])
				condBiases = append(condBiases, bp.BootBias)
				widths = append(widths, bp.Percentile95Upper-bp.Percentile95Lower)
				if bp.CorrectedOutOfBounds {
					oobCnt++
				}
				if bp.CorrectedNonPhysical {
					nonPhysCnt++
				}
			}
		}

		rawMean := stats.Mean(rawVals)
		rawBias := rawMean - trueV[i]
		rawSD := stats.SD(rawVals)
		rawMCSE := stats.SE(rawVals)
		rawRMSE := stats.RMSE(rawVals, trueV[i])

		corrMean := stats.Mean(corrVals)
		corrBias := corrMean - trueV[i]
		corrSD := stats.SD(corrVals)
		corrMCSE := stats.SE(corrVals)
		corrRMSE := stats.RMSE(corrVals, trueV[i])
		pairedDiff := stats.Mean(diffVals)
		pairedMCSE := stats.SE(diffVals)
		rmseRatio := 0.0
		if rawRMSE > 0 {
			rmseRatio = corrRMSE / rawRMSE
		}

		_, _, waldCnt, _, waldRate, waldLo, waldHi := stats.CoverageWithValidity(waldValid, waldCov)
		bValCnt, bInvCnt, bCovCnt, _, bRate, bLo, bHi := stats.CoverageWithValidity(bootValid, bootCov)

		meanWidth := stats.Mean(widths)
		medWidth := 0.0
		if len(widths) > 0 {
			medWidth = stats.Quantile(widths, 0.5)
		}

		resOut.Params[i] = BootParamSummary{
			ParamIdx:                   i,
			Name:                       kalman.ParamNames[i],
			TrueValue:                  trueV[i],
			TotalOuterReps:             nRep,
			ValidCorrectedOuterCount:   len(corrVals),
			RawMean:                    rawMean,
			RawBias:                    rawBias,
			RawBiasMCSE:                rawMCSE,
			RawEmpiricalSD:             rawSD,
			RawRMSE:                    rawRMSE,
			BootMeanOfConditionalBias:  stats.Mean(condBiases),
			BootCorrectedMean:          corrMean,
			BootCorrectedBias:          corrBias,
			BootCorrectedBiasMCSE:      corrMCSE,
			BootCorrectedEmpiricalSD:   corrSD,
			BootCorrectedRMSE:          corrRMSE,
			PairedBiasDiff:             pairedDiff,
			PairedBiasDiffMCSE:         pairedMCSE,
			AbsBiasReduction:           math.Abs(rawBias) - math.Abs(corrBias),
			RMSEDiff:                   corrRMSE - rawRMSE,
			RMSERatio:                  rmseRatio,
			OutOfBoundsCount:           oobCnt,
			OutOfBoundsPct:             100.0 * float64(oobCnt) / fRep,
			NonPhysicalCount:           nonPhysCnt,
			NonPhysicalPct:             100.0 * float64(nonPhysCnt) / fRep,
			WaldCoverage95:             waldRate,
			WaldWilsonLower:            waldLo,
			WaldWilsonUpper:            waldHi,
			WaldCoveredCount:           waldCnt,
			BootPercentileCoverage95:   bRate,
			BootPercentileWilsonLower:  bLo,
			BootPercentileWilsonUpper:  bHi,
			BootPercentileCoveredCount: bCovCnt,
			BootPercentileValidOuter:   bValCnt,
			BootPercentileInvalidOuter: bInvCnt,
			BootPercentileMeanWidth:    meanWidth,
			BootPercentileMedianWidth:  medWidth,
		}
	}

	resOut.ElapsedSeconds = time.Since(t0).Seconds()
	resOut.OuterReplications = outerReps
	return resOut
}
