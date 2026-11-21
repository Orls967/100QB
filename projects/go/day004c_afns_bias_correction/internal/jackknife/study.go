package jackknife

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

// OuterPath menyimpan satu lintasan panjang penuh bersama tebakan awal acak untuk studi berpasangan.
type OuterPath struct {
	RepIdx   int
	Yields   [][]float64
	InitSpec kalman.ModelSpec
}

// GenerateSharedPrefixOuterPaths membangkitkan R lintasan stasioner panjang maxSteps secara deterministik.
func GenerateSharedPrefixOuterPaths(trueSpec kalman.ModelSpec, taus []float64, maxSteps int, dt float64, nRep, nWorkers int, mainSeed, salt uint64) []OuterPath {
	return parallel.RunChunks(nRep, nWorkers, func(r int) OuterPath {
		rg := rng.New(rng.DeriveSeed(mainSeed, salt, uint64(r)))
		_, y := afns.SimulateStationaryPath(trueSpec.AFNS, taus, maxSteps, dt, trueSpec.SigmaEps, rg)
		initSpec := kalman.PerturbInitialGuess(trueSpec, rg)
		return OuterPath{RepIdx: r, Yields: y, InitSpec: initSpec}
	})
}

// SegmentDiagnostics merangkum statistik konvergensi dan alasan berhenti optimizer satu segmen waktu.
type SegmentDiagnostics struct {
	Segment                 string
	TotalReps               int
	FiniteCount             int
	ConvergedCount          int
	ConvergedPct            float64
	ScaledGradLe1e6Count    int
	ScaledGradLe1e6Pct      float64
	LLGeTrueCount           int
	BoundHitCount           int
	BoundHitPct             float64
	HessianRegCount         int
	GradTolCount            int
	BoundaryStationaryCount int
	LineSearchFailedCount   int
	MaxIterCount            int
	InvalidObjCount         int
	BestStage1Count         int
	BestStage2Count         int
}

func summarizeSegment(name string, reps []kalman.MLEResult) SegmentDiagnostics {
	n := len(reps)
	d := SegmentDiagnostics{Segment: name, TotalReps: n}
	for _, r := range reps {
		if r.Finite {
			d.FiniteCount++
		}
		if r.Converged {
			d.ConvergedCount++
		}
		if r.ScaledGradNormOpt <= 1e-6 {
			d.ScaledGradLe1e6Count++
		}
		if r.Finite && r.EstLogLik >= r.TrueLogLik-1e-8 {
			d.LLGeTrueCount++
		}
		if r.HitBound {
			d.BoundHitCount++
		}
		if r.HessianRegularized {
			d.HessianRegCount++
		}
		if r.BestStage == 1 {
			d.BestStage1Count++
		} else if r.BestStage == 2 {
			d.BestStage2Count++
		}
		switch r.TerminationReason {
		case optim.TermGradientTol:
			d.GradTolCount++
		case optim.TermBoundaryStationary:
			d.BoundaryStationaryCount++
		case optim.TermLineSearchFailed:
			d.LineSearchFailedCount++
		case optim.TermMaxIterations:
			d.MaxIterCount++
		case optim.TermInvalidObjective:
			d.InvalidObjCount++
		}
	}
	fn := math.Max(float64(n), 1.0)
	d.ConvergedPct = 100.0 * float64(d.ConvergedCount) / fn
	d.ScaledGradLe1e6Pct = 100.0 * float64(d.ScaledGradLe1e6Count) / fn
	d.BoundHitPct = 100.0 * float64(d.BoundHitCount) / fn
	return d
}

// PairedParamSummary membandingkan MLE mentah dan jackknife temporal secara berpasangan pada satu parameter.
type PairedParamSummary struct {
	ParamIdx                       int
	Name                           string
	TrueValue                      float64
	TotalReps                      int
	RawMean                        float64
	RawBias                        float64
	RawBiasMCSE                    float64
	RawEmpiricalSD                 float64
	RawRMSE                        float64
	RawMeanHessianSE               float64
	RawMedianHessianSE             float64
	RawMeanHessianSEToEmpSDRatio   float64
	RawMedianHessianSEToEmpSDRatio float64
	WaldCoverage95                 float64
	WaldWilsonLower                float64
	WaldWilsonUpper                float64
	WaldCoveredCount               int
	RawBoundHitCount               int
	RawBoundHitPct                 float64
	JackknifeMean                  float64
	JackknifeBias                  float64
	JackknifeBiasMCSE              float64
	JackknifeEmpiricalSD           float64
	JackknifeRMSE                  float64
	JackknifeOutOfBoundsCount      int
	JackknifeOutOfBoundsPct        float64
	JackknifeNonPhysicalCount      int
	JackknifeNonPhysicalPct        float64
	PairedBiasDiff                 float64 // Bias_JK - Bias_Raw = Mean(est_JK - est_Raw)
	PairedBiasDiffMCSE             float64
	AbsBiasReduction               float64 // |Bias_Raw| - |Bias_JK| (>0 berarti bias absolut mengecil)
	RMSEDiff                       float64 // RMSE_JK - RMSE_Raw
	RMSERatio                      float64 // RMSE_JK / RMSE_Raw
}

// RepJackknifeResult menyimpan seluruh hasil penuh, dua setengah sampel, dan jackknife satu replikasi.
type RepJackknifeResult struct {
	RepIdx             int
	Full               kalman.MLEResult
	Half1              kalman.MLEResult
	Half2              kalman.MLEResult
	JackknifeEst       [11]float64
	OutOfBoundsByParam [11]bool
	NonPhysicalByParam [11]bool
}

// HorizonJackknifeResult merangkum studi jackknife temporal untuk satu horizon kalender.
type HorizonJackknifeResult struct {
	Name           string
	HorizonYears   float64
	Dt             float64
	TSteps         int
	HalfTSteps     int
	NumReps        int
	NumWorkers     int
	ElapsedSeconds float64
	FullDiag       SegmentDiagnostics
	Half1Diag      SegmentDiagnostics
	Half2Diag      SegmentDiagnostics
	Params         [11]PairedParamSummary
	Replications   []RepJackknifeResult
}

// RunJackknifeHorizonStudy mengeksekusi studi koreksi jackknife temporal pada horizon H (T langkah)
// dengan membagi lintasan menjadi dua segmen berurutan [0:T/2] dan [T/2:T] secara paralel deterministik.
func RunJackknifeHorizonStudy(name string, trueSpec ModelSpecOrAlias, taus []float64, paths []OuterPath, horizonYears float64, dt float64, tSteps, nWorkers int) HorizonJackknifeResult {
	// catatan: ModelSpecOrAlias memastikan tipe bersih tanpa duplikasi
	return runHorizonStudyImpl(name, trueSpec, taus, paths, horizonYears, dt, tSteps, nWorkers)
}

// ModelSpecOrAlias adalah alias kalman.ModelSpec agar tanda tangan fungsi ringkas.
type ModelSpecOrAlias = kalman.ModelSpec

func runHorizonStudyImpl(name string, trueSpec kalman.ModelSpec, taus []float64, paths []OuterPath, horizonYears float64, dt float64, tSteps, nWorkers int) HorizonJackknifeResult {
	t0 := time.Now()
	nRep := len(paths)
	halfSteps := tSteps / 2

	reps := parallel.RunChunks(nRep, nWorkers, func(r int) RepJackknifeResult {
		p := paths[r]
		fullRes := kalman.EstimateMLE(p.InitSpec, trueSpec, taus, dt, p.Yields[:tSteps])
		h1Res := kalman.EstimateMLE(p.InitSpec, trueSpec, taus, dt, p.Yields[:halfSteps])
		h2Res := kalman.EstimateMLE(p.InitSpec, trueSpec, taus, dt, p.Yields[halfSteps:tSteps])
		jk := CombineVec(fullRes.EstParams, h1Res.EstParams, h2Res.EstParams)
		return RepJackknifeResult{
			RepIdx:             r,
			Full:               fullRes,
			Half1:              h1Res,
			Half2:              h2Res,
			JackknifeEst:       jk.EstParams,
			OutOfBoundsByParam: jk.OutOfBoundsByParam,
			NonPhysicalByParam: jk.NonPhysicalByParam,
		}
	})

	fullSlice := make([]kalman.MLEResult, nRep)
	h1Slice := make([]kalman.MLEResult, nRep)
	h2Slice := make([]kalman.MLEResult, nRep)
	for r, rep := range reps {
		fullSlice[r] = rep.Full
		h1Slice[r] = rep.Half1
		h2Slice[r] = rep.Half2
	}

	trueV := kalman.SpecToVec(trueSpec)
	var summaries [11]PairedParamSummary
	fRep := math.Max(float64(nRep), 1.0)

	for i := 0; i < 11; i++ {
		rawVals := make([]float64, nRep)
		jkVals := make([]float64, nRep)
		diffVals := make([]float64, nRep)
		hessVals := make([]float64, nRep)
		covCnt, hitRawCnt, oobCnt, nonPhysCnt := 0, 0, 0, 0
		for r, rep := range reps {
			rv := rep.Full.EstParams[i]
			jv := rep.JackknifeEst[i]
			se := rep.Full.HessianSE[i]
			rawVals[r] = rv
			jkVals[r] = jv
			diffVals[r] = jv - rv
			hessVals[r] = se
			if rep.Full.HitBoundByParam[i] {
				hitRawCnt++
			}
			if se > 0 && math.Abs(rv-trueV[i]) <= stats.Z95*se {
				covCnt++
			}
			if rep.OutOfBoundsByParam[i] {
				oobCnt++
			}
			if rep.NonPhysicalByParam[i] {
				nonPhysCnt++
			}
		}
		rawMean := stats.Mean(rawVals)
		rawBias := rawMean - trueV[i]
		rawSD := stats.SD(rawVals)
		rawMCSE := stats.SE(rawVals)
		rawRMSE := stats.RMSE(rawVals, trueV[i])
		meanHSE := stats.Mean(hessVals)
		medHSE := stats.Quantile(hessVals, 0.5)
		meanHRatio, medHRatio := 0.0, 0.0
		if rawSD > 0 {
			meanHRatio = meanHSE / rawSD
			medHRatio = medHSE / rawSD
		}
		cov, wLo, wHi := stats.WilsonCI(covCnt, nRep, stats.Z95)

		jkMean := stats.Mean(jkVals)
		jkBias := jkMean - trueV[i]
		jkSD := stats.SD(jkVals)
		jkMCSE := stats.SE(jkVals)
		jkRMSE := stats.RMSE(jkVals, trueV[i])
		pairedDiff := stats.Mean(diffVals)
		pairedMCSE := stats.SE(diffVals)
		rmseRatio := 0.0
		if rawRMSE > 0 {
			rmseRatio = jkRMSE / rawRMSE
		}

		summaries[i] = PairedParamSummary{
			ParamIdx:                       i,
			Name:                           kalman.ParamNames[i],
			TrueValue:                      trueV[i],
			TotalReps:                      nRep,
			RawMean:                        rawMean,
			RawBias:                        rawBias,
			RawBiasMCSE:                    rawMCSE,
			RawEmpiricalSD:                 rawSD,
			RawRMSE:                        rawRMSE,
			RawMeanHessianSE:               meanHSE,
			RawMedianHessianSE:             medHSE,
			RawMeanHessianSEToEmpSDRatio:   meanHRatio,
			RawMedianHessianSEToEmpSDRatio: medHRatio,
			WaldCoverage95:                 cov,
			WaldWilsonLower:                wLo,
			WaldWilsonUpper:                wHi,
			WaldCoveredCount:               covCnt,
			RawBoundHitCount:               hitRawCnt,
			RawBoundHitPct:                 100.0 * float64(hitRawCnt) / fRep,
			JackknifeMean:                  jkMean,
			JackknifeBias:                  jkBias,
			JackknifeBiasMCSE:              jkMCSE,
			JackknifeEmpiricalSD:           jkSD,
			JackknifeRMSE:                  jkRMSE,
			JackknifeOutOfBoundsCount:      oobCnt,
			JackknifeOutOfBoundsPct:        100.0 * float64(oobCnt) / fRep,
			JackknifeNonPhysicalCount:      nonPhysCnt,
			JackknifeNonPhysicalPct:        100.0 * float64(nonPhysCnt) / fRep,
			PairedBiasDiff:                 pairedDiff,
			PairedBiasDiffMCSE:             pairedMCSE,
			AbsBiasReduction:               math.Abs(rawBias) - math.Abs(jkBias),
			RMSEDiff:                       jkRMSE - rawRMSE,
			RMSERatio:                      rmseRatio,
		}
	}

	return HorizonJackknifeResult{
		Name:           name,
		HorizonYears:   horizonYears,
		Dt:             dt,
		TSteps:         tSteps,
		HalfTSteps:     halfSteps,
		NumReps:        nRep,
		NumWorkers:     nWorkers,
		ElapsedSeconds: time.Since(t0).Seconds(),
		FullDiag:       summarizeSegment("full", fullSlice),
		Half1Diag:      summarizeSegment("half1", h1Slice),
		Half2Diag:      summarizeSegment("half2", h2Slice),
		Params:         summaries,
		Replications:   reps,
	}
}
