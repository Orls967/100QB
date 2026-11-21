package bootstrap

import (
	"math"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/stats"
)

// MinValidFraction adalah ambang pra-run proporsi minimal refit bootstrap bernilai finite
// agar interval persentil 95% satu replikasi luar dinyatakan valid (B=100 -> >= 95 finite).
const MinValidFraction = 0.95

// MinRequiredFinite mengembalikan bilangan bulat minimum refit finite agar interval valid.
func MinRequiredFinite(numBoot int) int {
	if numBoot <= 0 {
		return 1
	}
	req := int(math.Ceil(MinValidFraction * float64(numBoot)))
	if req < 1 {
		req = 1
	}
	return req
}

// PerturbBootstrapGuess memberi gangguan deterministik kecil di sekitar parameter fit luar
// sambil menjaga titik awal strictly interior terhadap kotak pengaman [lo, hi].
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

// OuterParamBootstrap menyimpan metrik bootstrap parametrik satu parameter pada satu dataset luar.
type OuterParamBootstrap struct {
	OuterEst             float64
	TrueValue            float64
	BootMean             float64
	BootBias             float64 // mean(boot) - outerEst
	BootCorrected        float64 // outerEst - BootBias (tanpa clamping)
	BootSD               float64
	Percentile95Lower    float64
	Percentile95Upper    float64
	FiniteCount          int
	InvalidCount         int
	IntervalValid        bool
	Percentile95Covered  bool
	CorrectedOutOfBounds bool
	CorrectedNonPhysical bool
}

// EvalParamBootstrap menghitung bias bootstrap kondisional, koreksi bias tanpa clamping,
// dan interval persentil dua sisi 95% (kuantil 0.025 & 0.975) dari seluruh hasil refit finite.
func EvalParamBootstrap(paramIdx int, outerEst, trueVal float64, rawBootDraws []float64, totalBoot int) OuterParamBootstrap {
	finite := make([]float64, 0, len(rawBootDraws))
	for _, x := range rawBootDraws {
		if !math.IsNaN(x) && !math.IsInf(x, 0) {
			finite = append(finite, x)
		}
	}
	finCnt := len(finite)
	invCnt := totalBoot - finCnt
	valid := finCnt >= MinRequiredFinite(totalBoot) && finCnt >= 2

	if finCnt == 0 {
		return OuterParamBootstrap{
			OuterEst:            outerEst,
			TrueValue:           trueVal,
			BootMean:            math.NaN(),
			BootBias:            math.NaN(),
			BootCorrected:       math.NaN(),
			BootSD:              math.NaN(),
			Percentile95Lower:   math.NaN(),
			Percentile95Upper:   math.NaN(),
			FiniteCount:         0,
			InvalidCount:        invCnt,
			IntervalValid:       false,
			Percentile95Covered: false,
		}
	}

	bMean := stats.Mean(finite)
	bBias := bMean - outerEst
	bCorr := outerEst - bBias
	bSD := stats.SD(finite)
	qLo := stats.Quantile(finite, 0.025)
	qHi := stats.Quantile(finite, 0.975)

	outBnd := bCorr < kalman.LowerBoundsOriginal[paramIdx] || bCorr > kalman.UpperBoundsOriginal[paramIdx]
	nonPhys := kalman.IsPositiveParam[paramIdx] && bCorr <= 0.0

	covered := false
	if valid {
		covered = trueVal >= qLo && trueVal <= qHi
	}

	return OuterParamBootstrap{
		OuterEst:             outerEst,
		TrueValue:            trueVal,
		BootMean:             bMean,
		BootBias:             bBias,
		BootCorrected:        bCorr,
		BootSD:               bSD,
		Percentile95Lower:    qLo,
		Percentile95Upper:    qHi,
		FiniteCount:          finCnt,
		InvalidCount:         invCnt,
		IntervalValid:        valid,
		Percentile95Covered:  covered,
		CorrectedOutOfBounds: outBnd,
		CorrectedNonPhysical: nonPhys,
	}
}
