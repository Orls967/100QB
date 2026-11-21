package interval

import (
	"math"

	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/kalman"
	"github.com/Orls967/100QB/go/day004d_afns_interval_calibration/internal/stats"
)

const (
	MethodWald        = "wald"
	MethodPercentile  = "percentile"
	MethodBasic       = "basic"
	MethodStudentized = "studentized"
	MethodLogBasic    = "log_basic"

	// MinValidFraction adalah ambang pra-run proporsi minimal draw bootstrap valid (B=100 -> >= 95).
	MinValidFraction = 0.95

	ReasonNone                 = "none"
	ReasonOuterNonFinite       = "outer_non_finite"
	ReasonOuterSENonPositive   = "outer_se_non_positive"
	ReasonOuterNonPositiveLog  = "outer_non_positive_for_log"
	ReasonInsufficientFinite   = "insufficient_finite_draws"
	ReasonInsufficientValidSE  = "insufficient_valid_se_draws"
	ReasonInsufficientPositive = "insufficient_positive_draws"
)

// MinRequiredValid mengembalikan jumlah minimum draw bootstrap valid untuk ambang pra-run.
func MinRequiredValid(totalBoot int) int {
	if totalBoot <= 0 {
		return 1
	}
	req := int(math.Ceil(MinValidFraction * float64(totalBoot)))
	if req < 2 {
		req = 2
	}
	return req
}

// IntervalResult menyimpan satu selang kepercayaan untuk satu metode, satu parameter, dan satu replikasi luar.
type IntervalResult struct {
	Method          string
	Lower           float64
	Upper           float64
	Width           float64
	Valid           bool
	InvalidReason   string
	Covered         bool
	NonPhysical     bool // batas bawah <= 0 untuk parameter wajib positif
	OutOfBounds     bool // melanggar [LowerBoundsOriginal, UpperBoundsOriginal]
	ValidDrawsCount int
}

func isFinite(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

func finalizeInterval(method string, paramIdx int, lo, hi, trueVal float64, validDraws int) IntervalResult {
	if !isFinite(lo) || !isFinite(hi) {
		return IntervalResult{
			Method:          method,
			Lower:           math.NaN(),
			Upper:           math.NaN(),
			Width:           math.NaN(),
			Valid:           false,
			InvalidReason:   ReasonOuterNonFinite,
			Covered:         false,
			ValidDrawsCount: validDraws,
		}
	}
	nonPhys := kalman.IsPositiveParam[paramIdx] && lo <= 0.0
	oob := lo < kalman.LowerBoundsOriginal[paramIdx] || hi > kalman.UpperBoundsOriginal[paramIdx]
	return IntervalResult{
		Method:          method,
		Lower:           lo,
		Upper:           hi,
		Width:           hi - lo,
		Valid:           true,
		InvalidReason:   ReasonNone,
		Covered:         trueVal >= lo && trueVal <= hi,
		NonPhysical:     nonPhys,
		OutOfBounds:     oob,
		ValidDrawsCount: validDraws,
	}
}

func invalidInterval(method, reason string, validDraws int) IntervalResult {
	return IntervalResult{
		Method:          method,
		Lower:           math.NaN(),
		Upper:           math.NaN(),
		Width:           math.NaN(),
		Valid:           false,
		InvalidReason:   reason,
		Covered:         false,
		ValidDrawsCount: validDraws,
	}
}

// EvalWald membangun selang Wald 95% berbasis SE Hessian luar.
func EvalWald(paramIdx int, psiHat, seHat, trueVal float64, outerFinite bool) IntervalResult {
	if !outerFinite || !isFinite(psiHat) {
		return invalidInterval(MethodWald, ReasonOuterNonFinite, 0)
	}
	if !isFinite(seHat) || seHat <= 0 {
		return invalidInterval(MethodWald, ReasonOuterSENonPositive, 0)
	}
	lo := psiHat - stats.Z95*seHat
	hi := psiHat + stats.Z95*seHat
	return finalizeInterval(MethodWald, paramIdx, lo, hi, trueVal, 1)
}

// EvalPercentile membangun selang persentil bootstrap 95% pada skala asli tanpa pemangkasan.
func EvalPercentile(paramIdx int, trueVal float64, outerFinite bool, bootEsts []float64, totalBoot int) IntervalResult {
	if !outerFinite {
		return invalidInterval(MethodPercentile, ReasonOuterNonFinite, 0)
	}
	finite := make([]float64, 0, len(bootEsts))
	for _, x := range bootEsts {
		if isFinite(x) {
			finite = append(finite, x)
		}
	}
	if len(finite) < MinRequiredValid(totalBoot) {
		return invalidInterval(MethodPercentile, ReasonInsufficientFinite, len(finite))
	}
	lo := stats.Quantile(finite, 0.025)
	hi := stats.Quantile(finite, 0.975)
	return finalizeInterval(MethodPercentile, paramIdx, lo, hi, trueVal, len(finite))
}

// EvalBasic membangun selang dasar / reversed percentile [2*psiHat - q_0.975, 2*psiHat - q_0.025] tanpa pemangkasan.
func EvalBasic(paramIdx int, psiHat, trueVal float64, outerFinite bool, bootEsts []float64, totalBoot int) IntervalResult {
	if !outerFinite || !isFinite(psiHat) {
		return invalidInterval(MethodBasic, ReasonOuterNonFinite, 0)
	}
	finite := make([]float64, 0, len(bootEsts))
	for _, x := range bootEsts {
		if isFinite(x) {
			finite = append(finite, x)
		}
	}
	if len(finite) < MinRequiredValid(totalBoot) {
		return invalidInterval(MethodBasic, ReasonInsufficientFinite, len(finite))
	}
	qLo := stats.Quantile(finite, 0.025)
	qHi := stats.Quantile(finite, 0.975)
	return finalizeInterval(MethodBasic, paramIdx, 2.0*psiHat-qHi, 2.0*psiHat-qLo, trueVal, len(finite))
}

// EvalStudentized membangun selang bootstrap-t / studentized [psiHat - q_0.975(t*)*seHat, psiHat - q_0.025(t*)*seHat].
// Draw dengan psi* nonfinite atau SE* nonfinite / <= 0 dikeluarkan; bila jumlah draw valid < ambang, interval ditandai invalid.
func EvalStudentized(paramIdx int, psiHat, seHat, trueVal float64, outerFinite bool, bootEsts, bootSEs []float64, totalBoot int) IntervalResult {
	if !outerFinite || !isFinite(psiHat) {
		return invalidInterval(MethodStudentized, ReasonOuterNonFinite, 0)
	}
	if !isFinite(seHat) || seHat <= 0 {
		return invalidInterval(MethodStudentized, ReasonOuterSENonPositive, 0)
	}
	n := len(bootEsts)
	if len(bootSEs) < n {
		n = len(bootSEs)
	}
	tStars := make([]float64, 0, n)
	for b := 0; b < n; b++ {
		if isFinite(bootEsts[b]) && isFinite(bootSEs[b]) && bootSEs[b] > 0 {
			tStars = append(tStars, (bootEsts[b]-psiHat)/bootSEs[b])
		}
	}
	if len(tStars) < MinRequiredValid(totalBoot) {
		return invalidInterval(MethodStudentized, ReasonInsufficientValidSE, len(tStars))
	}
	qtLo := stats.Quantile(tStars, 0.025)
	qtHi := stats.Quantile(tStars, 0.975)
	lo := psiHat - qtHi*seHat
	hi := psiHat - qtLo*seHat
	return finalizeInterval(MethodStudentized, paramIdx, lo, hi, trueVal, len(tStars))
}

// EvalLogBasic membangun selang dasar pada skala log untuk parameter positif lalu mengeksponensialkan kembali.
func EvalLogBasic(paramIdx int, psiHat, trueVal float64, outerFinite bool, bootEsts []float64, totalBoot int) IntervalResult {
	if !outerFinite || !isFinite(psiHat) {
		return invalidInterval(MethodLogBasic, ReasonOuterNonFinite, 0)
	}
	if psiHat <= 0 {
		return invalidInterval(MethodLogBasic, ReasonOuterNonPositiveLog, 0)
	}
	logDraws := make([]float64, 0, len(bootEsts))
	for _, x := range bootEsts {
		if isFinite(x) && x > 0 {
			logDraws = append(logDraws, math.Log(x))
		}
	}
	if len(logDraws) < MinRequiredValid(totalBoot) {
		return invalidInterval(MethodLogBasic, ReasonInsufficientPositive, len(logDraws))
	}
	etaHat := math.Log(psiHat)
	qLo := stats.Quantile(logDraws, 0.025)
	qHi := stats.Quantile(logDraws, 0.975)
	lo := math.Exp(2.0*etaHat - qHi)
	hi := math.Exp(2.0*etaHat - qLo)
	return finalizeInterval(MethodLogBasic, paramIdx, lo, hi, trueVal, len(logDraws))
}
