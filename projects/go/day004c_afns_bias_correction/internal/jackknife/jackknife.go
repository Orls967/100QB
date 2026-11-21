package jackknife

import "github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/kalman"

// TargetParamIndices adalah indeks 4 parameter sasaran utama koreksi bias & inferensi:
// kp11 (1), kp22 (2), kp33 (3), theta1 (4).
var TargetParamIndices = [4]int{1, 2, 3, 4}

// CombineSplit menghitung estimator jackknife temporal dua-paruh:
// psi_jk = 2 * psi_full - 0.5 * (psi_half1 + psi_half2) tanpa pemangkasan ke batas optimizer.
func CombineSplit(full, half1, half2 float64) float64 {
	return 2.0*full - 0.5*(half1+half2)
}

// EvalResult menyimpan nilai jackknife mentah tanpa pemangkasan beserta penanda batas & fisik.
type EvalResult struct {
	EstParams          [11]float64
	OutOfBoundsByParam [11]bool // di luar kotak batas optimizer [LowerBoundsOriginal, UpperBoundsOriginal]
	NonPhysicalByParam [11]bool // <= 0 pada parameter yang wajib positif (lambda, kp_ii, sigma_i, sigma_eps)
}

// CombineVec menghitung estimator jackknife temporal untuk seluruh 11 parameter tanpa pemangkasan.
func CombineVec(full, half1, half2 [11]float64) EvalResult {
	var out EvalResult
	for i := 0; i < 11; i++ {
		val := CombineSplit(full[i], half1[i], half2[i])
		out.EstParams[i] = val
		if val < kalman.LowerBoundsOriginal[i] || val > kalman.UpperBoundsOriginal[i] {
			out.OutOfBoundsByParam[i] = true
		}
		if kalman.IsPositiveParam[i] && val <= 0.0 {
			out.NonPhysicalByParam[i] = true
		}
	}
	return out
}
