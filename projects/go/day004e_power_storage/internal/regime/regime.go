package regime

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

const (
	Normal     = 0
	Scarcity   = 1
	Oversupply = 2
	NumRegimes = 3
)

// Names adalah label rezim terurut [normal, scarcity, oversupply].
var Names = [NumRegimes]string{"normal", "scarcity", "oversupply"}

// DefaultTransitionMatrix mengembalikan matriks transisi Markov baris-ke-baris per jam.
func DefaultTransitionMatrix() [NumRegimes][NumRegimes]float64 {
	return [NumRegimes][NumRegimes]float64{
		{0.985, 0.010, 0.005},
		{0.250, 0.740, 0.010},
		{0.150, 0.050, 0.800},
	}
}

// ValidateTransitionMatrix memverifikasi tiap entri di [0, 1] dan jumlah tiap baris sama dengan 1.
func ValidateTransitionMatrix(P [NumRegimes][NumRegimes]float64) error {
	for i := 0; i < NumRegimes; i++ {
		sum := 0.0
		for j := 0; j < NumRegimes; j++ {
			v := P[i][j]
			if math.IsNaN(v) || math.IsInf(v, 0) || v < 0.0 || v > 1.0 {
				return fmt.Errorf("entri matriks transisi (%d,%d) di luar [0,1]: %v", i, j, v)
			}
			sum += v
		}
		if math.Abs(sum-1.0) > 1e-12 {
			return fmt.Errorf("jumlah baris %d matriks transisi = %.17g != 1", i, sum)
		}
	}
	return nil
}

// StationaryDistribution menyelesaikan pi * P = pi dengan syarat sum(pi) = 1 secara eksak (eliminasi Gauss).
func StationaryDistribution(P [NumRegimes][NumRegimes]float64) ([NumRegimes]float64, error) {
	if err := ValidateTransitionMatrix(P); err != nil {
		return [NumRegimes]float64{}, err
	}
	// Sistem A * pi = b: dua persamaan pertama dari (P^T - I)pi = 0 dan baris terakhir sum(pi) = 1.
	var a [NumRegimes][NumRegimes + 1]float64
	for eq := 0; eq < NumRegimes-1; eq++ {
		for j := 0; j < NumRegimes; j++ {
			val := P[j][eq]
			if j == eq {
				val -= 1.0
			}
			a[eq][j] = val
		}
		a[eq][NumRegimes] = 0.0
	}
	for j := 0; j < NumRegimes; j++ {
		a[NumRegimes-1][j] = 1.0
	}
	a[NumRegimes-1][NumRegimes] = 1.0

	// Eliminasi Gauss dengan pivoting parsial.
	for col := 0; col < NumRegimes; col++ {
		piv := col
		best := math.Abs(a[col][col])
		for row := col + 1; row < NumRegimes; row++ {
			if math.Abs(a[row][col]) > best {
				best = math.Abs(a[row][col])
				piv = row
			}
		}
		if best < 1e-14 {
			return [NumRegimes]float64{}, fmt.Errorf("matriks distribusi stasioner singular pada kolom %d", col)
		}
		if piv != col {
			a[col], a[piv] = a[piv], a[col]
		}
		for row := col + 1; row < NumRegimes; row++ {
			factor := a[row][col] / a[col][col]
			for k := col; k <= NumRegimes; k++ {
				a[row][k] -= factor * a[col][k]
			}
		}
	}

	var pi [NumRegimes]float64
	for row := NumRegimes - 1; row >= 0; row-- {
		sum := a[row][NumRegimes]
		for k := row + 1; k < NumRegimes; k++ {
			sum -= a[row][k] * pi[k]
		}
		pi[row] = sum / a[row][row]
		if pi[row] < -1e-14 {
			return [NumRegimes]float64{}, fmt.Errorf("distribusi stasioner menghasilkan komponen negatif: %v", pi[row])
		}
		if pi[row] < 0 {
			pi[row] = 0
		}
	}
	return pi, nil
}

// StationaryResidual mengembalikan max_j |(pi * P)_j - pi_j|.
func StationaryResidual(P [NumRegimes][NumRegimes]float64, pi [NumRegimes]float64) float64 {
	maxErr := 0.0
	for j := 0; j < NumRegimes; j++ {
		nextProb := 0.0
		for i := 0; i < NumRegimes; i++ {
			nextProb += pi[i] * P[i][j]
		}
		diff := math.Abs(nextProb - pi[j])
		if diff > maxErr {
			maxErr = diff
		}
	}
	return maxErr
}

// SampleNext mengambil sampel rezim berikutnya bersyarat pada rezim saat ini.
func SampleNext(P [NumRegimes][NumRegimes]float64, cur int, rg *rng.RNG) int {
	u := rg.Float64()
	cum := 0.0
	for j := 0; j < NumRegimes-1; j++ {
		cum += P[cur][j]
		if u < cum {
			return j
		}
	}
	return NumRegimes - 1
}

// SampleStationary mengambil sampel rezim awal dari distribusi stasioner pi.
func SampleStationary(pi [NumRegimes]float64, rg *rng.RNG) int {
	u := rg.Float64()
	cum := 0.0
	for j := 0; j < NumRegimes-1; j++ {
		cum += pi[j]
		if u < cum {
			return j
		}
	}
	return NumRegimes - 1
}
