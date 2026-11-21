package jumps

import (
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

// Spec menyimpan parameter compound Poisson untuk satu rezim tujuan.
type Spec struct {
	Lambda   float64
	MeanJump float64
	Sign     float64 // 0 untuk normal, +1 untuk scarcity, -1 untuk oversupply
}

// DefaultSpecs mengembalikan parameter jump compound Poisson per rezim tujuan.
func DefaultSpecs() [regime.NumRegimes]Spec {
	return [regime.NumRegimes]Spec{
		regime.Normal:     {Lambda: 0.002, MeanJump: 0.0, Sign: 0.0},
		regime.Scarcity:   {Lambda: 0.08, MeanJump: 120.0, Sign: 1.0},
		regime.Oversupply: {Lambda: 0.05, MeanJump: 60.0, Sign: -1.0},
	}
}

// TheoreticalMoments mengembalikan E[J] dan Var(J) bersyarat pada rezim tujuan.
func TheoreticalMoments(specs [regime.NumRegimes]Spec, targetRegime int) (mean, variance float64) {
	sp := specs[targetRegime]
	if sp.Sign == 0 || sp.MeanJump == 0 {
		return 0.0, 0.0
	}
	mean = sp.Sign * sp.Lambda * sp.MeanJump
	variance = 2.0 * sp.Lambda * sp.MeanJump * sp.MeanJump
	return mean, variance
}

// SampleJump membangkitkan cacah jump N ~ Poisson(lambda) dan nilai jump agregat J bersyarat rezim tujuan.
func SampleJump(specs [regime.NumRegimes]Spec, targetRegime int, rg *rng.RNG) (nJumps int, jumpVal float64) {
	sp := specs[targetRegime]
	nJumps = rg.Poisson(sp.Lambda)
	if nJumps == 0 || sp.Sign == 0 || sp.MeanJump <= 0 {
		return nJumps, 0.0
	}
	sumExp := 0.0
	for k := 0; k < nJumps; k++ {
		sumExp += rg.ExpFloat64(sp.MeanJump)
	}
	return nJumps, sp.Sign * sumExp
}
