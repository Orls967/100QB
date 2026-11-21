package prices

import (
	"math"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/jumps"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/parallel"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

const HoursPerDay = 24

// ModelSpec menyimpan seluruh parameter model harga listrik sintetis 3-rezim.
type ModelSpec struct {
	TransitionMatrix [regime.NumRegimes][regime.NumRegimes]float64
	Mu               [regime.NumRegimes]float64
	Phi              [regime.NumRegimes]float64
	Sigma            [regime.NumRegimes]float64
	Jumps            [regime.NumRegimes]jumps.Spec
	SeasonalAmp      float64
	SeasonalPeakHour float64
	BurnInDays       int
}

// DefaultSpec mengembalikan spesifikasi model pasar listrik sintetis sesuai tugas day004e.
func DefaultSpec() ModelSpec {
	return ModelSpec{
		TransitionMatrix: regime.DefaultTransitionMatrix(),
		Mu:               [regime.NumRegimes]float64{45.0, 90.0, -25.0},
		Phi:              [regime.NumRegimes]float64{0.82, 0.60, 0.75},
		Sigma:            [regime.NumRegimes]float64{18.0, 100.0, 40.0},
		Jumps:            jumps.DefaultSpecs(),
		SeasonalAmp:      15.0,
		SeasonalPeakHour: 19.0,
		BurnInDays:       60,
	}
}

func mod24(h int) int {
	m := h % HoursPerDay
	if m < 0 {
		m += HoursPerDay
	}
	return m
}

// SeasonalMean menghitung m(r, h) = mu[r] + 15*cos(2*pi*(h%24 - 19)/24).
func (s ModelSpec) SeasonalMean(r, h int) float64 {
	hm := float64(mod24(h))
	return s.Mu[r] + s.SeasonalAmp*math.Cos(2.0*math.Pi*(hm-s.SeasonalPeakHour)/float64(HoursPerDay))
}

// Step melakukan satu transisi harga per jam dari (rNow, h, pNow) tanpa pemotongan harga.
func (s ModelSpec) Step(rNow, h int, pNow float64, rg *rng.RNG) (rNext int, pNext float64, nJumps int, jumpVal float64) {
	hm := mod24(h)
	rNext = regime.SampleNext(s.TransitionMatrix, rNow, rg)
	nJumps, jumpVal = jumps.SampleJump(s.Jumps, rNext, rg)
	z := rg.NormFloat64()
	mNext := s.SeasonalMean(rNext, hm+1)
	mNow := s.SeasonalMean(rNext, hm)
	pNext = mNext + s.Phi[rNext]*(pNow-mNow) + s.Sigma[rNext]*z + jumpVal
	return rNext, pNext, nJumps, jumpVal
}

// ConditionalMomentsGivenNextRegime mengembalikan E[P_next | r_next, h, P_now] dan Var(P_next | r_next, h, P_now).
func (s ModelSpec) ConditionalMomentsGivenNextRegime(rNext, h int, pNow float64) (mean, variance float64) {
	hm := mod24(h)
	jm, jv := jumps.TheoreticalMoments(s.Jumps, rNext)
	mean = s.SeasonalMean(rNext, hm+1) + s.Phi[rNext]*(pNow-s.SeasonalMean(rNext, hm)) + jm
	variance = s.Sigma[rNext]*s.Sigma[rNext] + jv
	return mean, variance
}

// OneStepConditionalMoments mengembalikan E[P_next | r_now, h, P_now] dan Var(P_next | r_now, h, P_now)
// menggunakan Hukum Ekspektasi Total dan Hukum Varians Total atas r_next ~ TransitionMatrix[r_now, :].
func (s ModelSpec) OneStepConditionalMoments(rNow, h int, pNow float64) (mean, variance float64) {
	var condM, condV [regime.NumRegimes]float64
	mean = 0.0
	for rNext := 0; rNext < regime.NumRegimes; rNext++ {
		p := s.TransitionMatrix[rNow][rNext]
		cm, cv := s.ConditionalMomentsGivenNextRegime(rNext, h, pNow)
		condM[rNext] = cm
		condV[rNext] = cv
		mean += p * cm
	}
	var eVar, varE float64
	for rNext := 0; rNext < regime.NumRegimes; rNext++ {
		p := s.TransitionMatrix[rNow][rNext]
		eVar += p * condV[rNext]
		d := condM[rNext] - mean
		varE += p * d * d
	}
	return mean, eVar + varE
}

// DayPath menyimpan satu lintasan evaluasi 24 jam setelah burn-in 60 hari.
type DayPath struct {
	DayIdx   int
	Prices   [HoursPerDay]float64
	Regimes  [HoursPerDay]int
	NumJumps [HoursPerDay]int
	JumpVals [HoursPerDay]float64
}

// SimulateEvalDay membangkitkan satu hari evaluasi independen: mulai dari pi stasioner,
// menjalani burn-in 60 hari (1.440 jam), lalu mencatat jam evaluasi h=0..23.
func SimulateEvalDay(dayIdx int, spec ModelSpec, pi [regime.NumRegimes]float64, baseSeed, salt uint64) DayPath {
	rg := rng.New(rng.DeriveSeed(baseSeed, salt, uint64(dayIdx)))
	r := regime.SampleStationary(pi, rg)
	p := spec.SeasonalMean(r, 0)
	nj := 0
	jv := 0.0

	burnHours := spec.BurnInDays * HoursPerDay
	for step := 0; step < burnHours; step++ {
		r, p, nj, jv = spec.Step(r, step, p, rg)
	}

	var dp DayPath
	dp.DayIdx = dayIdx
	dp.Prices[0] = p
	dp.Regimes[0] = r
	dp.NumJumps[0] = nj
	dp.JumpVals[0] = jv
	for h := 0; h < HoursPerDay-1; h++ {
		r, p, nj, jv = spec.Step(r, h, p, rg)
		dp.Prices[h+1] = p
		dp.Regimes[h+1] = r
		dp.NumJumps[h+1] = nj
		dp.JumpVals[h+1] = jv
	}
	return dp
}

// SimulateEvalDays membangkitkan N hari evaluasi independen secara paralel deterministik.
func SimulateEvalDays(nDays, nWorkers int, spec ModelSpec, pi [regime.NumRegimes]float64, baseSeed, salt uint64) []DayPath {
	return parallel.RunChunks(nDays, nWorkers, func(dayIdx int) DayPath {
		return SimulateEvalDay(dayIdx, spec, pi, baseSeed, salt)
	})
}
