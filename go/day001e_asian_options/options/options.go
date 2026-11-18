package options

import (
	"errors"
	"math"
	"math/rand/v2"
	"runtime"
	"sync"
)

const (
	chunkSize = 4000
	ciZ95     = 1.96
)

// Params holds option and GBM simulation parameters.
type Params struct {
	S0     float64 `json:"s0"`
	K      float64 `json:"k"`
	R      float64 `json:"r"`
	Sigma  float64 `json:"sigma"`
	T      float64 `json:"t"`
	NSteps int     `json:"n_steps"`
}

// EstimatorResult holds price, standard error, and 95% confidence interval.
type EstimatorResult struct {
	Price  float64 `json:"price"`
	SE     float64 `json:"se"`
	CILow  float64 `json:"ci_low"`
	CIHigh float64 `json:"ci_high"`
}

// ControlVariateResult extends EstimatorResult with optimal beta and exact control price.
type ControlVariateResult struct {
	Price  float64 `json:"price"`
	SE     float64 `json:"se"`
	CILow  float64 `json:"ci_low"`
	CIHigh float64 `json:"ci_high"`
	Beta   float64 `json:"beta"`
	YExact float64 `json:"y_exact"`
}

// DeltaResult holds Delta estimate, standard error, and 95% confidence interval.
type DeltaResult struct {
	Delta  float64 `json:"delta"`
	SE     float64 `json:"se"`
	CILow  float64 `json:"ci_low"`
	CIHigh float64 `json:"ci_high"`
}

// NormCDF computes the standard normal cumulative distribution function.
func NormCDF(x float64) float64 {
	return 0.5 * (1.0 + math.Erf(x/math.Sqrt2))
}

// BSCall computes the closed-form Black-Scholes European call price.
func BSCall(p Params) float64 {
	if p.T <= 0 || p.Sigma <= 0 {
		return math.Max(p.S0-p.K*math.Exp(-p.R*p.T), 0)
	}
	volSqrtT := p.Sigma * math.Sqrt(p.T)
	d1 := (math.Log(p.S0/p.K) + (p.R+0.5*p.Sigma*p.Sigma)*p.T) / volSqrtT
	d2 := d1 - volSqrtT
	return p.S0*NormCDF(d1) - p.K*math.Exp(-p.R*p.T)*NormCDF(d2)
}

// BSPut computes the closed-form Black-Scholes European put price.
func BSPut(p Params) float64 {
	if p.T <= 0 || p.Sigma <= 0 {
		return math.Max(p.K*math.Exp(-p.R*p.T)-p.S0, 0)
	}
	volSqrtT := p.Sigma * math.Sqrt(p.T)
	d1 := (math.Log(p.S0/p.K) + (p.R+0.5*p.Sigma*p.Sigma)*p.T) / volSqrtT
	d2 := d1 - volSqrtT
	return p.K*math.Exp(-p.R*p.T)*NormCDF(-d2) - p.S0*NormCDF(-d1)
}

// GeometricAsianCall computes the discrete-monitoring Kemna-Vorst closed-form call price.
func GeometricAsianCall(p Params) (float64, error) {
	if p.NSteps < 1 {
		return 0, errors.New("n_steps must be >= 1")
	}
	n := float64(p.NSteps)
	muLnG := math.Log(p.S0) + (p.R-0.5*p.Sigma*p.Sigma)*p.T*(n+1.0)/(2.0*n)
	varLnG := (p.Sigma * p.Sigma * p.T) * (n + 1.0) * (2.0*n + 1.0) / (6.0 * n * n)
	sigmaLnG := math.Sqrt(varLnG)
	expectedG := math.Exp(muLnG + 0.5*varLnG)
	d1 := (muLnG - math.Log(p.K) + varLnG) / sigmaLnG
	d2 := d1 - sigmaLnG
	return math.Exp(-p.R*p.T) * (expectedG*NormCDF(d1) - p.K*NormCDF(d2)), nil
}

func resolveWorkers(workers []int) int {
	if len(workers) > 0 && workers[0] > 0 {
		return workers[0]
	}
	return runtime.GOMAXPROCS(0)
}

func chunkRNG(seed uint64, chunkIdx int) *rand.Rand {
	s1 := seed ^ (uint64(chunkIdx+1) * 0x9E3779B97F4A7C15)
	s2 := (seed + uint64(chunkIdx)*0xBF58476D1CE4E5B9) | 1
	return rand.New(rand.NewPCG(s1, s2))
}

func simulateOnePathAvgs(rng *rand.Rand, logS0, drift, vol float64, nSteps int) (sT, arithAvg, geomAvg float64) {
	logS := logS0
	sumS := 0.0
	sumLogS := 0.0
	for range nSteps {
		z := rng.NormFloat64()
		logS += drift + vol*z
		s := math.Exp(logS)
		sumS += s
		sumLogS += logS
	}
	invN := 1.0 / float64(nSteps)
	return math.Exp(logS), sumS * invN, math.Exp(sumLogS * invN)
}

func simulateAntitheticPairArith(rng *rand.Rand, logS0, drift, vol float64, nSteps int) (avgPos, avgNeg float64) {
	logSPos := logS0
	logSNeg := logS0
	sumPos := 0.0
	sumNeg := 0.0
	for range nSteps {
		inc := vol * rng.NormFloat64()
		logSPos += drift + inc
		logSNeg += drift - inc
		sumPos += math.Exp(logSPos)
		sumNeg += math.Exp(logSNeg)
	}
	invN := 1.0 / float64(nSteps)
	return sumPos * invN, sumNeg * invN
}

func summarizeMoments(count int, sum, sumSq float64) EstimatorResult {
	n := float64(count)
	mean := sum / n
	variance := math.Max((sumSq-n*mean*mean)/(n-1.0), 0.0)
	se := math.Sqrt(variance / n)
	return EstimatorResult{
		Price:  mean,
		SE:     se,
		CILow:  mean - ciZ95*se,
		CIHigh: mean + ciZ95*se,
	}
}

func validateMC(p Params, nPaths int) error {
	if p.NSteps < 1 {
		return errors.New("n_steps must be >= 1")
	}
	if nPaths < 2 {
		return errors.New("n_paths must be >= 2")
	}
	return nil
}

type momentAcc struct {
	sum   float64
	sumSq float64
}

func parallelMoments(nItems int, workers int, fn func(chunkIdx, len int) momentAcc) (float64, float64) {
	nChunks := (nItems + chunkSize - 1) / chunkSize
	results := make([]momentAcc, nChunks)
	w := min(workers, nChunks)
	var wg sync.WaitGroup
	wg.Add(w)
	for workerID := range w {
		go func(wid int) {
			defer wg.Done()
			for c := wid; c < nChunks; c += w {
				start := c * chunkSize
				length := min(nItems-start, chunkSize)
				results[c] = fn(c, length)
			}
		}(workerID)
	}
	wg.Wait()

	var totalSum, totalSq float64
	for _, r := range results {
		totalSum += r.sum
		totalSq += r.sumSq
	}
	return totalSum, totalSq
}

// PriceEuropeanMC estimates the European call price from simulated GBM paths.
func PriceEuropeanMC(p Params, nPaths int, seed uint64, workers ...int) (EstimatorResult, error) {
	if err := validateMC(p, nPaths); err != nil {
		return EstimatorResult{}, err
	}
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPaths, resolveWorkers(workers), func(c, length int) momentAcc {
		rng := chunkRNG(seed, c)
		var acc momentAcc
		for range length {
			sT, _, _ := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
			payoff := math.Max(sT-p.K, 0) * disc
			acc.sum += payoff
			acc.sumSq += payoff * payoff
		}
		return acc
	})
	return summarizeMoments(nPaths, sum, sumSq), nil
}

// PriceAsianNaive estimates the arithmetic Asian call price with naive Monte Carlo.
func PriceAsianNaive(p Params, nPaths int, seed uint64, workers ...int) (EstimatorResult, error) {
	if err := validateMC(p, nPaths); err != nil {
		return EstimatorResult{}, err
	}
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPaths, resolveWorkers(workers), func(c, length int) momentAcc {
		rng := chunkRNG(seed, c)
		var acc momentAcc
		for range length {
			_, arithAvg, _ := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
			payoff := math.Max(arithAvg-p.K, 0) * disc
			acc.sum += payoff
			acc.sumSq += payoff * payoff
		}
		return acc
	})
	return summarizeMoments(nPaths, sum, sumSq), nil
}

// PriceAsianAntithetic estimates the arithmetic Asian call price using antithetic pairs.
func PriceAsianAntithetic(p Params, nPaths int, seed uint64, workers ...int) (EstimatorResult, error) {
	if p.NSteps < 1 {
		return EstimatorResult{}, errors.New("n_steps must be >= 1")
	}
	if nPaths < 4 || nPaths%2 != 0 {
		return EstimatorResult{}, errors.New("n_paths must be even and >= 4")
	}
	nPairs := nPaths / 2
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPairs, resolveWorkers(workers), func(c, length int) momentAcc {
		rng := chunkRNG(seed, c)
		var acc momentAcc
		for range length {
			avgPos, avgNeg := simulateAntitheticPairArith(rng, logS0, drift, vol, p.NSteps)
			pairAvg := 0.5 * (math.Max(avgPos-p.K, 0) + math.Max(avgNeg-p.K, 0)) * disc
			acc.sum += pairAvg
			acc.sumSq += pairAvg * pairAvg
		}
		return acc
	})
	return summarizeMoments(nPairs, sum, sumSq), nil
}

// PriceAsianControlVariate estimates the arithmetic Asian call using the discrete geometric Asian call as control variate.
func PriceAsianControlVariate(p Params, nPaths int, seed uint64, workers ...int) (ControlVariateResult, error) {
	if err := validateMC(p, nPaths); err != nil {
		return ControlVariateResult{}, err
	}
	yExact, err := GeometricAsianCall(p)
	if err != nil {
		return ControlVariateResult{}, err
	}
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	type cvAcc struct {
		sx, sy, sxx, syy, sxy float64
	}
	nChunks := (nPaths + chunkSize - 1) / chunkSize
	results := make([]cvAcc, nChunks)
	w := min(resolveWorkers(workers), nChunks)
	var wg sync.WaitGroup
	wg.Add(w)
	for workerID := range w {
		go func(wid int) {
			defer wg.Done()
			for c := wid; c < nChunks; c += w {
				start := c * chunkSize
				length := min(nPaths-start, chunkSize)
				rng := chunkRNG(seed, c)
				var acc cvAcc
				for range length {
					_, arithAvg, geomAvg := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
					x := math.Max(arithAvg-p.K, 0) * disc
					y := math.Max(geomAvg-p.K, 0) * disc
					acc.sx += x
					acc.sy += y
					acc.sxx += x * x
					acc.syy += y * y
					acc.sxy += x * y
				}
				results[c] = acc
			}
		}(workerID)
	}
	wg.Wait()

	var sx, sy, sxx, syy, sxy float64
	for _, r := range results {
		sx += r.sx
		sy += r.sy
		sxx += r.sxx
		syy += r.syy
		sxy += r.sxy
	}

	n := float64(nPaths)
	mx := sx / n
	my := sy / n
	covXY := (sxy - n*mx*my) / (n - 1.0)
	varY := (syy - n*my*my) / (n - 1.0)
	if varY <= 0 {
		return ControlVariateResult{}, errors.New("var_y is zero; cannot compute control variate beta")
	}
	beta := covXY / varY
	zMean := mx - beta*(my-yExact)
	varX := (sxx - n*mx*mx) / (n - 1.0)
	varZ := math.Max(varX+beta*beta*varY-2.0*beta*covXY, 0.0)
	se := math.Sqrt(varZ / n)

	return ControlVariateResult{
		Price:  zMean,
		SE:     se,
		CILow:  zMean - ciZ95*se,
		CIHigh: zMean + ciZ95*se,
		Beta:   beta,
		YExact: yExact,
	}, nil
}

// DeltaAsianPathwise estimates the arithmetic Asian call Delta using the pathwise derivative.
func DeltaAsianPathwise(p Params, nPaths int, seed uint64, workers ...int) (DeltaResult, error) {
	if err := validateMC(p, nPaths); err != nil {
		return DeltaResult{}, err
	}
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPaths, resolveWorkers(workers), func(c, length int) momentAcc {
		rng := chunkRNG(seed, c)
		var acc momentAcc
		for range length {
			_, arithAvg, _ := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
			sample := 0.0
			if arithAvg > p.K {
				sample = disc * (arithAvg / p.S0)
			}
			acc.sum += sample
			acc.sumSq += sample * sample
		}
		return acc
	})
	res := summarizeMoments(nPaths, sum, sumSq)
	return DeltaResult{
		Delta:  res.Price,
		SE:     res.SE,
		CILow:  res.CILow,
		CIHigh: res.CIHigh,
	}, nil
}

// DeltaAsianBump estimates the arithmetic Asian call Delta via central difference with CRN.
func DeltaAsianBump(p Params, nPaths int, hRel float64, seed uint64, workers ...int) (DeltaResult, error) {
	if err := validateMC(p, nPaths); err != nil {
		return DeltaResult{}, err
	}
	if hRel <= 0 {
		return DeltaResult{}, errors.New("h_rel must be > 0")
	}
	h := hRel * p.S0
	scaleUp := (p.S0 + h) / p.S0
	scaleDn := (p.S0 - h) / p.S0
	dt := p.T / float64(p.NSteps)
	drift := (p.R - 0.5*p.Sigma*p.Sigma) * dt
	vol := p.Sigma * math.Sqrt(dt)
	disc := math.Exp(-p.R * p.T)
	logS0 := math.Log(p.S0)

	sum, sumSq := parallelMoments(nPaths, resolveWorkers(workers), func(c, length int) momentAcc {
		rng := chunkRNG(seed, c)
		var acc momentAcc
		for range length {
			_, arithAvg, _ := simulateOnePathAvgs(rng, logS0, drift, vol, p.NSteps)
			payUp := math.Max(arithAvg*scaleUp-p.K, 0)
			payDn := math.Max(arithAvg*scaleDn-p.K, 0)
			sample := disc * (payUp - payDn) / (2.0 * h)
			acc.sum += sample
			acc.sumSq += sample * sample
		}
		return acc
	})
	res := summarizeMoments(nPaths, sum, sumSq)
	return DeltaResult{
		Delta:  res.Price,
		SE:     res.SE,
		CILow:  res.CILow,
		CIHigh: res.CIHigh,
	}, nil
}

// FitLogLogSlope fits the OLS slope of log10(se) vs log10(n).
func FitLogLogSlope(ns []int, ses []float64) (float64, error) {
	if len(ns) != len(ses) || len(ns) < 2 {
		return 0, errors.New("ns and ses must have equal length >= 2")
	}
	n := float64(len(ns))
	xs := make([]float64, len(ns))
	ys := make([]float64, len(ses))
	var sumX, sumY float64
	for i := range ns {
		xs[i] = math.Log10(float64(ns[i]))
		ys[i] = math.Log10(ses[i])
		sumX += xs[i]
		sumY += ys[i]
	}
	mx := sumX / n
	my := sumY / n
	var num, den float64
	for i := range ns {
		dx := xs[i] - mx
		num += dx * (ys[i] - my)
		den += dx * dx
	}
	if den <= 0 {
		return 0, errors.New("zero variance in x")
	}
	return num / den, nil
}
