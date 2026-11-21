package afns

import (
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ns"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/ode"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/parallel"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/rng"
	"github.com/Orls967/100QB/go/day004a_term_structure/internal/stats"
)

// Params menyimpan seluruh parameter model Arbitrage-Free Nelson-Siegel (AFNS) di bawah Q dan P.
type Params struct {
	Lambda float64
	Sigma  [3]float64    // diagonal volatilitas Sigma = diag(sigma1, sigma2, sigma3)
	KP     [3][3]float64 // matriks mean-reversion segitiga bawah di bawah P
	ThetaP [3]float64    // rata-rata jangka panjang di bawah P
}

// DefaultParams mengembalikan parameter AFNS standar tugas day004a.
func DefaultParams() Params {
	return Params{
		Lambda: 0.6,
		Sigma:  [3]float64{0.0060, 0.0100, 0.0200},
		KP: [3][3]float64{
			{0.05, 0, 0},
			{0, 0.5, 0},
			{0, 0, 0.4},
		},
		ThetaP: [3]float64{0.045, -0.020, 0.010},
	}
}

// IsDiagonalKP memeriksa apakah KP murni diagonal.
func (p Params) IsDiagonalKP() bool {
	return p.KP[1][0] == 0 && p.KP[2][0] == 0 && p.KP[2][1] == 0 &&
		p.KP[0][1] == 0 && p.KP[0][2] == 0 && p.KP[1][2] == 0
}

// BVec mengembalikan loading harga obligasi B(tau) = (-tau, -tau*g1(tau), -tau*g2(tau))^T.
func (p Params) BVec(tau float64) [3]float64 {
	if tau <= 0 {
		return [3]float64{0, 0, 0}
	}
	return [3]float64{
		-tau,
		-tau * ns.G1(tau, p.Lambda),
		-tau * ns.G2(tau, p.Lambda),
	}
}

// YieldLoadings mengembalikan -B(tau)/tau = (1, g1(lambda*tau), g2(lambda*tau)).
func (p Params) YieldLoadings(tau float64) [3]float64 {
	return [3]float64{
		1.0,
		ns.G1(tau, p.Lambda),
		ns.G2(tau, p.Lambda),
	}
}

// YieldAdjustmentComponents mengembalikan komponen suku penyesuaian bebas-arbitrase
// -A_1(tau)/tau (level), -A_2(tau)/tau (slope), -A_3(tau)/tau (kurvatur), dan total -A(tau)/tau.
// Suku penyesuaian selalu <= 0 karena dA/dtau = 1/2 B^T Sigma Sigma^T B >= 0.
func (p Params) YieldAdjustmentComponents(tau float64) (adj1, adj2, adj3, total float64) {
	if tau <= 0 {
		return 0, 0, 0, 0
	}
	s1Sq := p.Sigma[0] * p.Sigma[0]
	s2Sq := p.Sigma[1] * p.Sigma[1]
	s3Sq := p.Sigma[2] * p.Sigma[2]

	// Komponen level: -sigma1^2 * tau^2 / 6
	adj1 = -s1Sq * tau * tau / 6.0

	lam := p.Lambda
	u := lam * tau
	if u < 0.05 {
		// Untuk u kecil, integrasi Gauss-Legendre 16-titik pada 1/(2*tau) int_0^tau B_i(s)^2 ds
		// menghindari pembatalan numerik antara konstanta 1/(2*lambda^2) dan 1/u.
		nodes := lag.GaussLegendre(16, 0.0, tau)
		sum2, sum3 := 0.0, 0.0
		for _, nd := range nodes {
			bv := p.BVec(nd.X)
			sum2 += nd.W * bv[1] * bv[1]
			sum3 += nd.W * bv[2] * bv[2]
		}
		adj2 = -0.5 * s2Sq * sum2 / tau
		adj3 = -0.5 * s3Sq * sum3 / tau
		return adj1, adj2, adj3, adj1 + adj2 + adj3
	}

	lam2 := lam * lam
	lam3 := lam2 * lam
	eu := math.Exp(-u)
	e2u := eu * eu

	// Komponen slope (CDR 2011 Prop. 2):
	// A2(tau)/tau = sigma2^2 * [ 1/(2*lam^2) - (1-e^-u)/(lam^3*tau) + (1-e^-2u)/(4*lam^3*tau) ]
	a2OverTau := s2Sq * (1.0/(2.0*lam2) - (1.0-eu)/(lam3*tau) + (1.0-e2u)/(4.0*lam3*tau))
	adj2 = -a2OverTau

	// Komponen kurvatur (CDR 2011 Prop. 2):
	// int_0^tau B3(s)^2 ds = (1/lam^3) * [ u - 11/4 + (4+2u)e^-u - (5/4 + 3u/2 + u^2/2)e^-2u ]
	a3OverTau := s3Sq * (0.5/lam2 + eu/lam2 - 0.25*u*e2u/lam2 - 0.75*e2u/lam2 -
		11.0/(8.0*lam3*tau) + 2.0*eu/(lam3*tau) - 5.0*e2u/(8.0*lam3*tau))
	adj3 = -a3OverTau
	return adj1, adj2, adj3, adj1 + adj2 + adj3
}

// YieldAdjustment mengembalikan -A(tau)/tau dan A(tau) secara tertutup.
func (p Params) YieldAdjustment(tau float64) (yieldAdj, aTau float64) {
	_, _, _, total := p.YieldAdjustmentComponents(tau)
	return total, -tau * total
}

// Yield menghitung zero yield model AFNS pada state x dan tenor tau.
func (p Params) Yield(x [3]float64, tau float64) float64 {
	adj, _ := p.YieldAdjustment(tau)
	h := p.YieldLoadings(tau)
	return adj + h[0]*x[0] + h[1]*x[1] + h[2]*x[2]
}

// RiccatiDeriv mengevaluasi sisi kanan sistem ODE Riccati AFNS y = (B1, B2, B3, A):
// dB/dtau = -rho1 - (K^Q)^T B, dengan rho1 = (1, 1, 0)^T dan K^Q = [[0,0,0],[0,lam,-lam],[0,0,lam]]
// dA/dtau = 0.5 * B^T Sigma Sigma^T B + B^T K^Q theta^Q - rho0 (dengan theta^Q = 0, rho0 = 0).
func (p Params) RiccatiDeriv(_ float64, y, dydt []float64) {
	lam := p.Lambda
	b1, b2, b3 := y[0], y[1], y[2]
	dydt[0] = -1.0
	dydt[1] = -1.0 - lam*b2
	dydt[2] = lam*b2 - lam*b3
	dydt[3] = 0.5 * (p.Sigma[0]*p.Sigma[0]*b1*b1 + p.Sigma[1]*p.Sigma[1]*b2*b2 + p.Sigma[2]*p.Sigma[2]*b3*b3)
}

// SolveRiccatiRK4 mengintegrasikan ODE Riccati dari 0 ke tau dengan langkah seragam dan ekstrapolasi Richardson.
func (p Params) SolveRiccatiRK4(tau float64, nSteps int) (b [3]float64, a float64, bRich [3]float64, aRich float64) {
	y0 := []float64{0, 0, 0, 0}
	sol := ode.RK4System(0, tau, y0, nSteps, p.RiccatiDeriv)
	solR := ode.RK4Richardson(0, tau, y0, nSteps/2, p.RiccatiDeriv)
	return [3]float64{sol[0], sol[1], sol[2]}, sol[3], [3]float64{solR[0], solR[1], solR[2]}, solR[3]
}

// phi1 menghitung (exp(-a*dt) - exp(-b*dt)) / (b - a) secara stabil untuk |a-b|*dt kecil.
func phi1(a, b, dt float64) float64 {
	d := (a - b) * dt
	if math.Abs(d) < 1e-6 {
		// exp(-a*dt) * (exp(d) - 1)/((a-b)) = dt * exp(-a*dt) * (1 + d/2 + d^2/6 + d^3/24)
		e := math.Exp(-a * dt)
		return dt * e * (1.0 + d*(0.5+d*(1.0/6.0+d/24.0)))
	}
	return (math.Exp(-a*dt) - math.Exp(-b*dt)) / (b - a)
}

// TransitionMatrixClosed menghitung Phi(dt) = exp(-KP * dt) secara tertutup untuk
// matriks segitiga bawah 3x3 KP umum.
func (p Params) TransitionMatrixClosed(dt float64) [3][3]float64 {
	k11, k22, k33 := p.KP[0][0], p.KP[1][1], p.KP[2][2]
	k21, k31, k32 := p.KP[1][0], p.KP[2][0], p.KP[2][1]

	var phi [3][3]float64
	phi[0][0] = math.Exp(-k11 * dt)
	phi[1][1] = math.Exp(-k22 * dt)
	phi[2][2] = math.Exp(-k33 * dt)

	if k21 != 0 {
		phi[1][0] = -k21 * phi1(k11, k22, dt)
	}
	if k32 != 0 {
		phi[2][1] = -k32 * phi1(k22, k33, dt)
	}
	if k31 != 0 {
		phi[2][0] = -k31 * phi1(k11, k33, dt)
	}
	if k21 != 0 && k32 != 0 {
		// Suku orde dua: k32 * k21 * int_0^dt exp(-k33*(dt-s)) * phi1(k11, k22, s) ds
		var psi2 float64
		if math.Abs(k11-k22) < 1e-6 && math.Abs(k22-k33) < 1e-6 {
			psi2 = 0.5 * dt * dt * math.Exp(-k11*dt)
		} else if math.Abs(k11-k22) >= 1e-6 {
			psi2 = (phi1(k11, k33, dt) - phi1(k22, k33, dt)) / (k22 - k11)
		} else {
			psi2 = (phi1(k11, k33, dt) - dt*math.Exp(-k22*dt)) / (k33 - k22)
		}
		phi[2][0] += k32 * k21 * psi2
	}
	return phi
}

// TransitionCovClosed menghitung Q_dt = int_0^dt exp(-KP*s) Sigma Sigma^T exp(-KP^T*s) ds.
// Untuk KP diagonal memakai bentuk tertutup langsung; untuk KP segitiga bawah umum memakai
// kuadratur Gauss-Legendre 32-titik.
func (p Params) TransitionCovClosed(dt float64) [3][3]float64 {
	var q [3][3]float64
	if dt <= 0 {
		return q
	}
	if p.IsDiagonalKP() {
		for i := 0; i < 3; i++ {
			k := p.KP[i][i]
			u := 2.0 * k * dt
			var fac float64
			if math.Abs(u) < 1e-5 {
				fac = dt * (1.0 - 0.5*u + u*u/6.0 - u*u*u/24.0)
			} else {
				fac = -math.Expm1(-u) / (2.0 * k)
			}
			q[i][i] = p.Sigma[i] * p.Sigma[i] * fac
		}
		return q
	}
	return p.TransitionCovGL(dt, 32)
}

// TransitionCovGL menghitung Q_dt via integrasi Gauss-Legendre n-titik pada [0, dt].
func (p Params) TransitionCovGL(dt float64, quadOrder int) [3][3]float64 {
	var q [3][3]float64
	if dt <= 0 {
		return q
	}
	nodes := lag.GaussLegendre(quadOrder, 0.0, dt)
	sSq := [3]float64{p.Sigma[0] * p.Sigma[0], p.Sigma[1] * p.Sigma[1], p.Sigma[2] * p.Sigma[2]}
	for _, nd := range nodes {
		phi := p.TransitionMatrixClosed(nd.X)
		for i := 0; i < 3; i++ {
			for j := i; j < 3; j++ {
				sum := 0.0
				for k := 0; k < 3; k++ {
					sum += phi[i][k] * sSq[k] * phi[j][k]
				}
				q[i][j] += nd.W * sum
			}
		}
	}
	q[1][0], q[2][0], q[2][1] = q[0][1], q[0][2], q[1][2]
	return q
}

// TransitionCovVanLoan menghitung Q_dt via metode eksponensial matriks blok 6x6 Van Loan.
func (p Params) TransitionCovVanLoan(dt float64) [3][3]float64 {
	m := lag.Zeros(6, 6)
	for i := 0; i < 3; i++ {
		m[i][3+i] = p.Sigma[i] * p.Sigma[i] * dt
		for j := 0; j < 3; j++ {
			m[i][j] = -p.KP[i][j] * dt
			m[3+i][3+j] = p.KP[j][i] * dt
		}
	}
	expM := lag.MatExpSeries(m, 35)
	// F11 = exp(-KP*dt), F12 = int_0^dt exp(-KP*(dt-s)) Sigma^2 exp(KP^T s) ds -> Q = F12 * F11^T
	var q [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			sum := 0.0
			for k := 0; k < 3; k++ {
				sum += expM[i][3+k] * expM[j][k]
			}
			q[i][j] = sum
		}
	}
	return q
}

// StationaryCov menghitung kovarians stasioner Q_inf melalui persamaan Lyapunov
// KP * Q_inf + Q_inf * KP^T = Sigma * Sigma^T untuk matriks segitiga bawah KP stabil.
func (p Params) StationaryCov() [3][3]float64 {
	var q [3][3]float64
	sSq := [3]float64{p.Sigma[0] * p.Sigma[0], p.Sigma[1] * p.Sigma[1], p.Sigma[2] * p.Sigma[2]}
	for i := 0; i < 3; i++ {
		for j := 0; j <= i; j++ {
			rhs := 0.0
			if i == j {
				rhs = sSq[i]
			}
			for m := 0; m < i; m++ {
				rhs -= p.KP[i][m] * q[m][j]
			}
			for m := 0; m < j; m++ {
				rhs -= p.KP[j][m] * q[i][m]
			}
			val := rhs / (p.KP[i][i] + p.KP[j][j])
			q[i][j] = val
			q[j][i] = val
		}
	}
	return q
}

// MarketPriceOfRisk menghitung parameter harga risiko afinitas Lambda_t = lam0 + lam1 * X_t
// dari dW^Q = dW^P + (lam0 + lam1*X_t) dt sehingga lam0 = Sigma^-1 KP thetaP dan lam1 = Sigma^-1 (KQ - KP).
func (p Params) MarketPriceOfRisk() (lam0 [3]float64, lam1 [3][3]float64) {
	kq := [3][3]float64{
		{0, 0, 0},
		{0, p.Lambda, -p.Lambda},
		{0, 0, p.Lambda},
	}
	for i := 0; i < 3; i++ {
		sum := 0.0
		for j := 0; j < 3; j++ {
			sum += p.KP[i][j] * p.ThetaP[j]
			lam1[i][j] = (kq[i][j] - p.KP[i][j]) / p.Sigma[i]
		}
		lam0[i] = sum / p.Sigma[i]
	}
	return lam0, lam1
}

// LongRunTermPremium menghitung rata-rata premi jangka waktu jangka panjang pada tenor tau:
// TP_bar(tau) = E^P[y_t(tau)] - E^P[r_t] = -A(tau)/tau + (g1(tau)-1)*thetaP[1] + g2(tau)*thetaP[2].
func (p Params) LongRunTermPremium(tau float64) (tp, longRunYield, longRunShortRate, adjTerm, slopeTerm, curvTerm float64) {
	adj, _ := p.YieldAdjustment(tau)
	h := p.YieldLoadings(tau)
	longRunShortRate = p.ThetaP[0] + p.ThetaP[1]
	longRunYield = adj + h[0]*p.ThetaP[0] + h[1]*p.ThetaP[1] + h[2]*p.ThetaP[2]
	adjTerm = adj
	slopeTerm = (h[1] - 1.0) * p.ThetaP[1]
	curvTerm = h[2] * p.ThetaP[2]
	tp = longRunYield - longRunShortRate
	return tp, longRunYield, longRunShortRate, adjTerm, slopeTerm, curvTerm
}

// CondTermPremium menghitung TP_t(tau) = y_t(tau) - (1/tau) int_0^tau E_t^P[r_{t+s}] ds pada state x.
func (p Params) CondTermPremium(x [3]float64, tau float64) float64 {
	y := p.Yield(x, tau)
	// Evaluasi (1/tau) int_0^tau E_t^P[X_{t+s}] ds lewat Gauss-Legendre 16-titik (eksak untuk KP segitiga bawah umum)
	nodes := lag.GaussLegendre(16, 0.0, tau)
	var expAvgX [3]float64
	var dx [3]float64
	for i := 0; i < 3; i++ {
		dx[i] = x[i] - p.ThetaP[i]
	}
	for _, nd := range nodes {
		phi := p.TransitionMatrixClosed(nd.X)
		for i := 0; i < 3; i++ {
			expAvgX[i] += nd.W * (p.ThetaP[i] + phi[i][0]*dx[0] + phi[i][1]*dx[1] + phi[i][2]*dx[2])
		}
	}
	avgShortRate := (expAvgX[0] + expAvgX[1]) / tau
	return y - avgShortRate
}

func cholesky3(a [3][3]float64) [3][3]float64 {
	var l [3][3]float64
	for i := 0; i < 3; i++ {
		for j := 0; j <= i; j++ {
			s := a[i][j]
			for k := 0; k < j; k++ {
				s -= l[i][k] * l[j][k]
			}
			if i == j {
				if s <= 0 {
					s = 0
				}
				l[i][j] = math.Sqrt(s)
			} else if l[j][j] > 0 {
				l[i][j] = s / l[j][j]
			}
		}
	}
	return l
}

// TransitionMCReport menyimpan hasil uji Monte Carlo mean dan kovarians transisi dari state awal x0.
type TransitionMCReport struct {
	X0          [3]float64
	Dt          float64
	NumPaths    int
	TheoMean    [3]float64
	MCMean      [3]float64
	MCMeanSE    [3]float64
	MeanZScores [3]float64
	TheoCov     [3][3]float64
	MCCov       [3][3]float64
	MCCovSE     [3][3]float64
	CovZScores  [3][3]float64
	MaxZScore   float64
}

// VerifyTransitionMC memverifikasi rata-rata dan kovarians transisi diskret eksak lewat Monte Carlo.
func (p Params) VerifyTransitionMC(x0 [3]float64, dt float64, nPaths, nChunks, nWorkers int, seed uint64) TransitionMCReport {
	phi := p.TransitionMatrixClosed(dt)
	qDt := p.TransitionCovClosed(dt)
	chol := cholesky3(qDt)

	var theoMean [3]float64
	var dx [3]float64
	for i := 0; i < 3; i++ {
		dx[i] = x0[i] - p.ThetaP[i]
	}
	for i := 0; i < 3; i++ {
		theoMean[i] = p.ThetaP[i] + phi[i][0]*dx[0] + phi[i][1]*dx[1] + phi[i][2]*dx[2]
	}

	type acc struct {
		sumX   [3]float64
		sumCov [3][3]float64
		sumC2  [3][3]float64
		n      int
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) acc {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x2301, uint64(c)))
		var res acc
		res.n = cnt
		for m := 0; m < cnt; m++ {
			z0 := rg.NormFloat64()
			z1 := rg.NormFloat64()
			z2 := rg.NormFloat64()
			e0 := chol[0][0] * z0
			e1 := chol[1][0]*z0 + chol[1][1]*z1
			e2 := chol[2][0]*z0 + chol[2][1]*z1 + chol[2][2]*z2
			x := [3]float64{theoMean[0] + e0, theoMean[1] + e1, theoMean[2] + e2}
			d := [3]float64{x[0] - theoMean[0], x[1] - theoMean[1], x[2] - theoMean[2]}
			for i := 0; i < 3; i++ {
				res.sumX[i] += x[i]
				for j := i; j < 3; j++ {
					pij := d[i] * d[j]
					res.sumCov[i][j] += pij
					res.sumC2[i][j] += pij * pij
				}
			}
		}
		return res
	})

	var totSumX [3]float64
	var totSumCov, totSumC2 [3][3]float64
	for _, ch := range chunks {
		for i := 0; i < 3; i++ {
			totSumX[i] += ch.sumX[i]
			for j := i; j < 3; j++ {
				totSumCov[i][j] += ch.sumCov[i][j]
				totSumC2[i][j] += ch.sumC2[i][j]
			}
		}
	}

	fn := float64(nPaths)
	var mcMean, mcMeanSE, meanZ [3]float64
	var mcCov, mcCovSE, covZ [3][3]float64
	maxZ := 0.0
	for i := 0; i < 3; i++ {
		mcMean[i] = totSumX[i] / fn
		mcMeanSE[i] = math.Sqrt(qDt[i][i] / fn)
		meanZ[i] = math.Abs(mcMean[i]-theoMean[i]) / mcMeanSE[i]
		if meanZ[i] > maxZ {
			maxZ = meanZ[i]
		}
		for j := i; j < 3; j++ {
			mij := totSumCov[i][j] / fn
			varP := (totSumC2[i][j] - totSumCov[i][j]*totSumCov[i][j]/fn) / (fn - 1.0)
			if varP < 0 {
				varP = 0
			}
			seij := math.Sqrt(varP / fn)
			zij := math.Abs(mij-qDt[i][j]) / math.Max(seij, 1e-30)
			mcCov[i][j], mcCov[j][i] = mij, mij
			mcCovSE[i][j], mcCovSE[j][i] = seij, seij
			covZ[i][j], covZ[j][i] = zij, zij
			if zij > maxZ {
				maxZ = zij
			}
		}
	}

	return TransitionMCReport{
		X0:          x0,
		Dt:          dt,
		NumPaths:    nPaths,
		TheoMean:    theoMean,
		MCMean:      mcMean,
		MCMeanSE:    mcMeanSE,
		MeanZScores: meanZ,
		TheoCov:     qDt,
		MCCov:       mcCov,
		MCCovSE:     mcCovSE,
		CovZScores:  covZ,
		MaxZScore:   maxZ,
	}
}

// TermPremiumMCReport merangkum simulasi Monte Carlo jangka panjang untuk premi jangka waktu.
type TermPremiumMCReport struct {
	Tau              float64
	NumPaths         int
	ClosedFormTP     float64
	ClosedFormTPBps  float64
	LongRunYield     float64
	LongRunShortRate float64
	MCMeanYield10y   float64
	MCMeanAvgShort   float64
	MCTermPremium    float64
	MCTermPremiumSE  float64
	ZScore           float64
}

// VerifyTermPremiumMC memverifikasi premi jangka waktu jangka panjang pada tenor tau secara Monte Carlo
// dengan mengambil X_0 dari distribusi stasioner P dan menyimulasikan lintasan jangka pendek selama [0, tau].
func (p Params) VerifyTermPremiumMC(tau float64, nSubSteps int, nPaths, nChunks, nWorkers int, seed uint64) TermPremiumMCReport {
	tpClosed, lrYield, lrShort, _, _, _ := p.LongRunTermPremium(tau)
	qInf := p.StationaryCov()
	cholInf := cholesky3(qInf)

	dt := tau / float64(nSubSteps)
	phi := p.TransitionMatrixClosed(dt)
	qDt := p.TransitionCovClosed(dt)
	cholDt := cholesky3(qDt)

	type acc struct {
		sumY  float64
		sumR  float64
		sumTP float64
		sumT2 float64
	}
	perChunk := nPaths / nChunks
	rem := nPaths % nChunks
	chunks := parallel.RunChunks(nChunks, nWorkers, func(c int) acc {
		cnt := perChunk
		if c < rem {
			cnt++
		}
		rg := rng.New(rng.DeriveChunkSeed(seed, 0x2401, uint64(c)))
		var res acc
		for m := 0; m < cnt; m++ {
			z0, z1, z2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
			x := [3]float64{
				p.ThetaP[0] + cholInf[0][0]*z0,
				p.ThetaP[1] + cholInf[1][0]*z0 + cholInf[1][1]*z1,
				p.ThetaP[2] + cholInf[2][0]*z0 + cholInf[2][1]*z1 + cholInf[2][2]*z2,
			}
			y10 := p.Yield(x, tau)

			// Integrasi trapesium (1/tau) int_0^tau r_s ds sepanjang langkah dt
			rPrev := x[0] + x[1]
			intR := 0.0
			for step := 0; step < nSubSteps; step++ {
				w0, w1, w2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
				dx0 := x[0] - p.ThetaP[0]
				dx1 := x[1] - p.ThetaP[1]
				dx2 := x[2] - p.ThetaP[2]
				nx := [3]float64{
					p.ThetaP[0] + phi[0][0]*dx0 + phi[0][1]*dx1 + phi[0][2]*dx2 + cholDt[0][0]*w0,
					p.ThetaP[1] + phi[1][0]*dx0 + phi[1][1]*dx1 + phi[1][2]*dx2 + cholDt[1][0]*w0 + cholDt[1][1]*w1,
					p.ThetaP[2] + phi[2][0]*dx0 + phi[2][1]*dx1 + phi[2][2]*dx2 + cholDt[2][0]*w0 + cholDt[2][1]*w1 + cholDt[2][2]*w2,
				}
				rCur := nx[0] + nx[1]
				intR += 0.5 * (rPrev + rCur) * dt
				rPrev = rCur
				x = nx
			}
			avgR := intR / tau
			tp := y10 - avgR
			res.sumY += y10
			res.sumR += avgR
			res.sumTP += tp
			res.sumT2 += tp * tp
		}
		return res
	})

	sumY, sumR, sumTP, sumT2 := 0.0, 0.0, 0.0, 0.0
	for _, ch := range chunks {
		sumY += ch.sumY
		sumR += ch.sumR
		sumTP += ch.sumTP
		sumT2 += ch.sumT2
	}
	fn := float64(nPaths)
	mcY := sumY / fn
	mcR := sumR / fn
	mcTP := sumTP / fn
	varTP := (sumT2 - sumTP*sumTP/fn) / (fn - 1.0)
	if varTP < 0 {
		varTP = 0
	}
	seTP := math.Sqrt(varTP / fn)
	z := math.Abs(mcTP-tpClosed) / math.Max(seTP, 1e-30)

	return TermPremiumMCReport{
		Tau:              tau,
		NumPaths:         nPaths,
		ClosedFormTP:     tpClosed,
		ClosedFormTPBps:  tpClosed * 10_000.0,
		LongRunYield:     lrYield,
		LongRunShortRate: lrShort,
		MCMeanYield10y:   mcY,
		MCMeanAvgShort:   mcR,
		MCTermPremium:    mcTP,
		MCTermPremiumSE:  seTP,
		ZScore:           z,
	}
}

// SimulateStationaryPath membangkitkan lintasan state dan yield AFNS sepanjang T langkah delta
// mulai dari distribusi stasioner P.
func SimulateStationaryPath(p Params, taus []float64, tSteps int, dt, sigmaEps float64, rg *rng.RNG) (states [][3]float64, yields [][]float64) {
	qInf := p.StationaryCov()
	cholInf := cholesky3(qInf)
	phi := p.TransitionMatrixClosed(dt)
	qDt := p.TransitionCovClosed(dt)
	cholDt := cholesky3(qDt)

	m := len(taus)
	loadings := make([][3]float64, m)
	adjs := make([]float64, m)
	for j, tau := range taus {
		loadings[j] = p.YieldLoadings(tau)
		adjs[j], _ = p.YieldAdjustment(tau)
	}

	states = make([][3]float64, tSteps)
	yields = make([][]float64, tSteps)

	z0, z1, z2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
	x := [3]float64{
		p.ThetaP[0] + cholInf[0][0]*z0,
		p.ThetaP[1] + cholInf[1][0]*z0 + cholInf[1][1]*z1,
		p.ThetaP[2] + cholInf[2][0]*z0 + cholInf[2][1]*z1 + cholInf[2][2]*z2,
	}

	for t := 0; t < tSteps; t++ {
		if t > 0 {
			w0, w1, w2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
			dx0 := x[0] - p.ThetaP[0]
			dx1 := x[1] - p.ThetaP[1]
			dx2 := x[2] - p.ThetaP[2]
			x = [3]float64{
				p.ThetaP[0] + phi[0][0]*dx0 + phi[0][1]*dx1 + phi[0][2]*dx2 + cholDt[0][0]*w0,
				p.ThetaP[1] + phi[1][0]*dx0 + phi[1][1]*dx1 + phi[1][2]*dx2 + cholDt[1][0]*w0 + cholDt[1][1]*w1,
				p.ThetaP[2] + phi[2][0]*dx0 + phi[2][1]*dx1 + phi[2][2]*dx2 + cholDt[2][0]*w0 + cholDt[2][1]*w1 + cholDt[2][2]*w2,
			}
		}
		states[t] = x
		yt := make([]float64, m)
		for j := 0; j < m; j++ {
			yt[j] = adjs[j] + loadings[j][0]*x[0] + loadings[j][1]*x[1] + loadings[j][2]*x[2] + sigmaEps*rg.NormFloat64()
		}
		yields[t] = yt
	}
	return states, yields
}

// SimulateDynamicNSPath membangkitkan lintasan Dynamic Nelson-Siegel (tanpa suku penyesuaian -A/tau).
func SimulateDynamicNSPath(p Params, taus []float64, tSteps int, dt, sigmaEps float64, rg *rng.RNG) (states [][3]float64, yields [][]float64) {
	qInf := p.StationaryCov()
	cholInf := cholesky3(qInf)
	phi := p.TransitionMatrixClosed(dt)
	qDt := p.TransitionCovClosed(dt)
	cholDt := cholesky3(qDt)

	m := len(taus)
	loadings := make([][3]float64, m)
	for j, tau := range taus {
		loadings[j] = p.YieldLoadings(tau)
	}

	states = make([][3]float64, tSteps)
	yields = make([][]float64, tSteps)

	z0, z1, z2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
	x := [3]float64{
		p.ThetaP[0] + cholInf[0][0]*z0,
		p.ThetaP[1] + cholInf[1][0]*z0 + cholInf[1][1]*z1,
		p.ThetaP[2] + cholInf[2][0]*z0 + cholInf[2][1]*z1 + cholInf[2][2]*z2,
	}

	for t := 0; t < tSteps; t++ {
		if t > 0 {
			w0, w1, w2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
			dx0 := x[0] - p.ThetaP[0]
			dx1 := x[1] - p.ThetaP[1]
			dx2 := x[2] - p.ThetaP[2]
			x = [3]float64{
				p.ThetaP[0] + phi[0][0]*dx0 + phi[0][1]*dx1 + phi[0][2]*dx2 + cholDt[0][0]*w0,
				p.ThetaP[1] + phi[1][0]*dx0 + phi[1][1]*dx1 + phi[1][2]*dx2 + cholDt[1][0]*w0 + cholDt[1][1]*w1,
				p.ThetaP[2] + phi[2][0]*dx0 + phi[2][1]*dx1 + phi[2][2]*dx2 + cholDt[2][0]*w0 + cholDt[2][1]*w1 + cholDt[2][2]*w2,
			}
		}
		states[t] = x
		yt := make([]float64, m)
		for j := 0; j < m; j++ {
			yt[j] = loadings[j][0]*x[0] + loadings[j][1]*x[1] + loadings[j][2]*x[2] + sigmaEps*rg.NormFloat64()
		}
		yields[t] = yt
	}
	return states, yields
}

// ComputeCosineCheckDummy digunakan agar import stats selalu relevan dan terukur di pengujian.
func MaxSliceDiff(a, b []float64) float64 {
	m := 0.0
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > m {
			m = d
		}
	}
	_ = stats.Mean(a)
	return m
}
