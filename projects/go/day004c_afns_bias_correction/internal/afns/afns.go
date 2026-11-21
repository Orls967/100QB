package afns

import (
	"math"

	"github.com/Orls967/100QB/go/day004c_afns_bias_correction/internal/rng"
)

// G1 menghitung pemuat faktor slope Nelson-Siegel g1(lam*tau) = (1 - exp(-lam*tau)) / (lam*tau)
// menggunakan ekspansi Taylor pada u kecil agar bebas pembatalan numerik.
func G1(tau, lam float64) float64 {
	u := lam * tau
	if math.Abs(u) < 1e-4 {
		return 1.0 - u*(0.5-u*(1.0/6.0-u*(1.0/24.0-u/120.0)))
	}
	return -math.Expm1(-u) / u
}

// G2 menghitung pemuat kurvatur Nelson-Siegel g2(lam*tau) = g1(u) - exp(-u)
// menggunakan ekspansi Taylor pada u kecil.
func G2(tau, lam float64) float64 {
	u := lam * tau
	if math.Abs(u) < 1e-4 {
		return u * (0.5 - u*(1.0/3.0-u*(1.0/8.0-u/30.0)))
	}
	return -math.Expm1(-u)/u - math.Exp(-u)
}

type glNode struct{ X, W float64 }

func gaussLegendre16(a, b float64) [16]glNode {
	var nodes [16]glNode
	const n = 16
	mid := 0.5 * (a + b)
	half := 0.5 * (b - a)
	for i := 0; i < n/2; i++ {
		z := math.Cos(math.Pi * (float64(i) + 0.75) / (float64(n) + 0.5))
		var pp float64
		for iter := 0; iter < 32; iter++ {
			p1, p2 := 1.0, 0.0
			for j := 1; j <= n; j++ {
				p3 := p2
				p2 = p1
				p1 = ((2.0*float64(j)-1.0)*z*p2 - float64(j-1)*p3) / float64(j)
			}
			pp = float64(n) * (z*p1 - p2) / (z*z - 1.0)
			dz := p1 / pp
			z -= dz
			if math.Abs(dz) < 1e-15 {
				break
			}
		}
		w := 2.0 / ((1.0 - z*z) * pp * pp)
		nodes[i] = glNode{X: mid - half*z, W: half * w}
		nodes[n-1-i] = glNode{X: mid + half*z, W: half * w}
	}
	return nodes
}

// Params menyimpan parameter AFNS di bawah ukuran Q dan P.
// Satuan: yield dan ThetaP dalam desimal tahunan, Lambda & KP dalam tahun^-1, Sigma dalam tahun^-1/2.
type Params struct {
	Lambda float64
	Sigma  [3]float64
	KP     [3][3]float64
	ThetaP [3]float64
}

// DefaultParams mengembalikan spesifikasi AFNS acuan yang sama dengan day004b.
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

// BVec mengembalikan loading harga obligasi B(tau) = (-tau, -tau*g1, -tau*g2)^T.
func (p Params) BVec(tau float64) [3]float64 {
	if tau <= 0 {
		return [3]float64{0, 0, 0}
	}
	return [3]float64{
		-tau,
		-tau * G1(tau, p.Lambda),
		-tau * G2(tau, p.Lambda),
	}
}

// YieldLoadings mengembalikan baris observasi -B(tau)/tau = (1, g1(lam*tau), g2(lam*tau)).
func (p Params) YieldLoadings(tau float64) [3]float64 {
	return [3]float64{1.0, G1(tau, p.Lambda), G2(tau, p.Lambda)}
}

// YieldAdjustmentComponents mengembalikan tiga suku konveksitas -A_i(tau)/tau dan totalnya.
func (p Params) YieldAdjustmentComponents(tau float64) (adj1, adj2, adj3, total float64) {
	if tau <= 0 {
		return 0, 0, 0, 0
	}
	s1Sq := p.Sigma[0] * p.Sigma[0]
	s2Sq := p.Sigma[1] * p.Sigma[1]
	s3Sq := p.Sigma[2] * p.Sigma[2]

	adj1 = -s1Sq * tau * tau / 6.0

	lam := p.Lambda
	u := lam * tau
	if u < 0.05 {
		nodes := gaussLegendre16(0.0, tau)
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

	a2OverTau := s2Sq * (1.0/(2.0*lam2) - (1.0-eu)/(lam3*tau) + (1.0-e2u)/(4.0*lam3*tau))
	adj2 = -a2OverTau

	a3OverTau := s3Sq * (0.5/lam2 + eu/lam2 - 0.25*u*e2u/lam2 - 0.75*e2u/lam2 -
		11.0/(8.0*lam3*tau) + 2.0*eu/(lam3*tau) - 5.0*e2u/(8.0*lam3*tau))
	adj3 = -a3OverTau
	return adj1, adj2, adj3, adj1 + adj2 + adj3
}

// YieldAdjustment mengembalikan (-A(tau)/tau, A(tau)) secara tertutup.
func (p Params) YieldAdjustment(tau float64) (yieldAdj, aTau float64) {
	_, _, _, total := p.YieldAdjustmentComponents(tau)
	return total, -tau * total
}

func phi1(a, b, dt float64) float64 {
	d := (a - b) * dt
	if math.Abs(d) < 1e-6 {
		e := math.Exp(-a * dt)
		return dt * e * (1.0 + d*(0.5+d*(1.0/6.0+d/24.0)))
	}
	return (math.Exp(-a*dt) - math.Exp(-b*dt)) / (b - a)
}

// TransitionMatrixClosed menghitung Phi(dt) = exp(-KP * dt).
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
	nodes := gaussLegendre16(0.0, dt)
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

// StationaryCov menghitung kovarians stasioner Q_inf lewat persamaan Lyapunov.
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

// SimulateStationaryPath mensimulasikan keadaan laten eksak dari distribusi stasioner P
// beserta observasi yield + noise iid sigma_eps.
func SimulateStationaryPath(p Params, taus []float64, tSteps int, dt, sigmaEps float64, rg *rng.RNG) (states [][3]float64, yields [][]float64) {
	states = make([][3]float64, tSteps)
	m := len(taus)
	yields = make([][]float64, tSteps)

	aVec := make([]float64, m)
	hMat := make([][3]float64, m)
	for j, tau := range taus {
		aVec[j], _ = p.YieldAdjustment(tau)
		hMat[j] = p.YieldLoadings(tau)
	}

	cholInf := cholesky3(p.StationaryCov())
	cholDt := cholesky3(p.TransitionCovClosed(dt))
	phi := p.TransitionMatrixClosed(dt)

	for t := 0; t < tSteps; t++ {
		if t == 0 {
			z0, z1, z2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
			states[0] = [3]float64{
				p.ThetaP[0] + cholInf[0][0]*z0,
				p.ThetaP[1] + cholInf[1][0]*z0 + cholInf[1][1]*z1,
				p.ThetaP[2] + cholInf[2][0]*z0 + cholInf[2][1]*z1 + cholInf[2][2]*z2,
			}
		} else {
			prev := states[t-1]
			dx := [3]float64{prev[0] - p.ThetaP[0], prev[1] - p.ThetaP[1], prev[2] - p.ThetaP[2]}
			z0, z1, z2 := rg.NormFloat64(), rg.NormFloat64(), rg.NormFloat64()
			states[t] = [3]float64{
				p.ThetaP[0] + phi[0][0]*dx[0] + phi[0][1]*dx[1] + phi[0][2]*dx[2] + cholDt[0][0]*z0,
				p.ThetaP[1] + phi[1][0]*dx[0] + phi[1][1]*dx[1] + phi[1][2]*dx[2] + cholDt[1][0]*z0 + cholDt[1][1]*z1,
				p.ThetaP[2] + phi[2][0]*dx[0] + phi[2][1]*dx[1] + phi[2][2]*dx[2] + cholDt[2][0]*z0 + cholDt[2][1]*z1 + cholDt[2][2]*z2,
			}
		}
		row := make([]float64, m)
		x := states[t]
		for j := 0; j < m; j++ {
			row[j] = aVec[j] + hMat[j][0]*x[0] + hMat[j][1]*x[1] + hMat[j][2]*x[2] + sigmaEps*rg.NormFloat64()
		}
		yields[t] = row
	}
	return states, yields
}
