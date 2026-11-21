package ns

import (
	"math"

	"github.com/Orls967/100QB/go/day004a_term_structure/internal/lag"
)

// Params menyimpan parameter kurva suku bunga Nelson-Siegel (dan perluasan Svensson).
type Params struct {
	L       float64
	S       float64
	C       float64
	Lambda  float64
	C2      float64
	Lambda2 float64
}

// DefaultNS mengembalikan parameter Nelson-Siegel standar tugas day004a.
func DefaultNS() Params {
	return Params{
		L:      0.045,
		S:      -0.025,
		C:      0.015,
		Lambda: 0.6,
	}
}

// DefaultNSS mengembalikan parameter Nelson-Siegel-Svensson dengan tambahan (C2, lambda2) = (0.01, 0.15).
func DefaultNSS() Params {
	p := DefaultNS()
	p.C2 = 0.01
	p.Lambda2 = 0.15
	return p
}

// G1 menghitung pemuatan slope g1(tau) = (1 - exp(-lambda*tau)) / (lambda*tau)
// secara stabil tanpa pembatalan numerik pada u = lambda*tau kecil.
func G1(tau, lambda float64) float64 {
	if tau <= 0 || lambda == 0 {
		return 1.0
	}
	u := lambda * tau
	if math.Abs(u) < 1e-3 {
		// Deret Taylor orde 7: 1 - u/2 + u^2/6 - u^3/24 + ...
		return 1.0 + u*(-0.5+u*(1.0/6.0+u*(-1.0/24.0+u*(1.0/120.0+u*(-1.0/720.0+u*(1.0/5040.0-u/40320.0))))))
	}
	return -math.Expm1(-u) / u
}

// G2 menghitung pemuatan kurvatur g2(tau) = g1(tau) - exp(-lambda*tau)
// secara stabil tanpa pembatalan numerik pada u = lambda*tau kecil.
func G2(tau, lambda float64) float64 {
	if tau <= 0 || lambda == 0 {
		return 0.0
	}
	u := lambda * tau
	if math.Abs(u) < 1e-3 {
		// Deret Taylor orde 7 bebas pembatalan: u/2 - u^2/3 + u^3/8 - u^4/30 + u^5/144 - u^6/840 + u^7/5760
		return u * (0.5 + u*(-1.0/3.0+u*(1.0/8.0+u*(-1.0/30.0+u*(1.0/144.0+u*(-1.0/840.0+u/5760.0))))))
	}
	return G1(tau, lambda) - math.Exp(-u)
}

// G1Naive menghitung (1 - exp(-u)) / u secara naif untuk mendemonstrasikan hilangnya digit pada tau kecil.
func G1Naive(tau, lambda float64) float64 {
	u := lambda * tau
	if u == 0 {
		return 1.0
	}
	return (1.0 - math.Exp(-u)) / u
}

// G2Naive menghitung g1_naive - exp(-u) secara naif untuk mendemonstrasikan hilangnya digit pada tau kecil.
func G2Naive(tau, lambda float64) float64 {
	u := lambda * tau
	if u == 0 {
		return 0.0
	}
	return (1.0-math.Exp(-u))/u - math.Exp(-u)
}

// Yield menghitung zero-coupon yield kontinu y(tau) menggunakan evaluasi stabil.
func (p Params) Yield(tau float64) float64 {
	if tau <= 0 {
		return p.L + p.S
	}
	y := p.L + p.S*G1(tau, p.Lambda) + p.C*G2(tau, p.Lambda)
	if p.C2 != 0 {
		y += p.C2 * G2(tau, p.Lambda2)
	}
	return y
}

// YieldNaive menghitung zero-coupon yield menggunakan rumus naif (rentan pembatalan pada tau kecil).
func (p Params) YieldNaive(tau float64) float64 {
	if tau <= 0 {
		return p.L + p.S
	}
	y := p.L + p.S*G1Naive(tau, p.Lambda) + p.C*G2Naive(tau, p.Lambda)
	if p.C2 != 0 {
		y += p.C2 * G2Naive(tau, p.Lambda2)
	}
	return y
}

// Forward menghitung suku bunga forward seketika f(tau) = d/dtau (tau * y(tau)).
func (p Params) Forward(tau float64) float64 {
	if tau <= 0 {
		return p.L + p.S
	}
	u := p.Lambda * tau
	eu := math.Exp(-u)
	f := p.L + p.S*eu + p.C*u*eu
	if p.C2 != 0 {
		u2 := p.Lambda2 * tau
		f += p.C2 * u2 * math.Exp(-u2)
	}
	return f
}

// Discount menghitung faktor diskonto P(tau) = exp(-tau * y(tau)).
func (p Params) Discount(tau float64) float64 {
	if tau <= 0 {
		return 1.0
	}
	return math.Exp(-tau * p.Yield(tau))
}

// ForwardIntegralAvg menghitung (1/tau) * integral_0^tau f(s) ds menggunakan kuadratur Gauss-Legendre.
func (p Params) ForwardIntegralAvg(tau float64, quadOrder int) float64 {
	if tau <= 0 {
		return p.L + p.S
	}
	nodes := lag.GaussLegendre(quadOrder, 0.0, tau)
	sum := 0.0
	for _, nd := range nodes {
		sum += nd.W * p.Forward(nd.X)
	}
	return sum / tau
}

// SmallTauRow menyimpan perbandingan evaluasi stabil vs naif pada tau sangat kecil.
type SmallTauRow struct {
	Tau         float64
	G1Stable    float64
	G1Naive     float64
	G1RelErr    float64
	G2Stable    float64
	G2Naive     float64
	G2RelErr    float64
	YieldStable float64
	YieldNaive  float64
	YieldAbsErr float64
}

// CompareSmallTauStability mengevaluasi pembatalan numerik g1, g2, dan y pada deretan tau kecil.
func (p Params) CompareSmallTauStability(taus []float64) []SmallTauRow {
	rows := make([]SmallTauRow, len(taus))
	for i, tau := range taus {
		g1s := G1(tau, p.Lambda)
		g1n := G1Naive(tau, p.Lambda)
		g2s := G2(tau, p.Lambda)
		g2n := G2Naive(tau, p.Lambda)
		ys := p.Yield(tau)
		yn := p.YieldNaive(tau)
		rows[i] = SmallTauRow{
			Tau:         tau,
			G1Stable:    g1s,
			G1Naive:     g1n,
			G1RelErr:    math.Abs(g1n-g1s) / math.Max(math.Abs(g1s), 1e-30),
			G2Stable:    g2s,
			G2Naive:     g2n,
			G2RelErr:    math.Abs(g2n-g2s) / math.Max(math.Abs(g2s), 1e-30),
			YieldStable: ys,
			YieldNaive:  yn,
			YieldAbsErr: math.Abs(yn - ys),
		}
	}
	return rows
}
