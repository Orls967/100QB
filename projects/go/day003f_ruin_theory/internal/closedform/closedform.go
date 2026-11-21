package closedform

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/special"
)

// Formula memuat representasi eksak psi(u) = C1*exp(-R1*u) + C2*exp(-R2*u)
// yang diturunkan langsung lewat transformasi Laplace rasional & pecahan parsial (LU).
type Formula struct {
	Dist     claims.Dist
	R1, R2   float64
	C1, C2   float64
	PsiZero  float64
	Residual float64
}

// Derive menurunkan parameter bentuk tertutup eksak psi(u) lewat akar persamaan
// penyebut transformasi Laplace Pollaczek-Khinchine dan sistem linear LU untuk residu.
func Derive(d claims.Dist) (Formula, error) {
	rho := claims.Rho
	switch d.Kind {
	case claims.KindExponential:
		// psi(u) = rho * exp(-theta / ((1+theta)*mu) * u)
		r := claims.Theta / ((1.0 + claims.Theta) * claims.Mu)
		return Formula{
			Dist:    d,
			R1:      r,
			C1:      rho,
			PsiZero: rho,
		}, nil

	case claims.KindErlang2:
		// Transformasi Laplace: psi_hat(s) = rho*(s + 3) / (s^2 + (4 - rho)*s + 4*(1 - rho))
		bQuad := 4.0 - rho
		cQuad := 4.0 * (1.0 - rho)
		disc := math.Sqrt(bQuad*bQuad - 4.0*cQuad)
		r1 := 0.5 * (bQuad - disc)
		r2 := 0.5 * (bQuad + disc)

		mat := [][]float64{
			{1.0, 1.0},
			{r2, r1},
		}
		rhs := []float64{rho, 3.0 * rho}
		sol, err := special.SolveLinearLU(mat, rhs)
		if err != nil {
			return Formula{}, err
		}
		return Formula{
			Dist:    d,
			R1:      r1,
			R2:      r2,
			C1:      sol[0],
			C2:      sol[1],
			PsiZero: sol[0] + sol[1],
		}, nil

	case claims.KindHyperexp2:
		// f_L(x) = w1*b1 e^{-b1 x} + w2*b2 e^{-b2 x} = p1 e^{-b1 x} + p2 e^{-b2 x}
		// => f_hat_L(s) = ((p1+p2) s + p1*b2 + p2*b1) / ((s+b1)(s+b2)) = (s + b1*b2) / D(s)
		// => 1 - f_hat_L(s) = s * (s + b1 + b2 - 1) / D(s)
		b1, b2 := d.HyperB1, d.HyperB2
		numZero := (b1 + b2) - 1.0

		bQuad := (b1 + b2) - rho
		cQuad := b1 * b2 * (1.0 - rho)
		disc := math.Sqrt(bQuad*bQuad - 4.0*cQuad)
		r1 := 0.5 * (bQuad - disc)
		r2 := 0.5 * (bQuad + disc)

		mat := [][]float64{
			{1.0, 1.0},
			{r2, r1},
		}
		rhs := []float64{rho, rho * numZero}
		sol, err := special.SolveLinearLU(mat, rhs)
		if err != nil {
			return Formula{}, err
		}
		return Formula{
			Dist:    d,
			R1:      r1,
			R2:      r2,
			C1:      sol[0],
			C2:      sol[1],
			PsiZero: sol[0] + sol[1],
		}, nil

	default:
		return Formula{}, fmt.Errorf("closed-form ruin formula does not exist for %s", d.Name)
	}
}

// Eval menghitung psi(u) secara eksak pada modal awal u >= 0.
func (f Formula) Eval(u float64) float64 {
	if u < 0 {
		return 1.0
	}
	v := f.C1 * math.Exp(-f.R1*u)
	if f.R2 > 0 {
		v += f.C2 * math.Exp(-f.R2*u)
	}
	return v
}
