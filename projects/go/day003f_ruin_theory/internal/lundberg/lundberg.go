package lundberg

import (
	"fmt"
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/claims"
	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/special"
)

// Result merangkum eksponen Lundberg R dan faktor depan asimtotik Cramér-Lundberg C.
type Result struct {
	R          float64
	C          float64
	MGFAtR     float64
	DMGFAtR    float64
	EqResidual float64
}

// Solve mencari eksponen penyesuaian R > 0 unik pemecah lambda*(M_X(R)-1) = c*R
// menggunakan metode Brent, lalu menghitung koefisien Cramér-Lundberg C = theta*mu/(M_X'(R) - (1+theta)*mu).
func Solve(d claims.Dist) (Result, error) {
	if d.HeavyTail {
		return Result{}, fmt.Errorf("distribution %s is heavy-tailed: Lundberg exponent does not exist", d.Code)
	}
	rMax := d.MaxMGFDomain()
	fn := func(r float64) float64 {
		mgf, _, ok := d.MGFAndDeriv(r)
		if !ok {
			return math.Inf(1)
		}
		return (mgf-1.0)/r - (1.0+claims.Theta)*claims.Mu
	}

	lo := 1e-12
	hi := rMax - 1e-8
	// Kencangkan batas atas agar bebas singularitas tepat di rMax.
	for fn(hi) < 0 && hi < rMax-1e-13 {
		hi = 0.5 * (hi + rMax)
	}
	r, err := special.BrentRoot(fn, lo, hi, 1e-14)
	if err != nil {
		return Result{}, err
	}
	mgf, dmgf, ok := d.MGFAndDeriv(r)
	if !ok {
		return Result{}, fmt.Errorf("MGF undefined at solved R=%g", r)
	}
	cPref := claims.Theta * claims.Mu / (dmgf - (1.0+claims.Theta)*claims.Mu)
	res := claims.Lambda*(mgf-1.0) - claims.Premium*r
	return Result{
		R:          r,
		C:          cPref,
		MGFAtR:     mgf,
		DMGFAtR:    dmgf,
		EqResidual: math.Abs(res),
	}, nil
}

// Bound mengembalikan batas atas Lundberg e^{-R u}.
func Bound(r, u float64) float64 {
	return math.Exp(-r * u)
}
