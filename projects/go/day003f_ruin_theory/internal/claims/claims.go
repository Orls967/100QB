package claims

import (
	"math"

	"github.com/Orls967/100QB/go/day003f_ruin_theory/internal/rng"
)

// Konstanta model risiko baku untuk seluruh sesi day003f.
const (
	Lambda  = 1.0
	Mu      = 1.0
	Theta   = 0.2
	Premium = (1.0 + Theta) * Lambda * Mu // 1.2
	Rho     = 1.0 / (1.0 + Theta)         // 5/6
)

// Kind mengidentifikasi kelima distribusi klaim pada model uji.
type Kind int

const (
	KindExponential Kind = iota
	KindErlang2
	KindHyperexp2
	KindPareto
	KindWeibull
)

// Dist mendeskripsikan satu distribusi klaim non-negatif ber-mean mu = 1.
type Dist struct {
	Kind      Kind
	Name      string
	Code      string
	HeavyTail bool

	// Parameter spesifik sesuai tipe:
	ExpRate float64 // A: 1.0; B: 2.0
	HyperP1 float64 // C: 0.5
	HyperB1 float64 // C: 2/3
	HyperB2 float64 // C: 2.0
	Alpha   float64 // D: 2.5
	Scale   float64 // D: 1.5; E: 0.5
	Shape   float64 // E: 0.5
}

// AllDists mengembalikan lima distribusi klaim terkalibrasi (mean = 1).
func AllDists() []Dist {
	return []Dist{
		Exponential(),
		Erlang2(),
		Hyperexp2(),
		Pareto(),
		Weibull(),
	}
}

// LightTailedDists mengembalikan klaim bereksponen Lundberg (A), (B), (C).
func LightTailedDists() []Dist {
	return []Dist{Exponential(), Erlang2(), Hyperexp2()}
}

// HeavyTailedDists mengembalikan klaim berekor subeksponensial (D) dan (E).
func HeavyTailedDists() []Dist {
	return []Dist{Pareto(), Weibull()}
}

// Exponential mengembalikan (A) Exp(mean = 1, rate = 1).
func Exponential() Dist {
	return Dist{Kind: KindExponential, Name: "Eksponensial(1)", Code: "A", ExpRate: 1.0}
}

// Erlang2 mengembalikan (B) Erlang(k = 2, rate = 2) dengan mean = 1.
func Erlang2() Dist {
	return Dist{Kind: KindErlang2, Name: "Erlang(2)", Code: "B", ExpRate: 2.0}
}

// Hyperexp2 mengembalikan (C) hipereksponensial 2-fase (p = 0.5, beta1 = 2/3, beta2 = 2),
// menghasilkan E[X] = 1, E[X^2] = 2.5, CV = sqrt(1.5) > 1.
func Hyperexp2() Dist {
	return Dist{
		Kind:    KindHyperexp2,
		Name:    "Hipereksponensial-2",
		Code:    "C",
		HyperP1: 0.5,
		HyperB1: 2.0 / 3.0,
		HyperB2: 2.0,
	}
}

// Pareto mengembalikan (D) Pareto Tipe II / Lomax (alpha = 2.5, scale = 1.5) dengan mean = 1.
func Pareto() Dist {
	return Dist{
		Kind:      KindPareto,
		Name:      "Pareto(alpha=2.5)",
		Code:      "D",
		HeavyTail: true,
		Alpha:     2.5,
		Scale:     1.5,
	}
}

// Weibull mengembalikan (E) Weibull (shape = 0.5, scale = 0.5) dengan mean = 1.
func Weibull() Dist {
	return Dist{
		Kind:      KindWeibull,
		Name:      "Weibull(k=0.5)",
		Code:      "E",
		HeavyTail: true,
		Shape:     0.5,
		Scale:     0.5,
	}
}

// Moments mengembalikan E[X], E[X^2], dan koefisien variasi CV.
func (d Dist) Moments() (mean, ex2, cv float64) {
	switch d.Kind {
	case KindExponential:
		return 1.0, 2.0, 1.0
	case KindErlang2:
		return 1.0, 1.5, math.Sqrt(0.5)
	case KindHyperexp2:
		p1, p2 := d.HyperP1, 1.0-d.HyperP1
		m1 := p1/d.HyperB1 + p2/d.HyperB2
		m2 := 2.0 * (p1/(d.HyperB1*d.HyperB1) + p2/(d.HyperB2*d.HyperB2))
		return m1, m2, math.Sqrt(m2-m1*m1) / m1
	case KindPareto:
		m1 := d.Scale / (d.Alpha - 1.0)
		m2 := 2.0 * d.Scale * d.Scale / ((d.Alpha - 1.0) * (d.Alpha - 2.0))
		return m1, m2, math.Sqrt(m2-m1*m1) / m1
	case KindWeibull:
		return 1.0, 6.0, math.Sqrt(5.0)
	default:
		return 0, 0, 0
	}
}

// Survival mengembalikan P(X > x) untuk x >= 0.
func (d Dist) Survival(x float64) float64 {
	if x <= 0 {
		return 1.0
	}
	switch d.Kind {
	case KindExponential:
		return math.Exp(-d.ExpRate * x)
	case KindErlang2:
		y := d.ExpRate * x
		return (1.0 + y) * math.Exp(-y)
	case KindHyperexp2:
		p2 := 1.0 - d.HyperP1
		return d.HyperP1*math.Exp(-d.HyperB1*x) + p2*math.Exp(-d.HyperB2*x)
	case KindPareto:
		return math.Pow(1.0+x/d.Scale, -d.Alpha)
	case KindWeibull:
		return math.Exp(-math.Sqrt(x / d.Scale))
	default:
		return 0
	}
}

// PDF mengembalikan f(x) untuk x > 0.
func (d Dist) PDF(x float64) float64 {
	if x < 0 {
		return 0
	}
	switch d.Kind {
	case KindExponential:
		return d.ExpRate * math.Exp(-d.ExpRate*x)
	case KindErlang2:
		return d.ExpRate * d.ExpRate * x * math.Exp(-d.ExpRate*x)
	case KindHyperexp2:
		p2 := 1.0 - d.HyperP1
		return d.HyperP1*d.HyperB1*math.Exp(-d.HyperB1*x) + p2*d.HyperB2*math.Exp(-d.HyperB2*x)
	case KindPareto:
		return (d.Alpha / d.Scale) * math.Pow(1.0+x/d.Scale, -(d.Alpha+1.0))
	case KindWeibull:
		if x == 0 {
			return math.Inf(1)
		}
		return (0.5 / math.Sqrt(d.Scale*x)) * math.Exp(-math.Sqrt(x/d.Scale))
	default:
		return 0
	}
}

// LadderSurvival mengembalikan fungsi survival distribusi tinggi tangga (ladder height)
// F_bar_L(x) = (1/mu) integral_x^inf F_bar(y) dy.
func (d Dist) LadderSurvival(x float64) float64 {
	if x <= 0 {
		return 1.0
	}
	switch d.Kind {
	case KindExponential:
		return math.Exp(-d.ExpRate * x)
	case KindErlang2:
		// E[X] = 1, integral_x^inf (1+2y)e^{-2y} dy = (1 + x)e^{-2x}
		return (1.0 + x) * math.Exp(-d.ExpRate*x)
	case KindHyperexp2:
		p2 := 1.0 - d.HyperP1
		w1 := d.HyperP1 / d.HyperB1
		w2 := p2 / d.HyperB2
		return w1*math.Exp(-d.HyperB1*x) + w2*math.Exp(-d.HyperB2*x)
	case KindPareto:
		// integral_x^inf (1 + y/1.5)^(-2.5) dy = (1 + x/1.5)^(-1.5)
		return math.Pow(1.0+x/d.Scale, -(d.Alpha - 1.0))
	case KindWeibull:
		// integral_x^inf e^{-sqrt(2y)} dy = (1 + sqrt(2x)) e^{-sqrt(2x)}
		s := math.Sqrt(x / d.Scale)
		return (1.0 + s) * math.Exp(-s)
	default:
		return 0
	}
}

// LadderCDF mengembalikan F_L(x) = 1 - F_bar_L(x).
func (d Dist) LadderCDF(x float64) float64 {
	if x <= 0 {
		return 0
	}
	return 1.0 - d.LadderSurvival(x)
}

// SubexpAsymptotic mengembalikan ekspansi asimtotik subeksponensial
// psi(u) ~ (1 / (theta * mu)) * integral_u^inf F_bar(x) dx = F_bar_L(u) / theta.
func (d Dist) SubexpAsymptotic(u float64) float64 {
	return d.LadderSurvival(u) / Theta
}

// MaxMGFDomain mengembalikan batas atas domain r di mana M_X(r) < inf (+Inf bila eksis? Tidak, untuk D & E adalah 0).
func (d Dist) MaxMGFDomain() float64 {
	switch d.Kind {
	case KindExponential:
		return d.ExpRate
	case KindErlang2:
		return d.ExpRate
	case KindHyperexp2:
		if d.HyperB1 < d.HyperB2 {
			return d.HyperB1
		}
		return d.HyperB2
	default:
		return 0
	}
}

// MGFAndDeriv mengembalikan M_X(r) dan dM_X(r)/dr untuk distribusi berekor ringan pada r < MaxMGFDomain().
func (d Dist) MGFAndDeriv(r float64) (mgf, dmgf float64, ok bool) {
	if d.HeavyTail || r >= d.MaxMGFDomain() {
		return math.Inf(1), math.Inf(1), false
	}
	switch d.Kind {
	case KindExponential:
		den := d.ExpRate - r
		mgf = d.ExpRate / den
		dmgf = d.ExpRate / (den * den)
		return mgf, dmgf, true
	case KindErlang2:
		ratio := d.ExpRate / (d.ExpRate - r)
		mgf = ratio * ratio
		dmgf = 2.0 * ratio * ratio / (d.ExpRate - r)
		return mgf, dmgf, true
	case KindHyperexp2:
		p2 := 1.0 - d.HyperP1
		d1 := d.HyperB1 - r
		d2 := d.HyperB2 - r
		mgf = d.HyperP1*(d.HyperB1/d1) + p2*(d.HyperB2/d2)
		dmgf = d.HyperP1*(d.HyperB1/(d1*d1)) + p2*(d.HyperB2/(d2*d2))
		return mgf, dmgf, true
	default:
		return math.NaN(), math.NaN(), false
	}
}

// Sample menghasilkan satu klaim X ~ F memakai RNG buatan sendiri.
func (d Dist) Sample(rg *rng.RNG) float64 {
	switch d.Kind {
	case KindExponential:
		return rg.ExpRate(d.ExpRate)
	case KindErlang2:
		return rg.Erlang(2, d.ExpRate)
	case KindHyperexp2:
		if rg.Float64() < d.HyperP1 {
			return rg.ExpRate(d.HyperB1)
		}
		return rg.ExpRate(d.HyperB2)
	case KindPareto:
		return rngPareto(rg, d.Alpha, d.Scale)
	case KindWeibull:
		return rg.Weibull(d.Shape, d.Scale)
	default:
		return 0
	}
}

func rngPareto(rg *rng.RNG, alpha, scale float64) float64 {
	return rg.ParetoLomax(alpha, scale)
}

// SampleLadder menghasilkan satu sampel tinggi tangga eksak Y ~ F_L untuk algoritma Asmussen-Kroese.
func (d Dist) SampleLadder(rg *rng.RNG) float64 {
	switch d.Kind {
	case KindExponential:
		return rg.ExpRate(d.ExpRate)
	case KindErlang2:
		// F_bar_L(x) = (1 + x) e^{-2x} = 0.5 e^{-2x} + 0.5 (1 + 2x)e^{-2x} (campuran Exp(2) & Erlang(2,2))
		if rg.Float64() < 0.5 {
			return rg.ExpRate(d.ExpRate)
		}
		return rg.Erlang(2, d.ExpRate)
	case KindHyperexp2:
		w1 := d.HyperP1 / d.HyperB1
		if rg.Float64() < w1 {
			return rg.ExpRate(d.HyperB1)
		}
		return rg.ExpRate(d.HyperB2)
	case KindPareto:
		// F_bar_L(x) = (1 + x/1.5)^(-1.5) -> Lomax dengan parameter alpha - 1 = 1.5
		return rg.ParetoLomax(d.Alpha-1.0, d.Scale)
	case KindWeibull:
		// F_bar_L(x) = P(0.5 * E_2^2 > x) dengan E_2 ~ Erlang(2, 1)
		t := rg.ExpFloat64() + rg.ExpFloat64()
		return d.Scale * t * t
	default:
		return 0
	}
}

// SampleTilted menghasilkan klaim di bawah ukuran Esscher-tilted eksponensial proporsional terhadap e^{r x} f(x).
func (d Dist) SampleTilted(rTilt float64, rg *rng.RNG) float64 {
	switch d.Kind {
	case KindExponential:
		return rg.ExpRate(d.ExpRate - rTilt)
	case KindErlang2:
		return rg.Erlang(2, d.ExpRate-rTilt)
	case KindHyperexp2:
		p2 := 1.0 - d.HyperP1
		t1 := d.HyperP1 * d.HyperB1 / (d.HyperB1 - rTilt)
		t2 := p2 * d.HyperB2 / (d.HyperB2 - rTilt)
		if rg.Float64() < t1/(t1+t2) {
			return rg.ExpRate(d.HyperB1 - rTilt)
		}
		return rg.ExpRate(d.HyperB2 - rTilt)
	default:
		return math.NaN()
	}
}
