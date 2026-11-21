package jumps

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day004e_power_storage/internal/regime"
	"github.com/Orls967/100QB/go/day004e_power_storage/internal/rng"
)

func TestCompoundPoissonMoments(t *testing.T) {
	specs := DefaultSpecs()

	// Cek formula eksak.
	mNorm, vNorm := TheoreticalMoments(specs, regime.Normal)
	if mNorm != 0 || vNorm != 0 {
		t.Fatalf("normal moments must be (0,0), got (%v,%v)", mNorm, vNorm)
	}
	mScar, vScar := TheoreticalMoments(specs, regime.Scarcity)
	if math.Abs(mScar-9.6) > 1e-12 || math.Abs(vScar-2304.0) > 1e-12 {
		t.Fatalf("scarcity moments mismatch: (%v,%v)", mScar, vScar)
	}
	mOver, vOver := TheoreticalMoments(specs, regime.Oversupply)
	if math.Abs(mOver-(-3.0)) > 1e-12 || math.Abs(vOver-360.0) > 1e-12 {
		t.Fatalf("oversupply moments mismatch: (%v,%v)", mOver, vOver)
	}

	// Verifikasi Monte Carlo pada N = 200.000 draw per rezim terhadap toleransi berbasis galat MC.
	const n = 200000
	for r := 0; r < regime.NumRegimes; r++ {
		rg := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x2100+uint64(r), 0))
		var sum, sum2 float64
		multiJumpCount := 0
		for i := 0; i < n; i++ {
			nj, jv := SampleJump(specs, r, rg)
			if nj > 1 {
				multiJumpCount++
			}
			sum += jv
			sum2 += jv * jv
		}
		theoM, theoV := TheoreticalMoments(specs, r)
		empM := sum / n
		empV := (sum2 - n*empM*empM) / (n - 1)
		if r == regime.Normal {
			if empM != 0 || empV != 0 {
				t.Fatalf("normal jump must always be 0, got mean=%v var=%v", empM, empV)
			}
			continue
		}
		if multiJumpCount == 0 {
			t.Fatalf("expected multi-jump events (N>1) for regime %s", regime.Names[r])
		}
		mcseMean := math.Sqrt(theoV / n)
		// E[(J - E[J])^4] jauh lebih besar daripada normal untuk compound Poisson dengan Exp marks;
		// estimasi SE varians langsung dari momen keempat sampel (delta method).
		var m4 float64
		rg2 := rng.New(rng.DeriveSeed(rng.VerifySeed, 0x2100+uint64(r), 0))
		for i := 0; i < n; i++ {
			_, jv := SampleJump(specs, r, rg2)
			d := jv - empM
			d2 := d * d
			m4 += d2 * d2
		}
		m4 /= n
		mcseVar := math.Sqrt(math.Max(0, (m4-empV*empV)/n))
		if math.Abs(empM-theoM) > 4.0*mcseMean {
			t.Fatalf("regime %s mean %.4f outside 4*MCSE (%.4f) of %.4f", regime.Names[r], empM, mcseMean, theoM)
		}
		if math.Abs(empV-theoV) > 4.0*mcseVar {
			t.Fatalf("regime %s var %.4f outside 4*MCSE (%.4f) of %.4f", regime.Names[r], empV, mcseVar, theoV)
		}
	}
}
