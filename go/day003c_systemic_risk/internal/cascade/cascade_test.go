package cascade_test

import (
	"math"
	"testing"

	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/cascade"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/graphs"
	"github.com/Orls967/100QB/go/day003c_systemic_risk/internal/rng"
)

func TestCascadeOrderInvariance(t *testing.T) {
	r := rng.New(20261121)
	n := 10000
	g := graphs.BuildUndirectedER(n, 4.0, r)
	perm := r.Perm(n)
	seeds := perm[:20]

	cntFIFO, maskFIFO := cascade.Simulate(g, seeds, 0.18, false)
	cntLIFO, maskLIFO := cascade.Simulate(g, seeds, 0.18, true)

	if cntFIFO != cntLIFO {
		t.Fatalf("Order invariance gagal: FIFO=%d vs LIFO=%d", cntFIFO, cntLIFO)
	}
	for i := 0; i < n; i++ {
		if maskFIFO[i] != maskLIFO[i] {
			t.Fatalf("Simpul %d berbeda antara FIFO dan LIFO", i)
		}
	}
}

func TestPhaseBoundariesAndMeanField(t *testing.T) {
	c1, c2, ok := cascade.PhaseBoundaries(0.18)
	if !ok || !(c1 > 0.5 && c1 < 2.0 && c2 > 5.0 && c2 < 10.0) {
		t.Fatalf("Batas fase phi=0.18 tidak wajar: c1=%f c2=%f ok=%v", c1, c2, ok)
	}
	if math.Abs(cascade.VulnerableBranching(c1, 0.18)) > 1e-9 || math.Abs(cascade.VulnerableBranching(c2, 0.18)) > 1e-9 {
		t.Fatalf("Akar bisection tidak nol")
	}

	// Di dalam jendela (c=4.0), S_teori besar; di luar jendela (c=10.0), S_teori dekat rho0
	_, sIn := cascade.MeanField(4.0, 0.001, 0.18)
	_, sOut := cascade.MeanField(10.0, 0.001, 0.18)
	if sIn < 0.5 {
		t.Fatalf("S_teori di c=4.0 terlalu kecil: %f", sIn)
	}
	if sOut > 0.01 {
		t.Fatalf("S_teori di c=10.0 terlalu besar: %f", sOut)
	}
}
