package battery

import (
	"math"
	"testing"
)

func TestBatteryTransitionsCashflowAndBounds(t *testing.T) {
	b := DefaultSpec()
	if math.Abs(SoCMWh(b.InitialSoCIdx)-2.0) > 1e-12 {
		t.Fatalf("unexpected initial SoC: %v", SoCMWh(b.InitialSoCIdx))
	}

	// Dari SoC 2.0 MWh (idx 20), semua 5 tindakan harus feasible.
	for _, a := range Actions {
		if !b.IsFeasible(20, a) {
			t.Fatalf("action %s should be feasible at socIdx=20", a.Name)
		}
	}

	// Dari SoC 3.2 MWh (idx 32): charge 2 MW (+18 -> 50) tidak feasible, charge 1 MW (+9 -> 41) tidak feasible.
	if b.IsFeasible(32, Actions[ActionCharge2MW]) || b.IsFeasible(32, Actions[ActionCharge1MW]) {
		t.Fatalf("charge actions should violate SoC<=4.0 at socIdx=32")
	}
	// Dari SoC 3.1 MWh (idx 31): charge 1 MW (+9 -> 40) tepat feasible, charge 2 MW tidak.
	if !b.IsFeasible(31, Actions[ActionCharge1MW]) || b.IsFeasible(31, Actions[ActionCharge2MW]) {
		t.Fatalf("boundary feasibility wrong at socIdx=31")
	}
	// Dari SoC 1.0 MWh (idx 10): discharge 0.9 MW (-10 -> 0) feasible, discharge 1.8 MW (-20 -> -10) tidak.
	if !b.IsFeasible(10, Actions[ActionDischarge09MW]) || b.IsFeasible(10, Actions[ActionDischarge18MW]) {
		t.Fatalf("boundary feasibility wrong at socIdx=10")
	}
	// Dari SoC 0.0 MWh (idx 0): semua discharge tidak feasible.
	if b.IsFeasible(0, Actions[ActionDischarge09MW]) || b.IsFeasible(0, Actions[ActionDischarge18MW]) {
		t.Fatalf("discharge should violate SoC>=0 at socIdx=0")
	}

	// Cek cashflow harga negatif & positif.
	// Harga = -50, charge 2 MW -> -(-50)*2 - 2*2 = +100 - 4 = +96 EUR.
	if cf := b.HourlyCashflow(-50.0, Actions[ActionCharge2MW]); math.Abs(cf-96.0) > 1e-12 {
		t.Fatalf("negative price charge cashflow mismatch: %.6f", cf)
	}
	// Harga = 250, discharge 1.8 MW -> 250*1.8 - 2*1.8 = 450 - 3.6 = 446.4 EUR.
	if cf := b.HourlyCashflow(250.0, Actions[ActionDischarge18MW]); math.Abs(cf-446.4) > 1e-12 {
		t.Fatalf("high price discharge cashflow mismatch: %.6f", cf)
	}
	if cf := b.HourlyCashflow(100.0, Actions[ActionIdle]); cf != 0.0 {
		t.Fatalf("idle cashflow must be 0, got %v", cf)
	}

	// Cek nilai terminal dan EFC.
	if tv := b.TerminalValue(20); math.Abs(tv-(45.0*0.90*2.0)) > 1e-12 {
		t.Fatalf("terminal value mismatch: %.6f", tv)
	}
	if efc := b.EquivalentFullCycles(8.0); math.Abs(efc-1.0) > 1e-12 {
		t.Fatalf("EFC for 8 MWh stored delta on 4 MWh battery must be 1.0, got %.6f", efc)
	}
}
