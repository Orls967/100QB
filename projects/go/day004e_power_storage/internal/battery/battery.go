package battery

import "math"

const (
	CapacityMWh     = 4.0
	InitialSoCMWh   = 2.0
	SoCStepMWh      = 0.1
	NumSoCStates    = 41
	InitialSoCIdx   = 20
	EtaC            = 0.90
	EtaD            = 0.90
	MaxChargeMW     = 2.0
	MaxDischargeMW  = 1.8
	DegradationCost = 2.0
	TerminalPrice   = 45.0
	NumActions      = 5
)

// ActionID mengindeks lima tindakan diskret baterai tanpa charge dan discharge bersamaan.
type ActionID int

const (
	ActionCharge2MW     ActionID = 0
	ActionCharge1MW     ActionID = 1
	ActionIdle          ActionID = 2
	ActionDischarge09MW ActionID = 3
	ActionDischarge18MW ActionID = 4
)

// Action mendefinisikan daya sisi-grid, perubahan indeks SoC eksak, dan keterangan tindakan.
type Action struct {
	ID            ActionID
	Name          string
	PowerMW       float64 // positif = charge dari grid, negatif = discharge ke grid
	DeltaSoCIdx   int     // kelipatan eksak grid 0.1 MWh
	GridThroughMW float64 // throughput sisi-grid per jam (MWh)
	StoredDelta   float64 // perubahan energi tersimpan di baterai (MWh)
}

// Actions adalah daftar tetap lima tindakan baterai sesuai spesifikasi day004e.
var Actions = [NumActions]Action{
	{
		ID:            ActionCharge2MW,
		Name:          "charge_2MW",
		PowerMW:       2.0,
		DeltaSoCIdx:   18, // +1.8 MWh
		GridThroughMW: 2.0,
		StoredDelta:   1.8,
	},
	{
		ID:            ActionCharge1MW,
		Name:          "charge_1MW",
		PowerMW:       1.0,
		DeltaSoCIdx:   9, // +0.9 MWh
		GridThroughMW: 1.0,
		StoredDelta:   0.9,
	},
	{
		ID:            ActionIdle,
		Name:          "idle",
		PowerMW:       0.0,
		DeltaSoCIdx:   0,
		GridThroughMW: 0.0,
		StoredDelta:   0.0,
	},
	{
		ID:            ActionDischarge09MW,
		Name:          "discharge_0.9MW",
		PowerMW:       -0.9,
		DeltaSoCIdx:   -10, // -1.0 MWh
		GridThroughMW: 0.9,
		StoredDelta:   -1.0,
	},
	{
		ID:            ActionDischarge18MW,
		Name:          "discharge_1.8MW",
		PowerMW:       -1.8,
		DeltaSoCIdx:   -20, // -2.0 MWh
		GridThroughMW: 1.8,
		StoredDelta:   -2.0,
	},
}

// Spec menampung parameter baterai agar dapat dipakai ulang serta diaudit.
type Spec struct {
	CapacityMWh     float64
	InitialSoCMWh   float64
	SoCStepMWh      float64
	NumSoCStates    int
	InitialSoCIdx   int
	EtaC            float64
	EtaD            float64
	MaxChargeMW     float64
	MaxDischargeMW  float64
	DegradationCost float64
	TerminalPrice   float64
}

// DefaultSpec mengembalikan spesifikasi baterai standar day004e.
func DefaultSpec() Spec {
	return Spec{
		CapacityMWh:     CapacityMWh,
		InitialSoCMWh:   InitialSoCMWh,
		SoCStepMWh:      SoCStepMWh,
		NumSoCStates:    NumSoCStates,
		InitialSoCIdx:   InitialSoCIdx,
		EtaC:            EtaC,
		EtaD:            EtaD,
		MaxChargeMW:     MaxChargeMW,
		MaxDischargeMW:  MaxDischargeMW,
		DegradationCost: DegradationCost,
		TerminalPrice:   TerminalPrice,
	}
}

// SoCMWh mengonversi indeks grid diskret [0..40] menjadi SoC dalam MWh.
func SoCMWh(socIdx int) float64 {
	return float64(socIdx) * SoCStepMWh
}

// IsFeasible memeriksa apakah tindakan sah pada socIdx tanpa pelanggaran batas SoC maupun batas daya.
func (s Spec) IsFeasible(socIdx int, a Action) bool {
	if socIdx < 0 || socIdx >= s.NumSoCStates {
		return false
	}
	if a.PowerMW > s.MaxChargeMW+1e-12 || (-a.PowerMW) > s.MaxDischargeMW+1e-12 {
		return false
	}
	nxt := socIdx + a.DeltaSoCIdx
	return nxt >= 0 && nxt < s.NumSoCStates
}

// NextSoCIdx mengembalikan indeks SoC berikutnya beserta status kelayakan.
func (s Spec) NextSoCIdx(socIdx int, a Action) (int, bool) {
	if !s.IsFeasible(socIdx, a) {
		return socIdx, false
	}
	return socIdx + a.DeltaSoCIdx, true
}

// HourlyCashflow menghitung arus kas satu jam dalam EUR:
// - charge c: -P_t*c - degradationCost*c
// - discharge d: P_t*d - degradationCost*d
// - idle: 0
func (s Spec) HourlyCashflow(price float64, a Action) float64 {
	if a.ID == ActionIdle {
		return 0.0
	}
	// PowerMW positif saat charge (-P*c), negatif saat discharge (+P*d = -P*(PowerMW)).
	return -price*a.PowerMW - s.DegradationCost*a.GridThroughMW
}

// TerminalValue menghitung nilai terminal konservatif energi tersisa: 45 * eta_d * soc (EUR).
func (s Spec) TerminalValue(socIdx int) float64 {
	return s.TerminalPrice * s.EtaD * SoCMWh(socIdx)
}

// EquivalentFullCycles menghitung EFC eksplisit: sum |Delta SoC_MWh| / (2 * E).
// Satu siklus penuh setara pengisian 0 -> E MWh lalu pengosongan E -> 0 MWh (total perubahan tersimpan 2*E).
func (s Spec) EquivalentFullCycles(totalAbsStoredDeltaMWh float64) float64 {
	if s.CapacityMWh <= 0 {
		return 0.0
	}
	return math.Abs(totalAbsStoredDeltaMWh) / (2.0 * s.CapacityMWh)
}
