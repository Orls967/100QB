package jsonout

import (
	"encoding/json"
	"math"
	"testing"
)

func TestMarshalPrecisionAndNaNRejection(t *testing.T) {
	obj := Obj{
		{Key: "pi", Val: math.Pi},
		{Key: "ok", Val: true},
	}
	b, err := Marshal(obj)
	if err != nil {
		t.Fatalf("Marshal error: %v", err)
	}
	var parsed map[string]any
	if err := json.Unmarshal(b, &parsed); err != nil {
		t.Fatalf("json.Unmarshal error: %v", err)
	}
	if _, err := Marshal(Obj{{Key: "bad", Val: math.NaN()}}); err == nil {
		t.Fatalf("expected error on NaN")
	}
}
