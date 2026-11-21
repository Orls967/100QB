package jsonout

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strings"
	"testing"
)

func TestMarshalPrecisionAndNaNRejection(t *testing.T) {
	obj := Obj{
		{Key: "pi", Val: math.Pi},
		{Key: "flags", Val: []bool{true, false}},
		{Key: "rows", Val: [][]float64{{1.25, 2.5}}},
	}
	b, err := Marshal(obj)
	if err != nil {
		t.Fatalf("Marshal error: %v", err)
	}
	if !strings.Contains(string(b), "3.141592653589793") {
		t.Fatalf("expected 17-digit pi precision")
	}
	var dst map[string]any
	if err := json.Unmarshal(b, &dst); err != nil {
		t.Fatalf("invalid JSON: %v", err)
	}
	if _, err := Marshal(Obj{{Key: "inf", Val: math.Inf(1)}}); err == nil {
		t.Fatalf("expected error for Inf")
	}
	if err := WriteFile(filepath.Join(t.TempDir(), "out.json"), obj); err != nil {
		t.Fatalf("WriteFile error: %v", err)
	}
}
