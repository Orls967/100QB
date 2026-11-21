package jsonout

import (
	"math"
	"strings"
	"testing"
)

func TestMarshalPrecisionAndNaNRejection(t *testing.T) {
	obj := Obj{
		{Key: "pi", Val: math.Pi},
		{Key: "arr", Val: []float64{1.2345678901234567, -50.5}},
		{Key: "flag", Val: true},
	}
	b, err := Marshal(obj)
	if err != nil {
		t.Fatalf("unexpected marshal error: %v", err)
	}
	s := string(b)
	if !strings.Contains(s, "3.1415926535897931") {
		t.Fatalf("expected 17-digit float64 precision in %s", s)
	}

	bad := Obj{{Key: "nan", Val: math.NaN()}}
	if _, err := Marshal(bad); err == nil {
		t.Fatalf("expected error for NaN float64")
	}
}
