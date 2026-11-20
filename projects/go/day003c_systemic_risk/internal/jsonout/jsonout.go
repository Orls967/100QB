package jsonout

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// Field merepresentasikan pasangan key-value terurut pada objek JSON.
type Field struct {
	Key string
	Val any
}

// Obj adalah slice Field terurut.
type Obj []Field

// WriteFile menulis nilai ke berkas JSON berindentasi dengan presisi 17 digit signifikan untuk float64.
func WriteFile(path string, v any) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	var sb strings.Builder
	writeVal(&sb, v, 0)
	sb.WriteByte('\n')
	return os.WriteFile(path, []byte(sb.String()), 0o644)
}

func indentStr(sb *strings.Builder, n int) {
	for i := 0; i < n; i++ {
		sb.WriteByte(' ')
	}
}

func writeVal(sb *strings.Builder, v any, indent int) {
	switch val := v.(type) {
	case nil:
		sb.WriteString("null")
	case bool:
		if val {
			sb.WriteString("true")
		} else {
			sb.WriteString("false")
		}
	case int:
		sb.WriteString(strconv.Itoa(val))
	case int64:
		sb.WriteString(strconv.FormatInt(val, 10))
	case float64:
		if math.IsNaN(val) || math.IsInf(val, 0) {
			sb.WriteString("null")
		} else {
			sb.WriteString(strconv.FormatFloat(val, 'g', 17, 64))
		}
	case string:
		sb.WriteString(strconv.Quote(val))
	case []float64:
		if len(val) == 0 {
			sb.WriteString("[]")
			return
		}
		sb.WriteString("[\n")
		for i, x := range val {
			indentStr(sb, indent+2)
			writeVal(sb, x, indent+2)
			if i+1 < len(val) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		indentStr(sb, indent)
		sb.WriteByte(']')
	case []int:
		if len(val) == 0 {
			sb.WriteString("[]")
			return
		}
		sb.WriteString("[\n")
		for i, x := range val {
			indentStr(sb, indent+2)
			writeVal(sb, x, indent+2)
			if i+1 < len(val) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		indentStr(sb, indent)
		sb.WriteByte(']')
	case []any:
		if len(val) == 0 {
			sb.WriteString("[]")
			return
		}
		sb.WriteString("[\n")
		for i, x := range val {
			indentStr(sb, indent+2)
			writeVal(sb, x, indent+2)
			if i+1 < len(val) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		indentStr(sb, indent)
		sb.WriteByte(']')
	case []Obj:
		if len(val) == 0 {
			sb.WriteString("[]")
			return
		}
		sb.WriteString("[\n")
		for i, x := range val {
			indentStr(sb, indent+2)
			writeVal(sb, x, indent+2)
			if i+1 < len(val) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		indentStr(sb, indent)
		sb.WriteByte(']')
	case Obj:
		if len(val) == 0 {
			sb.WriteString("{}")
			return
		}
		sb.WriteString("{\n")
		for i, f := range val {
			indentStr(sb, indent+2)
			sb.WriteString(strconv.Quote(f.Key))
			sb.WriteString(": ")
			writeVal(sb, f.Val, indent+2)
			if i+1 < len(val) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		indentStr(sb, indent)
		sb.WriteByte('}')
	default:
		sb.WriteString(strconv.Quote(fmt.Sprintf("%v", val)))
	}
}
