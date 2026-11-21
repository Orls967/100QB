package jsonout

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// Field adalah pasangan kunci-nilai berurut pada objek JSON.
type Field struct {
	Key string
	Val any
}

// Obj merepresentasikan objek JSON terurut.
type Obj []Field

func writeIndent(sb *strings.Builder, depth int) {
	for i := 0; i < depth; i++ {
		sb.WriteString("  ")
	}
}

func formatFloat(v float64) (string, error) {
	if math.IsNaN(v) || math.IsInf(v, 0) {
		return "", fmt.Errorf("non-finite float64 in JSON output: %v", v)
	}
	if v == 0 {
		return "0", nil
	}
	return strconv.FormatFloat(v, 'g', 17, 64), nil
}

func encodeValue(sb *strings.Builder, v any, depth int) error {
	switch x := v.(type) {
	case nil:
		sb.WriteString("null")
	case bool:
		if x {
			sb.WriteString("true")
		} else {
			sb.WriteString("false")
		}
	case []bool:
		sb.WriteByte('[')
		for i, item := range x {
			if i > 0 {
				sb.WriteString(", ")
			}
			if item {
				sb.WriteString("true")
			} else {
				sb.WriteString("false")
			}
		}
		sb.WriteByte(']')
	case int:
		sb.WriteString(strconv.Itoa(x))
	case int64:
		sb.WriteString(strconv.FormatInt(x, 10))
	case uint64:
		sb.WriteString(strconv.FormatUint(x, 10))
	case float64:
		s, err := formatFloat(x)
		if err != nil {
			return err
		}
		sb.WriteString(s)
	case string:
		sb.WriteString(strconv.Quote(x))
	case []string:
		sb.WriteByte('[')
		for i, item := range x {
			if i > 0 {
				sb.WriteString(", ")
			}
			sb.WriteString(strconv.Quote(item))
		}
		sb.WriteByte(']')
	case []float64:
		sb.WriteByte('[')
		for i, item := range x {
			if i > 0 {
				sb.WriteString(", ")
			}
			s, err := formatFloat(item)
			if err != nil {
				return err
			}
			sb.WriteString(s)
		}
		sb.WriteByte(']')
	case [][]float64:
		if len(x) == 0 {
			sb.WriteString("[]")
			return nil
		}
		sb.WriteString("[\n")
		for i, row := range x {
			writeIndent(sb, depth+1)
			if err := encodeValue(sb, row, depth+1); err != nil {
				return err
			}
			if i+1 < len(x) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		writeIndent(sb, depth)
		sb.WriteByte(']')
	case []int:
		sb.WriteByte('[')
		for i, item := range x {
			if i > 0 {
				sb.WriteString(", ")
			}
			sb.WriteString(strconv.Itoa(item))
		}
		sb.WriteByte(']')
	case []Obj:
		if len(x) == 0 {
			sb.WriteString("[]")
			return nil
		}
		sb.WriteString("[\n")
		for i, item := range x {
			writeIndent(sb, depth+1)
			if err := encodeObj(sb, item, depth+1); err != nil {
				return err
			}
			if i+1 < len(x) {
				sb.WriteByte(',')
			}
			sb.WriteByte('\n')
		}
		writeIndent(sb, depth)
		sb.WriteByte(']')
	case Obj:
		return encodeObj(sb, x, depth)
	default:
		return fmt.Errorf("unsupported JSON value type %T", v)
	}
	return nil
}

func encodeObj(sb *strings.Builder, obj Obj, depth int) error {
	if len(obj) == 0 {
		sb.WriteString("{}")
		return nil
	}
	sb.WriteString("{\n")
	for i, f := range obj {
		writeIndent(sb, depth+1)
		sb.WriteString(strconv.Quote(f.Key))
		sb.WriteString(": ")
		if err := encodeValue(sb, f.Val, depth+1); err != nil {
			return err
		}
		if i+1 < len(obj) {
			sb.WriteByte(',')
		}
		sb.WriteByte('\n')
	}
	writeIndent(sb, depth)
	sb.WriteByte('}')
	return nil
}

// Marshal mengkodekan Obj ke byte JSON berpresisi float64 17-digit.
func Marshal(root Obj) ([]byte, error) {
	var sb strings.Builder
	if err := encodeObj(&sb, root, 0); err != nil {
		return nil, err
	}
	sb.WriteByte('\n')
	return []byte(sb.String()), nil
}

// WriteFile menulis objek JSON berurut ke disk dengan presisi float64 17-digit.
func WriteFile(path string, root Obj) error {
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		return err
	}
	b, err := Marshal(root)
	if err != nil {
		return err
	}
	return os.WriteFile(path, b, 0o644)
}
