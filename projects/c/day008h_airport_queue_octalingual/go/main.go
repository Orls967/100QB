package main

import (
	"bufio"
	"bytes"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

func manifestScan(dir, prefix string, w *bufio.Writer) {
	ents, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}
	sort.Slice(ents, func(i, j int) bool { return ents[i].Name() < ents[j].Name() })
	for _, e := range ents {
		nm := e.Name()
		if strings.HasPrefix(nm, ".") {
			continue
		}
		p := filepath.Join(dir, nm)
		if e.IsDir() {
			manifestScan(p, prefix+nm+"/", w)
		} else if nm != "files_manifest.csv" {
			b, err := os.ReadFile(p)
			if err != nil {
				panic(err)
			}
			lc := bytes.Count(b, []byte{'\n'})
			fmt.Fprintf(w, "%s%s,%d,%d,%016x\n", prefix, nm, lc, len(b), fnv1aBytes(b))
		}
	}
}

func main() {
	outDir := "out/go"
	_ = os.MkdirAll(outDir+"/figures", 0o755)
	cfg := loadConfig("config.kv")
	if len(os.Args) > 1 && os.Args[1] == "s01" {
		writeS01ToS13(cfg, outDir)
		return
	}

	tot := 0.0
	for st := 1; st <= 13; st++ {
		fmt.Printf("[go] tahap %d/48 (s%02d)\n", st, st)
	}
	t0 := time.Now()
	writeS01ToS13(cfg, outDir)
	tot += time.Since(t0).Seconds()

	mainDays := int64(365)
	fmt.Println("[go] tahap 14/48 (s14)")
	t0 = time.Now()
	rpm := defaultRP(cfg, mainDays)
	rpm.auditOn = 1
	sMain := newSim(cfg, rpm)
	simRun(sMain)
	tot += time.Since(t0).Seconds()

	for st := 15; st <= 41; st++ {
		fmt.Printf("[go] tahap %d/48 (s%02d)\n", st, st)
	}
	t0 = time.Now()
	writeS14ToS45(cfg, outDir, sMain)
	tot += time.Since(t0).Seconds()

	fmt.Println("[go] tahap 42/48 (s42)")
	vals := make([]float64, 5)
	for r := 0; r < 5; r++ {
		sb := newSim(cfg, defaultRP(cfg, 3))
		tb := time.Now()
		simRun(sb)
		dt := time.Since(tb).Seconds()
		if dt < 1e-9 {
			dt = 1e-9
		}
		vals[r] = float64(sb.sysSumL) / dt
	}
	sort.Float64s(vals)
	med := vals[2]
	dev := make([]float64, 5)
	for r := 0; r < 5; r++ {
		dev[r] = math.Abs(vals[r] - med)
	}
	sort.Float64s(dev)
	f, w := openBW(outDir, "s42_bench.csv")
	fmt.Fprintf(w, "slot_pax_med,slot_pax_per_sec_med,slot_pax_per_sec_mad\n1,%.2f,%.2f\n", med, dev[2])
	closeBW(f, w)

	fmt.Println("[go] tahap 43/48 (s43)\n[go] tahap 44/48 (s44)\n[go] tahap 45/48 (s45)\n[go] tahap 46/48 (s46)")
	t0 = time.Now()
	writeS46SVGs(sMain, outDir)
	tot += time.Since(t0).Seconds()

	fmt.Println("[go] tahap 47/48 (s47)\n[go] tahap 48/48 (s48)")
	f, w = openBW(outDir, "s48_summary.json")
	fmt.Fprintf(w, "{\"lang\":\"go\",\"days\":%d,\"sw_days\":3,\"toolchain\":\"go1.27.4\",\"total_sec\":%.3f}\n", mainDays, tot)
	closeBW(f, w)

	f, w = openBW(outDir, "files_manifest.csv")
	fmt.Fprintln(w, "file,lines,bytes,fnv1a")
	manifestScan(outDir, "", w)
	closeBW(f, w)
}
