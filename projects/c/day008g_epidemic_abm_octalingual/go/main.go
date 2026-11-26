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
		writeS01ToS09(cfg, outDir)
		return
	}

	tot := 0.0
	for st := 1; st <= 9; st++ {
		fmt.Printf("[go] tahap %d/39 (s%02d)\n", st, st)
	}
	t0 := time.Now()
	writeS01ToS09(cfg, outDir)
	tot += time.Since(t0).Seconds()

	nMain, betaCal := int64(100000), int64(400)
	fmt.Println("[go] tahap 10/39 (s10)")
	t0 = time.Now()
	sNo := newSim(cfg, defaultRunParams(nMain, cfg.mainDays, betaCal, 0), true)
	simRun(sNo)
	tot += time.Since(t0).Seconds()

	fmt.Println("[go] tahap 11/39 (s11)")
	t0 = time.Now()
	sMain := newSim(cfg, defaultRunParams(nMain, cfg.mainDays, betaCal, 1), true)
	simRun(sMain)
	tot += time.Since(t0).Seconds()

	for st := 12; st <= 34; st++ {
		fmt.Printf("[go] tahap %d/39 (s%02d)\n", st, st)
	}
	t0 = time.Now()
	writeS10ToS37(sNo, sMain, cfg, outDir, betaCal)
	tot += time.Since(t0).Seconds()

	fmt.Println("[go] tahap 35/39 (s35)")
	vals := make([]float64, 5)
	for i := 0; i < 5; i++ {
		sb := newSim(cfg, defaultRunParams(2000, 60, betaCal, 1), false)
		tb := time.Now()
		simRun(sb)
		dt := time.Since(tb).Seconds()
		if dt < 1e-9 {
			dt = 1e-9
		}
		vals[i] = 120000.0 / dt
	}
	sort.Float64s(vals)
	med := vals[2]
	dev := make([]float64, 5)
	for i := 0; i < 5; i++ {
		dev[i] = math.Abs(vals[i] - med)
	}
	sort.Float64s(dev)
	f, w := openBW(outDir, "s35_bench.csv")
	fmt.Fprintf(w, "person_days,pd_per_sec_med,pd_per_sec_mad\n120000,%.2f,%.2f\n", med, dev[2])
	closeBW(f, w)

	fmt.Println("[go] tahap 36/39 (s36)\n[go] tahap 37/39 (s37)\n[go] tahap 38/39 (s38)")
	t0 = time.Now()
	writeS38SVGs(sNo, sMain, outDir)
	tot += time.Since(t0).Seconds()

	fmt.Println("[go] tahap 39/39 (s39)")
	f, w = openBW(outDir, "s39_summary.json")
	fmt.Fprintf(w, "{\"lang\":\"go\",\"N\":%d,\"days\":%d,\"N_sw\":1000,\"toolchain\":\"go1.27.4\",\"total_sec\":%.3f}\n", nMain, cfg.mainDays, tot)
	closeBW(f, w)

	f, w = openBW(outDir, "files_manifest.csv")
	fmt.Fprintln(w, "file,lines,bytes,fnv1a")
	manifestScan(outDir, "", w)
	closeBW(f, w)
}
