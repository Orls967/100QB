package main

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"time"
)

func manifestScan(dir, prefix string, lines *[]string) {
	ents, _ := os.ReadDir(dir)
	sort.Slice(ents, func(i, j int) bool { return ents[i].Name() < ents[j].Name() })
	for _, e := range ents {
		nm := e.Name()
		if nm[0] == '.' {
			continue
		}
		p := filepath.Join(dir, nm)
		if e.IsDir() {
			manifestScan(p, prefix+nm+"/", lines)
		} else if nm != "files_manifest.csv" {
			b, _ := os.ReadFile(p)
			lc := bytes.Count(b, []byte{'\n'})
			*lines = append(*lines, fmt.Sprintf("%s%s,%d,%d,%016x", prefix, nm, lc, len(b), fnv1aBytes(b)))
		}
	}
}

func main() {
	outDir := "out/go"
	os.MkdirAll(outDir+"/figures", 0755)
	cfg := loadConfig("config.kv")
	cmd := "all"
	if len(os.Args) > 1 {
		cmd = os.Args[1]
	}
	if cmd == "s01" {
		writeS01ToS09(cfg, outDir)
		return
	}
	tot := 0.0
	for st := 1; st <= 9; st++ {
		fmt.Printf("[go] tahap %d/40 (s%02d)\n", st, st)
	}
	t0 := time.Now()
	writeS01ToS09(cfg, outDir)
	tot += time.Since(t0).Seconds()
	fmt.Println("[go] tahap 10/40 (s10)")
	t0 = time.Now()
	sim := newSim(cfg, true, outDir)
	simRunDays(sim, cfg.numDays, 400000, 0, -1, 1)
	sim.closeLogs()
	tot += time.Since(t0).Seconds()
	for st := 11; st <= 38; st++ {
		fmt.Printf("[go] tahap %d/40 (s%02d)\n", st, st)
	}
	t0 = time.Now()
	writeS11ToS38(sim, cfg, outDir)
	tot += time.Since(t0).Seconds()
	fmt.Println("[go] tahap 39/40 (s39)")
	t0 = time.Now()
	writeS39SVGs(sim, outDir)
	tot += time.Since(t0).Seconds()
	fmt.Println("[go] tahap 40/40 (s40)")
	writeFile(outDir, "s40_summary.json", fmt.Sprintf("{\"lang\":\"go\",\"days\":5,\"ev_per_day\":400000,\"ev_total\":%d,\"toolchain\":\"go 1.27.4\",\"total_sec\":%.3f}\n", cfg.numDays*400000, tot))
	var mlines = []string{"file,lines,bytes,fnv1a"}
	manifestScan(outDir, "", &mlines)
	writeFile(outDir, "files_manifest.csv", fmt.Sprintf("%s\n", func() string {
		s := ""
		for i, l := range mlines {
			if i > 0 {
				s += "\n"
			}
			s += l
		}
		return s
	}()))
}
