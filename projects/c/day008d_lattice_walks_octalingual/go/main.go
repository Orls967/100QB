package main

import (
	"bufio"
	"fmt"
	"math"
	"math/bits"
	"os"
	"sort"
	"time"
)

const (
	masterSeed                           uint64 = 20261126
	wWalks, wRuin, benchSteps, tExactMax        = 200000, 20000, 10000000, 5000
)

var (
	tsMsd                   = [13]int{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096}
	tsRet                   = [4]int{10, 100, 1000, 10000}
	gMsd, gMsdSe            [3][13]float64
	gRetSim, gRetSe, gRetEx [3][4]float64
	gRuinP, gRuinPe         [6]float64
	gS031024                [3][]float64
	gS0410k                 [3][]uint8
	gStageSec               [12]float64
)

func sm64(z uint64) uint64 {
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
	z = (z ^ (z >> 27)) * 0x94D049BB133111EB
	return z ^ (z >> 31)
}

type Rng struct{ s [4]uint64 }

func newRng(seed uint64) Rng {
	var r Rng
	st := seed
	for i := 0; i < 4; i++ {
		st += 0x9E3779B97F4A7C15
		r.s[i] = sm64(st)
	}
	return r
}
func (r *Rng) u64() uint64 {
	res := bits.RotateLeft64(r.s[0]+r.s[3], 23) + r.s[0]
	t := r.s[1] << 17
	r.s[2] ^= r.s[0]
	r.s[3] ^= r.s[1]
	r.s[1] ^= r.s[2]
	r.s[0] ^= r.s[3]
	r.s[2] ^= t
	r.s[3] = bits.RotateLeft64(r.s[3], 45)
	return res
}
func (r *Rng) u01() float64 { return float64(r.u64()>>11) * math.Exp2(-53) }
func (r *Rng) rint(n int) int {
	k := int(math.Floor(r.u01() * float64(n)))
	if k >= n {
		return n - 1
	}
	return k
}
func (r *Rng) step(d int, p *[3]int) {
	k := r.rint(2 * d)
	p[k>>1] += 1 - 2*(k&1)
}
func taskSeed(d, st, idx int) uint64 {
	h := uint64(0x517CC1B727220A95)
	h = sm64(h ^ uint64(d))
	h = sm64(h ^ uint64(st))
	h = sm64(h ^ uint64(idx))
	return sm64(masterSeed ^ h)
}
func fnv1aI32(v []int32) uint64 {
	h := uint64(0xCBF29CE484222325)
	for _, val := range v {
		u := uint32(val)
		for b := 0; b < 4; b++ {
			h = (h ^ uint64((u>>(8*b))&0xFF)) * 0x100000001B3
		}
	}
	return h
}

func s01() {
	r := newRng(masterSeed)
	x0, x1, x2 := r.u64(), r.u64(), r.u64()
	ok := x0 == 0xA9920E0E7DAB3EF0 && x1 == 0x033B9CEB4B3BA7A0 && x2 == 0xE7473A6A39AFCEA6
	w := [5]float64{1.0, 1.0 / 3.0, 0, 0, 0}
	for n := 1; n < 4; n++ {
		nf := float64(n)
		w[n+1] = ((10*nf*nf+10*nf+3)*w[n] - nf*nf*w[n-1]) / (9 * (nf + 1) * (nf + 1))
	}
	eok := math.Abs(w[2]-15.0/81.0) < 1e-14 && math.Abs(w[3]-93.0/729.0) < 1e-14 && math.Abs(w[4]-639.0/6561.0) < 1e-14
	oi, ei := 0, 0
	if ok && eok {
		oi = 1
	}
	if eok {
		ei = 1
	}
	f, _ := os.Create("out/go/s01_selftest.txt")
	fmt.Fprintf(f, "ok,x0,x1,x2,enum_ok\n%d,%016x,%016x,%016x,%d\n", oi, x0, x1, x2, ei)
	f.Close()
	if oi == 0 {
		os.Exit(1)
	}
}

func s02() {
	f, _ := os.Create("out/go/s02_parity.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,walk,x0,x1,x2,fnv1a_hex")
	for d := 1; d <= 3; d++ {
		coords := make([]int32, 0, 8*d)
		var posAll [8][3]int32
		for w := 0; w < 8; w++ {
			r := newRng(taskSeed(d, 2, w))
			var p [3]int
			for s := 0; s < 200; s++ {
				r.step(d, &p)
			}
			for a := 0; a < 3; a++ {
				posAll[w][a] = int32(p[a])
			}
			for a := 0; a < d; a++ {
				coords = append(coords, int32(p[a]))
			}
		}
		h := fnv1aI32(coords)
		for w := 0; w < 8; w++ {
			fmt.Fprintf(bw, "%d,%d,%d,%d,%d,%016x\n", d, w, posAll[w][0], posAll[w][1], posAll[w][2], h)
		}
	}
	bw.Flush()
	f.Close()
}

func s03() {
	f, _ := os.Create("out/go/s03_msd.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,t,msd,se,exact")
	for d := 1; d <= 3; d++ {
		gS031024[d-1] = make([]float64, wWalks)
		var sum1, sum2 [13]float64
		for w := 0; w < wWalks; w++ {
			r := newRng(taskSeed(d, 3, w))
			var p [3]int
			ti := 0
			for t := 1; t <= 4096; t++ {
				r.step(d, &p)
				if t == tsMsd[ti] {
					r2 := float64(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
					sum1[ti] += r2
					sum2[ti] += r2 * r2
					if ti == 10 {
						gS031024[d-1][w] = r2
					}
					ti++
				}
			}
		}
		for ti := 0; ti < 13; ti++ {
			m := sum1[ti] / float64(wWalks)
			se := math.Sqrt(math.Max(0, sum2[ti]/float64(wWalks)-m*m) / float64(wWalks))
			gMsd[d-1][ti], gMsdSe[d-1][ti] = m, se
			fmt.Fprintf(bw, "%d,%d,%.8f,%.8f,%d\n", d, tsMsd[ti], m, se, tsMsd[ti])
		}
	}
	bw.Flush()
	f.Close()
}

func s04() {
	f, _ := os.Create("out/go/s04_return_sim.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,T,p_hat,se")
	for d := 1; d <= 3; d++ {
		gS0410k[d-1] = make([]uint8, wWalks)
		var hits [4]int
		for w := 0; w < wWalks; w++ {
			r := newRng(taskSeed(d, 4, w))
			var p [3]int
			retT := 0
			for t := 1; t <= 10000; t++ {
				r.step(d, &p)
				if (p[0] | p[1] | p[2]) == 0 {
					retT = t
					break
				}
			}
			if retT > 0 {
				gS0410k[d-1][w] = 1
				for ti := 0; ti < 4; ti++ {
					if retT <= tsRet[ti] {
						hits[ti]++
					}
				}
			}
		}
		for ti := 0; ti < 4; ti++ {
			ph := float64(hits[ti]) / float64(wWalks)
			se := math.Sqrt(math.Max(0, ph*(1-ph)) / float64(wWalks))
			gRetSim[d-1][ti], gRetSe[d-1][ti] = ph, se
			fmt.Fprintf(bw, "%d,%d,%.8f,%.8f\n", d, tsRet[ti], ph, se)
		}
	}
	bw.Flush()
	f.Close()
}

func s05() {
	u1, u2, w3, u3, fv := make([]float64, tExactMax+1), make([]float64, tExactMax+1), make([]float64, tExactMax+1), make([]float64, tExactMax+1), make([]float64, tExactMax+1)
	u1[0], u2[0], w3[0], u3[0], w3[1] = 1, 1, 1, 1, 1.0/3.0
	for n := 1; n < tExactMax; n++ {
		nf := float64(n)
		w3[n+1] = ((10*nf*nf+10*nf+3)*w3[n] - nf*nf*w3[n-1]) / (9 * (nf + 1) * (nf + 1))
	}
	for n := 1; n <= tExactMax; n++ {
		nf := float64(n)
		u1[n] = u1[n-1] * (2*nf - 1) / (2 * nf)
		u2[n] = u1[n] * u1[n]
		u3[n] = u1[n] * w3[n]
	}
	us := [3][]float64{u1, u2, u3}
	f, _ := os.Create("out/go/s05_return_exact.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,T,F_exact")
	for d := 1; d <= 3; d++ {
		u := us[d-1]
		acc := 0.0
		for n := 1; n <= tExactMax; n++ {
			s := u[n]
			for k := 1; k < n; k++ {
				s -= fv[k] * u[n-k]
			}
			fv[n] = s
			acc += s
			for ti := 0; ti < 4; ti++ {
				if 2*n == tsRet[ti] {
					gRetEx[d-1][ti] = acc
					fmt.Fprintf(bw, "%d,%d,%.10f\n", d, tsRet[ti], acc)
				}
			}
		}
	}
	bw.Flush()
	f.Close()
}

func s06() {
	f, _ := os.Create("out/go/s06_return_compare.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,T,p_hat,se,F_exact,z_score,p_inf_est")
	p3Inf := gRetSim[2][3] + (gRetSim[2][3]-gRetSim[2][2])/(math.Sqrt(10)-1)
	for d := 1; d <= 3; d++ {
		for ti := 0; ti < 4; ti++ {
			z, pinf := 0.0, 0.0
			if gRetSe[d-1][ti] > 1e-12 {
				z = (gRetSim[d-1][ti] - gRetEx[d-1][ti]) / gRetSe[d-1][ti]
			}
			if d == 3 && ti == 3 {
				pinf = p3Inf
			}
			fmt.Fprintf(bw, "%d,%d,%.8f,%.8f,%.10f,%.4f,%.8f\n", d, tsRet[ti], gRetSim[d-1][ti], gRetSe[d-1][ti], gRetEx[d-1][ti], z, pinf)
		}
	}
	bw.Flush()
	f.Close()
}

func s07() {
	ns, xs := [6]int{10, 10, 10, 50, 50, 50}, [6]int{2, 5, 8, 10, 25, 40}
	f, _ := os.Create("out/go/s07_ruin.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z")
	for ci := 0; ci < 6; ci++ {
		N, x0, wins, st1, st2 := ns[ci], xs[ci], 0, 0.0, 0.0
		for w := 0; w < wRuin; w++ {
			r := newRng(taskSeed(1, 7, ci*100000+w))
			x, steps := x0, 0
			for x > 0 && x < N {
				x += 1 - 2*r.rint(2)
				steps++
			}
			if x == N {
				wins++
			}
			st1 += float64(steps)
			st2 += float64(steps) * float64(steps)
		}
		ph := float64(wins) / float64(wRuin)
		pse := math.Sqrt(math.Max(0, ph*(1-ph)) / float64(wRuin))
		pe := float64(x0) / float64(N)
		th := st1 / float64(wRuin)
		tse := math.Sqrt(math.Max(0, st2/float64(wRuin)-th*th) / float64(wRuin))
		te := float64(x0 * (N - x0))
		gRuinP[ci], gRuinPe[ci] = ph, pe
		fmt.Fprintf(bw, "%d,%d,%.8f,%.8f,%.8f,%.4f,%.4f,%.4f,%.4f,%.4f\n", N, x0, ph, pse, pe, (ph-pe)/pse, th, tse, te, (th-te)/tse)
	}
	bw.Flush()
	f.Close()
}

func s08() {
	f, _ := os.Create("out/go/s08_endpoint.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count")
	var hist [40]int
	nEp := 20000
	for d := 1; d <= 3; d++ {
		m1, m2, m4 := 0.0, 0.0, 0.0
		for w := 0; w < nEp; w++ {
			r := newRng(taskSeed(d, 8, w))
			var p [3]int
			for s := 0; s < 1024; s++ {
				r.step(d, &p)
			}
			x := float64(p[0])
			m1 += x
			m2 += x * x
			m4 += x * x * x * x
			if d == 1 {
				b := int(math.Floor((x + 128.0) / 6.4))
				if b < 0 {
					b = 0
				} else if b > 39 {
					b = 39
				}
				hist[b]++
			}
		}
		m1 /= float64(nEp)
		v := m2/float64(nEp) - m1*m1
		fmt.Fprintf(bw, "%d,%.6f,%.6f,%.6f,%.6f,-1,0,0,%d\n", d, v, 1024.0/float64(d), (m4/float64(nEp))/(v*v), 3.0+(float64(d)-3.0)/1024.0, nEp)
	}
	for b := 0; b < 40; b++ {
		fmt.Fprintf(bw, "1,0,0,0,0,%d,%.2f,%.2f,%d\n", b, -128.0+float64(b)*6.4, -128.0+float64(b+1)*6.4, hist[b])
	}
	bw.Flush()
	f.Close()
}

func s09() {
	f, _ := os.Create("out/go/s09_ci.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "metric,d,estimate,boot_se,ci_lo,ci_hi")
	boot := make([]float64, 200)
	for d := 1; d <= 3; d++ {
		r1 := newRng(taskSeed(d, 9, 0))
		for b := 0; b < 200; b++ {
			s := 0.0
			for i := 0; i < wWalks; i++ {
				s += gS031024[d-1][r1.rint(wWalks)]
			}
			boot[b] = s / float64(wWalks)
		}
		sort.Float64s(boot)
		bm, bv := 0.0, 0.0
		for _, v := range boot {
			bm += v
		}
		bm /= 200.0
		for _, v := range boot {
			bv += (v - bm) * (v - bm)
		}
		fmt.Fprintf(bw, "msd_1024,%d,%.6f,%.6f,%.6f,%.6f\n", d, gMsd[d-1][10], math.Sqrt(bv/199.0), boot[5], boot[194])
		r2 := newRng(taskSeed(d, 9, 1))
		for b := 0; b < 200; b++ {
			s := 0
			for i := 0; i < wWalks; i++ {
				s += int(gS0410k[d-1][r2.rint(wWalks)])
			}
			boot[b] = float64(s) / float64(wWalks)
		}
		sort.Float64s(boot)
		bm, bv = 0.0, 0.0
		for _, v := range boot {
			bm += v
		}
		bm /= 200.0
		for _, v := range boot {
			bv += (v - bm) * (v - bm)
		}
		fmt.Fprintf(bw, "return_10000,%d,%.6f,%.6f,%.6f,%.6f\n", d, gRetSim[d-1][3], math.Sqrt(bv/199.0), boot[5], boot[194])
	}
	bw.Flush()
	f.Close()
}

func s10() {
	f, _ := os.Create("out/go/s10_fit.csv")
	fmt.Fprintln(f, "d,slope,slope_se,intercept")
	for d := 1; d <= 3; d++ {
		mx, my, num, den, rss := 0.0, 0.0, 0.0, 0.0, 0.0
		for i := 0; i < 13; i++ {
			mx += math.Log(float64(tsMsd[i]))
			my += math.Log(gMsd[d-1][i])
		}
		mx /= 13.0
		my /= 13.0
		for i := 0; i < 13; i++ {
			dx := math.Log(float64(tsMsd[i])) - mx
			num += dx * (math.Log(gMsd[d-1][i]) - my)
			den += dx * dx
		}
		sl := num / den
		ic := my - sl*mx
		for i := 0; i < 13; i++ {
			e := math.Log(gMsd[d-1][i]) - (sl*math.Log(float64(tsMsd[i])) + ic)
			rss += e * e
		}
		fmt.Fprintf(f, "%d,%.6f,%.6f,%.6f\n", d, sl, math.Sqrt((rss/11.0)/den), ic)
	}
	f.Close()
}

func s11() {
	f, _ := os.Create("out/go/s11_bench.csv")
	fmt.Fprintln(f, "d,steps_per_rep,reps,ns_med,ns_mad")
	for d := 1; d <= 3; d++ {
		vals, dev := make([]float64, 7), make([]float64, 7)
		for rep := 0; rep < 7; rep++ {
			r := newRng(taskSeed(d, 11, rep))
			var p [3]int
			t0 := time.Now()
			for s := 0; s < benchSteps; s++ {
				r.step(d, &p)
			}
			vals[rep] = float64(time.Since(t0).Nanoseconds()) / float64(benchSteps)
			_ = p[0]
		}
		sort.Float64s(vals)
		med := vals[3]
		for rep := 0; rep < 7; rep++ {
			dev[rep] = math.Abs(vals[rep] - med)
		}
		sort.Float64s(dev)
		fmt.Fprintf(f, "%d,%d,7,%.4f,%.4f\n", d, benchSteps, med, dev[3])
	}
	f.Close()
}

func s12() {
	cols := [3]string{"#1f77b4", "#ff7f0e", "#2ca02c"}
	f1, _ := os.Create("out/go/figures/msd.svg")
	fmt.Fprint(f1, `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>`)
	for d := 0; d < 3; d++ {
		fmt.Fprintf(f1, `<polyline fill="none" stroke="%s" stroke-width="1.6" points="`, cols[d])
		for i := 0; i < 13; i++ {
			fmt.Fprintf(f1, "%.1f,%.1f ", 60.0+500.0*math.Log10(float64(tsMsd[i]))/math.Log10(4096), 310.0-280.0*math.Log10(math.Max(1, gMsd[d][i]))/math.Log10(4096))
		}
		fmt.Fprint(f1, `"/>`)
	}
	fmt.Fprintln(f1, `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go MSD Log-Log (y=t)</text></svg>`)
	f1.Close()

	f2, _ := os.Create("out/go/figures/return_prob.svg")
	fmt.Fprint(f2, `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>`)
	for d := 0; d < 3; d++ {
		fmt.Fprintf(f2, `<polyline fill="none" stroke="%s" stroke-width="1.6" points="`, cols[d])
		for i := 0; i < 4; i++ {
			fmt.Fprintf(f2, "%.1f,%.1f ", 60.0+500.0*(float64(i)/3.0), 310.0-280.0*gRetEx[d][i])
		}
		fmt.Fprint(f2, `"/>`)
		for i := 0; i < 4; i++ {
			fmt.Fprintf(f2, `<circle cx="%.1f" cy="%.1f" r="4" fill="%s"/>`, 60.0+500.0*(float64(i)/3.0), 310.0-280.0*gRetSim[d][i], cols[d])
		}
	}
	fmt.Fprintln(f2, `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go Return Prob F_d(T) Exact vs Sim</text></svg>`)
	f2.Close()

	f3, _ := os.Create("out/go/figures/ruin.svg")
	fmt.Fprint(f3, `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>`)
	for i := 0; i < 6; i++ {
		fmt.Fprintf(f3, `<circle cx="%.1f" cy="%.1f" r="5" fill="#1f77b4"/>`, 60.0+500.0*gRuinPe[i], 310.0-280.0*gRuinP[i])
	}
	fmt.Fprintln(f3, `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go Gambler's Ruin Exact vs Sim</text></svg>`)
	f3.Close()

	tot := 0.0
	for _, s := range gStageSec {
		tot += s
	}
	fj, _ := os.Create("out/go/s12_summary.json")
	fmt.Fprintf(fj, "{\"lang\":\"go\",\"W\":%d,\"ruin_W\":%d,\"bench_steps\":%d,\"toolchain\":\"go 1.27.4\",\"total_sec\":%.3f}\n", wWalks, wRuin, benchSteps, tot)
	fj.Close()
}

func main() {
	os.MkdirAll("out/go/figures", 0755)
	cmd := "all"
	if len(os.Args) > 1 {
		cmd = os.Args[1]
	}
	stages := []struct {
		name string
		fn   func()
	}{
		{"s01", s01}, {"s02", s02}, {"s03", s03}, {"s04", s04}, {"s05", s05}, {"s06", s06},
		{"s07", s07}, {"s08", s08}, {"s09", s09}, {"s10", s10}, {"s11", s11}, {"s12", s12},
	}
	if cmd == "all" {
		for i, st := range stages {
			fmt.Printf("[go] tahap %d/12 (%s)\n", i+1, st.name)
			t0 := time.Now()
			st.fn()
			gStageSec[i] = time.Since(t0).Seconds()
		}
		return
	}
	for _, st := range stages {
		if cmd == st.name {
			st.fn()
			return
		}
	}
	os.Exit(1)
}
