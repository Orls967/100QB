package main

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"sort"
	"time"
)

var piVal float64 = 3.14159265358979323846
var nsAll = [6]int{15, 31, 63, 127, 255, 511}
var mnames = [7]string{"jacobi", "gs", "rbgs", "sor", "cg", "pcg_ssor", "mg"}

//go:noinline
func fm(a, b float64) float64 { return a * b }

func sinSpec(x float64) float64 {
	sgn := 1.0
	if x < 0 {
		x = -x
		sgn = -sgn
	}
	for x > piVal {
		x = x - piVal
		sgn = -sgn
	}
	if x > 0.5*piVal {
		x = piVal - x
	}
	x2 := x * x
	c := [10]float64{
		1.0, -1.0 / 6.0, 1.0 / 120.0, -1.0 / 5040.0, 1.0 / 362880.0, -1.0 / 39916800.0,
		1.0 / 6227020800.0, -1.0 / 1307674368000.0, 1.0 / 355687428096000.0, -1.0 / 121645100408832000.0,
	}
	p := c[9]
	for k := 8; k >= 0; k-- {
		p = fm(p, x2) + c[k]
	}
	return sgn * (x * p)
}
func cosSpec(x float64) float64 { return sinSpec(0.5*piVal - x) }
func omegaOpt(N int) float64 {
	m := map[int]float64{
		15: 1.673513677715992, 31: 1.8214651907890225, 63: 1.906454701582762,
		127: 1.952093233850055, 255: 1.9757544535797149, 511: 1.9878030696593354,
	}
	return m[N]
}
func fnv1aF64(v []float64) uint64 {
	h := uint64(0xCBF29CE484222325)
	for _, x := range v {
		u := math.Float64bits(x)
		for b := 0; b < 8; b++ {
			h = (h ^ ((u >> (8 * b)) & 0xFF)) * 0x100000001B3
		}
	}
	return h
}
func nbr4(N int, u []float64, i, j int) (w, e, s, n float64) {
	k := i*N + j
	if j > 0 {
		w = u[k-1]
	}
	if j+1 < N {
		e = u[k+1]
	}
	if i > 0 {
		s = u[k-N]
	}
	if i+1 < N {
		n = u[k+N]
	}
	return
}
func applyA(N int, u, out []float64) {
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			k := i*N + j
			w, e, s, n := nbr4(N, u, i, j)
			out[k] = (((float64(4.0*u[k]) - w) - e) - s) - n
		}
	}
}
func dotSeq(a, b []float64) float64 {
	s := 0.0
	for i := range a {
		s = s + fm(a[i], b[i])
	}
	return s
}
func norm2Seq(a []float64) float64 { return math.Sqrt(dotSeq(a, a)) }
func maxAbsDiff(a, b []float64) float64 {
	m := 0.0
	for i := range a {
		if d := math.Abs(a[i] - b[i]); d > m {
			m = d
		}
	}
	return m
}
func buildProb(N int) ([]float64, []float64) {
	h := 1.0 / float64(N+1)
	h2 := h * h
	sx, sy := make([]float64, N), make([]float64, N)
	for i := 0; i < N; i++ {
		sx[i] = sinSpec(piVal * float64(i+1) * h)
		sy[i] = sinSpec(2.0 * piVal * float64(i+1) * h)
	}
	b, ue := make([]float64, N*N), make([]float64, N*N)
	c5 := float64(5.0*piVal) * piVal
	for i := 0; i < N; i++ {
		x := float64(i+1) * h
		xx := x * (1.0 - x)
		for j := 0; j < N; j++ {
			y := float64(j+1) * h
			yy := y * (1.0 - y)
			sxy := sx[i] * sy[j]
			ue[i*N+j] = sxy + fm(xx, yy)
			b[i*N+j] = h2 * ((fm(c5, sxy) + 2.0*xx) + 2.0*yy)
		}
	}
	return b, ue
}
func computeRes(N int, u, b, r []float64) float64 {
	applyA(N, u, r)
	for i := range r {
		r[i] = b[i] - r[i]
	}
	return norm2Seq(r)
}
func rbgsSweep(N int, u, b []float64) {
	for col := 0; col < 2; col++ {
		for i := 0; i < N; i++ {
			for j := 0; j < N; j++ {
				if ((i + j) & 1) != col {
					continue
				}
				k := i*N + j
				w, e, s, n := nbr4(N, u, i, j)
				u[k] = 0.25 * ((((b[k] + w) + e) + s) + n)
			}
		}
	}
}
func ssorPre(N int, r, y, z []float64) {
	for i := 0; i < N; i++ {
		for j := 0; j < N; j++ {
			k := i*N + j
			yw, ys := 0.0, 0.0
			if j > 0 {
				yw = y[k-1]
			}
			if i > 0 {
				ys = y[k-N]
			}
			y[k] = 0.25 * ((r[k] + yw) + ys)
		}
	}
	for i := N - 1; i >= 0; i-- {
		for j := N - 1; j >= 0; j-- {
			k := i*N + j
			ze, zn := 0.0, 0.0
			if j+1 < N {
				ze = z[k+1]
			}
			if i+1 < N {
				zn = z[k+N]
			}
			z[k] = y[k] + float64(0.25*(ze+zn))
		}
	}
}
func restrictFW(Nf int, rf []float64) []float64 {
	Nc := (Nf - 1) / 2
	rc := make([]float64, Nc*Nc)
	for ic := 0; ic < Nc; ic++ {
		fi := 2*ic + 1
		for jc := 0; jc < Nc; jc++ {
			fj := 2*jc + 1
			k := fi*Nf + fj
			edge := ((rf[k-1] + rf[k+1]) + rf[k-Nf]) + rf[k+Nf]
			corn := ((rf[k-Nf-1] + rf[k-Nf+1]) + rf[k+Nf-1]) + rf[k+Nf+1]
			rc[ic*Nc+jc] = 4.0 * ((float64(0.25*rf[k]) + float64(0.125*edge)) + float64(0.0625*corn))
		}
	}
	return rc
}
func prolongAdd(Nf int, uf, ec []float64) {
	Nc := (Nf - 1) / 2
	for ic := 0; ic < Nc; ic++ {
		fi := 2*ic + 1
		for jc := 0; jc < Nc; jc++ {
			fj := 2*jc + 1
			k := fi*Nf + fj
			v := ec[ic*Nc+jc]
			uf[k] += v
			for _, d := range [4]int{-1, 1, -Nf, Nf} {
				uf[k+d] += 0.5 * v
			}
			for _, d := range [4]int{-Nf - 1, -Nf + 1, Nf - 1, Nf + 1} {
				uf[k+d] += 0.25 * v
			}
		}
	}
}
func solveDense3(b, x []float64) {
	var M [9][10]float64
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			r := i*3 + j
			M[r][r], M[r][9] = 4.0, b[r]
			if j > 0 {
				M[r][r-1] = -1
			}
			if j < 2 {
				M[r][r+1] = -1
			}
			if i > 0 {
				M[r][r-3] = -1
			}
			if i < 2 {
				M[r][r+3] = -1
			}
		}
	}
	for c := 0; c < 9; c++ {
		piv := c
		for r := c + 1; r < 9; r++ {
			if math.Abs(M[r][c]) > math.Abs(M[piv][c]) {
				piv = r
			}
		}
		M[c], M[piv] = M[piv], M[c]
		for r := c + 1; r < 9; r++ {
			f := M[r][c] / M[c][c]
			for k := c; k < 10; k++ {
				M[r][k] = M[r][k] - fm(f, M[c][k])
			}
		}
	}
	for r := 8; r >= 0; r-- {
		s := M[r][9]
		for k := r + 1; k < 9; k++ {
			s = s - fm(M[r][k], x[k])
		}
		x[r] = s / M[r][r]
	}
}
func vcycle(N int, u, b []float64) {
	if N == 3 {
		solveDense3(b, u)
		return
	}
	rbgsSweep(N, u, b)
	rbgsSweep(N, u, b)
	r := make([]float64, N*N)
	computeRes(N, u, b, r)
	rc := restrictFW(N, r)
	Nc := (N - 1) / 2
	ec := make([]float64, Nc*Nc)
	vcycle(Nc, ec, rc)
	prolongAdd(N, u, ec)
	rbgsSweep(N, u, b)
}
func fmgSolve(N int, u, b []float64) {
	if N == 3 {
		solveDense3(b, u)
		return
	}
	bc := restrictFW(N, b)
	Nc := (N - 1) / 2
	uc := make([]float64, Nc*Nc)
	fmgSolve(Nc, uc, bc)
	clear(u)
	prolongAdd(N, u, uc)
	vcycle(N, u, b)
}

type SolveRes struct {
	iters               int
	relRes, errInf, rho float64
	histHash, solHash   uint64
}

var (
	gHist63  [7][]float64
	gRes     [7][6]SolveRes
	gDiscErr [6]float64
	gStageS  [19]float64
)

func runSolver(meth, N int, save63 bool) SolveRes {
	n2 := N * N
	b, ue := buildProb(N)
	u, r, p, ap, z := make([]float64, n2), make([]float64, n2), make([]float64, n2), make([]float64, n2), make([]float64, n2)
	nb, om := norm2Seq(b), omegaOpt(N)
	maxIt := [7]int{100000, 50000, 50000, 10000, 5000, 3000, 100}[meth]
	hist := make([]float64, 0, 256)
	rr, rz := 0.0, 0.0
	if meth == 4 || meth == 5 {
		copy(r, b)
		rr = dotSeq(r, r)
	}
	if meth == 4 {
		copy(p, r)
	}
	if meth == 5 {
		ssorPre(N, r, ap, z)
		copy(p, z)
		rz = dotSeq(r, z)
	}
	for it := 1; it <= maxIt; it++ {
		rel := 0.0
		if meth == 0 {
			for i := 0; i < N; i++ {
				for j := 0; j < N; j++ {
					k := i*N + j
					w, e, s, n := nbr4(N, u, i, j)
					z[k] = 0.25 * ((((b[k] + w) + e) + s) + n)
				}
			}
			u, z = z, u
			rel = computeRes(N, u, b, r) / nb
		} else if meth == 1 || meth == 3 {
			for i := 0; i < N; i++ {
				for j := 0; j < N; j++ {
					k := i*N + j
					w, e, s, n := nbr4(N, u, i, j)
					gs := 0.25 * ((((b[k] + w) + e) + s) + n)
					if meth == 1 {
						u[k] = gs
					} else {
						u[k] = fm(1.0-om, u[k]) + fm(om, gs)
					}
				}
			}
			rel = computeRes(N, u, b, r) / nb
		} else if meth == 2 {
			rbgsSweep(N, u, b)
			rel = computeRes(N, u, b, r) / nb
		} else if meth == 4 || meth == 5 {
			applyA(N, p, ap)
			num := rr
			if meth == 5 {
				num = rz
			}
			alpha := num / dotSeq(p, ap)
			for k := 0; k < n2; k++ {
				u[k] = u[k] + fm(alpha, p[k])
				r[k] = r[k] - fm(alpha, ap[k])
			}
			if meth == 4 {
				rrNew := dotSeq(r, r)
				rel = math.Sqrt(rrNew) / nb
				hist = append(hist, rel)
				if rel <= 1e-8 {
					break
				}
				beta := rrNew / rr
				for k := 0; k < n2; k++ {
					p[k] = r[k] + fm(beta, p[k])
				}
				rr = rrNew
			} else {
				rel = norm2Seq(r) / nb
				hist = append(hist, rel)
				if rel <= 1e-8 {
					break
				}
				ssorPre(N, r, ap, z)
				rzNew := dotSeq(r, z)
				beta := rzNew / rz
				for k := 0; k < n2; k++ {
					p[k] = z[k] + fm(beta, p[k])
				}
				rz = rzNew
			}
			continue
		} else {
			vcycle(N, u, b)
			rel = computeRes(N, u, b, r) / nb
		}
		hist = append(hist, rel)
		if rel <= 1e-8 {
			break
		}
	}
	used := len(hist)
	rho := 0.0
	if used >= 51 {
		rho = math.Pow(hist[used-1]/hist[used-51], 0.02)
	} else if used >= 2 {
		rho = math.Pow(hist[used-1]/hist[0], 1.0/float64(used-1))
	}
	if save63 && N == 63 {
		gHist63[meth] = hist
	}
	return SolveRes{used, hist[used-1], maxAbsDiff(u, ue), rho, fnv1aF64(hist), fnv1aF64(u)}
}

func s01() {
	u := []float64{1, 2, 3, 4, 5, 6, 7, 8, 9}
	au := make([]float64, 9)
	applyA(3, u, au)
	h1, h2, h3 := fnv1aF64(u), fnv1aF64(au), fnv1aF64([]float64{1.0, -2.5, 3.25})
	ok := 0
	if maxAbsDiff(au, []float64{-2, -1, 4, 3, 0, 7, 16, 11, 22}) == 0 && h1 == 0xc56de10c6671e6c2 && h2 == 0x59ed72ec82ee959c && h3 == 0xafc5dc1ca48a2366 {
		ok = 1
	}
	f, _ := os.Create("out/go/s01_selftest.txt")
	fmt.Fprintf(f, "ok,h_u3,h_Au3,h_test\n%d,%016x,%016x,%016x\n", ok, h1, h2, h3)
	f.Close()
	if ok == 0 {
		os.Exit(1)
	}
}
func s02() {
	f, _ := os.Create("out/go/s02_problem.csv")
	fmt.Fprintln(f, "N,norm_b,h_b,h_ue")
	for _, N := range nsAll {
		b, ue := buildProb(N)
		fmt.Fprintf(f, "%d,%.17g,%016x,%016x\n", N, norm2Seq(b), fnv1aF64(b), fnv1aF64(ue))
	}
	f.Close()
}
func s03() {
	f, _ := os.Create("out/go/s03_eigen.csv")
	fmt.Fprintln(f, "N,j,k,lambda_exact,rel_err_inf")
	modes := [3][2]int{{1, 1}, {1, 2}, {3, 5}}
	for _, N := range nsAll[:3] {
		h := 1.0 / float64(N+1)
		v, av := make([]float64, N*N), make([]float64, N*N)
		for _, md := range modes {
			lam := 4.0 - 2.0*cosSpec(float64(md[0])*piVal*h) - 2.0*cosSpec(float64(md[1])*piVal*h)
			for r := 0; r < N; r++ {
				for c := 0; c < N; c++ {
					v[r*N+c] = sinSpec(float64(md[0])*piVal*float64(r+1)*h) * sinSpec(float64(md[1])*piVal*float64(c+1)*h)
				}
			}
			applyA(N, v, av)
			me := 0.0
			for q := range v {
				if d := math.Abs(av[q] - float64(lam*v[q])); d > me {
					me = d
				}
			}
			fmt.Fprintf(f, "%d,%d,%d,%.17g,%.3e\n", N, md[0], md[1], lam, me/lam)
		}
	}
	f.Close()
}
func runM(m int, fn string, nmax int) {
	f, _ := os.Create(fn)
	fmt.Fprintln(f, "N,iters,rel_res,err_inf,rho_emp")
	for i := 0; i < nmax; i++ {
		r := runSolver(m, nsAll[i], nsAll[i] == 63)
		gRes[m][i] = r
		fmt.Fprintf(f, "%d,%d,%.10e,%.10e,%.10f\n", nsAll[i], r.iters, r.relRes, r.errInf, r.rho)
	}
	f.Close()
}
func s11() {
	f, _ := os.Create("out/go/s11_fmg.csv")
	fmt.Fprintln(f, "N,fmg_err_inf,disc_err_inf,ratio")
	for i, N := range nsAll {
		b, ue := buildProb(N)
		u := make([]float64, N*N)
		fmgSolve(N, u, b)
		fe := maxAbsDiff(u, ue)
		gDiscErr[i] = gRes[6][i].errInf
		fmt.Fprintf(f, "%d,%.10e,%.10e,%.6f\n", N, fe, gDiscErr[i], fe/gDiscErr[i])
	}
	f.Close()
}
func s12() {
	f, _ := os.Create("out/go/s12_order.csv")
	fmt.Fprintln(f, "N,h,err_inf,slope")
	mx, my, num, den := 0.0, 0.0, 0.0, 0.0
	for i, N := range nsAll {
		mx += math.Log(1.0 / float64(N+1))
		my += math.Log(gDiscErr[i])
	}
	mx /= 6
	my /= 6
	for i, N := range nsAll {
		dx := math.Log(1.0/float64(N+1)) - mx
		num += dx * (math.Log(gDiscErr[i]) - my)
		den += dx * dx
	}
	for i, N := range nsAll {
		fmt.Fprintf(f, "%d,%.10e,%.10e,%.6f\n", N, 1.0/float64(N+1), gDiscErr[i], num/den)
	}
	f.Close()
}
func s13() {
	f, _ := os.Create("out/go/s13_scaling.csv")
	fmt.Fprintln(f, "method,N,iters,slope,slope_se")
	for m := 0; m < 7; m++ {
		cnt := 6
		if m <= 2 {
			cnt = 4
		}
		mx, my, num, den, rss := 0.0, 0.0, 0.0, 0.0, 0.0
		for i := 0; i < cnt; i++ {
			mx += math.Log(float64(nsAll[i]))
			my += math.Log(float64(gRes[m][i].iters))
		}
		mx /= float64(cnt)
		my /= float64(cnt)
		for i := 0; i < cnt; i++ {
			dx := math.Log(float64(nsAll[i])) - mx
			num += dx * (math.Log(float64(gRes[m][i].iters)) - my)
			den += dx * dx
		}
		sl := num / den
		ic := my - sl*mx
		for i := 0; i < cnt; i++ {
			e := math.Log(float64(gRes[m][i].iters)) - (sl*math.Log(float64(nsAll[i])) + ic)
			rss += e * e
		}
		for i := 0; i < cnt; i++ {
			fmt.Fprintf(f, "%s,%d,%d,%.6f,%.6f\n", mnames[m], nsAll[i], gRes[m][i].iters, sl, math.Sqrt((rss/float64(cnt-2))/den))
		}
	}
	f.Close()
}
func s14() {
	f, _ := os.Create("out/go/s14_history.csv")
	bw := bufio.NewWriter(f)
	fmt.Fprintln(bw, "method,iter,rel_res")
	for m := 0; m < 7; m++ {
		ln := len(gHist63[m])
		step := max(1, ln/200)
		for k := 0; k < ln; k += step {
			fmt.Fprintf(bw, "%s,%d,%.10e\n", mnames[m], k+1, gHist63[m][k])
		}
		if (ln-1)%step != 0 {
			fmt.Fprintf(bw, "%s,%d,%.10e\n", mnames[m], ln, gHist63[m][ln-1])
		}
	}
	bw.Flush()
	f.Close()
}
func s15() {
	f, _ := os.Create("out/go/s15_mg_factor.csv")
	fmt.Fprintln(f, "N,iters,mg_factor_last5")
	for _, N := range []int{63, 127} {
		b, _ := buildProb(N)
		u, r := make([]float64, N*N), make([]float64, N*N)
		nb := norm2Seq(b)
		var h []float64
		for it := 1; it <= 20; it++ {
			vcycle(N, u, b)
			rel := computeRes(N, u, b, r) / nb
			h = append(h, rel)
			if rel <= 1e-8 {
				break
			}
		}
		it := len(h)
		fmt.Fprintf(f, "%d,%d,%.8f\n", N, it, math.Pow(h[it-1]/h[it-6], 0.2))
	}
	f.Close()
}
func s16() {
	flops := [7]int{9, 9, 9, 12, 21, 35, 36}
	f, _ := os.Create("out/go/s16_cost.csv")
	fmt.Fprintln(f, "method,N,flops_per_unk_iter,ns_med,ns_mad")
	N := 63
	n2 := N * N
	b, _ := buildProb(N)
	u, z := make([]float64, n2), make([]float64, n2)
	for m := 0; m < 7; m++ {
		iters := 80
		if m == 6 {
			iters = 10
		}
		vals, dev := make([]float64, 7), make([]float64, 7)
		for rp := 0; rp < 7; rp++ {
			clear(u)
			t0 := time.Now()
			for it := 0; it < iters; it++ {
				if m == 6 {
					vcycle(N, u, b)
				} else if m == 2 {
					rbgsSweep(N, u, b)
				} else {
					applyA(N, u, z)
					for k := 0; k < n2; k++ {
						u[k] += 0.25 * (b[k] - z[k])
					}
				}
			}
			vals[rp] = float64(time.Since(t0).Nanoseconds()) / float64(iters*n2)
		}
		sort.Float64s(vals)
		med := vals[3]
		for rp := 0; rp < 7; rp++ {
			dev[rp] = math.Abs(vals[rp] - med)
		}
		sort.Float64s(dev)
		fmt.Fprintf(f, "%s,%d,%d,%.4f,%.4f\n", mnames[m], N, flops[m], med, dev[3])
	}
	f.Close()
}
func s17() {
	f, _ := os.Create("out/go/s17_parity.csv")
	fmt.Fprintln(f, "method,N,iters,hist_fnv1a,sol_fnv1a")
	for m := 0; m < 7; m++ {
		for ni := 0; ni < 2; ni++ {
			fmt.Fprintf(f, "%s,%d,%d,%016x,%016x\n", mnames[m], nsAll[ni], gRes[m][ni].iters, gRes[m][ni].histHash, gRes[m][ni].solHash)
		}
	}
	f.Close()
}
func s18() {
	cols := [7]string{"#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"}
	hdr := `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>`
	f1, _ := os.Create("out/go/figures/history.svg")
	fmt.Fprint(f1, hdr)
	for m := 0; m < 7; m++ {
		fmt.Fprintf(f1, `<polyline fill="none" stroke="%s" stroke-width="1.4" points="`, cols[m])
		ln := len(gHist63[m])
		step := max(1, ln/150)
		for k := 0; k < ln; k += step {
			fmt.Fprintf(f1, "%.1f,%.1f ", 60.0+500.0*math.Log10(float64(k+1))/math.Log10(float64(len(gHist63[0]))), 310.0-280.0*math.Min(1, math.Max(0, -math.Log10(math.Max(1e-12, gHist63[m][k]))/9.0)))
		}
		fmt.Fprint(f1, `"/>`)
	}
	fmt.Fprintln(f1, `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go Residual History (N=63)</text></svg>`)
	f1.Close()
	f2, _ := os.Create("out/go/figures/scaling.svg")
	fmt.Fprint(f2, hdr)
	for m := 0; m < 7; m++ {
		cnt := 6
		if m <= 2 {
			cnt = 4
		}
		fmt.Fprintf(f2, `<polyline fill="none" stroke="%s" stroke-width="1.5" points="`, cols[m])
		for i := 0; i < cnt; i++ {
			fmt.Fprintf(f2, "%.1f,%.1f ", 60.0+500.0*(math.Log10(float64(nsAll[i]))-math.Log10(15))/(math.Log10(511)-math.Log10(15)), 310.0-280.0*math.Log10(float64(gRes[m][i].iters))/5.0)
		}
		fmt.Fprint(f2, `"/>`)
	}
	fmt.Fprintln(f2, `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go Iterations vs N</text></svg>`)
	f2.Close()
	f3, _ := os.Create("out/go/figures/order.svg")
	fmt.Fprint(f3, hdr, `<polyline fill="none" stroke="#1f77b4" stroke-width="1.8" points="`)
	for i := 0; i < 6; i++ {
		fmt.Fprintf(f3, "%.1f,%.1f ", 60.0+500.0*(float64(i)/5.0), 310.0-280.0*(math.Log10(gDiscErr[i])+5.0)/4.0)
	}
	fmt.Fprintln(f3, `"/><text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Go Discretization Error vs h (Order ~2)</text></svg>`)
	f3.Close()
}
func s19() {
	tot := 0.0
	for _, s := range gStageS {
		tot += s
	}
	f, _ := os.Create("out/go/s19_summary.json")
	fmt.Fprintf(f, "{\"lang\":\"go\",\"N_max\":511,\"N_stat_max\":127,\"toolchain\":\"go 1.27.4\",\"total_sec\":%.3f}\n", tot)
	f.Close()
}
func main() {
	os.MkdirAll("out/go/figures", 0755)
	cmd := "all"
	if len(os.Args) > 1 {
		cmd = os.Args[1]
	}
	st := []struct {
		nm string
		fn func()
	}{
		{"s01", s01}, {"s02", s02}, {"s03", s03}, {"s04", func() { runM(0, "out/go/s04_jacobi.csv", 4) }},
		{"s05", func() { runM(1, "out/go/s05_gs.csv", 4) }}, {"s06", func() { runM(2, "out/go/s06_rbgs.csv", 4) }},
		{"s07", func() { runM(3, "out/go/s07_sor.csv", 6) }}, {"s08", func() { runM(4, "out/go/s08_cg.csv", 6) }},
		{"s09", func() { runM(5, "out/go/s09_pcg_ssor.csv", 6) }}, {"s10", func() { runM(6, "out/go/s10_mg.csv", 6) }},
		{"s11", s11}, {"s12", s12}, {"s13", s13}, {"s14", s14}, {"s15", s15}, {"s16", s16}, {"s17", s17}, {"s18", s18}, {"s19", s19},
	}
	if cmd == "all" {
		for i, x := range st {
			fmt.Printf("[go] tahap %d/19 (%s)\n", i+1, x.nm)
			t0 := time.Now()
			x.fn()
			gStageS[i] = time.Since(t0).Seconds()
		}
		return
	}
	for _, x := range st {
		if cmd == x.nm {
			x.fn()
			return
		}
	}
	os.Exit(1)
}
