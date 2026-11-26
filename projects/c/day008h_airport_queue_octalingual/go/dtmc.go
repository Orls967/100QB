package main

import "math"

type DtmcRes struct {
	pBlock, L, Lq, U, W, Wq, lamEff float64
}

func binomProbs(m int, mu float64) [8]float64 {
	var out [8]float64
	out[0] = 1.0
	for s := 0; s < m; s++ {
		for k := m; k >= 1; k-- {
			out[k] = out[k]*(1.0-mu) + out[k-1]*mu
		}
		out[0] *= 1.0 - mu
	}
	return out
}

func solveDtmcSingle(c, K int, lam, mu float64) DtmcRes {
	var P [64][64]float64
	for n := 0; n <= K; n++ {
		m := n
		if c < m {
			m = c
		}
		bp := binomProbs(m, mu)
		for d := 0; d <= m; d++ {
			if n < K {
				P[n][n-d] += (1.0 - lam) * bp[d]
				P[n][n+1-d] += lam * bp[d]
			} else {
				P[n][n-d] += bp[d]
			}
		}
	}
	var pi [64]float64
	for i := 0; i <= K; i++ {
		pi[i] = 1.0 / float64(K+1)
	}
	for it := 0; it < 25000; it++ {
		maxd, sum := 0.0, 0.0
		for j := 0; j <= K; j++ {
			nv := 0.0
			for i := 0; i <= K; i++ {
				nv += pi[i] * P[i][j]
			}
			if df := math.Abs(nv - pi[j]); df > maxd {
				maxd = df
			}
			pi[j] = nv
			sum += nv
		}
		for j := 0; j <= K; j++ {
			pi[j] /= sum
		}
		if maxd < 1e-14 {
			break
		}
	}
	var r DtmcRes
	r.pBlock = pi[K]
	r.lamEff = lam * (1.0 - pi[K])
	for n := 0; n <= K; n++ {
		m := n
		if c < m {
			m = c
		}
		r.L += float64(n) * pi[n]
		r.Lq += float64(n-m) * pi[n]
		r.U += (float64(m) / float64(c)) * pi[n]
	}
	if r.lamEff > 0 {
		r.W = r.L / r.lamEff
		r.Wq = r.Lq / r.lamEff
	}
	return r
}

func solveGeoClosedForm(lam, mu float64) DtmcRes {
	var r DtmcRes
	r.lamEff = lam
	r.U = lam / mu
	r.L = (lam * (1.0 - lam)) / (mu - lam)
	r.W = (1.0 - lam) / (mu - lam)
	r.Wq = r.W - 1.0/mu
	r.Lq = lam * r.Wq
	return r
}

func solveDtmcTandem(K1, K2 int, lam, mu1, mu2 float64) (DtmcRes, DtmcRes, float64) {
	S := (K1 + 1) * (K2 + 1)
	var P [64][64]float64
	for n1 := 0; n1 <= K1; n1++ {
		for n2 := 0; n2 <= K2; n2++ {
			sIdx := n1*(K2+1) + n2
			for a := 0; a <= 1; a++ {
				pa := 1.0 - lam
				if a != 0 {
					pa = lam
				}
				d1Max := 0
				if n1 > 0 {
					d1Max = 1
				}
				for d1 := 0; d1 <= d1Max; d1++ {
					pd1 := 1.0
					if n1 > 0 {
						if d1 != 0 {
							pd1 = mu1
						} else {
							pd1 = 1.0 - mu1
						}
					}
					d2Max := 0
					if n2 > 0 {
						d2Max = 1
					}
					for d2 := 0; d2 <= d2Max; d2++ {
						pd2 := 1.0
						if n2 > 0 {
							if d2 != 0 {
								pd2 = mu2
							} else {
								pd2 = 1.0 - mu2
							}
						}
						nn1 := n1 - d1
						if n1 < K1 {
							nn1 = n1 + a - d1
						}
						n2Mid := n2 - d2
						nn2 := n2Mid
						if n2Mid < K2 {
							nn2 = n2Mid + d1
						}
						P[sIdx][nn1*(K2+1)+nn2] += pa * pd1 * pd2
					}
				}
			}
		}
	}
	var pi [64]float64
	for i := 0; i < S; i++ {
		pi[i] = 1.0 / float64(S)
	}
	for it := 0; it < 25000; it++ {
		maxd, sum := 0.0, 0.0
		for j := 0; j < S; j++ {
			nv := 0.0
			for i := 0; i < S; i++ {
				nv += pi[i] * P[i][j]
			}
			if df := math.Abs(nv - pi[j]); df > maxd {
				maxd = df
			}
			pi[j] = nv
			sum += nv
		}
		for j := 0; j < S; j++ {
			pi[j] /= sum
		}
		if maxd < 1e-14 {
			break
		}
	}
	var r1, r2 DtmcRes
	blk1, blk2 := 0.0, 0.0
	for n1 := 0; n1 <= K1; n1++ {
		for n2 := 0; n2 <= K2; n2++ {
			p := pi[n1*(K2+1)+n2]
			r1.L += float64(n1) * p
			r2.L += float64(n2) * p
			if n1 > 0 {
				r1.U += p
				r1.Lq += float64(n1-1) * p
			}
			if n2 > 0 {
				r2.U += p
				r2.Lq += float64(n2-1) * p
			}
			if n1 == K1 {
				blk1 += p * lam
			}
			if n1 > 0 && n2 == K2 {
				blk2 += p * mu1 * (1.0 - mu2)
			}
		}
	}
	r1.pBlock = blk1
	r1.lamEff = lam - blk1
	if r1.lamEff > 0 {
		r1.W = r1.L / r1.lamEff
		r1.Wq = r1.Lq / r1.lamEff
	}
	r2.pBlock = blk2
	r2.lamEff = r1.lamEff - blk2
	if r2.lamEff > 0 {
		r2.W = r2.L / r2.lamEff
		r2.Wq = r2.Lq / r2.lamEff
	}
	den := r1.lamEff
	if den <= 0 {
		den = 1.0
	}
	return r1, r2, blk2 / den
}
