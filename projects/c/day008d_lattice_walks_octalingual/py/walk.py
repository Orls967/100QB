from __future__ import annotations
import math, os, struct, sys, time
from pathlib import Path

MASTER_SEED = 20261126
MASK64 = 0xFFFFFFFFFFFFFFFF
W_WALKS = 4000
W_RUIN = 2000
BENCH_STEPS = 200000
T_EXACT_MAX = 5000
TS_MSD = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
TS_RET = [10, 100, 1000, 10000]
TWO_M53 = 2.0 ** -53

def sm64(z: int) -> int:
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & MASK64
    return (z ^ (z >> 31)) & MASK64

class Rng:
    def __init__(self, seed: int):
        st = seed & MASK64
        self.s = [0, 0, 0, 0]
        for i in range(4):
            st = (st + 0x9E3779B97F4A7C15) & MASK64
            self.s[i] = sm64(st)
    def u64(self) -> int:
        s = self.s
        x = (s[0] + s[3]) & MASK64
        res = ((((x << 23) | (x >> 41)) & MASK64) + s[0]) & MASK64
        t = (s[1] << 17) & MASK64
        s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3]; s[2] ^= t
        s[3] = ((s[3] << 45) | (s[3] >> 19)) & MASK64
        return res
    def u01(self) -> float: return (self.u64() >> 11) * TWO_M53
    def rint(self, n: int) -> int:
        k = int(math.floor(self.u01() * n))
        return n - 1 if k >= n else k

def task_seed(d: int, st: int, idx: int) -> int:
    h = 0x517CC1B727220A95
    for v in (d, st, idx): h = sm64(h ^ (v & MASK64))
    return sm64(MASTER_SEED ^ h)

def fnv1a_i32(coords: list[int]) -> int:
    h = 0xCBF29CE484222325
    for b in struct.pack(f"<{len(coords)}i", *coords):
        h = ((h ^ b) * 0x100000001B3) & MASK64
    return h

Path("out/py/figures").mkdir(parents=True, exist_ok=True)
G_MSD = [[0.0]*13 for _ in range(3)]
G_RET_SIM = [[0.0]*4 for _ in range(3)]
G_RET_SE = [[0.0]*4 for _ in range(3)]
G_RET_EX = [[0.0]*4 for _ in range(3)]
G_RUIN_P, G_RUIN_PE = [0.0]*6, [0.0]*6
G_S03_1024 = [[0.0]*W_WALKS for _ in range(3)]
G_S04_10K = [[0]*W_WALKS for _ in range(3)]
STAGE_SEC = [0.0]*12

def s01():
    r = Rng(MASTER_SEED)
    x0, x1, x2 = r.u64(), r.u64(), r.u64()
    ok = (x0 == 0xA9920E0E7DAB3EF0 and x1 == 0x033B9CEB4B3BA7A0 and x2 == 0xE7473A6A39AFCEA6)
    w = [1.0, 1.0/3.0, 0.0, 0.0, 0.0]
    for n in range(1, 4): w[n+1] = ((10.0*n*n + 10.0*n + 3.0)*w[n] - n*n*w[n-1]) / (9.0*(n+1)*(n+1))
    eok = abs(w[2]-15.0/81.0)<1e-14 and abs(w[3]-93.0/729.0)<1e-14 and abs(w[4]-639.0/6561.0)<1e-14
    Path("out/py/s01_selftest.txt").write_text(f"ok,x0,x1,x2,enum_ok\n{int(ok and eok)},{x0:016x},{x1:016x},{x2:016x},{int(eok)}\n")
    if not (ok and eok): sys.exit(1)

def s02():
    lines = ["d,walk,x0,x1,x2,fnv1a_hex"]
    for d in (1, 2, 3):
        coords, pos_all = [], []
        for w in range(8):
            r = Rng(task_seed(d, 2, w)); p = [0, 0, 0]
            for _ in range(200): k = r.rint(2*d); p[k >> 1] += -1 if (k & 1) else 1
            pos_all.append(p); coords.extend(p[:d])
        h = fnv1a_i32(coords)
        for w in range(8): lines.append(f"{d},{w},{pos_all[w][0]},{pos_all[w][1]},{pos_all[w][2]},{h:016x}")
    Path("out/py/s02_parity.csv").write_text("\n".join(lines) + "\n")

def s03():
    lines = ["d,t,msd,se,exact"]
    for d in (1, 2, 3):
        s1, s2 = [0.0]*13, [0.0]*13
        for w in range(W_WALKS):
            r = Rng(task_seed(d, 3, w)); p = [0, 0, 0]; ti = 0
            for t in range(1, 4097):
                k = r.rint(2*d); p[k >> 1] += -1 if (k & 1) else 1
                if t == TS_MSD[ti]:
                    r2 = float(p[0]*p[0] + p[1]*p[1] + p[2]*p[2])
                    s1[ti] += r2; s2[ti] += r2*r2
                    if ti == 10: G_S03_1024[d-1][w] = r2
                    ti += 1
        for ti in range(13):
            m = s1[ti]/W_WALKS; se = math.sqrt(max(0.0, s2[ti]/W_WALKS - m*m)/W_WALKS)
            G_MSD[d-1][ti] = m
            lines.append(f"{d},{TS_MSD[ti]},{m:.8f},{se:.8f},{TS_MSD[ti]}")
    Path("out/py/s03_msd.csv").write_text("\n".join(lines) + "\n")

def s04():
    lines = ["d,T,p_hat,se"]
    for d in (1, 2, 3):
        hits = [0, 0, 0, 0]
        for w in range(W_WALKS):
            r = Rng(task_seed(d, 4, w)); p = [0, 0, 0]; ret_t = 0
            for t in range(1, 10001):
                k = r.rint(2*d); p[k >> 1] += -1 if (k & 1) else 1
                if (p[0] | p[1] | p[2]) == 0: ret_t = t; break
            if ret_t > 0:
                G_S04_10K[d-1][w] = 1
                for ti in range(4):
                    if ret_t <= TS_RET[ti]: hits[ti] += 1
        for ti in range(4):
            ph = hits[ti]/W_WALKS; se = math.sqrt(max(0.0, ph*(1.0-ph))/W_WALKS)
            G_RET_SIM[d-1][ti] = ph; G_RET_SE[d-1][ti] = se
            lines.append(f"{d},{TS_RET[ti]},{ph:.8f},{se:.8f}")
    Path("out/py/s04_return_sim.csv").write_text("\n".join(lines) + "\n")

def s05():
    u1, u2, w3, u3, fv = [1.0]*(T_EXACT_MAX+1), [1.0]*(T_EXACT_MAX+1), [1.0]*(T_EXACT_MAX+1), [1.0]*(T_EXACT_MAX+1), [0.0]*(T_EXACT_MAX+1)
    w3[1] = 1.0/3.0
    for n in range(1, T_EXACT_MAX): w3[n+1] = ((10.0*n*n+10.0*n+3.0)*w3[n] - n*n*w3[n-1])/(9.0*(n+1)*(n+1))
    for n in range(1, T_EXACT_MAX+1):
        u1[n] = u1[n-1]*(2.0*n-1.0)/(2.0*n); u2[n] = u1[n]*u1[n]; u3[n] = u1[n]*w3[n]
    lines = ["d,T,F_exact"]
    for d, u in enumerate((u1, u2, u3), 1):
        acc = 0.0
        for n in range(1, T_EXACT_MAX+1):
            s = u[n]
            for k in range(1, n): s -= fv[k]*u[n-k]
            fv[n] = s; acc += s
            for ti in range(4):
                if 2*n == TS_RET[ti]: G_RET_EX[d-1][ti] = acc; lines.append(f"{d},{TS_RET[ti]},{acc:.10f}")
    Path("out/py/s05_return_exact.csv").write_text("\n".join(lines) + "\n")

def s06():
    lines = ["d,T,p_hat,se,F_exact,z_score,p_inf_est"]
    p3_inf = G_RET_SIM[2][3] + (G_RET_SIM[2][3] - G_RET_SIM[2][2])/(math.sqrt(10.0)-1.0)
    for d in (1, 2, 3):
        for ti in range(4):
            z = (G_RET_SIM[d-1][ti]-G_RET_EX[d-1][ti])/G_RET_SE[d-1][ti] if G_RET_SE[d-1][ti] > 1e-12 else 0.0
            lines.append(f"{d},{TS_RET[ti]},{G_RET_SIM[d-1][ti]:.8f},{G_RET_SE[d-1][ti]:.8f},{G_RET_EX[d-1][ti]:.10f},{z:.4f},{p3_inf if (d==3 and ti==3) else 0.0:.8f}")
    Path("out/py/s06_return_compare.csv").write_text("\n".join(lines) + "\n")

def s07():
    ns, xs = [10,10,10,50,50,50], [2,5,8,10,25,40]
    lines = ["N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z"]
    for ci, (n, x0) in enumerate(zip(ns, xs)):
        wins, st1, st2 = 0, 0.0, 0.0
        for w in range(W_RUIN):
            r = Rng(task_seed(1, 7, ci*100000+w)); x, steps = x0, 0
            while 0 < x < n: x += 1 if r.rint(2)==0 else -1; steps += 1
            if x == n: wins += 1
            st1 += steps; st2 += steps*steps
        ph = wins/W_RUIN; pse = math.sqrt(max(0.0, ph*(1.0-ph))/W_RUIN); pe = x0/n
        th = st1/W_RUIN; tse = math.sqrt(max(0.0, st2/W_RUIN - th*th)/W_RUIN); te = float(x0*(n-x0))
        G_RUIN_P[ci], G_RUIN_PE[ci] = ph, pe
        lines.append(f"{n},{x0},{ph:.8f},{pse:.8f},{pe:.8f},{(ph-pe)/pse:.4f},{th:.4f},{tse:.4f},{te:.4f},{(th-te)/tse:.4f}")
    Path("out/py/s07_ruin.csv").write_text("\n".join(lines) + "\n")

def s08():
    lines = ["d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count"]
    hist, n_ep = [0]*40, 4000
    for d in (1, 2, 3):
        m1 = m2 = m4 = 0.0
        for w in range(n_ep):
            r = Rng(task_seed(d, 8, w)); p = [0, 0, 0]
            for _ in range(1024): k = r.rint(2*d); p[k >> 1] += -1 if (k & 1) else 1
            x = float(p[0]); m1 += x; m2 += x*x; m4 += x*x*x*x
            if d == 1: hist[min(39, max(0, int(math.floor((x+128.0)/6.4))))] += 1
        m1 /= n_ep; v = m2/n_ep - m1*m1
        lines.append(f"{d},{v:.6f},{1024.0/d:.6f},{(m4/n_ep)/(v*v):.6f},{3.0+(d-3.0)/1024.0:.6f},-1,0,0,{n_ep}")
    for b in range(40): lines.append(f"1,0,0,0,0,{b},{-128.0+b*6.4:.2f},{-128.0+(b+1)*6.4:.2f},{hist[b]}")
    Path("out/py/s08_endpoint.csv").write_text("\n".join(lines) + "\n")

def s09():
    lines = ["metric,d,estimate,boot_se,ci_lo,ci_hi"]
    for d in (1, 2, 3):
        r1 = Rng(task_seed(d, 9, 0))
        boot = sorted(sum(G_S03_1024[d-1][r1.rint(W_WALKS)] for _ in range(W_WALKS))/W_WALKS for _ in range(200))
        bm = sum(boot)/200.0; bv = sum((v-bm)**2 for v in boot)/199.0
        lines.append(f"msd_1024,{d},{G_MSD[d-1][10]:.6f},{math.sqrt(bv):.6f},{boot[5]:.6f},{boot[194]:.6f}")
        r2 = Rng(task_seed(d, 9, 1))
        boot2 = sorted(sum(G_S04_10K[d-1][r2.rint(W_WALKS)] for _ in range(W_WALKS))/W_WALKS for _ in range(200))
        bm2 = sum(boot2)/200.0; bv2 = sum((v-bm2)**2 for v in boot2)/199.0
        lines.append(f"return_10000,{d},{G_RET_SIM[d-1][3]:.6f},{math.sqrt(bv2):.6f},{boot2[5]:.6f},{boot2[194]:.6f}")
    Path("out/py/s09_ci.csv").write_text("\n".join(lines) + "\n")

def s10():
    lines = ["d,slope,slope_se,intercept"]
    xs = [math.log(t) for t in TS_MSD]; mx = sum(xs)/13.0; den = sum((x-mx)**2 for x in xs)
    for d in (1, 2, 3):
        ys = [math.log(m) for m in G_MSD[d-1]]; my = sum(ys)/13.0
        sl = sum((xs[i]-mx)*(ys[i]-my) for i in range(13))/den; ic = my - sl*mx
        rss = sum((ys[i]-(sl*xs[i]+ic))**2 for i in range(13))
        lines.append(f"{d},{sl:.6f},{math.sqrt((rss/11.0)/den):.6f},{ic:.6f}")
    Path("out/py/s10_fit.csv").write_text("\n".join(lines) + "\n")

def s11():
    lines = ["d,steps_per_rep,reps,ns_med,ns_mad"]
    for d in (1, 2, 3):
        vals = []
        for rep in range(7):
            r = Rng(task_seed(d, 11, rep)); p = [0, 0, 0]
            t0 = time.perf_counter()
            for _ in range(BENCH_STEPS): k = r.rint(2*d); p[k >> 1] += -1 if (k & 1) else 1
            vals.append((time.perf_counter()-t0)*1e9/BENCH_STEPS)
        vals.sort(); med = vals[3]; dev = sorted(abs(v-med) for v in vals)
        lines.append(f"{d},{BENCH_STEPS},7,{med:.4f},{dev[3]:.4f}")
    Path("out/py/s11_bench.csv").write_text("\n".join(lines) + "\n")

def s12():
    cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    s1 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>'
    for d in range(3):
        pts = " ".join(f"{60.0+500.0*math.log10(TS_MSD[i])/math.log10(4096):.1f},{310.0-280.0*math.log10(max(1.0,G_MSD[d][i]))/math.log10(4096):.1f}" for i in range(13))
        s1 += f'<polyline fill="none" stroke="{cols[d]}" stroke-width="1.6" points="{pts}"/>'
    Path("out/py/figures/msd.svg").write_text(s1 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python MSD Log-Log (y=t)</text></svg>\n')
    s2 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>'
    for d in range(3):
        pts = " ".join(f"{60.0+500.0*(i/3.0):.1f},{310.0-280.0*G_RET_EX[d][i]:.1f}" for i in range(4))
        s2 += f'<polyline fill="none" stroke="{cols[d]}" stroke-width="1.6" points="{pts}"/>' + "".join(f'<circle cx="{60.0+500.0*(i/3.0):.1f}" cy="{310.0-280.0*G_RET_SIM[d][i]:.1f}" r="4" fill="{cols[d]}"/>' for i in range(4))
    Path("out/py/figures/return_prob.svg").write_text(s2 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python Return Prob F_d(T) Exact vs Sim</text></svg>\n')
    s3 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>'
    s3 += "".join(f'<circle cx="{60.0+500.0*G_RUIN_PE[i]:.1f}" cy="{310.0-280.0*G_RUIN_P[i]:.1f}" r="5" fill="#1f77b4"/>' for i in range(6))
    Path("out/py/figures/ruin.svg").write_text(s3 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python Gambler\'s Ruin Exact vs Sim</text></svg>\n')
    Path("out/py/s12_summary.json").write_text(f'{{"lang":"py","W":{W_WALKS},"ruin_W":{W_RUIN},"bench_steps":{BENCH_STEPS},"toolchain":"python 3.11.18","total_sec":{sum(STAGE_SEC):.3f}}}\n')

STAGES = [("s01",s01),("s02",s02),("s03",s03),("s04",s04),("s05",s05),("s06",s06),("s07",s07),("s08",s08),("s09",s09),("s10",s10),("s11",s11),("s12",s12)]
cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
if cmd == "all":
    for i, (nm, fn) in enumerate(STAGES):
        print(f"[py] tahap {i+1}/12 ({nm})")
        t0 = time.perf_counter(); fn(); STAGE_SEC[i] = time.perf_counter() - t0
else:
    dict(STAGES)[cmd]()
