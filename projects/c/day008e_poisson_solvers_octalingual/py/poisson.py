from __future__ import annotations
import math, struct, sys, time
from pathlib import Path

PI_VAL = 3.14159265358979323846
NS_ALL = [15, 31, 63, 127]
MNAMES = ["jacobi", "gs", "rbgs", "sor", "cg", "pcg_ssor", "mg"]
SIN_C = [1.0, -1.0/6.0, 1.0/120.0, -1.0/5040.0, 1.0/362880.0, -1.0/39916800.0, 1.0/6227020800.0, -1.0/1307674368000.0, 1.0/355687428096000.0, -1.0/121645100408832000.0]

def sin_spec(x: float) -> float:
    sgn = 1.0
    if x < 0.0: x = -x; sgn = -sgn
    if x > PI_VAL: x -= PI_VAL; sgn = -sgn
    if x > 0.5 * PI_VAL: x = PI_VAL - x
    x2 = x * x; p = SIN_C[9]
    for k in range(8, -1, -1): p = p * x2 + SIN_C[k]
    return sgn * (x * p)

def cos_spec(x: float) -> float: return sin_spec(0.5 * PI_VAL - x)
def omega_opt(n: int) -> float: return {15: 1.673513677715992, 31: 1.8214651907890225, 63: 1.906454701582762, 127: 1.952093233850055}[n]

def fnv1a_f64(v: list[float]) -> int:
    h, mask = 0xCBF29CE484222325, 0xFFFFFFFFFFFFFFFF
    for b in struct.pack(f"<{len(v)}d", *v): h = ((h ^ b) * 0x100000001B3) & mask
    return h

def apply_a(n: int, u: list[float], out: list[float]):
    for i in range(n):
        r = i * n
        for j in range(n):
            k = r + j
            w = u[k - 1] if j > 0 else 0.0; e = u[k + 1] if j + 1 < n else 0.0
            s = u[k - n] if i > 0 else 0.0; nn = u[k + n] if i + 1 < n else 0.0
            out[k] = (((4.0 * u[k] - w) - e) - s) - nn

def dot_seq(a: list[float], b: list[float]) -> float:
    s = 0.0
    for i in range(len(a)): s += a[i] * b[i]
    return s

def norm2_seq(a: list[float]) -> float: return math.sqrt(dot_seq(a, a))
def max_abs_diff(a: list[float], b: list[float]) -> float: return max(abs(a[i] - b[i]) for i in range(len(a)))

def build_prob(n: int):
    h = 1.0 / (n + 1); h2 = h * h
    sx = [sin_spec(PI_VAL * (i + 1) * h) for i in range(n)]
    sy = [sin_spec(2.0 * PI_VAL * (i + 1) * h) for i in range(n)]
    b, ue = [0.0] * (n * n), [0.0] * (n * n)
    for i in range(n):
        x = (i + 1) * h; xx = x * (1.0 - x)
        for j in range(n):
            y = (j + 1) * h; yy = y * (1.0 - y); sxy = sx[i] * sy[j]
            ue[i * n + j] = sxy + xx * yy
            b[i * n + j] = h2 * (5.0 * PI_VAL * PI_VAL * sxy + 2.0 * xx + 2.0 * yy)
    return b, ue

def compute_res(n: int, u: list[float], b: list[float], r: list[float]) -> float:
    apply_a(n, u, r)
    for i in range(n * n): r[i] = b[i] - r[i]
    return norm2_seq(r)

def rbgs_sweep(n: int, u: list[float], b: list[float]):
    for col in (0, 1):
        for i in range(n):
            row = i * n
            for j in range(n):
                if ((i + j) & 1) != col: continue
                k = row + j
                w = u[k - 1] if j > 0 else 0.0; e = u[k + 1] if j + 1 < n else 0.0
                s = u[k - n] if i > 0 else 0.0; nn = u[k + n] if i + 1 < n else 0.0
                u[k] = 0.25 * ((((b[k] + w) + e) + s) + nn)

def restrict_fw(nf: int, rf: list[float]) -> list[float]:
    nc = (nf - 1) // 2; rc = [0.0] * (nc * nc)
    for ic in range(nc):
        fi = 2 * ic + 1
        for jc in range(nc):
            fj = 2 * jc + 1; k = fi * nf + fj
            edge = ((rf[k - 1] + rf[k + 1]) + rf[k - nf]) + rf[k + nf]
            corn = ((rf[k - nf - 1] + rf[k - nf + 1]) + rf[k + nf - 1]) + rf[k + nf + 1]
            rc[ic * nc + jc] = 4.0 * (0.25 * rf[k] + 0.125 * edge + 0.0625 * corn)
    return rc

def prolong_add(nf: int, uf: list[float], ec: list[float]):
    nc = (nf - 1) // 2
    for ic in range(nc):
        fi = 2 * ic + 1
        for jc in range(nc):
            fj = 2 * jc + 1; k = fi * nf + fj; v = ec[ic * nc + jc]
            uf[k] += v
            uf[k - 1] += 0.5 * v; uf[k + 1] += 0.5 * v; uf[k - nf] += 0.5 * v; uf[k + nf] += 0.5 * v
            uf[k - nf - 1] += 0.25 * v; uf[k - nf + 1] += 0.25 * v; uf[k + nf - 1] += 0.25 * v; uf[k + nf + 1] += 0.25 * v

def solve_dense3(b: list[float], x: list[float]):
    m = [[0.0] * 10 for _ in range(9)]
    for i in range(3):
        for j in range(3):
            r = i * 3 + j; m[r][r] = 4.0
            if j > 0: m[r][r - 1] = -1.0
            if j < 2: m[r][r + 1] = -1.0
            if i > 0: m[r][r - 3] = -1.0
            if i < 2: m[r][r + 3] = -1.0
            m[r][9] = b[r]
    for c in range(9):
        piv = max(range(c, 9), key=lambda r: abs(m[r][c]))
        if piv != c: m[c], m[piv] = m[piv], m[c]
        for r in range(c + 1, 9):
            f = m[r][c] / m[c][c]
            for k in range(c, 10): m[r][k] -= f * m[c][k]
    for r in range(8, -1, -1):
        s = m[r][9]
        for k in range(r + 1, 9): s -= m[r][k] * x[k]
        x[r] = s / m[r][r]

def vcycle(n: int, u: list[float], b: list[float]):
    if n == 3: solve_dense3(b, u); return
    rbgs_sweep(n, u, b); rbgs_sweep(n, u, b)
    r = [0.0] * (n * n); compute_res(n, u, b, r)
    rc = restrict_fw(n, r); nc = (n - 1) // 2; ec = [0.0] * (nc * nc)
    vcycle(nc, ec, rc); prolong_add(n, u, ec); rbgs_sweep(n, u, b)

def fmg_solve(n: int, u: list[float], b: list[float]):
    if n == 3: solve_dense3(b, u); return
    bc = restrict_fw(n, b); nc = (n - 1) // 2; uc = [0.0] * (nc * nc)
    fmg_solve(nc, uc, bc)
    for i in range(n * n): u[i] = 0.0
    prolong_add(n, u, uc); vcycle(n, u, b)

Path("out/py/figures").mkdir(parents=True, exist_ok=True)
G_RES = [[None]*4 for _ in range(7)]
G_HIST63 = [[] for _ in range(7)]
G_DISC_ERR = [0.0]*4
STAGE_S = [0.0]*19

def run_solver(meth: int, n: int, save63: bool):
    n2 = n * n; b, ue = build_prob(n)
    u, r, p, ap, z = [0.0]*n2, [0.0]*n2, [0.0]*n2, [0.0]*n2, [0.0]*n2
    nb = norm2_seq(b); om = omega_opt(n)
    max_it = [100000, 50000, 50000, 10000, 5000, 3000, 100][meth]
    hist, rr, rz = [], 0.0, 0.0
    if meth in (4, 5): r = b[:]; rr = dot_seq(r, r)
    if meth == 4: p = r[:]
    if meth == 5:
        for i in range(n):
            for j in range(n): k = i * n + j; z[k] = 0.25 * ((r[k] + (z[k - 1] if j > 0 else 0.0)) + (z[k - n] if i > 0 else 0.0))
        for i in range(n - 1, -1, -1):
            for j in range(n - 1, -1, -1): k = i * n + j; p[k] = z[k] + 0.25 * ((p[k + 1] if j + 1 < n else 0.0) + (p[k + n] if i + 1 < n else 0.0))
        z = p[:]; rz = dot_seq(r, z)
    for _ in range(max_it):
        if meth == 0:
            for i in range(n):
                row = i * n
                for j in range(n):
                    k = row + j; w = u[k - 1] if j > 0 else 0.0; e = u[k + 1] if j + 1 < n else 0.0
                    s = u[k - n] if i > 0 else 0.0; nn = u[k + n] if i + 1 < n else 0.0
                    z[k] = 0.25 * ((((b[k] + w) + e) + s) + nn)
            u, z = z, u; rel = compute_res(n, u, b, r) / nb
        elif meth in (1, 3):
            for i in range(n):
                row = i * n
                for j in range(n):
                    k = row + j; w = u[k - 1] if j > 0 else 0.0; e = u[k + 1] if j + 1 < n else 0.0
                    s = u[k - n] if i > 0 else 0.0; nn = u[k + n] if i + 1 < n else 0.0
                    gs = 0.25 * ((((b[k] + w) + e) + s) + nn)
                    u[k] = gs if meth == 1 else ((1.0 - om) * u[k] + om * gs)
            rel = compute_res(n, u, b, r) / nb
        elif meth == 2:
            rbgs_sweep(n, u, b); rel = compute_res(n, u, b, r) / nb
        elif meth == 4:
            apply_a(n, p, ap); alpha = rr / dot_seq(p, ap)
            for k in range(n2): u[k] += alpha * p[k]; r[k] -= alpha * ap[k]
            rr_new = dot_seq(r, r); rel = math.sqrt(rr_new) / nb; hist.append(rel)
            if rel <= 1e-8: break
            beta = rr_new / rr
            for k in range(n2): p[k] = r[k] + beta * p[k]
            rr = rr_new; continue
        elif meth == 5:
            apply_a(n, p, ap); alpha = rz / dot_seq(p, ap)
            for k in range(n2): u[k] += alpha * p[k]; r[k] -= alpha * ap[k]
            rel = norm2_seq(r) / nb; hist.append(rel)
            if rel <= 1e-8: break
            for i in range(n):
                for j in range(n): k = i * n + j; ap[k] = 0.25 * ((r[k] + (ap[k - 1] if j > 0 else 0.0)) + (ap[k - n] if i > 0 else 0.0))
            for i in range(n - 1, -1, -1):
                for j in range(n - 1, -1, -1): k = i * n + j; z[k] = ap[k] + 0.25 * ((z[k + 1] if j + 1 < n else 0.0) + (z[k + n] if i + 1 < n else 0.0))
            rz_new = dot_seq(r, z); beta = rz_new / rz
            for k in range(n2): p[k] = z[k] + beta * p[k]
            rz = rz_new; continue
        else:
            vcycle(n, u, b); rel = compute_res(n, u, b, r) / nb
        hist.append(rel)
        if rel <= 1e-8: break
    used = len(hist)
    rho = (hist[-1] / hist[-51]) ** 0.02 if used >= 51 else ((hist[-1] / hist[0]) ** (1.0 / (used - 1)) if used >= 2 else 0.0)
    if save63 and n == 63: G_HIST63[meth] = hist
    return {"iters": used, "rel_res": hist[-1], "err_inf": max_abs_diff(u, ue), "rho_emp": rho, "hist_hash": fnv1a_f64(hist), "sol_hash": fnv1a_f64(u)}

def s01():
    u = [float(k) for k in range(1, 10)]; au = [0.0]*9; apply_a(3, u, au)
    h1, h2, h3 = fnv1a_f64(u), fnv1a_f64(au), fnv1a_f64([1.0, -2.5, 3.25])
    ok = max_abs_diff(au, [-2,-1,4,3,0,7,16,11,22]) == 0.0 and h1 == 0xc56de10c6671e6c2 and h2 == 0x59ed72ec82ee959c and h3 == 0xafc5dc1ca48a2366
    Path("out/py/s01_selftest.txt").write_text(f"ok,h_u3,h_Au3,h_test\n{int(ok)},{h1:016x},{h2:016x},{h3:016x}\n")
    if not ok: sys.exit(1)

def s02():
    lines = ["N,norm_b,h_b,h_ue"]
    for n in NS_ALL: b, ue = build_prob(n); lines.append(f"{n},{norm2_seq(b):.17g},{fnv1a_f64(b):016x},{fnv1a_f64(ue):016x}")
    Path("out/py/s02_problem.csv").write_text("\n".join(lines) + "\n")

def s03():
    lines = ["N,j,k,lambda_exact,rel_err_inf"]
    for n in NS_ALL[:3]:
        h = 1.0 / (n + 1); v, av = [0.0]*(n*n), [0.0]*(n*n)
        for j, k in ((1,1), (1,2), (3,5)):
            lam = 4.0 - 2.0 * cos_spec(j * PI_VAL * h) - 2.0 * cos_spec(k * PI_VAL * h)
            for r in range(n):
                for c in range(n): v[r*n+c] = sin_spec(j * PI_VAL * (r+1) * h) * sin_spec(k * PI_VAL * (c+1) * h)
            apply_a(n, v, av); me = max(abs(av[q] - lam * v[q]) for q in range(n*n))
            lines.append(f"{n},{j},{k},{lam:.17g},{me/lam:.3e}")
    Path("out/py/s03_eigen.csv").write_text("\n".join(lines) + "\n")

def run_m(m: int, fn: str, nmax: int):
    lines = ["N,iters,rel_res,err_inf,rho_emp"]
    for i in range(nmax): r = run_solver(m, NS_ALL[i], NS_ALL[i] == 63); G_RES[m][i] = r; lines.append(f"{NS_ALL[i]},{r['iters']},{r['rel_res']:.10e},{r['err_inf']:.10e},{r['rho_emp']:.10f}")
    Path(fn).write_text("\n".join(lines) + "\n")

def s04(): run_m(0, "out/py/s04_jacobi.csv", 3)
def s05(): run_m(1, "out/py/s05_gs.csv", 3)
def s06(): run_m(2, "out/py/s06_rbgs.csv", 3)
def s07(): run_m(3, "out/py/s07_sor.csv", 4)
def s08(): run_m(4, "out/py/s08_cg.csv", 4)
def s09(): run_m(5, "out/py/s09_pcg_ssor.csv", 4)
def s10(): run_m(6, "out/py/s10_mg.csv", 4)
def s11():
    lines = ["N,fmg_err_inf,disc_err_inf,ratio"]
    for i, n in enumerate(NS_ALL):
        b, ue = build_prob(n); u = [0.0]*(n*n); fmg_solve(n, u, b); fe = max_abs_diff(u, ue); G_DISC_ERR[i] = G_RES[6][i]["err_inf"]
        lines.append(f"{n},{fe:.10e},{G_DISC_ERR[i]:.10e},{fe/G_DISC_ERR[i]:.6f}")
    Path("out/py/s11_fmg.csv").write_text("\n".join(lines) + "\n")
def s12():
    xs = [math.log(1.0/(n+1)) for n in NS_ALL]; ys = [math.log(e) for e in G_DISC_ERR]
    mx, my = sum(xs)/4.0, sum(ys)/4.0; sl = sum((xs[i]-mx)*(ys[i]-my) for i in range(4))/sum((x-mx)**2 for x in xs)
    Path("out/py/s12_order.csv").write_text("\n".join(["N,h,err_inf,slope"] + [f"{NS_ALL[i]},{1.0/(NS_ALL[i]+1):.10e},{G_DISC_ERR[i]:.10e},{sl:.6f}" for i in range(4)]) + "\n")
def s13():
    lines = ["method,N,iters,slope,slope_se"]
    for m in range(7):
        cnt = 3 if m <= 2 else 4
        xs = [math.log(NS_ALL[i]) for i in range(cnt)]; ys = [math.log(G_RES[m][i]["iters"]) for i in range(cnt)]
        mx, my = sum(xs)/cnt, sum(ys)/cnt; den = sum((x-mx)**2 for x in xs)
        sl = sum((xs[i]-mx)*(ys[i]-my) for i in range(cnt))/den; ic = my - sl*mx
        rss = sum((ys[i]-(sl*xs[i]+ic))**2 for i in range(cnt)); se = math.sqrt((rss/(cnt-2))/den)
        for i in range(cnt): lines.append(f"{MNAMES[m]},{NS_ALL[i]},{G_RES[m][i]['iters']},{sl:.6f},{se:.6f}")
    Path("out/py/s13_scaling.csv").write_text("\n".join(lines) + "\n")
def s14():
    lines = ["method,iter,rel_res"]
    for m in range(7):
        ln = len(G_HIST63[m]); step = ln // 200 if ln > 200 else 1
        for k in range(0, ln, step): lines.append(f"{MNAMES[m]},{k+1},{G_HIST63[m][k]:.10e}")
        if (ln - 1) % step != 0: lines.append(f"{MNAMES[m]},{ln},{G_HIST63[m][-1]:.10e}")
    Path("out/py/s14_history.csv").write_text("\n".join(lines) + "\n")
def s15():
    lines = ["N,iters,mg_factor_last5"]
    for n in (63, 127):
        b, _ = build_prob(n); u, r = [0.0]*(n*n), [0.0]*(n*n); nb = norm2_seq(b); h = []
        for _ in range(20): vcycle(n, u, b); rel = compute_res(n, u, b, r)/nb; h.append(rel); (lambda: None)() if rel > 1e-8 else None
        for idx, val in enumerate(h):
            if val <= 1e-8: h = h[:idx+1]; break
        lines.append(f"{n},{len(h)},{(h[-1]/h[-6])**0.2:.8f}")
    Path("out/py/s15_mg_factor.csv").write_text("\n".join(lines) + "\n")
def s16():
    flops = [9,9,9,12,21,35,36]; lines = ["method,N,flops_per_unk_iter,ns_med,ns_mad"]
    n = 63; n2 = n*n; b, _ = build_prob(n); u, z = [0.0]*n2, [0.0]*n2
    for m in range(7):
        iters = 3 if m == 6 else 10; vals = []
        for _ in range(7):
            u = [0.0]*n2; t0 = time.perf_counter()
            for _ in range(iters):
                if m == 6: vcycle(n, u, b)
                elif m == 2: rbgs_sweep(n, u, b)
                else:
                    apply_a(n, u, z)
                    for k in range(n2): u[k] += 0.25 * (b[k] - z[k])
            vals.append((time.perf_counter()-t0)*1e9/(iters*n2))
        vals.sort(); med = vals[3]; dev = sorted(abs(v-med) for v in vals)
        lines.append(f"{MNAMES[m]},{n},{flops[m]},{med:.4f},{dev[3]:.4f}")
    Path("out/py/s16_cost.csv").write_text("\n".join(lines) + "\n")
def s17():
    lines = ["method,N,iters,hist_fnv1a,sol_fnv1a"]
    for m in range(7):
        for ni in (0, 1): r = G_RES[m][ni]; lines.append(f"{MNAMES[m]},{NS_ALL[ni]},{r['iters']},{r['hist_hash']:016x},{r['sol_hash']:016x}")
    Path("out/py/s17_parity.csv").write_text("\n".join(lines) + "\n")
def s18():
    cols = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2"]
    s1 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>'
    for m in range(7):
        ln = len(G_HIST63[m]); step = ln // 150 if ln > 150 else 1
        pts = " ".join(f"{60.0+500.0*math.log10(k+1)/math.log10(len(G_HIST63[0])):.1f},{310.0-280.0*min(1.0,max(0.0,-math.log10(max(1e-12,G_HIST63[m][k]))/9.0)):.1f}" for k in range(0, ln, step))
        s1 += f'<polyline fill="none" stroke="{cols[m]}" stroke-width="1.4" points="{pts}"/>'
    Path("out/py/figures/history.svg").write_text(s1 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python Residual History (N=63)</text></svg>\n')
    s2 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>'
    for m in range(7):
        cnt = 3 if m <= 2 else 4
        pts = " ".join(f"{60.0+500.0*(math.log10(NS_ALL[i])-math.log10(15))/(math.log10(127)-math.log10(15)):.1f},{310.0-280.0*math.log10(G_RES[m][i]['iters'])/5.0:.1f}" for i in range(cnt))
        s2 += f'<polyline fill="none" stroke="{cols[m]}" stroke-width="1.5" points="{pts}"/>'
    Path("out/py/figures/scaling.svg").write_text(s2 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python Iterations vs N</text></svg>\n')
    pts3 = " ".join(f"{60.0+500.0*(i/3.0):.1f},{310.0-280.0*(math.log10(G_DISC_ERR[i])+5.0)/4.0:.1f}" for i in range(4))
    Path("out/py/figures/order.svg").write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><polyline fill="none" stroke="#1f77b4" stroke-width="1.8" points="{pts3}"/><text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">Python Discretization Error vs h</text></svg>\n')
def s19():
    Path("out/py/s19_summary.json").write_text(f'{{"lang":"py","N_max":127,"N_stat_max":63,"toolchain":"python 3.11.18","total_sec":{sum(STAGE_S):.3f}}}\n')

STAGES = [(f"s{i:02d}", globals()[f"s{i:02d}"]) for i in range(1, 20)]
cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
if cmd == "all":
    for i, (nm, fn) in enumerate(STAGES):
        print(f"[py] tahap {i+1}/19 ({nm})"); t0 = time.perf_counter(); fn(); STAGE_S[i] = time.perf_counter() - t0
else: dict(STAGES)[cmd]()
