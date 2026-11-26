import { mkdirSync, writeFileSync } from "node:fs";

const PI_VAL = 3.14159265358979323846;
const NS_ALL = [15, 31, 63, 127, 255];
const MNAMES = ["jacobi", "gs", "rbgs", "sor", "cg", "pcg_ssor", "mg"];
const OMEGA: Record<number, number> = {
  15: 1.673513677715992, 31: 1.8214651907890225, 63: 1.906454701582762,
  127: 1.952093233850055, 255: 1.9757544535797149,
};
const C_SIN = [
  1.0, -1.0 / 6.0, 1.0 / 120.0, -1.0 / 5040.0, 1.0 / 362880.0, -1.0 / 39916800.0,
  1.0 / 6227020800.0, -1.0 / 1307674368000.0, 1.0 / 355687428096000.0, -1.0 / 121645100408832000.0,
];

function sinSpec(x: number): number {
  let sgn = 1.0;
  if (x < 0) { x = -x; sgn = -sgn; }
  if (x > PI_VAL) { x = x - PI_VAL; sgn = -sgn; }
  if (x > 0.5 * PI_VAL) x = PI_VAL - x;
  const x2 = x * x;
  let p = C_SIN[9];
  for (let k = 8; k >= 0; k--) p = p * x2 + C_SIN[k];
  return sgn * (x * p);
}
const cosSpec = (x: number): number => sinSpec(0.5 * PI_VAL - x);

function fnv1aF64(v: Float64Array | number[]): string {
  const arr = v instanceof Float64Array ? v : Float64Array.from(v);
  const dv = new DataView(arr.buffer, arr.byteOffset, arr.byteLength);
  let h = 0xcbf29ce484222325n;
  const p = 0x100000001b3n;
  for (let i = 0; i < arr.length; i++) {
    let u = dv.getBigUint64(i * 8, true);
    for (let b = 0; b < 8; b++) {
      h = BigInt.asUintN(64, (h ^ (u & 0xffn)) * p);
      u >>= 8n;
    }
  }
  return h.toString(16).padStart(16, "0");
}

function applyA(N: number, u: Float64Array, out: Float64Array): void {
  for (let i = 0; i < N; i++) {
    const r = i * N;
    for (let j = 0; j < N; j++) {
      const k = r + j;
      const w = j > 0 ? u[k - 1] : 0.0, e = j + 1 < N ? u[k + 1] : 0.0;
      const s = i > 0 ? u[k - N] : 0.0, n = i + 1 < N ? u[k + N] : 0.0;
      out[k] = (((4.0 * u[k] - w) - e) - s) - n;
    }
  }
}

function dotSeq(a: Float64Array, b: Float64Array): number {
  let s = 0.0;
  for (let i = 0; i < a.length; i++) s += a[i] * b[i];
  return s;
}
const norm2Seq = (a: Float64Array): number => Math.sqrt(dotSeq(a, a));
function maxAbsDiff(a: Float64Array, b: Float64Array): number {
  let m = 0.0;
  for (let i = 0; i < a.length; i++) {
    const d = Math.abs(a[i] - b[i]);
    if (d > m) m = d;
  }
  return m;
}

function buildProb(N: number): [Float64Array, Float64Array] {
  const h = 1.0 / (N + 1), h2 = h * h;
  const sx = new Float64Array(N), sy = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    sx[i] = sinSpec(PI_VAL * (i + 1) * h);
    sy[i] = sinSpec(2.0 * PI_VAL * (i + 1) * h);
  }
  const b = new Float64Array(N * N), ue = new Float64Array(N * N);
  const c5 = 5.0 * PI_VAL * PI_VAL;
  for (let i = 0; i < N; i++) {
    const x = (i + 1) * h, xx = x * (1.0 - x);
    for (let j = 0; j < N; j++) {
      const y = (j + 1) * h, yy = y * (1.0 - y), sxy = sx[i] * sy[j];
      ue[i * N + j] = sxy + xx * yy;
      b[i * N + j] = h2 * ((c5 * sxy + 2.0 * xx) + 2.0 * yy);
    }
  }
  return [b, ue];
}

function computeRes(N: number, u: Float64Array, b: Float64Array, r: Float64Array): number {
  applyA(N, u, r);
  for (let i = 0; i < r.length; i++) r[i] = b[i] - r[i];
  return norm2Seq(r);
}

function rbgsSweep(N: number, u: Float64Array, b: Float64Array): void {
  for (let col = 0; col < 2; col++) {
    for (let i = 0; i < N; i++) {
      const row = i * N;
      for (let j = 0; j < N; j++) {
        if (((i + j) & 1) !== col) continue;
        const k = row + j;
        const w = j > 0 ? u[k - 1] : 0.0, e = j + 1 < N ? u[k + 1] : 0.0;
        const s = i > 0 ? u[k - N] : 0.0, n = i + 1 < N ? u[k + N] : 0.0;
        u[k] = 0.25 * ((((b[k] + w) + e) + s) + n);
      }
    }
  }
}

function ssorPre(N: number, r: Float64Array, y: Float64Array, z: Float64Array): void {
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      const k = i * N + j;
      const yw = j > 0 ? y[k - 1] : 0.0, ys = i > 0 ? y[k - N] : 0.0;
      y[k] = 0.25 * ((r[k] + yw) + ys);
    }
  }
  for (let i = N - 1; i >= 0; i--) {
    for (let j = N - 1; j >= 0; j--) {
      const k = i * N + j;
      const ze = j + 1 < N ? z[k + 1] : 0.0, zn = i + 1 < N ? z[k + N] : 0.0;
      z[k] = y[k] + 0.25 * (ze + zn);
    }
  }
}

function restrictFW(Nf: number, rf: Float64Array): Float64Array {
  const Nc = (Nf - 1) >> 1;
  const rc = new Float64Array(Nc * Nc);
  for (let ic = 0; ic < Nc; ic++) {
    const fi = 2 * ic + 1;
    for (let jc = 0; jc < Nc; jc++) {
      const fj = 2 * jc + 1, k = fi * Nf + fj;
      const edge = ((rf[k - 1] + rf[k + 1]) + rf[k - Nf]) + rf[k + Nf];
      const corn = ((rf[k - Nf - 1] + rf[k - Nf + 1]) + rf[k + Nf - 1]) + rf[k + Nf + 1];
      rc[ic * Nc + jc] = 4.0 * ((0.25 * rf[k] + 0.125 * edge) + 0.0625 * corn);
    }
  }
  return rc;
}

function prolongAdd(Nf: number, uf: Float64Array, ec: Float64Array): void {
  const Nc = (Nf - 1) >> 1;
  for (let ic = 0; ic < Nc; ic++) {
    const fi = 2 * ic + 1;
    for (let jc = 0; jc < Nc; jc++) {
      const fj = 2 * jc + 1, k = fi * Nf + fj, v = ec[ic * Nc + jc];
      uf[k] += v;
      uf[k - 1] += 0.5 * v; uf[k + 1] += 0.5 * v; uf[k - Nf] += 0.5 * v; uf[k + Nf] += 0.5 * v;
      uf[k - Nf - 1] += 0.25 * v; uf[k - Nf + 1] += 0.25 * v; uf[k + Nf - 1] += 0.25 * v; uf[k + Nf + 1] += 0.25 * v;
    }
  }
}

function solveDense3(b: Float64Array, x: Float64Array): void {
  const M = Array.from({ length: 9 }, () => new Float64Array(10));
  for (let i = 0; i < 3; i++) {
    for (let j = 0; j < 3; j++) {
      const r = i * 3 + j;
      M[r][r] = 4.0; M[r][9] = b[r];
      if (j > 0) M[r][r - 1] = -1;
      if (j < 2) M[r][r + 1] = -1;
      if (i > 0) M[r][r - 3] = -1;
      if (i < 2) M[r][r + 3] = -1;
    }
  }
  for (let c = 0; c < 9; c++) {
    let piv = c;
    for (let r = c + 1; r < 9; r++) if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) piv = r;
    if (piv !== c) [M[c], M[piv]] = [M[piv], M[c]];
    for (let r = c + 1; r < 9; r++) {
      const f = M[r][c] / M[c][c];
      for (let k = c; k < 10; k++) M[r][k] -= f * M[c][k];
    }
  }
  for (let r = 8; r >= 0; r--) {
    let s = M[r][9];
    for (let k = r + 1; k < 9; k++) s -= M[r][k] * x[k];
    x[r] = s / M[r][r];
  }
}

function vcycle(N: number, u: Float64Array, b: Float64Array): void {
  if (N === 3) { solveDense3(b, u); return; }
  rbgsSweep(N, u, b);
  rbgsSweep(N, u, b);
  const r = new Float64Array(N * N);
  computeRes(N, u, b, r);
  const rc = restrictFW(N, r);
  const Nc = (N - 1) >> 1;
  const ec = new Float64Array(Nc * Nc);
  vcycle(Nc, ec, rc);
  prolongAdd(N, u, ec);
  rbgsSweep(N, u, b);
}

function fmgSolve(N: number, u: Float64Array, b: Float64Array): void {
  if (N === 3) { solveDense3(b, u); return; }
  const bc = restrictFW(N, b);
  const Nc = (N - 1) >> 1;
  const uc = new Float64Array(Nc * Nc);
  fmgSolve(Nc, uc, bc);
  u.fill(0);
  prolongAdd(N, u, uc);
  vcycle(N, u, b);
}

interface SolveRes {
  iters: number; relRes: number; errInf: number; rho: number;
  histHash: string; solHash: string;
}

const gHist63: number[][] = Array.from({ length: 7 }, () => []);
const gRes: SolveRes[][] = Array.from({ length: 7 }, () => []);
const gDiscErr = new Float64Array(5);
const gStageS = new Float64Array(19);

function runSolver(meth: number, N: number, save63: boolean): SolveRes {
  const n2 = N * N;
  const [b, ue] = buildProb(N);
  let u = new Float64Array(n2), r = new Float64Array(n2), p = new Float64Array(n2);
  const ap = new Float64Array(n2);
  let z = new Float64Array(n2);
  const nb = norm2Seq(b), om = OMEGA[N];
  const maxIt = [100000, 50000, 50000, 10000, 5000, 3000, 100][meth];
  const hist: number[] = [];
  let rr = 0.0, rz = 0.0;
  if (meth === 4 || meth === 5) { r.set(b); rr = dotSeq(r, r); }
  if (meth === 4) p.set(r);
  if (meth === 5) { ssorPre(N, r, ap, z); p.set(z); rz = dotSeq(r, z); }
  for (let it = 1; it <= maxIt; it++) {
    let rel = 0.0;
    if (meth === 0) {
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
          const k = i * N + j;
          const w = j > 0 ? u[k - 1] : 0.0, e = j + 1 < N ? u[k + 1] : 0.0;
          const s = i > 0 ? u[k - N] : 0.0, n = i + 1 < N ? u[k + N] : 0.0;
          z[k] = 0.25 * ((((b[k] + w) + e) + s) + n);
        }
      }
      [u, z] = [z, u];
      rel = computeRes(N, u, b, r) / nb;
    } else if (meth === 1 || meth === 3) {
      for (let i = 0; i < N; i++) {
        for (let j = 0; j < N; j++) {
          const k = i * N + j;
          const w = j > 0 ? u[k - 1] : 0.0, e = j + 1 < N ? u[k + 1] : 0.0;
          const s = i > 0 ? u[k - N] : 0.0, n = i + 1 < N ? u[k + N] : 0.0;
          const gs = 0.25 * ((((b[k] + w) + e) + s) + n);
          u[k] = meth === 1 ? gs : (1.0 - om) * u[k] + om * gs;
        }
      }
      rel = computeRes(N, u, b, r) / nb;
    } else if (meth === 2) {
      rbgsSweep(N, u, b);
      rel = computeRes(N, u, b, r) / nb;
    } else if (meth === 4 || meth === 5) {
      applyA(N, p, ap);
      const alpha = (meth === 4 ? rr : rz) / dotSeq(p, ap);
      for (let k = 0; k < n2; k++) { u[k] += alpha * p[k]; r[k] -= alpha * ap[k]; }
      if (meth === 4) {
        const rrNew = dotSeq(r, r);
        rel = Math.sqrt(rrNew) / nb;
        hist.push(rel);
        if (rel <= 1e-8) break;
        const beta = rrNew / rr;
        for (let k = 0; k < n2; k++) p[k] = r[k] + beta * p[k];
        rr = rrNew;
      } else {
        rel = norm2Seq(r) / nb;
        hist.push(rel);
        if (rel <= 1e-8) break;
        ssorPre(N, r, ap, z);
        const rzNew = dotSeq(r, z);
        const beta = rzNew / rz;
        for (let k = 0; k < n2; k++) p[k] = z[k] + beta * p[k];
        rz = rzNew;
      }
      continue;
    } else {
      vcycle(N, u, b);
      rel = computeRes(N, u, b, r) / nb;
    }
    hist.push(rel);
    if (rel <= 1e-8) break;
  }
  const used = hist.length;
  let rho = 0.0;
  if (used >= 51) rho = Math.pow(hist[used - 1] / hist[used - 51], 0.02);
  else if (used >= 2) rho = Math.pow(hist[used - 1] / hist[0], 1.0 / (used - 1));
  if (save63 && N === 63) gHist63[meth] = hist;
  return { iters: used, relRes: hist[used - 1], errInf: maxAbsDiff(u, ue), rho, histHash: fnv1aF64(hist), solHash: fnv1aF64(u) };
}

function s01(): void {
  const u = Float64Array.from([1, 2, 3, 4, 5, 6, 7, 8, 9]), au = new Float64Array(9);
  applyA(3, u, au);
  const h1 = fnv1aF64(u), h2 = fnv1aF64(au), h3 = fnv1aF64([1.0, -2.5, 3.25]);
  const ok = maxAbsDiff(au, Float64Array.from([-2, -1, 4, 3, 0, 7, 16, 11, 22])) === 0 &&
    h1 === "c56de10c6671e6c2" && h2 === "59ed72ec82ee959c" && h3 === "afc5dc1ca48a2366" ? 1 : 0;
  writeFileSync("out/ts/s01_selftest.txt", `ok,h_u3,h_Au3,h_test\n${ok},${h1},${h2},${h3}\n`);
  if (!ok) process.exit(1);
}
function s02(): void {
  const lines = ["N,norm_b,h_b,h_ue"];
  for (const N of NS_ALL) {
    const [b, ue] = buildProb(N);
    lines.push(`${N},${norm2Seq(b).toPrecision(17)},${fnv1aF64(b)},${fnv1aF64(ue)}`);
  }
  writeFileSync("out/ts/s02_problem.csv", lines.join("\n") + "\n");
}
function s03(): void {
  const lines = ["N,j,k,lambda_exact,rel_err_inf"];
  const modes = [[1, 1], [1, 2], [3, 5]];
  for (const N of NS_ALL.slice(0, 3)) {
    const h = 1.0 / (N + 1), v = new Float64Array(N * N), av = new Float64Array(N * N);
    for (const [mj, mk] of modes) {
      const lam = 4.0 - 2.0 * cosSpec(mj * PI_VAL * h) - 2.0 * cosSpec(mk * PI_VAL * h);
      for (let r = 0; r < N; r++) for (let c = 0; c < N; c++)
        v[r * N + c] = sinSpec(mj * PI_VAL * (r + 1) * h) * sinSpec(mk * PI_VAL * (c + 1) * h);
      applyA(N, v, av);
      let me = 0.0;
      for (let q = 0; q < N * N; q++) { const d = Math.abs(av[q] - lam * v[q]); if (d > me) me = d; }
      lines.push(`${N},${mj},${mk},${lam.toPrecision(17)},${(me / lam).toExponential(3)}`);
    }
  }
  writeFileSync("out/ts/s03_eigen.csv", lines.join("\n") + "\n");
}
function runM(m: number, fn: string, nmax: number): void {
  const lines = ["N,iters,rel_res,err_inf,rho_emp"];
  for (let i = 0; i < nmax; i++) {
    const r = runSolver(m, NS_ALL[i], NS_ALL[i] === 63);
    gRes[m][i] = r;
    lines.push(`${NS_ALL[i]},${r.iters},${r.relRes.toExponential(10)},${r.errInf.toExponential(10)},${r.rho.toFixed(10)}`);
  }
  writeFileSync(fn, lines.join("\n") + "\n");
}
function s11(): void {
  const lines = ["N,fmg_err_inf,disc_err_inf,ratio"];
  for (let i = 0; i < NS_ALL.length; i++) {
    const N = NS_ALL[i], [b, ue] = buildProb(N), u = new Float64Array(N * N);
    fmgSolve(N, u, b);
    const fe = maxAbsDiff(u, ue);
    gDiscErr[i] = gRes[6][i].errInf;
    lines.push(`${N},${fe.toExponential(10)},${gDiscErr[i].toExponential(10)},${(fe / gDiscErr[i]).toFixed(6)}`);
  }
  writeFileSync("out/ts/s11_fmg.csv", lines.join("\n") + "\n");
}
function s12(): void {
  const xs = NS_ALL.map((N) => Math.log(1.0 / (N + 1))), ys = Array.from(gDiscErr, (e) => Math.log(e));
  const mx = xs.reduce((a, b) => a + b, 0) / xs.length, my = ys.reduce((a, b) => a + b, 0) / ys.length;
  let num = 0.0, den = 0.0;
  for (let i = 0; i < xs.length; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) ** 2; }
  const sl = num / den;
  const lines = ["N,h,err_inf,slope"];
  for (let i = 0; i < NS_ALL.length; i++)
    lines.push(`${NS_ALL[i]},${(1.0 / (NS_ALL[i] + 1)).toExponential(10)},${gDiscErr[i].toExponential(10)},${sl.toFixed(6)}`);
  writeFileSync("out/ts/s12_order.csv", lines.join("\n") + "\n");
}
function s13(): void {
  const lines = ["method,N,iters,slope,slope_se"];
  for (let m = 0; m < 7; m++) {
    const cnt = m <= 2 ? 3 : 5;
    const xs = Array.from({ length: cnt }, (_, i) => Math.log(NS_ALL[i]));
    const ys = Array.from({ length: cnt }, (_, i) => Math.log(gRes[m][i].iters));
    const mx = xs.reduce((a, b) => a + b, 0) / cnt, my = ys.reduce((a, b) => a + b, 0) / cnt;
    let num = 0.0, den = 0.0, rss = 0.0;
    for (let i = 0; i < cnt; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) ** 2; }
    const sl = num / den, ic = my - sl * mx;
    for (let i = 0; i < cnt; i++) rss += (ys[i] - (sl * xs[i] + ic)) ** 2;
    const se = Math.sqrt((rss / (cnt - 2)) / den);
    for (let i = 0; i < cnt; i++) lines.push(`${MNAMES[m]},${NS_ALL[i]},${gRes[m][i].iters},${sl.toFixed(6)},${se.toFixed(6)}`);
  }
  writeFileSync("out/ts/s13_scaling.csv", lines.join("\n") + "\n");
}
function s14(): void {
  const lines = ["method,iter,rel_res"];
  for (let m = 0; m < 7; m++) {
    const h = gHist63[m], step = Math.max(1, Math.floor(h.length / 200));
    for (let k = 0; k < h.length; k += step) lines.push(`${MNAMES[m]},${k + 1},${h[k].toExponential(10)}`);
    if ((h.length - 1) % step !== 0) lines.push(`${MNAMES[m]},${h.length},${h[h.length - 1].toExponential(10)}`);
  }
  writeFileSync("out/ts/s14_history.csv", lines.join("\n") + "\n");
}
function s15(): void {
  const lines = ["N,iters,mg_factor_last5"];
  for (const N of [63, 127]) {
    const [b] = buildProb(N), u = new Float64Array(N * N), r = new Float64Array(N * N), nb = norm2Seq(b);
    const h: number[] = [];
    for (let it = 1; it <= 20; it++) {
      vcycle(N, u, b);
      const rel = computeRes(N, u, b, r) / nb;
      h.push(rel);
      if (rel <= 1e-8) break;
    }
    const it = h.length;
    lines.push(`${N},${it},${Math.pow(h[it - 1] / h[it - 6], 0.2).toFixed(8)}`);
  }
  writeFileSync("out/ts/s15_mg_factor.csv", lines.join("\n") + "\n");
}
function s16(): void {
  const flops = [9, 9, 9, 12, 21, 35, 36];
  const lines = ["method,N,flops_per_unk_iter,ns_med,ns_mad"];
  const N = 63, n2 = N * N, [b] = buildProb(N), u = new Float64Array(n2), z = new Float64Array(n2);
  for (let m = 0; m < 7; m++) {
    const iters = m === 6 ? 10 : 80, vals: number[] = [];
    for (let rp = 0; rp < 7; rp++) {
      u.fill(0);
      const t0 = process.hrtime.bigint();
      for (let it = 0; it < iters; it++) {
        if (m === 6) vcycle(N, u, b);
        else if (m === 2) rbgsSweep(N, u, b);
        else { applyA(N, u, z); for (let k = 0; k < n2; k++) u[k] += 0.25 * (b[k] - z[k]); }
      }
      vals.push(Number(process.hrtime.bigint() - t0) / (iters * n2));
    }
    vals.sort((a, b) => a - b);
    const med = vals[3], dev = vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b);
    lines.push(`${MNAMES[m]},${N},${flops[m]},${med.toFixed(4)},${dev[3].toFixed(4)}`);
  }
  writeFileSync("out/ts/s16_cost.csv", lines.join("\n") + "\n");
}
function s17(): void {
  const lines = ["method,N,iters,hist_fnv1a,sol_fnv1a"];
  for (let m = 0; m < 7; m++) for (let ni = 0; ni < 2; ni++) {
    const r = gRes[m][ni];
    lines.push(`${MNAMES[m]},${NS_ALL[ni]},${r.iters},${r.histHash},${r.solHash}`);
  }
  writeFileSync("out/ts/s17_parity.csv", lines.join("\n") + "\n");
}
function s18(): void {
  const cols = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2"];
  const hdr = `<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>`;
  let s1 = hdr;
  for (let m = 0; m < 7; m++) {
    const h = gHist63[m], step = Math.max(1, Math.floor(h.length / 150)), pts: string[] = [];
    for (let k = 0; k < h.length; k += step) {
      const px = 60.0 + 500.0 * Math.log10(k + 1) / Math.log10(gHist63[0].length);
      const py = 310.0 - 280.0 * Math.min(1, Math.max(0, -Math.log10(Math.max(1e-12, h[k])) / 9.0));
      pts.push(`${px.toFixed(1)},${py.toFixed(1)}`);
    }
    s1 += `<polyline fill="none" stroke="${cols[m]}" stroke-width="1.4" points="${pts.join(" ")}"/>`;
  }
  writeFileSync("out/ts/figures/history.svg", s1 + `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript Residual History (N=63)</text></svg>\n`);
  let s2 = hdr;
  for (let m = 0; m < 7; m++) {
    const cnt = m <= 2 ? 3 : 5, pts: string[] = [];
    for (let i = 0; i < cnt; i++) {
      const px = 60.0 + 500.0 * (Math.log10(NS_ALL[i]) - Math.log10(15)) / (Math.log10(255) - Math.log10(15));
      const py = 310.0 - 280.0 * Math.log10(gRes[m][i].iters) / 5.0;
      pts.push(`${px.toFixed(1)},${py.toFixed(1)}`);
    }
    s2 += `<polyline fill="none" stroke="${cols[m]}" stroke-width="1.5" points="${pts.join(" ")}"/>`;
  }
  writeFileSync("out/ts/figures/scaling.svg", s2 + `<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript Iterations vs N</text></svg>\n`);
  const pts3 = NS_ALL.map((_, i) => `${(60 + 500 * i / 4).toFixed(1)},${(310 - 280 * (Math.log10(gDiscErr[i]) + 5.0) / 4.0).toFixed(1)}`).join(" ");
  writeFileSync("out/ts/figures/order.svg", hdr + `<polyline fill="none" stroke="#1f77b4" stroke-width="1.8" points="${pts3}"/><text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript Discretization Error vs h (Order ~2)</text></svg>\n`);
}
function s19(): void {
  const tot = gStageS.reduce((a, b) => a + b, 0);
  writeFileSync("out/ts/s19_summary.json", `{"lang":"ts","N_max":255,"N_stat_max":63,"toolchain":"node v24.21.0 (--experimental-strip-types)","total_sec":${tot.toFixed(3)}}\n`);
}

mkdirSync("out/ts/figures", { recursive: true });
const cmd = process.argv[2] ?? "all";
const st: [string, () => void][] = [
  ["s01", s01], ["s02", s02], ["s03", s03],
  ["s04", () => runM(0, "out/ts/s04_jacobi.csv", 3)],
  ["s05", () => runM(1, "out/ts/s05_gs.csv", 3)],
  ["s06", () => runM(2, "out/ts/s06_rbgs.csv", 3)],
  ["s07", () => runM(3, "out/ts/s07_sor.csv", 5)],
  ["s08", () => runM(4, "out/ts/s08_cg.csv", 5)],
  ["s09", () => runM(5, "out/ts/s09_pcg_ssor.csv", 5)],
  ["s10", () => runM(6, "out/ts/s10_mg.csv", 5)],
  ["s11", s11], ["s12", s12], ["s13", s13], ["s14", s14], ["s15", s15], ["s16", s16], ["s17", s17], ["s18", s18], ["s19", s19],
];
if (cmd === "all") {
  for (let i = 0; i < st.length; i++) {
    console.log(`[ts] tahap ${i + 1}/19 (${st[i][0]})`);
    const t0 = performance.now();
    st[i][1]();
    gStageS[i] = (performance.now() - t0) / 1000.0;
  }
} else {
  const f = st.find(([nm]) => nm === cmd);
  if (!f) process.exit(1);
  f[1]();
}
