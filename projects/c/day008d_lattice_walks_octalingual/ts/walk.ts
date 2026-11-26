import { mkdirSync, writeFileSync } from "node:fs";

const MASTER_SEED = 20261126n;
const MASK64 = 0xFFFFFFFFFFFFFFFFn;
const W_WALKS = 20000;
const W_RUIN = 5000;
const BENCH_STEPS = 2000000;
const T_EXACT_MAX = 5000;
const TS_MSD = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096];
const TS_RET = [10, 100, 1000, 10000];
const TWO_M53 = 2 ** -53;

function sm64(z: bigint): bigint {
  z = ((z ^ (z >> 30n)) * 0xBF58476D1CE4E5B9n) & MASK64;
  z = ((z ^ (z >> 27n)) * 0x94D049BB133111EBn) & MASK64;
  return (z ^ (z >> 31n)) & MASK64;
}

class Rng {
  lo = new Uint32Array(4);
  hi = new Uint32Array(4);
  constructor(seed: bigint) {
    let st = seed & MASK64;
    for (let i = 0; i < 4; i++) {
      st = (st + 0x9E3779B97F4A7C15n) & MASK64;
      const v = sm64(st);
      this.lo[i] = Number(v & 0xFFFFFFFFn);
      this.hi[i] = Number((v >> 32n) & 0xFFFFFFFFn);
    }
  }
  u64Parts(): [number, number] {
    const l0 = this.lo[0], h0 = this.hi[0], l3 = this.lo[3], h3 = this.hi[3];
    const sumLo = (l0 + l3) >>> 0;
    const carry = l0 + l3 >= 0x100000000 ? 1 : 0;
    const sumHi = (h0 + h3 + carry) >>> 0;
    const rotLo = ((sumLo << 23) | (sumHi >>> 9)) >>> 0;
    const rotHi = ((sumHi << 23) | (sumLo >>> 9)) >>> 0;
    const resLo = (rotLo + l0) >>> 0;
    const c2 = rotLo + l0 >= 0x100000000 ? 1 : 0;
    const resHi = (rotHi + h0 + c2) >>> 0;
    const tLo = (this.lo[1] << 17) >>> 0;
    const tHi = ((this.hi[1] << 17) | (this.lo[1] >>> 15)) >>> 0;
    this.lo[2] ^= l0; this.hi[2] ^= h0;
    this.lo[3] ^= this.lo[1]; this.hi[3] ^= this.hi[1];
    this.lo[1] ^= this.lo[2]; this.hi[1] ^= this.hi[2];
    this.lo[0] ^= this.lo[3]; this.hi[0] ^= this.hi[3];
    this.lo[2] ^= tLo; this.hi[2] ^= tHi;
    const s3l = this.lo[3], s3h = this.hi[3];
    this.lo[3] = ((s3h << 13) | (s3l >>> 19)) >>> 0;
    this.hi[3] = ((s3l << 13) | (s3h >>> 19)) >>> 0;
    return [resLo, resHi];
  }
  u01(): number {
    const [lo, hi] = this.u64Parts();
    const top53 = hi * 2097152 + (lo >>> 11);
    return top53 * TWO_M53;
  }
  rint(n: number): number {
    const k = Math.floor(this.u01() * n);
    return k >= n ? n - 1 : k;
  }
}

function taskSeed(d: number, st: number, idx: number): bigint {
  let h = 0x517CC1B727220A95n;
  for (const v of [d, st, idx]) h = sm64(h ^ (BigInt(v) & MASK64));
  return sm64(MASTER_SEED ^ h);
}

function fnv1aI32(coords: number[]): bigint {
  let h = 0xCBF29CE484222325n;
  for (const val of coords) {
    const u = val >>> 0;
    for (let b = 0; b < 4; b++) {
      h = ((h ^ BigInt((u >>> (8 * b)) & 0xFF)) * 0x100000001B3n) & MASK64;
    }
  }
  return h;
}

mkdirSync("out/ts/figures", { recursive: true });
const gMsd = Array.from({ length: 3 }, () => new Float64Array(13));
const gRetSim = Array.from({ length: 3 }, () => new Float64Array(4));
const gRetSe = Array.from({ length: 3 }, () => new Float64Array(4));
const gRetEx = Array.from({ length: 3 }, () => new Float64Array(4));
const gRuinP = new Float64Array(6), gRuinPe = new Float64Array(6);
const gS03_1024 = Array.from({ length: 3 }, () => new Float64Array(W_WALKS));
const gS04_10k = Array.from({ length: 3 }, () => new Uint8Array(W_WALKS));
const stageSec = new Float64Array(12);

function hex64(lo: number, hi: number): string {
  return hi.toString(16).padStart(8, "0") + lo.toString(16).padStart(8, "0");
}

function s01() {
  const r = new Rng(MASTER_SEED);
  const [l0, h0] = r.u64Parts(), [l1, h1] = r.u64Parts(), [l2, h2] = r.u64Parts();
  const x0 = hex64(l0, h0), x1 = hex64(l1, h1), x2 = hex64(l2, h2);
  const ok = x0 === "a9920e0e7dab3ef0" && x1 === "033b9ceb4b3ba7a0" && x2 === "e7473a6a39afcea6";
  const w = [1.0, 1.0 / 3.0, 0, 0, 0];
  for (let n = 1; n < 4; n++) w[n + 1] = ((10 * n * n + 10 * n + 3) * w[n] - n * n * w[n - 1]) / (9 * (n + 1) * (n + 1));
  const eok = Math.abs(w[2] - 15 / 81) < 1e-14 && Math.abs(w[3] - 93 / 729) < 1e-14 && Math.abs(w[4] - 639 / 6561) < 1e-14;
  writeFileSync("out/ts/s01_selftest.txt", `ok,x0,x1,x2,enum_ok\n${ok && eok ? 1 : 0},${x0},${x1},${x2},${eok ? 1 : 0}\n`);
  if (!(ok && eok)) process.exit(1);
}

function s02() {
  const lines = ["d,walk,x0,x1,x2,fnv1a_hex"];
  for (let d = 1; d <= 3; d++) {
    const coords: number[] = [], posAll: number[][] = [];
    for (let w = 0; w < 8; w++) {
      const r = new Rng(taskSeed(d, 2, w));
      const p = [0, 0, 0];
      for (let s = 0; s < 200; s++) { const k = r.rint(2 * d); p[k >> 1] += (k & 1) ? -1 : 1; }
      posAll.push(p);
      for (let a = 0; a < d; a++) coords.push(p[a]);
    }
    const h = fnv1aI32(coords).toString(16).padStart(16, "0");
    for (let w = 0; w < 8; w++) lines.push(`${d},${w},${posAll[w][0]},${posAll[w][1]},${posAll[w][2]},${h}`);
  }
  writeFileSync("out/ts/s02_parity.csv", lines.join("\n") + "\n");
}

function s03() {
  const lines = ["d,t,msd,se,exact"];
  for (let d = 1; d <= 3; d++) {
    const s1 = new Float64Array(13), s2 = new Float64Array(13);
    for (let w = 0; w < W_WALKS; w++) {
      const r = new Rng(taskSeed(d, 3, w));
      let p0 = 0, p1 = 0, p2 = 0, ti = 0;
      for (let t = 1; t <= 4096; t++) {
        const k = r.rint(2 * d), delta = (k & 1) ? -1 : 1, ax = k >> 1;
        if (ax === 0) p0 += delta; else if (ax === 1) p1 += delta; else p2 += delta;
        if (t === TS_MSD[ti]) {
          const r2 = p0 * p0 + p1 * p1 + p2 * p2;
          s1[ti] += r2; s2[ti] += r2 * r2;
          if (ti === 10) gS03_1024[d - 1][w] = r2;
          ti++;
        }
      }
    }
    for (let ti = 0; ti < 13; ti++) {
      const m = s1[ti] / W_WALKS, se = Math.sqrt(Math.max(0, s2[ti] / W_WALKS - m * m) / W_WALKS);
      gMsd[d - 1][ti] = m;
      lines.push(`${d},${TS_MSD[ti]},${m.toFixed(8)},${se.toFixed(8)},${TS_MSD[ti]}`);
    }
  }
  writeFileSync("out/ts/s03_msd.csv", lines.join("\n") + "\n");
}

function s04() {
  const lines = ["d,T,p_hat,se"];
  for (let d = 1; d <= 3; d++) {
    const hits = [0, 0, 0, 0];
    for (let w = 0; w < W_WALKS; w++) {
      const r = new Rng(taskSeed(d, 4, w));
      let p0 = 0, p1 = 0, p2 = 0, retT = 0;
      for (let t = 1; t <= 10000; t++) {
        const k = r.rint(2 * d), delta = (k & 1) ? -1 : 1, ax = k >> 1;
        if (ax === 0) p0 += delta; else if (ax === 1) p1 += delta; else p2 += delta;
        if ((p0 | p1 | p2) === 0) { retT = t; break; }
      }
      if (retT > 0) {
        gS04_10k[d - 1][w] = 1;
        for (let ti = 0; ti < 4; ti++) if (retT <= TS_RET[ti]) hits[ti]++;
      }
    }
    for (let ti = 0; ti < 4; ti++) {
      const ph = hits[ti] / W_WALKS, se = Math.sqrt(Math.max(0, ph * (1 - ph)) / W_WALKS);
      gRetSim[d - 1][ti] = ph; gRetSe[d - 1][ti] = se;
      lines.push(`${d},${TS_RET[ti]},${ph.toFixed(8)},${se.toFixed(8)}`);
    }
  }
  writeFileSync("out/ts/s04_return_sim.csv", lines.join("\n") + "\n");
}

function s05() {
  const u1 = new Float64Array(T_EXACT_MAX + 1).fill(1), u2 = new Float64Array(T_EXACT_MAX + 1).fill(1);
  const w3 = new Float64Array(T_EXACT_MAX + 1).fill(1), u3 = new Float64Array(T_EXACT_MAX + 1).fill(1), fv = new Float64Array(T_EXACT_MAX + 1);
  w3[1] = 1 / 3;
  for (let n = 1; n < T_EXACT_MAX; n++) w3[n + 1] = ((10 * n * n + 10 * n + 3) * w3[n] - n * n * w3[n - 1]) / (9 * (n + 1) * (n + 1));
  for (let n = 1; n <= T_EXACT_MAX; n++) { u1[n] = u1[n - 1] * (2 * n - 1) / (2 * n); u2[n] = u1[n] * u1[n]; u3[n] = u1[n] * w3[n]; }
  const us = [u1, u2, u3], lines = ["d,T,F_exact"];
  for (let d = 1; d <= 3; d++) {
    const u = us[d - 1]; let acc = 0;
    for (let n = 1; n <= T_EXACT_MAX; n++) {
      let s = u[n]; for (let k = 1; k < n; k++) s -= fv[k] * u[n - k];
      fv[n] = s; acc += s;
      for (let ti = 0; ti < 4; ti++) if (2 * n === TS_RET[ti]) { gRetEx[d - 1][ti] = acc; lines.push(`${d},${TS_RET[ti]},${acc.toFixed(10)}`); }
    }
  }
  writeFileSync("out/ts/s05_return_exact.csv", lines.join("\n") + "\n");
}

function s06() {
  const lines = ["d,T,p_hat,se,F_exact,z_score,p_inf_est"];
  const p3Inf = gRetSim[2][3] + (gRetSim[2][3] - gRetSim[2][2]) / (Math.sqrt(10) - 1);
  for (let d = 1; d <= 3; d++)
    for (let ti = 0; ti < 4; ti++) {
      const z = gRetSe[d - 1][ti] > 1e-12 ? (gRetSim[d - 1][ti] - gRetEx[d - 1][ti]) / gRetSe[d - 1][ti] : 0;
      lines.push(`${d},${TS_RET[ti]},${gRetSim[d - 1][ti].toFixed(8)},${gRetSe[d - 1][ti].toFixed(8)},${gRetEx[d - 1][ti].toFixed(10)},${z.toFixed(4)},${(d === 3 && ti === 3 ? p3Inf : 0).toFixed(8)}`);
    }
  writeFileSync("out/ts/s06_return_compare.csv", lines.join("\n") + "\n");
}

function s07() {
  const ns = [10, 10, 10, 50, 50, 50], xs = [2, 5, 8, 10, 25, 40], lines = ["N,x,p_hat,p_se,p_exact,p_z,t_hat,t_se,t_exact,t_z"];
  for (let ci = 0; ci < 6; ci++) {
    const N = ns[ci], x0 = xs[ci]; let wins = 0, st1 = 0, st2 = 0;
    for (let w = 0; w < W_RUIN; w++) {
      const r = new Rng(taskSeed(1, 7, ci * 100000 + w)); let x = x0, steps = 0;
      while (x > 0 && x < N) { x += r.rint(2) === 0 ? 1 : -1; steps++; }
      if (x === N) wins++; st1 += steps; st2 += steps * steps;
    }
    const ph = wins / W_RUIN, pse = Math.sqrt(Math.max(0, ph * (1 - ph)) / W_RUIN), pe = x0 / N;
    const th = st1 / W_RUIN, tse = Math.sqrt(Math.max(0, st2 / W_RUIN - th * th) / W_RUIN), te = x0 * (N - x0);
    gRuinP[ci] = ph; gRuinPe[ci] = pe;
    lines.push(`${N},${x0},${ph.toFixed(8)},${pse.toFixed(8)},${pe.toFixed(8)},${((ph - pe) / pse).toFixed(4)},${th.toFixed(4)},${tse.toFixed(4)},${te.toFixed(4)},${((th - te) / tse).toFixed(4)}`);
  }
  writeFileSync("out/ts/s07_ruin.csv", lines.join("\n") + "\n");
}

function s08() {
  const lines = ["d,var_x0,var_exact,kurt_x0,kurt_exact,bin_idx,bin_lo,bin_hi,count"], hist = new Int32Array(40), nEp = 5000;
  for (let d = 1; d <= 3; d++) {
    let m1 = 0, m2 = 0, m4 = 0;
    for (let w = 0; w < nEp; w++) {
      const r = new Rng(taskSeed(d, 8, w)); let p0 = 0;
      for (let s = 0; s < 1024; s++) { const k = r.rint(2 * d); if ((k >> 1) === 0) p0 += (k & 1) ? -1 : 1; }
      m1 += p0; m2 += p0 * p0; m4 += p0 * p0 * p0 * p0;
      if (d === 1) hist[Math.min(39, Math.max(0, Math.floor((p0 + 128) / 6.4)))]++;
    }
    m1 /= nEp; const v = m2 / nEp - m1 * m1;
    lines.push(`${d},${v.toFixed(6)},${(1024 / d).toFixed(6)},${((m4 / nEp) / (v * v)).toFixed(6)},${(3 + (d - 3) / 1024).toFixed(6)},-1,0,0,${nEp}`);
  }
  for (let b = 0; b < 40; b++) lines.push(`1,0,0,0,0,${b},${(-128 + b * 6.4).toFixed(2)},${(-128 + (b + 1) * 6.4).toFixed(2)},${hist[b]}`);
  writeFileSync("out/ts/s08_endpoint.csv", lines.join("\n") + "\n");
}

function s09() {
  const lines = ["metric,d,estimate,boot_se,ci_lo,ci_hi"];
  for (let d = 1; d <= 3; d++) {
    const r1 = new Rng(taskSeed(d, 9, 0)), boot: number[] = [];
    for (let b = 0; b < 200; b++) { let s = 0; for (let i = 0; i < W_WALKS; i++) s += gS03_1024[d - 1][r1.rint(W_WALKS)]; boot.push(s / W_WALKS); }
    boot.sort((a, b) => a - b); const bm = boot.reduce((a, b) => a + b, 0) / 200;
    const bv = boot.reduce((a, v) => a + (v - bm) ** 2, 0) / 199;
    lines.push(`msd_1024,${d},${gMsd[d - 1][10].toFixed(6)},${Math.sqrt(bv).toFixed(6)},${boot[5].toFixed(6)},${boot[194].toFixed(6)}`);

    const r2 = new Rng(taskSeed(d, 9, 1)), b2: number[] = [];
    for (let b = 0; b < 200; b++) { let s = 0; for (let i = 0; i < W_WALKS; i++) s += gS04_10k[d - 1][r2.rint(W_WALKS)]; b2.push(s / W_WALKS); }
    b2.sort((a, b) => a - b); const bm2 = b2.reduce((a, b) => a + b, 0) / 200;
    const bv2 = b2.reduce((a, v) => a + (v - bm2) ** 2, 0) / 199;
    lines.push(`return_10000,${d},${gRetSim[d - 1][3].toFixed(6)},${Math.sqrt(bv2).toFixed(6)},${b2[5].toFixed(6)},${b2[194].toFixed(6)}`);
  }
  writeFileSync("out/ts/s09_ci.csv", lines.join("\n") + "\n");
}

function s10() {
  const lines = ["d,slope,slope_se,intercept"], xs = TS_MSD.map((t) => Math.log(t));
  const mx = xs.reduce((a, b) => a + b, 0) / 13, den = xs.reduce((a, x) => a + (x - mx) ** 2, 0);
  for (let d = 1; d <= 3; d++) {
    const ys = Array.from(gMsd[d - 1], (m) => Math.log(m)), my = ys.reduce((a, b) => a + b, 0) / 13;
    const sl = xs.reduce((a, x, i) => a + (x - mx) * (ys[i] - my), 0) / den, ic = my - sl * mx;
    const rss = ys.reduce((a, y, i) => a + (y - (sl * xs[i] + ic)) ** 2, 0);
    lines.push(`${d},${sl.toFixed(6)},${Math.sqrt((rss / 11) / den).toFixed(6)},${ic.toFixed(6)}`);
  }
  writeFileSync("out/ts/s10_fit.csv", lines.join("\n") + "\n");
}

function s11() {
  const lines = ["d,steps_per_rep,reps,ns_med,ns_mad"];
  for (let d = 1; d <= 3; d++) {
    const vals: number[] = [];
    for (let rep = 0; rep < 7; rep++) {
      const r = new Rng(taskSeed(d, 11, rep)); let p0 = 0;
      const t0 = performance.now();
      for (let s = 0; s < BENCH_STEPS; s++) { const k = r.rint(2 * d); p0 += (k & 1) ? -1 : 1; }
      vals.push((performance.now() - t0) * 1e6 / BENCH_STEPS + (p0 & 0));
    }
    vals.sort((a, b) => a - b); const med = vals[3], dev = vals.map((v) => Math.abs(v - med)).sort((a, b) => a - b);
    lines.push(`${d},${BENCH_STEPS},7,${med.toFixed(4)},${dev[3].toFixed(4)}`);
  }
  writeFileSync("out/ts/s11_bench.csv", lines.join("\n") + "\n");
}

function s12() {
  const cols = ["#1f77b4", "#ff7f0e", "#2ca02c"];
  let s1 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>';
  for (let d = 0; d < 3; d++) {
    const pts = TS_MSD.map((t, i) => `${(60 + 500 * Math.log10(t) / Math.log10(4096)).toFixed(1)},${(310 - 280 * Math.log10(Math.max(1, gMsd[d][i])) / Math.log10(4096)).toFixed(1)}`).join(" ");
    s1 += `<polyline fill="none" stroke="${cols[d]}" stroke-width="1.6" points="${pts}"/>`;
  }
  writeFileSync("out/ts/figures/msd.svg", s1 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript MSD Log-Log (y=t)</text></svg>\n');
  let s2 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/>';
  for (let d = 0; d < 3; d++) {
    const pts = [0, 1, 2, 3].map((i) => `${(60 + 500 * (i / 3)).toFixed(1)},${(310 - 280 * gRetEx[d][i]).toFixed(1)}`).join(" ");
    s2 += `<polyline fill="none" stroke="${cols[d]}" stroke-width="1.6" points="${pts}"/>` + [0, 1, 2, 3].map((i) => `<circle cx="${(60 + 500 * (i / 3)).toFixed(1)}" cy="${(310 - 280 * gRetSim[d][i]).toFixed(1)}" r="4" fill="${cols[d]}"/>`).join("");
  }
  writeFileSync("out/ts/figures/return_prob.svg", s2 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript Return Prob F_d(T) Exact vs Sim</text></svg>\n');
  let s3 = '<svg xmlns="http://www.w3.org/2000/svg" width="600" height="360"><rect width="600" height="360" fill="#fff"/><rect x="60" y="30" width="500" height="280" fill="none" stroke="#ccc"/><line x1="60" y1="310" x2="560" y2="30" stroke="#000" stroke-dasharray="4,3"/>';
  for (let i = 0; i < 6; i++) s3 += `<circle cx="${(60 + 500 * gRuinPe[i]).toFixed(1)}" cy="${(310 - 280 * gRuinP[i]).toFixed(1)}" r="5" fill="#1f77b4"/>`;
  writeFileSync("out/ts/figures/ruin.svg", s3 + '<text x="300" y="20" text-anchor="middle" font-family="sans-serif" font-size="13">TypeScript Gambler\'s Ruin Exact vs Sim</text></svg>\n');
  const tot = stageSec.reduce((a, b) => a + b, 0);
  writeFileSync("out/ts/s12_summary.json", `{"lang":"ts","W":${W_WALKS},"ruin_W":${W_RUIN},"bench_steps":${BENCH_STEPS},"toolchain":"node v24.21.0 (--experimental-strip-types)","total_sec":${tot.toFixed(3)}}\n`);
}

const STAGES: [string, () => void][] = [["s01",s01],["s02",s02],["s03",s03],["s04",s04],["s05",s05],["s06",s06],["s07",s07],["s08",s08],["s09",s09],["s10",s10],["s11",s11],["s12",s12]];
const cmd = process.argv[2] ?? "all";
if (cmd === "all") {
  STAGES.forEach(([nm, fn], i) => { console.log(`[ts] tahap ${i + 1}/12 (${nm})`); const t0 = performance.now(); fn(); stageSec[i] = (performance.now() - t0) / 1000; });
} else {
  const f = STAGES.find(([nm]) => nm === cmd); if (!f) process.exit(1); f[1]();
}
