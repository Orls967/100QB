export interface DtmcRes {
  pBlock: number;
  l: number;
  lQ: number;
  u: number;
  w: number;
  wQ: number;
  lamEff: number;
}

function emptyDtmc(): DtmcRes {
  return { pBlock: 0, l: 0, lQ: 0, u: 0, w: 0, wQ: 0, lamEff: 0 };
}

export function binomProbs(m: number, mu: number): number[] {
  const out = new Array<number>(8).fill(0.0);
  out[0] = 1.0;
  for (let i = 0; i < m; i++) {
    for (let k = m; k >= 1; k--) {
      out[k] = out[k] * (1.0 - mu) + out[k - 1] * mu;
    }
    out[0] *= 1.0 - mu;
  }
  return out;
}

export function solveDtmcSingle(c: number, kCap: number, lam: number, mu: number): DtmcRes {
  const p: number[][] = Array.from({ length: 64 }, () => new Array<number>(64).fill(0.0));
  for (let n = 0; n <= kCap; n++) {
    const m = Math.min(n, c);
    const bp = binomProbs(m, mu);
    for (let d = 0; d <= m; d++) {
      if (n < kCap) {
        p[n][n - d] += (1.0 - lam) * bp[d];
        p[n][n + 1 - d] += lam * bp[d];
      } else {
        p[n][n - d] += bp[d];
      }
    }
  }
  const pi = new Array<number>(64).fill(0.0);
  for (let i = 0; i <= kCap; i++) pi[i] = 1.0 / (kCap + 1);
  for (let it = 0; it < 25000; it++) {
    let maxd = 0.0;
    let sum = 0.0;
    for (let j = 0; j <= kCap; j++) {
      let nv = 0.0;
      for (let i = 0; i <= kCap; i++) nv += pi[i] * p[i][j];
      const df = Math.abs(nv - pi[j]);
      if (df > maxd) maxd = df;
      pi[j] = nv;
      sum += nv;
    }
    for (let j = 0; j <= kCap; j++) pi[j] /= sum;
    if (maxd < 1e-14) break;
  }
  const r = emptyDtmc();
  r.pBlock = pi[kCap];
  r.lamEff = lam * (1.0 - pi[kCap]);
  for (let n = 0; n <= kCap; n++) {
    const m = Math.min(n, c);
    r.l += n * pi[n];
    r.lQ += (n - m) * pi[n];
    r.u += (m / c) * pi[n];
  }
  r.w = r.lamEff > 0.0 ? r.l / r.lamEff : 0.0;
  r.wQ = r.lamEff > 0.0 ? r.lQ / r.lamEff : 0.0;
  return r;
}

export function solveGeoClosedForm(lam: number, mu: number): DtmcRes {
  const l = (lam * (1.0 - lam)) / (mu - lam);
  const w = (1.0 - lam) / (mu - lam);
  const wq = w - 1.0 / mu;
  return { pBlock: 0.0, l, lQ: lam * wq, u: lam / mu, w, wQ: wq, lamEff: lam };
}

export function solveDtmcTandem(k1: number, k2: number, lam: number, mu1: number, mu2: number): [DtmcRes, DtmcRes, number] {
  const sTot = (k1 + 1) * (k2 + 1);
  const p: number[][] = Array.from({ length: 64 }, () => new Array<number>(64).fill(0.0));
  for (let n1 = 0; n1 <= k1; n1++) {
    for (let n2 = 0; n2 <= k2; n2++) {
      const sIdx = n1 * (k2 + 1) + n2;
      for (let a = 0; a <= 1; a++) {
        const pa = a !== 0 ? lam : 1.0 - lam;
        const d1Max = n1 > 0 ? 1 : 0;
        for (let d1 = 0; d1 <= d1Max; d1++) {
          const pd1 = n1 === 0 ? 1.0 : d1 !== 0 ? mu1 : 1.0 - mu1;
          const d2Max = n2 > 0 ? 1 : 0;
          for (let d2 = 0; d2 <= d2Max; d2++) {
            const pd2 = n2 === 0 ? 1.0 : d2 !== 0 ? mu2 : 1.0 - mu2;
            const nn1 = n1 < k1 ? n1 + a - d1 : n1 - d1;
            const n2Mid = n2 - d2;
            const nn2 = n2Mid < k2 ? n2Mid + d1 : n2Mid;
            p[sIdx][nn1 * (k2 + 1) + nn2] += pa * pd1 * pd2;
          }
        }
      }
    }
  }
  const pi = new Array<number>(64).fill(0.0);
  for (let i = 0; i < sTot; i++) pi[i] = 1.0 / sTot;
  for (let it = 0; it < 25000; it++) {
    let maxd = 0.0;
    let sum = 0.0;
    for (let j = 0; j < sTot; j++) {
      let nv = 0.0;
      for (let i = 0; i < sTot; i++) nv += pi[i] * p[i][j];
      const df = Math.abs(nv - pi[j]);
      if (df > maxd) maxd = df;
      pi[j] = nv;
      sum += nv;
    }
    for (let j = 0; j < sTot; j++) pi[j] /= sum;
    if (maxd < 1e-14) break;
  }
  const r1 = emptyDtmc();
  const r2 = emptyDtmc();
  let blk1 = 0.0;
  let blk2 = 0.0;
  for (let n1 = 0; n1 <= k1; n1++) {
    for (let n2 = 0; n2 <= k2; n2++) {
      const pv = pi[n1 * (k2 + 1) + n2];
      r1.l += n1 * pv;
      r2.l += n2 * pv;
      if (n1 > 0) {
        r1.u += pv;
        r1.lQ += (n1 - 1) * pv;
      }
      if (n2 > 0) {
        r2.u += pv;
        r2.lQ += (n2 - 1) * pv;
      }
      if (n1 === k1) blk1 += pv * lam;
      if (n1 > 0 && n2 === k2) blk2 += pv * mu1 * (1.0 - mu2);
    }
  }
  r1.pBlock = blk1;
  r1.lamEff = lam - blk1;
  r1.w = r1.lamEff > 0.0 ? r1.l / r1.lamEff : 0.0;
  r1.wQ = r1.lamEff > 0.0 ? r1.lQ / r1.lamEff : 0.0;
  r2.pBlock = blk2;
  r2.lamEff = r1.lamEff - blk2;
  r2.w = r2.lamEff > 0.0 ? r2.l / r2.lamEff : 0.0;
  r2.wQ = r2.lamEff > 0.0 ? r2.lQ / r2.lamEff : 0.0;
  return [r1, r2, blk2 / (r1.lamEff > 0.0 ? r1.lamEff : 1.0)];
}
