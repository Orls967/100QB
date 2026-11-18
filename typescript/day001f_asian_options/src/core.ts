import { SeededRng, deriveChunkSeed } from "./rng.js";
import { normCdf } from "./stats.js";

export const CHUNK_SIZE = 4_000;

export interface OptionParams {
  s0: number;
  k: number;
  r: number;
  sigma: number;
  t: number;
  n_steps: number;
}

export type ChunkTask = "european" | "naive" | "antithetic" | "cv" | "delta_pw" | "delta_bump";

export function bsCall(p: OptionParams): number {
  if (p.t <= 0 || p.sigma <= 0) {
    return Math.max(p.s0 - p.k * Math.exp(-p.r * p.t), 0);
  }
  const volSqrtT = p.sigma * Math.sqrt(p.t);
  const d1 = (Math.log(p.s0 / p.k) + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / volSqrtT;
  const d2 = d1 - volSqrtT;
  return p.s0 * normCdf(d1) - p.k * Math.exp(-p.r * p.t) * normCdf(d2);
}

export function bsPut(p: OptionParams): number {
  if (p.t <= 0 || p.sigma <= 0) {
    return Math.max(p.k * Math.exp(-p.r * p.t) - p.s0, 0);
  }
  const volSqrtT = p.sigma * Math.sqrt(p.t);
  const d1 = (Math.log(p.s0 / p.k) + (p.r + 0.5 * p.sigma * p.sigma) * p.t) / volSqrtT;
  const d2 = d1 - volSqrtT;
  return p.k * Math.exp(-p.r * p.t) * normCdf(-d2) - p.s0 * normCdf(-d1);
}

export function geometricAsianCall(p: OptionParams): number {
  if (p.n_steps < 1) {
    throw new Error("n_steps must be >= 1");
  }
  const n = p.n_steps;
  const muLnG = Math.log(p.s0) + ((p.r - 0.5 * p.sigma * p.sigma) * p.t * (n + 1)) / (2 * n);
  const varLnG = (p.sigma * p.sigma * p.t * (n + 1) * (2 * n + 1)) / (6 * n * n);
  const sigmaLnG = Math.sqrt(varLnG);
  const expectedG = Math.exp(muLnG + 0.5 * varLnG);
  const d1 = (muLnG - Math.log(p.k) + varLnG) / sigmaLnG;
  const d2 = d1 - sigmaLnG;
  return Math.exp(-p.r * p.t) * (expectedG * normCdf(d1) - p.k * normCdf(d2));
}

export function runChunk(
  task: ChunkTask,
  p: OptionParams,
  chunkIdx: number,
  chunkLen: number,
  seed: number,
  hRel = 0.01,
): Float64Array {
  const rng = new SeededRng(deriveChunkSeed(seed, chunkIdx));
  const dt = p.t / p.n_steps;
  const drift = (p.r - 0.5 * p.sigma * p.sigma) * dt;
  const vol = p.sigma * Math.sqrt(dt);
  const disc = Math.exp(-p.r * p.t);
  const logS0 = Math.log(p.s0);
  const invN = 1.0 / p.n_steps;

  if (task === "cv") {
    let sx = 0.0;
    let sy = 0.0;
    let sxx = 0.0;
    let syy = 0.0;
    let sxy = 0.0;
    for (let i = 0; i < chunkLen; i += 1) {
      let logS = logS0;
      let sumS = 0.0;
      let sumLogS = 0.0;
      for (let step = 0; step < p.n_steps; step += 1) {
        logS += drift + vol * rng.nextNormal();
        sumS += Math.exp(logS);
        sumLogS += logS;
      }
      const arithAvg = sumS * invN;
      const geomAvg = Math.exp(sumLogS * invN);
      const x = Math.max(arithAvg - p.k, 0.0) * disc;
      const y = Math.max(geomAvg - p.k, 0.0) * disc;
      sx += x;
      sy += y;
      sxx += x * x;
      syy += y * y;
      sxy += x * y;
    }
    return new Float64Array([sx, sy, sxx, syy, sxy]);
  }

  let sum = 0.0;
  let sumSq = 0.0;

  if (task === "antithetic") {
    for (let i = 0; i < chunkLen; i += 1) {
      let logSPos = logS0;
      let logSNeg = logS0;
      let sumPos = 0.0;
      let sumNeg = 0.0;
      for (let step = 0; step < p.n_steps; step += 1) {
        const inc = vol * rng.nextNormal();
        logSPos += drift + inc;
        logSNeg += drift - inc;
        sumPos += Math.exp(logSPos);
        sumNeg += Math.exp(logSNeg);
      }
      const pairAvg =
        0.5 * (Math.max(sumPos * invN - p.k, 0.0) + Math.max(sumNeg * invN - p.k, 0.0)) * disc;
      sum += pairAvg;
      sumSq += pairAvg * pairAvg;
    }
    return new Float64Array([sum, sumSq]);
  }

  const h = hRel * p.s0;
  const scaleUp = (p.s0 + h) / p.s0;
  const scaleDn = (p.s0 - h) / p.s0;

  for (let i = 0; i < chunkLen; i += 1) {
    let logS = logS0;
    let sumS = 0.0;
    for (let step = 0; step < p.n_steps; step += 1) {
      logS += drift + vol * rng.nextNormal();
      sumS += Math.exp(logS);
    }
    let val = 0.0;
    if (task === "european") {
      val = Math.max(Math.exp(logS) - p.k, 0.0) * disc;
    } else if (task === "naive") {
      val = Math.max(sumS * invN - p.k, 0.0) * disc;
    } else if (task === "delta_pw") {
      const arithAvg = sumS * invN;
      val = arithAvg > p.k ? disc * (arithAvg / p.s0) : 0.0;
    } else {
      const arithAvg = sumS * invN;
      const payUp = Math.max(arithAvg * scaleUp - p.k, 0.0);
      const payDn = Math.max(arithAvg * scaleDn - p.k, 0.0);
      val = (disc * (payUp - payDn)) / (2.0 * h);
    }
    sum += val;
    sumSq += val * val;
  }

  return new Float64Array([sum, sumSq]);
}
