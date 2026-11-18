import { existsSync } from "node:fs";
import os from "node:os";
import { Worker } from "node:worker_threads";
import {
  CHUNK_SIZE,
  type ChunkTask,
  type OptionParams,
  bsCall,
  bsPut,
  geometricAsianCall,
  runChunk,
} from "./core.js";
import { type EstimatorResult, fitLogLogSlope, summarizeMoments } from "./stats.js";
import type { WorkerRequest, WorkerResponse } from "./worker.js";

export {
  bsCall,
  bsPut,
  geometricAsianCall,
  fitLogLogSlope,
  type OptionParams,
  type EstimatorResult,
};

export interface ControlVariateResult extends EstimatorResult {
  beta: number;
  y_exact: number;
}

export interface DeltaResult {
  delta: number;
  se: number;
  ci_low: number;
  ci_high: number;
}

function validateMC(p: OptionParams, nPaths: number): void {
  if (p.n_steps < 1) {
    throw new Error("n_steps must be >= 1");
  }
  if (nPaths < 2) {
    throw new Error("n_paths must be >= 2");
  }
}

function resolveWorkerUrl(): URL {
  const sameDirJs = new URL("./worker.js", import.meta.url);
  if (existsSync(sameDirJs)) {
    return sameDirJs;
  }
  return new URL("../dist/worker.js", import.meta.url);
}

async function runParallelChunks(
  task: ChunkTask,
  params: OptionParams,
  nItems: number,
  seed: number,
  hRel = 0.01,
  workers?: number,
): Promise<Float64Array[]> {
  const nChunks = Math.ceil(nItems / CHUNK_SIZE);
  const w = Math.min(workers && workers > 0 ? workers : os.availableParallelism(), nChunks);

  if (w <= 1) {
    const out: Float64Array[] = new Array(nChunks);
    for (let c = 0; c < nChunks; c += 1) {
      const start = c * CHUNK_SIZE;
      const len = Math.min(nItems - start, CHUNK_SIZE);
      out[c] = runChunk(task, params, c, len, seed, hRel);
    }
    return out;
  }

  const buckets: Array<Array<{ idx: number; len: number }>> = Array.from({ length: w }, () => []);
  for (let c = 0; c < nChunks; c += 1) {
    const start = c * CHUNK_SIZE;
    const len = Math.min(nItems - start, CHUNK_SIZE);
    buckets[c % w]?.push({ idx: c, len });
  }

  const workerUrl = resolveWorkerUrl();
  const out: Float64Array[] = new Array(nChunks);

  await Promise.all(
    buckets.map(
      (chunks) =>
        new Promise<void>((resolve, reject) => {
          const worker = new Worker(workerUrl);
          worker.once("message", (resp: WorkerResponse) => {
            for (const item of resp.results) {
              out[item.idx] = Float64Array.from(item.acc);
            }
            void worker.terminate();
            resolve();
          });
          worker.once("error", (err) => {
            void worker.terminate();
            reject(err);
          });
          const req: WorkerRequest = { task, params, seed, hRel, chunks };
          worker.postMessage(req);
        }),
    ),
  );

  return out;
}

export async function priceEuropeanMC(
  p: OptionParams,
  nPaths: number,
  seed = 42,
  workers?: number,
): Promise<EstimatorResult> {
  validateMC(p, nPaths);
  const chunks = await runParallelChunks("european", p, nPaths, seed, 0.01, workers);
  let sum = 0.0;
  let sumSq = 0.0;
  for (const c of chunks) {
    sum += c[0] ?? 0;
    sumSq += c[1] ?? 0;
  }
  return summarizeMoments(nPaths, sum, sumSq);
}

export async function priceAsianNaive(
  p: OptionParams,
  nPaths: number,
  seed = 42,
  workers?: number,
): Promise<EstimatorResult> {
  validateMC(p, nPaths);
  const chunks = await runParallelChunks("naive", p, nPaths, seed, 0.01, workers);
  let sum = 0.0;
  let sumSq = 0.0;
  for (const c of chunks) {
    sum += c[0] ?? 0;
    sumSq += c[1] ?? 0;
  }
  return summarizeMoments(nPaths, sum, sumSq);
}

export async function priceAsianAntithetic(
  p: OptionParams,
  nPaths: number,
  seed = 42,
  workers?: number,
): Promise<EstimatorResult> {
  if (p.n_steps < 1) {
    throw new Error("n_steps must be >= 1");
  }
  if (nPaths < 4 || nPaths % 2 !== 0) {
    throw new Error("n_paths must be even and >= 4");
  }
  const nPairs = nPaths / 2;
  const chunks = await runParallelChunks("antithetic", p, nPairs, seed, 0.01, workers);
  let sum = 0.0;
  let sumSq = 0.0;
  for (const c of chunks) {
    sum += c[0] ?? 0;
    sumSq += c[1] ?? 0;
  }
  return summarizeMoments(nPairs, sum, sumSq);
}

export async function priceAsianControlVariate(
  p: OptionParams,
  nPaths: number,
  seed = 42,
  workers?: number,
): Promise<ControlVariateResult> {
  validateMC(p, nPaths);
  const yExact = geometricAsianCall(p);
  const chunks = await runParallelChunks("cv", p, nPaths, seed, 0.01, workers);
  let sx = 0.0;
  let sy = 0.0;
  let sxx = 0.0;
  let syy = 0.0;
  let sxy = 0.0;
  for (const c of chunks) {
    sx += c[0] ?? 0;
    sy += c[1] ?? 0;
    sxx += c[2] ?? 0;
    syy += c[3] ?? 0;
    sxy += c[4] ?? 0;
  }
  const n = nPaths;
  const mx = sx / n;
  const my = sy / n;
  const covXY = (sxy - n * mx * my) / (n - 1);
  const varY = (syy - n * my * my) / (n - 1);
  if (varY <= 0.0) {
    throw new Error("var_y is zero; cannot compute control variate beta");
  }
  const beta = covXY / varY;
  const zMean = mx - beta * (my - yExact);
  const varX = (sxx - n * mx * mx) / (n - 1);
  const varZ = Math.max(varX + beta * beta * varY - 2.0 * beta * covXY, 0.0);
  const se = Math.sqrt(varZ / n);
  return {
    price: zMean,
    se,
    ci_low: zMean - 1.96 * se,
    ci_high: zMean + 1.96 * se,
    beta,
    y_exact: yExact,
  };
}

export async function deltaAsianPathwise(
  p: OptionParams,
  nPaths: number,
  seed = 42,
  workers?: number,
): Promise<DeltaResult> {
  validateMC(p, nPaths);
  const chunks = await runParallelChunks("delta_pw", p, nPaths, seed, 0.01, workers);
  let sum = 0.0;
  let sumSq = 0.0;
  for (const c of chunks) {
    sum += c[0] ?? 0;
    sumSq += c[1] ?? 0;
  }
  const res = summarizeMoments(nPaths, sum, sumSq);
  return { delta: res.price, se: res.se, ci_low: res.ci_low, ci_high: res.ci_high };
}

export async function deltaAsianBump(
  p: OptionParams,
  nPaths: number,
  hRel = 0.01,
  seed = 42,
  workers?: number,
): Promise<DeltaResult> {
  validateMC(p, nPaths);
  if (hRel <= 0) {
    throw new Error("h_rel must be > 0");
  }
  const chunks = await runParallelChunks("delta_bump", p, nPaths, seed, hRel, workers);
  let sum = 0.0;
  let sumSq = 0.0;
  for (const c of chunks) {
    sum += c[0] ?? 0;
    sumSq += c[1] ?? 0;
  }
  const res = summarizeMoments(nPaths, sum, sumSq);
  return { delta: res.price, se: res.se, ci_low: res.ci_low, ci_high: res.ci_high };
}
