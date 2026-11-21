import { fitGpdMle } from "./evt.js";
import { gpdPotEs, gpdPotVar } from "./distributions.js";
import { stdTEs, stdTQuantile } from "./distributions.js";
import {
  evaluateRollingBacktestSeries,
  evaluateRegimeChangeSeries,
  kupiecFromCounts,
  christoffersenFromTransitions,
  simulateMarkovChainCounts,
  simulateZ2Single,
  type RollingBacktestSeriesResult,
  type RegimeDetectionMetrics,
} from "./backtest.js";
import {
  fitGarchGauss,
  fitGarchT,
  simulateGarchT,
  TRUE_GARCH_PARAMS,
  type GarchGaussFitResult,
  type GarchFitResult,
} from "./garch.js";
import { estimateRiskMethodsOnSample, type SampleRiskEstimates } from "./risk.js";
import { deriveChunkSeed, Xoshiro256PlusPlus } from "./rng.js";
import { sampleAcf, sampleExcessKurtosis, variance } from "./stats.js";

export type ParallelTaskSpec =
  | {
      kind: "pbh_sample";
      nu: number;
    }
  | {
      kind: "garch_fit_rep";
      nObs: number;
      includeRisk: boolean;
    }
  | {
      kind: "kupiec_h0";
      tLen: number;
      pLevel: number;
      sigLevel: number;
    }
  | {
      kind: "chris_mc";
      tLen: number;
      pLevel: number;
      p11: number | null; // null => iid H0; number => Markov H1
      sigLevel: number;
    }
  | {
      kind: "z2_mc";
      tLen: number;
      pLevels: number[];
      mode: "true" | "gaussian" | "scale90" | "scale80";
    }
  | {
      kind: "rolling_bt";
      totalDays: number;
      windowSize: number;
      stepDays: number;
    }
  | {
      kind: "regime_bt";
      z2CritT1000: number;
      z2CritT3000: number;
      totalDays: number;
      breakDay: number;
      windowSize: number;
      stepDays: number;
    };

export interface GarchRepOutput {
  empVar: number;
  zExcessKurt: number;
  eps2Acf1: number;
  eps2Acf2: number;
  tFit: {
    params: GarchFitResult["params"];
    seInfo: GarchFitResult["seInfo"];
    seSandwich: GarchFitResult["seSandwich"];
    converged: boolean;
    hitBoundary: boolean;
  };
  gFit: {
    params: GarchGaussFitResult["params"];
    converged: boolean;
    hitBoundary: boolean;
  };
  risk?: SampleRiskEstimates;
}

export type ChunkOutput =
  | { kind: "pbh_sample"; samples: Float64Array }
  | { kind: "garch_fit_rep"; items: GarchRepOutput[] }
  | { kind: "kupiec_h0"; rejExact: number; rejAsymp: number; total: number }
  | { kind: "chris_mc"; rejInd: number; rejCc: number; rejKupiecExact: number; total: number }
  | { kind: "z2_mc"; z2Flat: Float64Array; nSeries: number; nLevels: number }
  | { kind: "rolling_bt"; items: RollingBacktestSeriesResult[] }
  | { kind: "regime_bt"; items: RegimeDetectionMetrics[] };

export function executeChunk(
  task: ParallelTaskSpec,
  mainSeed: number,
  salt: number,
  chunkIdx: number,
  len: number,
): ChunkOutput {
  const seed = deriveChunkSeed(mainSeed, salt, chunkIdx);
  const rng = new Xoshiro256PlusPlus(seed);

  switch (task.kind) {
    case "pbh_sample": {
      const samples = new Float64Array(len);
      for (let i = 0; i < len; i += 1) {
        samples[i] = stdTQuantile(rng.nextUniformOpen(), task.nu);
      }
      return { kind: "pbh_sample", samples };
    }

    case "garch_fit_rep": {
      const items: GarchRepOutput[] = new Array(len);
      for (let i = 0; i < len; i += 1) {
        const sim = simulateGarchT(task.nObs, TRUE_GARCH_PARAMS, rng);
        const eps = new Float64Array(task.nObs);
        const eps2 = new Float64Array(task.nObs);
        for (let t = 0; t < task.nObs; t += 1) {
          const e = (sim.returns[t] ?? 0.0) - TRUE_GARCH_PARAMS.mu;
          eps[t] = e;
          eps2[t] = e * e;
        }
        const empVar = variance(eps, 0);
        const zExcessKurt = sampleExcessKurtosis(sim.z);
        const eps2Acf1 = sampleAcf(eps2, 1);
        const eps2Acf2 = sampleAcf(eps2, 2);

        const fitT = fitGarchT(sim.returns, undefined, !task.includeRisk);
        const fitG = fitGarchGauss(sim.returns);

        const entry: GarchRepOutput = {
          empVar,
          zExcessKurt,
          eps2Acf1,
          eps2Acf2,
          tFit: {
            params: fitT.params,
            seInfo: fitT.seInfo,
            seSandwich: fitT.seSandwich,
            converged: fitT.converged,
            hitBoundary: fitT.hitBoundary,
          },
          gFit: {
            params: fitG.params,
            converged: fitG.converged,
            hitBoundary: fitG.hitBoundary,
          },
        };
        if (task.includeRisk) {
          entry.risk = estimateRiskMethodsOnSample(
            fitT.params,
            fitT.residualsHat,
            Math.sqrt(fitT.sigma2NextHat),
            Math.sqrt(sim.sigma2NextTrue),
            TRUE_GARCH_PARAMS,
          );
        }
        items[i] = entry;
      }
      return { kind: "garch_fit_rep", items };
    }

    case "kupiec_h0": {
      const alpha = 1.0 - task.pLevel;
      let rejExact = 0;
      let rejAsymp = 0;
      for (let i = 0; i < len; i += 1) {
        let x = 0;
        for (let t = 0; t < task.tLen; t += 1) {
          if (rng.nextUniform() < alpha) x += 1;
        }
        const k = kupiecFromCounts(x, task.tLen, alpha);
        if (k.pValueExact < task.sigLevel) rejExact += 1;
        if (k.pValueAsymp < task.sigLevel) rejAsymp += 1;
      }
      return { kind: "kupiec_h0", rejExact, rejAsymp, total: len };
    }

    case "chris_mc": {
      const alpha = 1.0 - task.pLevel;
      let rejInd = 0;
      let rejCc = 0;
      let rejKupiecExact = 0;
      for (let i = 0; i < len; i += 1) {
        let n00 = 0;
        let n01 = 0;
        let n10 = 0;
        let n11 = 0;
        let xTotal = 0;
        if (task.p11 === null) {
          let prev = rng.nextUniform() < alpha ? 1 : 0;
          xTotal = prev;
          for (let t = 1; t < task.tLen; t += 1) {
            const cur = rng.nextUniform() < alpha ? 1 : 0;
            xTotal += cur;
            if (prev === 0 && cur === 0) n00 += 1;
            else if (prev === 0 && cur === 1) n01 += 1;
            else if (prev === 1 && cur === 0) n10 += 1;
            else n11 += 1;
            prev = cur;
          }
        } else {
          const mc = simulateMarkovChainCounts(task.tLen, alpha, task.p11, rng);
          n00 = mc.n00;
          n01 = mc.n01;
          n10 = mc.n10;
          n11 = mc.n11;
          xTotal = mc.xTotal;
        }
        const ch = christoffersenFromTransitions(n00, n01, n10, n11, xTotal, task.tLen, alpha);
        const ku = kupiecFromCounts(xTotal, task.tLen, alpha);
        if (ch.pInd < task.sigLevel) rejInd += 1;
        if (ch.pCc < task.sigLevel) rejCc += 1;
        if (ku.pValueExact < task.sigLevel) rejKupiecExact += 1;
      }
      return { kind: "chris_mc", rejInd, rejCc, rejKupiecExact, total: len };
    }

    case "z2_mc": {
      const nLevels = task.pLevels.length;
      const z2Flat = new Float64Array(len * nLevels);
      for (let i = 0; i < len; i += 1) {
        const vals = simulateZ2Single(task.tLen, task.pLevels, rng, task.mode);
        for (let j = 0; j < nLevels; j += 1) {
          z2Flat[i * nLevels + j] = vals[j] ?? 0.0;
        }
      }
      return { kind: "z2_mc", z2Flat, nSeries: len, nLevels };
    }

    case "rolling_bt": {
      const items: RollingBacktestSeriesResult[] = new Array(len);
      for (let i = 0; i < len; i += 1) {
        items[i] = evaluateRollingBacktestSeries(
          rng,
          task.totalDays,
          task.windowSize,
          task.stepDays,
        );
      }
      return { kind: "rolling_bt", items };
    }

    case "regime_bt": {
      const items: RegimeDetectionMetrics[] = new Array(len);
      for (let i = 0; i < len; i += 1) {
        items[i] = evaluateRegimeChangeSeries(
          rng,
          task.z2CritT1000,
          task.z2CritT3000,
          task.totalDays,
          task.breakDay,
          task.windowSize,
          task.stepDays,
        );
      }
      return { kind: "regime_bt", items };
    }
  }
}

// Unused helper imports silenced if not directly referenced
void fitGpdMle;
void gpdPotVar;
void gpdPotEs;
void stdTEs;
