/**
 * Full experiment runner for day003d:
 *   Part 1: Exact distributions (Student-t, ES, GPD) + Pickands-Balkema-de Haan & Hill (N=2,000,000)
 *   Part 2: GARCH(1,1)-t simulation, analytical score, deterministic reconstruction, R=200 t-MLE vs Gauss QMLE
 *   Part 3: Filtered EVT, Conditional VaR/ES comparison across 4 methods (R=200) & threshold sensitivity
 *   Part 4: Backtesting — Kupiec (100k), Christoffersen (20k), Acerbi-Szekely Z2, Rolling (40x6000d), Regime Change (40x6000d)
 */

import { existsSync, mkdirSync, writeFileSync } from "node:fs";
import path from "node:path";
import { performance } from "node:perf_hooks";
import {
  gpdCdf,
  gpdPotEs,
  gpdPotVar,
  gpdQuantile,
  locationScaleTCdf,
  locationScaleTQuantile,
  stdTCdf,
  stdTEs,
  stdTEsNumerical,
  stdTQuantile,
  studentTCdf,
  studentTCdfNumerical,
  studentTQuantile,
} from "./distributions.js";
import { hillEstimator, potFromSample } from "./evt.js";
import {
  exactGarchSqAcf,
  fitGarchT,
  garchTLogLikAndScore,
  reconstructSigma2,
  scoreVsCentralFd,
  simulateGarchT,
  theoreticalGarchMoments,
  TRUE_GARCH_PARAMS,
  type GarchTParams,
} from "./garch.js";
import {
  kupiecExactAnalyticalSize,
  type BacktestMethod,
  type RegimeDetectionMetrics,
  type RollingBacktestSeriesResult,
} from "./backtest.js";
import { serializeJson } from "./json.js";
import { nelderMead } from "./optim.js";
import { type GarchRepOutput } from "./core_tasks.js";
import { runParallelTask } from "./parallel.js";
import { type MethodVarEs, type RiskMethodName } from "./risk.js";
import { Xoshiro256PlusPlus } from "./rng.js";
import { adaptiveSimpson, normalQuantile } from "./special.js";
import {
  batchMeansSe,
  computeHistogram,
  mean,
  median,
  quantile,
  sampleAcf,
  sampleExcessKurtosis,
  sd,
  se,
  summarizeSeries,
  variance,
} from "./stats.js";
import { renderGroupedBarSvg, renderHistogramSvg, renderLineChartSvg } from "./svg.js";

const MAIN_SEED = 20261121;

function findRepoRoot(): string {
  let dir = process.cwd();
  while (true) {
    if (existsSync(path.join(dir, "src", "qb100"))) {
      return dir;
    }
    const parent = path.dirname(dir);
    if (parent === dir) {
      throw new Error("Repository root not found");
    }
    dir = parent;
  }
}

function summarizeReplications(values: number[], trueVal: number) {
  const s = summarizeSeries(values);
  const bias = s.mean - trueVal;
  let sqErr = 0.0;
  for (const v of values) {
    const d = v - trueVal;
    sqErr += d * d;
  }
  const rmse = Math.sqrt(sqErr / values.length);
  return {
    trueVal,
    mean: s.mean,
    sd: s.sd,
    seRep: s.se,
    bias,
    biasSeRatio: bias / s.se,
    rmse,
    q025: s.q025,
    median: s.median,
    q975: s.q975,
  };
}

function computeRelBiasRmse(ests: number[], trues: number[]) {
  const n = ests.length;
  const relErrs = new Array<number>(n);
  const sqRelErrs = new Array<number>(n);
  for (let i = 0; i < n; i += 1) {
    const re = ((ests[i] ?? 0.0) - (trues[i] ?? 1.0)) / (trues[i] ?? 1.0);
    relErrs[i] = re;
    sqRelErrs[i] = re * re;
  }
  const relBias = mean(relErrs);
  const relBiasSe = se(relErrs);
  const mseRel = mean(sqRelErrs);
  const relRmse = Math.sqrt(mseRel);
  const seMse = se(sqRelErrs);
  const relRmseSe = relRmse > 0.0 ? seMse / (2.0 * relRmse) : 0.0;
  return { relBias, relBiasSe, relRmse, relRmseSe, sqRelErrs };
}

async function main(): Promise<void> {
  const tStartAll = performance.now();
  const root = findRepoRoot();
  const outDir = path.join(root, "outputs", "day003d");
  const figDir = path.join(root, "outputs", "figures");
  mkdirSync(outDir, { recursive: true });
  mkdirSync(figDir, { recursive: true });

  // ============================================================================
  // BAGIAN 1: Exact Distributions, ES, GPD & Pickands-Balkema-de Haan (N=2M)
  // ============================================================================
  const tPart1Start = performance.now();

  // 1.1 Student-t CDF & Quantile roundtrip
  let maxCdfNumErr = 0.0;
  let maxQuantileRoundtripErr = 0.0;
  const cdfChecks: Array<{ nu: number; x: number; cdf: number; cdfNum: number; absErr: number }> =
    [];
  for (const nu of [3.5, 5, 10, 30]) {
    for (const x of [-4.0, -1.5, 0.5, 2.2, 4.5]) {
      const c = studentTCdf(x, nu);
      const cn = studentTCdfNumerical(x, nu, 1e-13);
      const err = Math.abs(c - cn);
      if (err > maxCdfNumErr) maxCdfNumErr = err;
      if (x === 2.2) cdfChecks.push({ nu, x, cdf: c, cdfNum: cn, absErr: err });
    }
    for (const p of [0.001, 0.01, 0.05, 0.5, 0.95, 0.975, 0.99, 0.999]) {
      const q = studentTQuantile(p, nu);
      const errRaw = Math.abs(studentTCdf(q, nu) - p);
      if (errRaw > maxQuantileRoundtripErr) maxQuantileRoundtripErr = errRaw;
      const qs = stdTQuantile(p, nu);
      const errStd = Math.abs(stdTCdf(qs, nu) - p);
      if (errStd > maxQuantileRoundtripErr) maxQuantileRoundtripErr = errStd;
    }
  }
  const lsRoundtripErr = Math.abs(
    locationScaleTCdf(locationScaleTQuantile(0.99, 0.0003, 0.01, 5, true), 0.0003, 0.01, 5, true) -
      0.99,
  );

  // 1.2 Exact standardized Student-t ES vs adaptive quadrature
  const esNus = [3.5, 5, 10, 30];
  const esPs = [0.95, 0.975, 0.99, 0.999];
  const esGrid: Array<{
    nu: number;
    p: number;
    esAnalytic: number;
    esNumerical: number;
    absErr: number;
  }> = [];
  let maxEsErr = 0.0;
  for (const nu of esNus) {
    for (const p of esPs) {
      const esAnalytic = stdTEs(p, nu);
      const esNumerical = stdTEsNumerical(p, nu, 1e-11);
      const absErr = Math.abs(esAnalytic - esNumerical);
      if (absErr > maxEsErr) maxEsErr = absErr;
      esGrid.push({ nu, p, esAnalytic, esNumerical, absErr });
    }
  }

  // 1.3 GPD numerical verification
  const gpdChecks: Array<{
    xi: number;
    beta: number;
    u: number;
    cdfRoundtripErr: number;
    esFormula: number;
    esNumerical: number;
    esAbsErr: number;
  }> = [];
  for (const xi of [-0.1, 0.0, 0.2, 0.4]) {
    const beta = 0.65;
    const u = 1.2;
    const pExceed = 0.1;
    const q = gpdQuantile(0.9, xi, beta);
    const cdfRoundtripErr = Math.abs(gpdCdf(q, xi, beta) - 0.9);
    const var99 = gpdPotVar(u, pExceed, xi, beta, 0.99);
    const esFormula = gpdPotEs(u, var99, xi, beta);
    const condSurv = (1.0 - 0.99) / pExceed;
    const esNumerical =
      u +
      adaptiveSimpson(
        (v) => {
          if (v <= 0.0) return 0.0;
          const v5 = v ** 5;
          const tp = condSurv * v5 * v;
          if (tp <= 1e-300) return 0.0;
          const yv =
            Math.abs(xi) < 1e-8 ? -beta * Math.log(tp) : (beta / xi) * (Math.pow(tp, -xi) - 1.0);
          return 6.0 * v5 * yv;
        },
        0.0,
        1.0,
        1e-11,
        16,
      );
    gpdChecks.push({
      xi,
      beta,
      u,
      cdfRoundtripErr,
      esFormula,
      esNumerical,
      esAbsErr: Math.abs(esFormula - esNumerical),
    });
  }

  // 1.4 Pickands-Balkema-de Haan: N = 2,000,000 exact inversion draws from std t(5)
  const nPbh = 2_000_000;
  const pbhParallel = await runParallelTask(
    { kind: "pbh_sample", nu: 5.0 },
    nPbh,
    200_000,
    MAIN_SEED,
    11,
  );
  const workerCount = pbhParallel.workerCount;
  const pbhSamples = new Float64Array(nPbh);
  let pos = 0;
  for (const out of pbhParallel.outputs) {
    if (out.kind === "pbh_sample") {
      pbhSamples.set(out.samples, pos);
      pos += out.samples.length;
    }
  }
  pbhSamples.sort();

  const trueXiPbh = 1.0 / 5.0; // 0.2
  const qGridPbh = [0.9, 0.95, 0.99, 0.999];
  const pbhResults: Array<{
    q: number;
    u: number;
    kExceed: number;
    xiPot: number;
    seXiPot: number;
    biasPot: number;
    betaPot: number;
    xiHill: number;
    seXiHill: number;
    biasHill: number;
    var999Pot: number;
    var999True: number;
    es999Pot: number;
    es999True: number;
  }> = [];

  // Descending positive tail above u_0.9 > 0 for fast Hill evaluation
  const pbhDesc = new Float64Array(nPbh);
  for (let i = 0; i < nPbh; i += 1) {
    pbhDesc[i] = pbhSamples[nPbh - 1 - i] ?? 0.0;
  }

  const var999TrueT5 = stdTQuantile(0.999, 5.0);
  const es999TrueT5 = stdTEs(0.999, 5.0);

  for (const q of qGridPbh) {
    const pot = potFromSample(pbhSamples, q, [0.999], [0.999], true);
    const hill = hillEstimator(pbhDesc, pot.nExceed, true);
    pbhResults.push({
      q,
      u: pot.u,
      kExceed: pot.nExceed,
      xiPot: pot.gpd.xi,
      seXiPot: pot.gpd.seXi,
      biasPot: pot.gpd.xi - trueXiPbh,
      betaPot: pot.gpd.beta,
      xiHill: hill.xiHat,
      seXiHill: hill.se,
      biasHill: hill.xiHat - trueXiPbh,
      var999Pot: pot.vars["0.999"] ?? NaN,
      var999True: var999TrueT5,
      es999Pot: pot.ess["0.999"] ?? NaN,
      es999True: es999TrueT5,
    });
  }

  // Build SVG 1 (t(5) vs GPD q=0.90 tail quantiles & ES) and SVG 2 (xi_POT & xi_Hill vs q)
  const pot90ForPlot = potFromSample(
    pbhSamples,
    0.9,
    [0.95, 0.975, 0.99, 0.995, 0.999],
    [0.95, 0.975, 0.99, 0.995, 0.999],
    true,
  );
  const tailPs = [0.95, 0.975, 0.99, 0.995, 0.999];
  const tailLogInvP = tailPs.map((p) => -Math.log10(1.0 - p));
  writeFileSync(
    path.join(figDir, "day003d_t-vs-gpd-tail.svg"),
    renderLineChartSvg({
      title: "Ekor Kuantil & ES: t_std(5) Eksak vs Aproksimasi GPD POT (q=0.90)",
      subtitle: "N = 2.000.000 sampel inversi t_std(5) | sumbu-x = -log10(1 - p)",
      xLabel: "-log10(1 - p)  [1.30 = 95%, 2.00 = 99%, 3.00 = 99.9%]",
      yLabel: "Nilai Risiko Ekor (Skala Baku)",
      series: [
        {
          label: "t(5) VaR Eksak",
          color: "#1d4ed8",
          x: tailLogInvP,
          y: tailPs.map((p) => stdTQuantile(p, 5)),
        },
        {
          label: "GPD POT VaR (u=q90)",
          color: "#0ea5e9",
          dashed: true,
          x: tailLogInvP,
          y: tailPs.map((p) => pot90ForPlot.vars[p.toString()] ?? 0),
        },
        {
          label: "t(5) ES Eksak",
          color: "#b91c1c",
          x: tailLogInvP,
          y: tailPs.map((p) => stdTEs(p, 5)),
        },
        {
          label: "GPD POT ES (u=q90)",
          color: "#f97316",
          dashed: true,
          x: tailLogInvP,
          y: tailPs.map((p) => pot90ForPlot.ess[p.toString()] ?? 0),
        },
      ],
    }),
  );

  writeFileSync(
    path.join(figDir, "day003d_xi-and-hill-vs-q.svg"),
    renderLineChartSvg({
      title: "Konvergensi Indeks Ekor Pickands–Balkema–de Haan: POT vs Hill",
      subtitle: "N = 2.000.000 sampel t_std(5) pada ambang kuantil q in {0.90, 0.95, 0.99, 0.999}",
      xLabel: "Kuantil Ambang q",
      yLabel: "Estimasi Indeks Ekor xi_hat",
      hLine: { y: trueXiPbh, label: "Target 1/nu = 0.200", color: "#111827" },
      series: [
        {
          label: "POT GPD MLE (± 1.96 SE)",
          color: "#2563eb",
          x: pbhResults.map((r) => r.q),
          y: pbhResults.map((r) => r.xiPot),
          errorBars: pbhResults.map((r) => 1.96 * r.seXiPot),
        },
        {
          label: "Estimator Hill (± 1.96 SE)",
          color: "#dc2626",
          x: pbhResults.map((r) => r.q),
          y: pbhResults.map((r) => r.xiHill),
          errorBars: pbhResults.map((r) => 1.96 * r.seXiHill),
        },
      ],
    }),
  );

  const part1TimeSec = (performance.now() - tPart1Start) / 1000.0;

  // ============================================================================
  // BAGIAN 2: GARCH(1,1)-t Simulation, Score Verification & R=200 Study
  // ============================================================================
  const tPart2Start = performance.now();

  const theoTrue = theoreticalGarchMoments(TRUE_GARCH_PARAMS);
  // Long stationary trajectory N=200,000 for batch-means moment verification
  const rngLong = new Xoshiro256PlusPlus(BigInt(MAIN_SEED) ^ 0x2222n);
  const longSim = simulateGarchT(200_000, TRUE_GARCH_PARAMS, rngLong);
  const longEps2 = new Float64Array(200_000);
  for (let t = 0; t < 200_000; t += 1) {
    const e = (longSim.returns[t] ?? 0.0) - TRUE_GARCH_PARAMS.mu;
    longEps2[t] = e * e;
  }
  const longVar = mean(longEps2);
  const longVarBatchSe = batchMeansSe(longEps2, 50);
  const longZKurt = sampleExcessKurtosis(longSim.z);
  const longAcf1 = sampleAcf(longEps2, 1);
  const longAcf2 = sampleAcf(longEps2, 2);

  // Finite-4th-moment ACF anchor: t(12), alpha=0.10, beta=0.80 => D = 0.1625 > 0
  const anchorParams: GarchTParams = { mu: 0.0003, omega: 1e-5, alpha: 0.1, beta: 0.8, nu: 12.0 };
  const anchorTheo = exactGarchSqAcf(
    anchorParams.alpha,
    anchorParams.beta,
    3.0 + 6.0 / (12.0 - 4.0),
  );
  const rngAnchor = new Xoshiro256PlusPlus(BigInt(MAIN_SEED) ^ 0x3333n);
  const anchorVars: number[] = [];
  const anchorKurts: number[] = [];
  const anchorAcf1: number[] = [];
  const anchorAcf2: number[] = [];
  for (let rep = 0; rep < 60; rep += 1) {
    const s = simulateGarchT(5000, anchorParams, rngAnchor);
    const eps = new Float64Array(5000);
    const e2 = new Float64Array(5000);
    for (let t = 0; t < 5000; t += 1) {
      const e = (s.returns[t] ?? 0.0) - anchorParams.mu;
      eps[t] = e;
      e2[t] = e * e;
    }
    anchorVars.push(variance(eps, 0));
    anchorKurts.push(sampleExcessKurtosis(s.z));
    anchorAcf1.push(sampleAcf(e2, 1));
    anchorAcf2.push(sampleAcf(e2, 2));
  }

  // Deterministic sigma2 reconstruction + analytical score vs FD check + Nelder-Mead check on a n=2500 sample
  const rngCheck = new Xoshiro256PlusPlus(BigInt(MAIN_SEED) ^ 0x4444n);
  const simCheck = simulateGarchT(2500, TRUE_GARCH_PARAMS, rngCheck);
  const recCheck = reconstructSigma2(simCheck.returns, TRUE_GARCH_PARAMS);
  let maxSigma2ReconDiff = 0.0;
  for (let t = 0; t < 2500; t += 1) {
    const d = Math.abs((recCheck.sigma2[t] ?? 0.0) - (simCheck.sigma2[t] ?? 0.0));
    if (d > maxSigma2ReconDiff) maxSigma2ReconDiff = d;
  }
  const scoreCheckTrue = scoreVsCentralFd(simCheck.returns, TRUE_GARCH_PARAMS);
  const scoreCheckAlt = scoreVsCentralFd(simCheck.returns, {
    mu: -0.0001,
    omega: 4e-6,
    alpha: 0.11,
    beta: 0.85,
    nu: 7.0,
  });
  const fitCheck = fitGarchT(simCheck.returns, undefined, true);
  const nmCheck = nelderMead(
    (theta) => {
      const p: GarchTParams = {
        mu: theta[0] ?? 0,
        omega: theta[1] ?? 1e-6,
        alpha: theta[2] ?? 0.08,
        beta: theta[3] ?? 0.9,
        nu: theta[4] ?? 5,
      };
      if (p.omega <= 0 || p.alpha < 0 || p.beta < 0 || p.alpha + p.beta >= 0.9999 || p.nu <= 2.01) {
        return 1e12;
      }
      return -garchTLogLikAndScore(simCheck.returns, p).logLik / 2500;
    },
    [
      fitCheck.params.mu,
      fitCheck.params.omega,
      fitCheck.params.alpha,
      fitCheck.params.beta,
      fitCheck.params.nu,
    ],
    500,
    1e-11,
  );
  const bfgsVsNmLogLikDiff = Math.abs(fitCheck.logLik / 2500 + nmCheck.value);

  // R = 200 replication study (n = 2500) with both Info & Sandwich SEs for Bagian 2
  const nRepGarch = 200;
  const garchSeRes = await runParallelTask(
    { kind: "garch_fit_rep", nObs: 2500, includeRisk: false },
    nRepGarch,
    10,
    MAIN_SEED,
    29,
  );
  const seRepItems: GarchRepOutput[] = [];
  for (const out of garchSeRes.outputs) {
    if (out.kind === "garch_fit_rep") {
      seRepItems.push(...out.items);
    }
  }

  const paramKeys = ["mu", "omega", "alpha", "beta", "nu"] as const;
  const tMleSummary: Record<string, unknown> = {};
  const zCrit95 = normalQuantile(0.975);
  let tConvCount = 0;
  let tBoundaryHits = 0;

  for (const item of seRepItems) {
    if (item.tFit.converged) tConvCount += 1;
    if (item.tFit.hitBoundary) tBoundaryHits += 1;
  }

  for (const k of paramKeys) {
    const vals = seRepItems.map((it) => it.tFit.params[k]);
    const tv = TRUE_GARCH_PARAMS[k];
    const summ = summarizeReplications(vals, tv);
    let covInfo = 0;
    let covSand = 0;
    let validInfo = 0;
    let validSand = 0;
    const meanSeInfoArr: number[] = [];
    const meanSeSandArr: number[] = [];
    for (const it of seRepItems) {
      const est = it.tFit.params[k];
      const sI = it.tFit.seInfo[k];
      const sS = it.tFit.seSandwich[k];
      if (Number.isFinite(sI) && sI > 0) {
        validInfo += 1;
        meanSeInfoArr.push(sI);
        if (Math.abs(est - tv) <= zCrit95 * sI) covInfo += 1;
      }
      if (Number.isFinite(sS) && sS > 0) {
        validSand += 1;
        meanSeSandArr.push(sS);
        if (Math.abs(est - tv) <= zCrit95 * sS) covSand += 1;
      }
    }
    tMleSummary[k] = {
      ...summ,
      meanSeInfo: mean(meanSeInfoArr),
      meanSeSandwich: mean(meanSeSandArr),
      covInfo95: validInfo > 0 ? covInfo / validInfo : NaN,
      covSandwich95: validSand > 0 ? covSand / validSand : NaN,
      validInfoCount: validInfo,
      validSandwichCount: validSand,
    };
  }

  // Gaussian QMLE vs t-MLE comparison on the exact same seRepItems
  const qmleKeys = ["mu", "omega", "alpha", "beta"] as const;
  const qmleSummary: Record<string, unknown> = {};
  let gConvCount = 0;
  let gBoundaryHits = 0;
  for (const item of seRepItems) {
    if (item.gFit.converged) gConvCount += 1;
    if (item.gFit.hitBoundary) gBoundaryHits += 1;
  }
  for (const k of qmleKeys) {
    const valsG = seRepItems.map((it) => it.gFit.params[k]);
    const valsT = seRepItems.map((it) => it.tFit.params[k]);
    const tv = TRUE_GARCH_PARAMS[k];
    const sg = summarizeReplications(valsG, tv);
    const st = summarizeReplications(valsT, tv);
    qmleSummary[k] = {
      ...sg,
      relEfficiencyMseRatio: (sg.rmse * sg.rmse) / (st.rmse * st.rmse),
    };
  }

  // Rep moments across 200 series (n=2500)
  const repEmpVar = seRepItems.map((it) => it.empVar);
  const repZKurt = seRepItems.map((it) => it.zExcessKurt);
  const repAcf1 = seRepItems.map((it) => it.eps2Acf1);
  const repAcf2 = seRepItems.map((it) => it.eps2Acf2);

  // SVG 3: GARCH parameter estimate histograms (alpha_hat t-MLE vs QMLE)
  const tAlphas = seRepItems.map((it) => it.tFit.params.alpha);
  const gAlphas = seRepItems.map((it) => it.gFit.params.alpha);
  const minA = Math.min(...tAlphas, ...gAlphas, 0.03);
  const maxA = Math.max(...tAlphas, ...gAlphas, 0.14);
  const histT = computeHistogram(tAlphas, 18, minA, maxA);
  const histG = computeHistogram(gAlphas, 18, minA, maxA);
  writeFileSync(
    path.join(figDir, "day003d_garch-mle-dist.svg"),
    renderHistogramSvg({
      title: "Distribusi Estimasi alpha pada R = 200 Replikasi GARCH(1,1)-t (n = 2500)",
      subtitle: "Perbandingan Efisiensi t-MLE (biru) vs Gaussian QMLE (merah)",
      xLabel: "Estimasi alpha_hat",
      edges: histT.edges,
      trueValue: TRUE_GARCH_PARAMS.alpha,
      series: [
        { label: `t-MLE (SD=${sd(tAlphas).toFixed(4)})`, density: histT.density, color: "#2563eb" },
        {
          label: `Gauss QMLE (SD=${sd(gAlphas).toFixed(4)})`,
          density: histG.density,
          color: "#dc2626",
        },
      ],
    }),
  );

  const part2TimeSec = (performance.now() - tPart2Start) / 1000.0;

  // ============================================================================
  // BAGIAN 3: Filtered EVT, Conditional VaR & ES (R=200)
  // ============================================================================
  const tPart3Start = performance.now();

  const garchRepRes = await runParallelTask(
    { kind: "garch_fit_rep", nObs: 2500, includeRisk: true },
    nRepGarch,
    10,
    MAIN_SEED,
    23,
  );
  const repItems: GarchRepOutput[] = [];
  for (const out of garchRepRes.outputs) {
    if (out.kind === "garch_fit_rep") {
      repItems.push(...out.items);
    }
  }

  const riskMethods: RiskMethodName[] = ["gaussian", "t_parametric", "fhs", "fhs_evt"];
  const metricKeys: Array<keyof MethodVarEs> = ["var95", "var99", "var999", "es975", "es99"];

  const riskSummary: Record<string, Record<string, unknown>> = {};
  for (const m of riskMethods) {
    riskSummary[m] = {};
    for (const mk of metricKeys) {
      const ests = repItems.map((it) => it.risk?.methods[m][mk] ?? 0.0);
      const trues = repItems.map((it) => it.risk?.trueValues[mk] ?? 1.0);
      riskSummary[m][mk] = computeRelBiasRmse(ests, trues);
    }
  }

  // Qualitative paired test: MSE(FHS) - MSE(EVT) and MSE(FHS) - MSE(t_parametric) on VaR_0.999
  const trues999 = repItems.map((it) => it.risk?.trueValues.var999 ?? 1.0);
  const fhs999 = repItems.map((it) => it.risk?.methods.fhs.var999 ?? 0.0);
  const evt999 = repItems.map((it) => it.risk?.methods.fhs_evt.var999 ?? 0.0);
  const t999 = repItems.map((it) => it.risk?.methods.t_parametric.var999 ?? 0.0);
  const gauss999 = repItems.map((it) => it.risk?.methods.gaussian.var999 ?? 0.0);

  const sqFhs = fhs999.map((v, i) => ((v - (trues999[i] ?? 1)) / (trues999[i] ?? 1)) ** 2);
  const sqEvt = evt999.map((v, i) => ((v - (trues999[i] ?? 1)) / (trues999[i] ?? 1)) ** 2);
  const sqT = t999.map((v, i) => ((v - (trues999[i] ?? 1)) / (trues999[i] ?? 1)) ** 2);
  const diffMseFhsMinusEvt = sqFhs.map((v, i) => v - (sqEvt[i] ?? 0));
  const diffMseFhsMinusT = sqFhs.map((v, i) => v - (sqT[i] ?? 0));
  const gaussRelBias999 = computeRelBiasRmse(gauss999, trues999);

  // Threshold sensitivity q in {0.85, 0.90, 0.95} & xi_hat vs 1/nu_hat
  const thresholdSens: Array<Record<string, unknown>> = [];
  for (let ti = 0; ti < 3; ti += 1) {
    const qVal = repItems[0]?.risk?.thresholdSensitivity[ti]?.q ?? 0.9;
    const xis = repItems.map((it) => it.risk?.thresholdSensitivity[ti]?.xiHat ?? 0.0);
    const betas = repItems.map((it) => it.risk?.thresholdSensitivity[ti]?.betaHat ?? 0.0);
    const v999 = repItems.map((it) => it.risk?.thresholdSensitivity[ti]?.var999 ?? 0.0);
    const e99 = repItems.map((it) => it.risk?.thresholdSensitivity[ti]?.es99 ?? 0.0);
    const t99 = repItems.map((it) => it.risk?.trueValues.es99 ?? 1.0);
    const rV999 = computeRelBiasRmse(v999, trues999);
    const rE99 = computeRelBiasRmse(e99, t99);
    thresholdSens.push({
      q: qVal,
      xi: summarizeReplications(xis, 0.2),
      meanBeta: mean(betas),
      var999: {
        relBias: rV999.relBias,
        relBiasSe: rV999.relBiasSe,
        relRmse: rV999.relRmse,
        relRmseSe: rV999.relRmseSe,
      },
      es99: {
        relBias: rE99.relBias,
        relBiasSe: rE99.relBiasSe,
        relRmse: rE99.relRmse,
        relRmseSe: rE99.relRmseSe,
      },
    });
  }
  const xiHatVsInvNu = {
    xiHat90: summarizeSeries(repItems.map((it) => it.risk?.evtXiHat ?? 0.0)),
    invNuHat: summarizeSeries(repItems.map((it) => it.risk?.invNuHat ?? 0.0)),
    diffXiMinusInvNu: summarizeSeries(
      repItems.map((it) => (it.risk?.evtXiHat ?? 0.0) - (it.risk?.invNuHat ?? 0.0)),
    ),
  };

  // Strip sqRelErrs before putting riskSummary into results.json
  const riskSummaryClean: Record<string, Record<string, unknown>> = {};
  for (const m of riskMethods) {
    riskSummaryClean[m] = {};
    for (const mk of metricKeys) {
      const r = riskSummary[m]?.[mk] as {
        relBias: number;
        relBiasSe: number;
        relRmse: number;
        relRmseSe: number;
      };
      riskSummaryClean[m][mk] = {
        relBias: r.relBias,
        relBiasSe: r.relBiasSe,
        relRmse: r.relRmse,
        relRmseSe: r.relRmseSe,
      };
    }
  }

  // SVG 4: VaR 0.999 & ES 0.99 relative bias & RMSE per method
  writeFileSync(
    path.join(figDir, "day003d_var-bias-rmse.svg"),
    renderGroupedBarSvg({
      title: "Akurasi VaR 99.9% & ES 99% Bersyarat per Metode (R = 200, n = 2500)",
      subtitle: "Bias relatif vs RMSE relatif terhadap nilai benar GARCH-t(5)",
      yLabel: "Kesalahan Relatif",
      groups: [
        {
          category: "Gaussian",
          values: [
            {
              label: "Bias Rel VaR 99.9%",
              value: (riskSummaryClean.gaussian?.var999 as { relBias: number }).relBias,
              color: "#dc2626",
            },
            {
              label: "RMSE Rel VaR 99.9%",
              value: (riskSummaryClean.gaussian?.var999 as { relRmse: number }).relRmse,
              color: "#2563eb",
            },
            {
              label: "RMSE Rel ES 99%",
              value: (riskSummaryClean.gaussian?.es99 as { relRmse: number }).relRmse,
              color: "#059669",
            },
          ],
        },
        {
          category: "t-Parametrik",
          values: [
            {
              label: "Bias Rel VaR 99.9%",
              value: (riskSummaryClean.t_parametric?.var999 as { relBias: number }).relBias,
              color: "#dc2626",
            },
            {
              label: "RMSE Rel VaR 99.9%",
              value: (riskSummaryClean.t_parametric?.var999 as { relRmse: number }).relRmse,
              color: "#2563eb",
            },
            {
              label: "RMSE Rel ES 99%",
              value: (riskSummaryClean.t_parametric?.es99 as { relRmse: number }).relRmse,
              color: "#059669",
            },
          ],
        },
        {
          category: "FHS Empiris",
          values: [
            {
              label: "Bias Rel VaR 99.9%",
              value: (riskSummaryClean.fhs?.var999 as { relBias: number }).relBias,
              color: "#dc2626",
            },
            {
              label: "RMSE Rel VaR 99.9%",
              value: (riskSummaryClean.fhs?.var999 as { relRmse: number }).relRmse,
              color: "#2563eb",
            },
            {
              label: "RMSE Rel ES 99%",
              value: (riskSummaryClean.fhs?.es99 as { relRmse: number }).relRmse,
              color: "#059669",
            },
          ],
        },
        {
          category: "FHS-EVT (q90)",
          values: [
            {
              label: "Bias Rel VaR 99.9%",
              value: (riskSummaryClean.fhs_evt?.var999 as { relBias: number }).relBias,
              color: "#dc2626",
            },
            {
              label: "RMSE Rel VaR 99.9%",
              value: (riskSummaryClean.fhs_evt?.var999 as { relRmse: number }).relRmse,
              color: "#2563eb",
            },
            {
              label: "RMSE Rel ES 99%",
              value: (riskSummaryClean.fhs_evt?.es99 as { relRmse: number }).relRmse,
              color: "#059669",
            },
          ],
        },
      ],
    }),
  );

  // SVG 5: Q-Q plot of GARCH-t standardized residuals (200 evenly spaced quantiles) vs t_std(nu_hat) & Normal
  const resSorted = Array.from(fitCheck.residualsHat).sort((a, b) => a - b);
  const qqGridPs = Array.from({ length: 180 }, (_, idx) => 0.002 + (0.996 * idx) / 179);
  const qqEmp = qqGridPs.map((p) => quantile(resSorted, p, true));
  const qqT = qqGridPs.map((p) => stdTQuantile(p, fitCheck.params.nu));
  const qqN = qqGridPs.map((p) => normalQuantile(p));
  writeFileSync(
    path.join(figDir, "day003d_residuals-qq.svg"),
    renderLineChartSvg({
      title: "Plot Q-Q Residu Terstandar GARCH(1,1)-t (n = 2500)",
      subtitle: `nu_hat = ${fitCheck.params.nu.toFixed(3)} — perbandingan garis referensi y = x terhadap t_std(nu_hat) & Normal`,
      xLabel: "Kuantil Teoretis Model",
      yLabel: "Kuantil Empiris z_hat",
      series: [
        {
          label: `vs t_std(nu_hat=${fitCheck.params.nu.toFixed(2)})`,
          color: "#2563eb",
          x: qqT,
          y: qqEmp,
        },
        { label: "vs N(0, 1) (ekor menyimpang)", color: "#dc2626", x: qqN, y: qqEmp },
        { label: "Garis Identitas y = x", color: "#111827", dashed: true, x: [-6, 6], y: [-6, 6] },
      ],
    }),
  );

  const part3TimeSec = (performance.now() - tPart3Start) / 1000.0;

  // ============================================================================
  // BAGIAN 4: Backtesting (Kupiec 100k, Christoffersen, Acerbi-Szekely Z2, Rolling, Regime)
  // ============================================================================
  const tPart4Start = performance.now();

  // 4.1 Kupiec exact & asymp size on 100,000 H0 series (T=1000) for p in {0.95, 0.99}
  const kupiecResults: Record<string, unknown> = {};
  for (const pLevel of [0.95, 0.99]) {
    const nSeries = 100_000;
    const tLen = 1000;
    const alpha = 1.0 - pLevel;
    const anal = kupiecExactAnalyticalSize(tLen, alpha, 0.05);
    const res = await runParallelTask(
      { kind: "kupiec_h0", tLen, pLevel, sigLevel: 0.05 },
      nSeries,
      5000,
      MAIN_SEED,
      pLevel === 0.95 ? 41 : 42,
    );
    let rejExact = 0;
    let rejAsymp = 0;
    for (const out of res.outputs) {
      if (out.kind === "kupiec_h0") {
        rejExact += out.rejExact;
        rejAsymp += out.rejAsymp;
      }
    }
    const rateExact = rejExact / nSeries;
    const rateAsymp = rejAsymp / nSeries;
    const binomSe = Math.sqrt((0.05 * 0.95) / nSeries);
    kupiecResults[pLevel.toString()] = {
      nSeries,
      tLen,
      alpha,
      rateExact,
      rateAsymp,
      binomSeFrom05: binomSe,
      zExactVs05: (rateExact - 0.05) / binomSe,
      zAsympVs05: (rateAsymp - 0.05) / binomSe,
      analyticalDiscreteSizeExact: anal.exactSize,
      analyticalDiscreteSizeAsymp: anal.asympSize,
      zExactVsDiscreteTheory:
        (rateExact - anal.exactSize) / Math.sqrt((anal.exactSize * (1 - anal.exactSize)) / nSeries),
    };
  }

  // 4.2 Christoffersen size (iid H0) & power (Markov P(1|1)=0.30) on R=20,000, T=1000
  const chrisResults: Record<string, unknown> = {};
  for (const pLevel of [0.95, 0.99]) {
    const nSeries = 20_000;
    const tLen = 1000;
    const h0 = await runParallelTask(
      { kind: "chris_mc", tLen, pLevel, p11: null, sigLevel: 0.05 },
      nSeries,
      2000,
      MAIN_SEED,
      pLevel === 0.95 ? 43 : 44,
    );
    const h1 = await runParallelTask(
      { kind: "chris_mc", tLen, pLevel, p11: 0.3, sigLevel: 0.05 },
      nSeries,
      2000,
      MAIN_SEED,
      pLevel === 0.95 ? 45 : 46,
    );
    const agg = (outs: typeof h0.outputs) => {
      let rInd = 0;
      let rCc = 0;
      let rKu = 0;
      for (const o of outs) {
        if (o.kind === "chris_mc") {
          rInd += o.rejInd;
          rCc += o.rejCc;
          rKu += o.rejKupiecExact;
        }
      }
      return { indRate: rInd / nSeries, ccRate: rCc / nSeries, kupiecRate: rKu / nSeries };
    };
    chrisResults[pLevel.toString()] = {
      nSeries,
      tLen,
      h0Size: agg(h0.outputs),
      h1PowerMarkov03: agg(h1.outputs),
    };
  }

  // 4.3 Acerbi-Szekely Z2 null distribution (simulated from true GARCH-t(5) model!) + size & power
  const z2Levels = [0.975, 0.99];
  const extractZ2 = (outs: Awaited<ReturnType<typeof runParallelTask>>["outputs"]) => {
    const vals: number[][] = [[], []];
    for (const o of outs) {
      if (o.kind === "z2_mc") {
        for (let i = 0; i < o.nSeries; i += 1) {
          vals[0]?.push(o.z2Flat[i * 2] ?? 0.0);
          vals[1]?.push(o.z2Flat[i * 2 + 1] ?? 0.0);
        }
      }
    }
    return vals;
  };

  const nZ2Null = 30_000;
  const nZ2Eval = 10_000;
  const z2Null1000 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 1000, pLevels: z2Levels, mode: "true" },
        nZ2Null,
        1500,
        MAIN_SEED,
        47,
      )
    ).outputs,
  );
  const z2Null3000 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 3000, pLevels: z2Levels, mode: "true" },
        10_000,
        500,
        MAIN_SEED,
        48,
      )
    ).outputs,
  );
  const z2Null5000 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 5000, pLevels: z2Levels, mode: "true" },
        10_000,
        500,
        MAIN_SEED,
        49,
      )
    ).outputs,
  );

  const crit1000_975 = quantile(z2Null1000[0] ?? [], 0.05);
  const crit1000_99 = quantile(z2Null1000[1] ?? [], 0.05);
  const crit3000_99 = quantile(z2Null3000[1] ?? [], 0.05);
  const crit5000_975 = quantile(z2Null5000[0] ?? [], 0.05);
  const crit5000_99 = quantile(z2Null5000[1] ?? [], 0.05);

  const z2FreshH0 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 1000, pLevels: z2Levels, mode: "true" },
        nZ2Eval,
        1000,
        MAIN_SEED,
        50,
      )
    ).outputs,
  );
  const z2Gauss = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 1000, pLevels: z2Levels, mode: "gaussian" },
        nZ2Eval,
        1000,
        MAIN_SEED,
        51,
      )
    ).outputs,
  );
  const z2Scale90 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 1000, pLevels: z2Levels, mode: "scale90" },
        nZ2Eval,
        1000,
        MAIN_SEED,
        52,
      )
    ).outputs,
  );
  const z2Scale80 = extractZ2(
    (
      await runParallelTask(
        { kind: "z2_mc", tLen: 1000, pLevels: z2Levels, mode: "scale80" },
        nZ2Eval,
        1000,
        MAIN_SEED,
        53,
      )
    ).outputs,
  );

  const summarizeZ2 = (arr: number[], crit: number) => {
    let rej = 0;
    for (const v of arr) {
      if (v < crit) rej += 1;
    }
    return {
      mean: mean(arr),
      se: se(arr),
      q05: quantile(arr, 0.05),
      rejectRateAtCrit05: rej / arr.length,
    };
  };

  const z2Summary = {
    nNullSimT1000: nZ2Null,
    nEvalSim: nZ2Eval,
    crit05: {
      t1000_975: crit1000_975,
      t1000_99: crit1000_99,
      t3000_99: crit3000_99,
      t5000_975: crit5000_975,
      t5000_99: crit5000_99,
    },
    nullDist1000: {
      p975: summarizeSeries(z2Null1000[0] ?? []),
      p99: summarizeSeries(z2Null1000[1] ?? []),
    },
    freshH0Size: {
      p975: summarizeZ2(z2FreshH0[0] ?? [], crit1000_975),
      p99: summarizeZ2(z2FreshH0[1] ?? [], crit1000_99),
    },
    powerGaussianModel: {
      p975: summarizeZ2(z2Gauss[0] ?? [], crit1000_975),
      p99: summarizeZ2(z2Gauss[1] ?? [], crit1000_99),
    },
    powerScale90Es: {
      p975: summarizeZ2(z2Scale90[0] ?? [], crit1000_975),
      p99: summarizeZ2(z2Scale90[1] ?? [], crit1000_99),
    },
    powerScale80Es: {
      p975: summarizeZ2(z2Scale80[0] ?? [], crit1000_975),
      p99: summarizeZ2(z2Scale80[1] ?? [], crit1000_99),
    },
  };

  // 4.4 Rolling backtest: R = 40 independent 6000-day series, W=1000, step=20
  const nRollingSeries = 40;
  const rollingParallel = await runParallelTask(
    { kind: "rolling_bt", totalDays: 6000, windowSize: 1000, stepDays: 20 },
    nRollingSeries,
    2,
    MAIN_SEED,
    61,
  );
  const rollingItems: RollingBacktestSeriesResult[] = [];
  for (const out of rollingParallel.outputs) {
    if (out.kind === "rolling_bt") rollingItems.push(...out.items);
  }

  const btMethods: BacktestMethod[] = [
    "true_control",
    "gaussian_garch",
    "t_garch",
    "fhs",
    "fhs_evt",
  ];
  const rollingSummary: Record<string, unknown> = {};
  for (const m of btMethods) {
    const r99 = rollingItems.map((it) => it.methods[m].rate99);
    const r975 = rollingItems.map((it) => it.methods[m].rate975);
    const kpAsymp = rollingItems.map((it) => it.methods[m].kupiec99AsympP);
    const kpExact = rollingItems.map((it) => it.methods[m].kupiec99ExactP);
    const chInd = rollingItems.map((it) => it.methods[m].chris99IndP);
    const chCc = rollingItems.map((it) => it.methods[m].chris99CcP);
    const z99 = rollingItems.map((it) => it.methods[m].z2_99);
    const z975 = rollingItems.map((it) => it.methods[m].z2_975);

    const fracBelow05 = (arr: number[]): number => arr.filter((v) => v < 0.05).length / arr.length;
    const z2Rej99 = z99.filter((v) => v < crit5000_99).length / z99.length;
    const z2Rej975 = z975.filter((v) => v < crit5000_975).length / z975.length;

    rollingSummary[m] = {
      rate99Mean: mean(r99),
      rate99SeAcrossSeries: se(r99),
      rate99PooledBinomSe: Math.sqrt((0.01 * 0.99) / (nRollingSeries * 5000)),
      rate975Mean: mean(r975),
      rate975SeAcrossSeries: se(r975),
      kupiec99ExactMedianP: median(kpExact),
      kupiec99ExactRejectFrac: fracBelow05(kpExact),
      kupiec99AsympMedianP: median(kpAsymp),
      kupiec99AsympRejectFrac: fracBelow05(kpAsymp),
      chris99IndMedianP: median(chInd),
      chris99IndRejectFrac: fracBelow05(chInd),
      chris99CcMedianP: median(chCc),
      chris99CcRejectFrac: fracBelow05(chCc),
      z2_99Mean: mean(z99),
      z2_99Se: se(z99),
      z2_99RejectFracAtSimCrit05: z2Rej99,
      z2_975Mean: mean(z975),
      z2_975Se: se(z975),
      z2_975RejectFracAtSimCrit05: z2Rej975,
    };
  }

  // 4.5 Regime change at day 3000: alpha 0.08 -> 0.15, beta 0.90 -> 0.83 (R=40)
  const regimeParallel = await runParallelTask(
    {
      kind: "regime_bt",
      z2CritT1000: crit1000_99,
      z2CritT3000: crit3000_99,
      totalDays: 6000,
      breakDay: 3000,
      windowSize: 1000,
      stepDays: 20,
    },
    nRollingSeries,
    2,
    MAIN_SEED,
    71,
  );
  const regimeItems: RegimeDetectionMetrics[] = [];
  for (const out of regimeParallel.outputs) {
    if (out.kind === "regime_bt") regimeItems.push(...out.items);
  }

  const summarizeDelay = (vals: Array<number | null>) => {
    const detected = vals.filter((v): v is number => v !== null);
    return {
      detectedCount: detected.length,
      detectedFrac: detected.length / vals.length,
      meanDays: detected.length > 0 ? mean(detected) : null,
      medianDays: detected.length > 0 ? median(detected) : null,
    };
  };

  const nPostBins = regimeItems[0]?.fixedBinRates99.length ?? 30;
  const fixedBinTrajectory = Array.from({ length: nPostBins }, (_, b) =>
    mean(regimeItems.map((it) => it.fixedBinRates99[b] ?? 0.0)),
  );
  const rollingBinTrajectory = Array.from({ length: nPostBins }, (_, b) =>
    mean(regimeItems.map((it) => it.rollingBinRates99[b] ?? 0.0)),
  );
  const rollingAlphaTrajectory = Array.from({ length: nPostBins }, (_, b) =>
    mean(regimeItems.map((it) => it.rollingAlphaTrack[b] ?? 0.0)),
  );

  const regimeSummary = {
    nSeries: nRollingSeries,
    breakDay: 3000,
    postBreakTrueParams: { alpha: 0.15, beta: 0.83 },
    fixedRate99PostMean: mean(regimeItems.map((it) => it.fixedRate99Post)),
    fixedRate99PostSe: se(regimeItems.map((it) => it.fixedRate99Post)),
    rollingRate99PostMean: mean(regimeItems.map((it) => it.rollingRate99Post)),
    rollingRate99PostSe: se(regimeItems.map((it) => it.rollingRate99Post)),
    fixedPostZ2Mean: mean(regimeItems.map((it) => it.fixedPostZ2)),
    rollingPostZ2Mean: mean(regimeItems.map((it) => it.rollingPostZ2)),
    fixedDelay: {
      trail250Kupiec: summarizeDelay(regimeItems.map((it) => it.fixedDetectTrail250Kupiec)),
      trail250Cc: summarizeDelay(regimeItems.map((it) => it.fixedDetectTrail250Cc)),
      cumKupiec: summarizeDelay(regimeItems.map((it) => it.fixedDetectCumKupiec)),
      cumZ2: summarizeDelay(regimeItems.map((it) => it.fixedDetectCumZ2)),
      finalKupiecRejectedFrac:
        regimeItems.filter((it) => it.fixedFinalKupiecRejected).length / regimeItems.length,
      finalCcRejectedFrac:
        regimeItems.filter((it) => it.fixedFinalCcRejected).length / regimeItems.length,
      finalZ2RejectedFrac:
        regimeItems.filter((it) => it.fixedFinalZ2Rejected).length / regimeItems.length,
    },
    rollingDelay: {
      trail250Kupiec: summarizeDelay(regimeItems.map((it) => it.rollingDetectTrail250Kupiec)),
      trail250Cc: summarizeDelay(regimeItems.map((it) => it.rollingDetectTrail250Cc)),
      cumKupiec: summarizeDelay(regimeItems.map((it) => it.rollingDetectCumKupiec)),
      cumZ2: summarizeDelay(regimeItems.map((it) => it.rollingDetectCumZ2)),
      finalKupiecRejectedFrac:
        regimeItems.filter((it) => it.rollingFinalKupiecRejected).length / regimeItems.length,
      finalCcRejectedFrac:
        regimeItems.filter((it) => it.rollingFinalCcRejected).length / regimeItems.length,
      finalZ2RejectedFrac:
        regimeItems.filter((it) => it.rollingFinalZ2Rejected).length / regimeItems.length,
    },
    fixedBinRates99Per100d: fixedBinTrajectory,
    rollingBinRates99Per100d: rollingBinTrajectory,
    rollingAlphaTrackPer100d: rollingAlphaTrajectory,
  };

  // SVG 6: Regime change violation trajectory + alpha_hat adaptation
  const binDays = Array.from({ length: nPostBins }, (_, b) => (b + 0.5) * 100);
  writeFileSync(
    path.join(figDir, "day003d_regime-change-violations.svg"),
    renderLineChartSvg({
      title: "Perubahan Rezim Hari 3000 (alpha: 0.08 -> 0.15): Laju Pelanggaran 99% Pasca-Break",
      subtitle:
        "Rata-rata R = 40 seri per blok 100 hari: Estimasi Tetap vs Jendela Bergulir W = 1000",
      xLabel: "Hari Setelah Perubahan Rezim (t - 3000)",
      yLabel: "Tingkat Pelanggaran VaR 99% (Target Nominal = 0.01)",
      hLine: { y: 0.01, label: "Target Nominal 1%", color: "#111827" },
      series: [
        {
          label: "Estimasi Tetap (Dibekukan di t=3000)",
          color: "#dc2626",
          x: binDays,
          y: fixedBinTrajectory,
        },
        {
          label: "Jendela Bergulir (Re-est tiap 20d)",
          color: "#2563eb",
          x: binDays,
          y: rollingBinTrajectory,
        },
        {
          label: "Rata-rata alpha_hat / 10 (skala visual)",
          color: "#059669",
          dashed: true,
          x: binDays,
          y: rollingAlphaTrajectory.map((a) => a / 10.0),
        },
      ],
    }),
  );

  // SVG 7: Backtest size & power comparison
  const ku95 = kupiecResults["0.95"] as { rateExact: number; rateAsymp: number };
  const ku99 = kupiecResults["0.99"] as { rateExact: number; rateAsymp: number };
  const ch95 = chrisResults["0.95"] as {
    h0Size: { indRate: number; ccRate: number };
    h1PowerMarkov03: { indRate: number; ccRate: number };
  };
  const ch99 = chrisResults["0.99"] as {
    h0Size: { indRate: number; ccRate: number };
    h1PowerMarkov03: { indRate: number; ccRate: number };
  };
  writeFileSync(
    path.join(figDir, "day003d_backtest-size-power.svg"),
    renderGroupedBarSvg({
      title: "Ukuran (H0) dan Daya (H1) Uji Backtesting pada Taraf Signifikansi 5%",
      subtitle:
        "Kupiec (100k), Christoffersen H0 & Markov P(1|1)=0.3 (20k), Acerbi-Szekely Z2 (10k)",
      yLabel: "Fraksi Penolakan H0",
      hLine: { y: 0.05, label: "Taraf Nominal 5%", color: "#111827" },
      groups: [
        {
          category: "Kupiec H0 (Size)",
          values: [
            { label: "p = 0.95 / 0.975", value: ku95.rateExact, color: "#2563eb" },
            { label: "p = 0.99", value: ku99.rateExact, color: "#f97316" },
          ],
        },
        {
          category: "Chris Ind H0 (Size)",
          values: [
            { label: "p = 0.95 / 0.975", value: ch95.h0Size.indRate, color: "#2563eb" },
            { label: "p = 0.99", value: ch99.h0Size.indRate, color: "#f97316" },
          ],
        },
        {
          category: "Chris CC H1 (Cluster)",
          values: [
            { label: "p = 0.95 / 0.975", value: ch95.h1PowerMarkov03.ccRate, color: "#2563eb" },
            { label: "p = 0.99", value: ch99.h1PowerMarkov03.ccRate, color: "#f97316" },
          ],
        },
        {
          category: "Z2 H0 (Size)",
          values: [
            {
              label: "p = 0.95 / 0.975",
              value: z2Summary.freshH0Size.p975.rejectRateAtCrit05,
              color: "#2563eb",
            },
            {
              label: "p = 0.99",
              value: z2Summary.freshH0Size.p99.rejectRateAtCrit05,
              color: "#f97316",
            },
          ],
        },
        {
          category: "Z2 H1 (Gauss)",
          values: [
            {
              label: "p = 0.95 / 0.975",
              value: z2Summary.powerGaussianModel.p975.rejectRateAtCrit05,
              color: "#2563eb",
            },
            {
              label: "p = 0.99",
              value: z2Summary.powerGaussianModel.p99.rejectRateAtCrit05,
              color: "#f97316",
            },
          ],
        },
        {
          category: "Z2 H1 (0.8*ES)",
          values: [
            {
              label: "p = 0.95 / 0.975",
              value: z2Summary.powerScale80Es.p975.rejectRateAtCrit05,
              color: "#2563eb",
            },
            {
              label: "p = 0.99",
              value: z2Summary.powerScale80Es.p99.rejectRateAtCrit05,
              color: "#f97316",
            },
          ],
        },
      ],
    }),
  );

  const part4TimeSec = (performance.now() - tPart4Start) / 1000.0;
  const totalTimeSec = (performance.now() - tStartAll) / 1000.0;

  const results = {
    meta: {
      mainSeed: MAIN_SEED,
      workerCount,
      trueModel: TRUE_GARCH_PARAMS,
    },
    timingsSec: {
      part1: part1TimeSec,
      part2: part2TimeSec,
      part3: part3TimeSec,
      part4: part4TimeSec,
      total: totalTimeSec,
    },
    part1: {
      maxCdfNumErr,
      maxQuantileRoundtripErr,
      locationScaleRoundtripErr: lsRoundtripErr,
      cdfChecks,
      maxEsQuadratureErr: maxEsErr,
      esGrid,
      gpdChecks,
      pbh: {
        nSamples: nPbh,
        trueXi: trueXiPbh,
        byThresholdQ: pbhResults,
      },
    },
    part2: {
      theoreticalTrueModel: theoTrue,
      longSeriesMoments200k: {
        empVar: longVar,
        batchMeansSe50: longVarBatchSe,
        zExcessKurt: longZKurt,
        eps2Acf1: longAcf1,
        eps2Acf2: longAcf2,
      },
      repMoments200x2500: {
        empVar: summarizeReplications(repEmpVar, theoTrue.unconditionalVariance),
        zExcessKurt: summarizeReplications(repZKurt, theoTrue.innovationExcessKurtosisZ),
        eps2Acf1VsFormal: summarizeReplications(repAcf1, theoTrue.formalAcfLag1),
        eps2Acf2VsFormal: summarizeReplications(repAcf2, theoTrue.formalAcfLag2),
      },
      finiteFourthMomentAnchor: {
        params: anchorParams,
        theoD: anchorTheo.D,
        theoRho1: anchorTheo.rho1,
        theoRho2: anchorTheo.rho2,
        nRep: 60,
        nObs: 5000,
        empVar: summarizeReplications(
          anchorVars,
          anchorParams.omega / (1 - anchorParams.alpha - anchorParams.beta),
        ),
        zExcessKurt: summarizeReplications(anchorKurts, 6.0 / (12.0 - 4.0)),
        eps2Acf1: summarizeReplications(anchorAcf1, anchorTheo.rho1),
        eps2Acf2: summarizeReplications(anchorAcf2, anchorTheo.rho2),
      },
      deterministicReconstructionMaxDiff: maxSigma2ReconDiff,
      scoreVsCentralFd: {
        truePointMaxRelErr: scoreCheckTrue.maxRelErr,
        truePointRelErrs: scoreCheckTrue.relErrs,
        altPointMaxRelErr: scoreCheckAlt.maxRelErr,
        altPointRelErrs: scoreCheckAlt.relErrs,
      },
      bfgsVsNelderMeadNllDiff: bfgsVsNmLogLikDiff,
      replicationStudy: {
        nRep: nRepGarch,
        nObs: 2500,
        tMleConvergedCount: tConvCount,
        tMleActiveBoundaryHits: tBoundaryHits,
        gQmleConvergedCount: gConvCount,
        gQmleActiveBoundaryHits: gBoundaryHits,
        tMle: tMleSummary,
        gaussQmle: qmleSummary,
      },
    },
    part3: {
      nRep: nRepGarch,
      nObs: 2500,
      methods: riskSummaryClean,
      qualitativeChecks: {
        gaussianVar999RelBias: gaussRelBias999.relBias,
        gaussianVar999RelBiasSe: gaussRelBias999.relBiasSe,
        gaussianVar999ZScore: gaussRelBias999.relBias / gaussRelBias999.relBiasSe,
        mseFhsMinusEvt999Mean: mean(diffMseFhsMinusEvt),
        mseFhsMinusEvt999Se: se(diffMseFhsMinusEvt),
        mseFhsMinusT999Mean: mean(diffMseFhsMinusT),
        mseFhsMinusT999Se: se(diffMseFhsMinusT),
        evtLowerRmseThanFhs999:
          (riskSummaryClean.fhs_evt?.var999 as { relRmse: number }).relRmse <
          (riskSummaryClean.fhs?.var999 as { relRmse: number }).relRmse,
        tLowerRmseThanFhs999:
          (riskSummaryClean.t_parametric?.var999 as { relRmse: number }).relRmse <
          (riskSummaryClean.fhs?.var999 as { relRmse: number }).relRmse,
      },
      thresholdSensitivity: thresholdSens,
      xiHatVsInvNuHat: xiHatVsInvNu,
    },
    part4: {
      kupiecH0Size: kupiecResults,
      christoffersenSizeAndPower: chrisResults,
      acerbiSzekelyZ2: z2Summary,
      rollingBacktest: rollingSummary,
      regimeChange: regimeSummary,
    },
  };

  writeFileSync(path.join(outDir, "results.json"), serializeJson(results));

  console.log(`Workers: ${workerCount} | Total runtime: ${totalTimeSec.toFixed(2)}s`);
  console.log(
    `Part 1 (${part1TimeSec.toFixed(2)}s): CDF err=${maxCdfNumErr.toExponential(2)}, ES err=${maxEsErr.toExponential(2)}, xi_POT(0.999)=${pbhResults[3]?.xiPot.toFixed(4)}`,
  );
  console.log(
    `Part 2 (${part2TimeSec.toFixed(2)}s): Recon=${maxSigma2ReconDiff.toExponential(2)}, ScoreErr=${scoreCheckTrue.maxRelErr.toExponential(2)}, Boundaries=${tBoundaryHits}`,
  );
  console.log(
    `Part 3 (${part3TimeSec.toFixed(2)}s): Gauss VaR999 bias=${gaussRelBias999.relBias.toFixed(4)}, EVT RMSE=${(riskSummaryClean.fhs_evt?.var999 as { relRmse: number }).relRmse.toFixed(4)}, FHS RMSE=${(riskSummaryClean.fhs?.var999 as { relRmse: number }).relRmse.toFixed(4)}`,
  );
  console.log(
    `Part 4 (${part4TimeSec.toFixed(2)}s): TrueCtrl rate99=${(rollingSummary.true_control as { rate99Mean: number }).rate99Mean.toFixed(5)}`,
  );
}

await main();
