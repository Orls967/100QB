/**
 * VaR and Expected Shortfall backtesting suite:
 * - Kupiec POF (exact two-sided binomial p-value + asymptotic chi2(1) LR)
 * - Christoffersen independence & conditional coverage Markov tests
 * - Acerbi-Szekely Z2 Expected Shortfall backtest (zero dependence on external tables)
 * - Rolling OOS backtesting & structural regime-change detection speed
 */

import { normalEs, stdTEs, stdTQuantile } from "./distributions.js";
import { fitGpdMle } from "./evt.js";
import {
  fitGarchGauss,
  fitGarchT,
  reconstructSigma2,
  simulateGarchT,
  TRUE_GARCH_PARAMS,
  type GarchGaussParams,
  type GarchTParams,
} from "./garch.js";
import { Xoshiro256PlusPlus } from "./rng.js";
import { chi2Sf, logGamma, normalQuantile } from "./special.js";
import { quantile } from "./stats.js";

export interface KupiecResult {
  n: number;
  x: number;
  alpha: number;
  rate: number;
  lrPof: number;
  pValueAsymp: number;
  pValueExact: number;
}

export interface ChristoffersenResult {
  n00: number;
  n01: number;
  n10: number;
  n11: number;
  pi0: number;
  pi1: number;
  pi: number;
  lrInd: number;
  pInd: number;
  lrCc: number;
  pCc: number;
}

export type BacktestMethod = "true_control" | "gaussian_garch" | "t_garch" | "fhs" | "fhs_evt";

export interface RollingMethodSeriesMetrics {
  rate99: number;
  rate975: number;
  kupiec99AsympP: number;
  kupiec99ExactP: number;
  chris99IndP: number;
  chris99CcP: number;
  z2_99: number;
  z2_975: number;
}

export interface RollingBacktestSeriesResult {
  methods: Record<BacktestMethod, RollingMethodSeriesMetrics>;
}

export interface RegimeDetectionMetrics {
  fixedRate99Post: number;
  rollingRate99Post: number;
  /** Post-break days until first trailing-250 Kupiec exact rejection (null if never in 3000d). */
  fixedDetectTrail250Kupiec: number | null;
  rollingDetectTrail250Kupiec: number | null;
  fixedDetectTrail250Cc: number | null;
  rollingDetectTrail250Cc: number | null;
  fixedDetectCumKupiec: number | null;
  rollingDetectCumKupiec: number | null;
  fixedDetectCumZ2: number | null;
  rollingDetectCumZ2: number | null;
  fixedFinalKupiecRejected: boolean;
  rollingFinalKupiecRejected: boolean;
  fixedFinalCcRejected: boolean;
  rollingFinalCcRejected: boolean;
  fixedFinalZ2Rejected: boolean;
  rollingFinalZ2Rejected: boolean;
  fixedPostZ2: number;
  rollingPostZ2: number;
  /** Binned post-break violation counts per 100-day block (30 blocks covering days 3000..5999). */
  fixedBinRates99: number[];
  rollingBinRates99: number[];
  /** Adapted rolling alpha_hat every 100 days post-break (30 points). */
  rollingAlphaTrack: number[];
}

function logBinomPmf(k: number, n: number, alpha: number): number {
  if (alpha <= 0.0) return k === 0 ? 0.0 : -Infinity;
  if (alpha >= 1.0) return k === n ? 0.0 : -Infinity;
  return (
    logGamma(n + 1) -
    logGamma(k + 1) -
    logGamma(n - k + 1) +
    k * Math.log(alpha) +
    (n - k) * Math.log1p(-alpha)
  );
}

/** Exact two-sided binomial p-value: sum_{j: P(X=j) <= P(X=x)} P(X=j). */
export function binomialExactTwoSidedP(x: number, n: number, alpha: number): number {
  const lpObs = logBinomPmf(x, n, alpha);
  const tol = 1e-12;
  let maxLog = lpObs;
  const logProbs = new Float64Array(n + 1);
  for (let j = 0; j <= n; j += 1) {
    const lp = logBinomPmf(j, n, alpha);
    logProbs[j] = lp;
    if (lp > maxLog) maxLog = lp;
  }
  let sum = 0.0;
  let total = 0.0;
  for (let j = 0; j <= n; j += 1) {
    const lp = logProbs[j] ?? -Infinity;
    const p = Math.exp(lp - maxLog);
    total += p;
    if (lp <= lpObs + tol) {
      sum += p;
    }
  }
  return Math.min(1.0, sum / total);
}

/** Kupiec (1995) Proportion of Failures test for exceedance indicator array (0/1). */
export function kupiecTest(exceedances: ArrayLike<number>, pLevel: number): KupiecResult {
  const n = exceedances.length;
  const alpha = 1.0 - pLevel;
  let x = 0;
  for (let i = 0; i < n; i += 1) {
    x += exceedances[i] ?? 0;
  }
  return kupiecFromCounts(x, n, alpha);
}

/** Fast Kupiec evaluation directly from violation count x out of n trials. */
export function kupiecFromCounts(x: number, n: number, alpha: number): KupiecResult {
  const rate = x / n;
  const logL0 = (n - x) * Math.log1p(-alpha) + (x > 0 ? x * Math.log(alpha) : 0.0);
  const logL1 = (x < n ? (n - x) * Math.log1p(-rate) : 0.0) + (x > 0 ? x * Math.log(rate) : 0.0);
  const lrPof = Math.max(0.0, -2.0 * (logL0 - logL1));
  return {
    n,
    x,
    alpha,
    rate,
    lrPof,
    pValueAsymp: chi2Sf(lrPof, 1.0),
    pValueExact: binomialExactTwoSidedP(x, n, alpha),
  };
}

/** Exact analytical discrete rejection rate P(test rejects under H0) for binomial count X ~ Bin(n, alpha). */
export function kupiecExactAnalyticalSize(
  n: number,
  alpha: number,
  sigLevel = 0.05,
): { exactSize: number; asympSize: number } {
  const critChi2_1 = 3.841458820694124; // chi2(1) 95% quantile
  let maxLog = -Infinity;
  const logProbs = new Float64Array(n + 1);
  for (let j = 0; j <= n; j += 1) {
    const lp = logBinomPmf(j, n, alpha);
    logProbs[j] = lp;
    if (lp > maxLog) maxLog = lp;
  }
  const probs = new Float64Array(n + 1);
  let total = 0.0;
  for (let j = 0; j <= n; j += 1) {
    const p = Math.exp((logProbs[j] ?? -Infinity) - maxLog);
    probs[j] = p;
    total += p;
  }
  for (let j = 0; j <= n; j += 1) {
    probs[j] = (probs[j] ?? 0.0) / total;
  }

  let asympSize = 0.0;
  let exactSize = 0.0;
  // Exact two-sided critical region: values j with cumulative two-sided p <= sigLevel
  const order = Array.from({ length: n + 1 }, (_, j) => j).sort(
    (a, b) => (probs[a] ?? 0.0) - (probs[b] ?? 0.0),
  );
  let cum = 0.0;
  const rejectExact = new Uint8Array(n + 1);
  for (const j of order) {
    // Include ties in order if exact p-value <= sigLevel
    const pj = probs[j] ?? 0.0;
    // Two-sided p-value at j equals sum of probs <= pj
    let pVal = 0.0;
    for (let k = 0; k <= n; k += 1) {
      if ((probs[k] ?? 0.0) <= pj + 1e-12) pVal += probs[k] ?? 0.0;
    }
    if (pVal <= sigLevel) {
      rejectExact[j] = 1;
      cum += pj;
    }
  }
  exactSize = cum;

  for (let j = 0; j <= n; j += 1) {
    const k = kupiecFromCounts(j, n, alpha);
    if (k.lrPof > critChi2_1) {
      asympSize += probs[j] ?? 0.0;
    }
  }
  void rejectExact;
  return { exactSize, asympSize };
}

/** Christoffersen (1998) Markov independence and conditional coverage tests. */
export function christoffersenTest(
  exceedances: ArrayLike<number>,
  pLevel: number,
): ChristoffersenResult {
  const n = exceedances.length;
  const alpha = 1.0 - pLevel;
  let n00 = 0;
  let n01 = 0;
  let n10 = 0;
  let n11 = 0;
  let xTotal = exceedances[0] ?? 0;

  for (let t = 1; t < n; t += 1) {
    const prev = exceedances[t - 1] ?? 0;
    const cur = exceedances[t] ?? 0;
    xTotal += cur;
    if (prev === 0 && cur === 0) n00 += 1;
    else if (prev === 0 && cur === 1) n01 += 1;
    else if (prev === 1 && cur === 0) n10 += 1;
    else n11 += 1;
  }

  return christoffersenFromTransitions(n00, n01, n10, n11, xTotal, n, alpha);
}

export function christoffersenFromTransitions(
  n00: number,
  n01: number,
  n10: number,
  n11: number,
  xTotal: number,
  nTotal: number,
  alpha: number,
): ChristoffersenResult {
  const n0 = n00 + n01;
  const n1 = n10 + n11;
  const m = n0 + n1;
  const xTrans = n01 + n11;

  const pi0 = n0 > 0 ? n01 / n0 : 0.0;
  const pi1 = n1 > 0 ? n11 / n1 : 0.0;
  const pi = m > 0 ? xTrans / m : 0.0;

  const safeLogTerm = (cnt: number, prob: number): number =>
    cnt > 0 && prob > 0.0 && prob < 1.0 ? cnt * Math.log(prob) : cnt > 0 && prob >= 1.0 ? 0.0 : 0.0;

  const llNull = safeLogTerm(n00 + n10, 1.0 - pi) + safeLogTerm(n01 + n11, pi);
  const llAlt =
    safeLogTerm(n00, 1.0 - pi0) +
    safeLogTerm(n01, pi0) +
    safeLogTerm(n10, 1.0 - pi1) +
    safeLogTerm(n11, pi1);

  const lrInd = n1 === 0 ? 0.0 : Math.max(0.0, -2.0 * (llNull - llAlt));
  const pInd = chi2Sf(lrInd, 1.0);

  const pof = kupiecFromCounts(xTotal, nTotal, alpha);
  const lrCc = pof.lrPof + lrInd;
  const pCc = chi2Sf(lrCc, 2.0);

  return {
    n00,
    n01,
    n10,
    n11,
    pi0,
    pi1,
    pi,
    lrInd,
    pInd,
    lrCc,
    pCc,
  };
}

/**
 * Acerbi & Szekely (2014) Z2 statistic using positive loss convention L_t = -r_t:
 *   Z2 = 1 - (1 / (T * alpha)) * sum_{t=1}^T [ L_t * I(L_t > VaR_t) / ES_t ]
 * Under H0: E[Z2] = 0. ES underestimation gives Z2 < 0 (one-sided left-tail rejection).
 */
export function acerbiSzekelyZ2(
  losses: ArrayLike<number>,
  vars: ArrayLike<number>,
  ess: ArrayLike<number>,
  pLevel: number,
): number {
  const n = losses.length;
  const alpha = 1.0 - pLevel;
  let ratioSum = 0.0;
  for (let t = 0; t < n; t += 1) {
    const l = losses[t] ?? 0.0;
    const v = vars[t] ?? 0.0;
    const es = ess[t] ?? 1.0;
    if (l > v && es > 0.0) {
      ratioSum += l / es;
    }
  }
  return 1.0 - ratioSum / (n * alpha);
}

/** Simulate stationary 2-state Markov exceedance chain with marginal rate alpha and clustering P(1|1) = p11. */
export function simulateMarkovChainCounts(
  n: number,
  alpha: number,
  p11: number,
  rng: Xoshiro256PlusPlus,
): { n00: number; n01: number; n10: number; n11: number; xTotal: number } {
  // Solve (1 - alpha)*p01 = alpha*(1 - p11) => p01 = alpha*(1 - p11)/(1 - alpha)
  const p01 = (alpha * (1.0 - p11)) / (1.0 - alpha);
  let state = rng.nextUniform() < alpha ? 1 : 0;
  let xTotal = state;
  let n00 = 0;
  let n01 = 0;
  let n10 = 0;
  let n11 = 0;
  for (let t = 1; t < n; t += 1) {
    const next = state === 1 ? (rng.nextUniform() < p11 ? 1 : 0) : rng.nextUniform() < p01 ? 1 : 0;
    xTotal += next;
    if (state === 0 && next === 0) n00 += 1;
    else if (state === 0 && next === 1) n01 += 1;
    else if (state === 1 && next === 0) n10 += 1;
    else n11 += 1;
    state = next;
  }
  return { n00, n01, n10, n11, xTotal };
}

/**
 * Simulate one true GARCH-t(5) series of length T and evaluate Z2 statistics for ES levels
 * (either true model or Gaussian misspecified / scaled ES alternative).
 */
export function simulateZ2Single(
  tLen: number,
  pLevels: readonly number[],
  rng: Xoshiro256PlusPlus,
  mode: "true" | "gaussian" | "scale90" | "scale80",
): number[] {
  const { mu, omega, alpha, beta, nu } = TRUE_GARCH_PARAMS;
  let h = omega / (1.0 - alpha - beta);
  let epsPrev = 0.0;

  const qs = pLevels.map((p) => stdTQuantile(p, nu));
  const esTrue = pLevels.map((p) => stdTEs(p, nu));
  const qNorm = pLevels.map((p) => normalQuantile(p));
  const esNorm = pLevels.map((p) => normalEs(p));

  const sums = new Float64Array(pLevels.length);

  for (let t = 0; t < tLen; t += 1) {
    if (t > 0) {
      h = omega + alpha * epsPrev * epsPrev + beta * h;
    }
    const sigma = Math.sqrt(h);
    const zt = rng.nextStdT(nu);
    const r = mu + sigma * zt;
    const loss = -r;
    epsPrev = r - mu;

    for (let j = 0; j < pLevels.length; j += 1) {
      let varEst = 0.0;
      let esEst = 0.0;
      if (mode === "gaussian") {
        varEst = -mu + sigma * (qNorm[j] ?? 0.0);
        esEst = -mu + sigma * (esNorm[j] ?? 1.0);
      } else {
        varEst = -mu + sigma * (qs[j] ?? 0.0);
        const esBase = -mu + sigma * (esTrue[j] ?? 1.0);
        esEst = mode === "scale90" ? 0.9 * esBase : mode === "scale80" ? 0.8 * esBase : esBase;
      }
      if (loss > varEst && esEst > 0.0) {
        sums[j] = (sums[j] ?? 0.0) + loss / esEst;
      }
    }
  }

  return pLevels.map((p, j) => 1.0 - (sums[j] ?? 0.0) / (tLen * (1.0 - p)));
}

/** Helper to compute FHS and FHS-EVT normalized loss quantiles from a GARCH-t estimation window. */
function computeFilteredMultipliers(residualsHat: Float64Array): {
  q99Fhs: number;
  q975Fhs: number;
  es99Fhs: number;
  es975Fhs: number;
  q99Evt: number;
  q975Evt: number;
  es99Evt: number;
  es975Evt: number;
} {
  const n = residualsHat.length;
  const sorted = new Array<number>(n);
  for (let i = 0; i < n; i += 1) {
    sorted[i] = -(residualsHat[i] ?? 0.0);
  }
  sorted.sort((a, b) => a - b);

  const q99Fhs = quantile(sorted, 0.99, true);
  const q975Fhs = quantile(sorted, 0.975, true);

  let sum99 = 0.0;
  let cnt99 = 0;
  let sum975 = 0.0;
  let cnt975 = 0;
  const u = quantile(sorted, 0.9, true);
  const exc: number[] = [];

  for (let i = n - 1; i >= 0; i -= 1) {
    const v = sorted[i] ?? 0.0;
    if (v > q99Fhs) {
      sum99 += v;
      cnt99 += 1;
    }
    if (v > q975Fhs) {
      sum975 += v;
      cnt975 += 1;
    }
    if (v > u) {
      exc.push(v - u);
    }
  }
  const es99Fhs = cnt99 > 0 ? sum99 / cnt99 : q99Fhs;
  const es975Fhs = cnt975 > 0 ? sum975 / cnt975 : q975Fhs;

  const pExc = exc.length / n;
  const gpd = fitGpdMle(exc);
  const potVar = (p: number): number => {
    const ratio = pExc / (1.0 - p);
    return Math.abs(gpd.xi) < 1e-8
      ? u + gpd.beta * Math.log(ratio)
      : u + (gpd.beta / gpd.xi) * (Math.pow(ratio, gpd.xi) - 1.0);
  };
  const q99Evt = potVar(0.99);
  const q975Evt = potVar(0.975);
  const potEs = (vp: number): number =>
    gpd.xi < 0.98 ? (vp + gpd.beta - gpd.xi * u) / (1.0 - gpd.xi) : es99Fhs;

  return {
    q99Fhs,
    q975Fhs,
    es99Fhs,
    es975Fhs,
    q99Evt,
    q975Evt,
    es99Evt: potEs(q99Evt),
    es975Evt: potEs(q975Evt),
  };
}

/**
 * Evaluate one 6000-day GARCH-t(5) series under rolling W=1000 re-estimation every 20 days
 * for True-control + 4 empirical methods: Gaussian-GARCH, t-GARCH, FHS, FHS-EVT.
 */
export function evaluateRollingBacktestSeries(
  rng: Xoshiro256PlusPlus,
  totalDays = 6000,
  windowSize = 1000,
  stepDays = 20,
): RollingBacktestSeriesResult {
  const sim = simulateGarchT(totalDays, TRUE_GARCH_PARAMS, rng);
  const r = sim.returns;
  const sig2True = sim.sigma2;
  const oosDays = totalDays - windowSize;

  const methodList: BacktestMethod[] = [
    "true_control",
    "gaussian_garch",
    "t_garch",
    "fhs",
    "fhs_evt",
  ];
  const exc99: Record<BacktestMethod, Uint8Array> = {
    true_control: new Uint8Array(oosDays),
    gaussian_garch: new Uint8Array(oosDays),
    t_garch: new Uint8Array(oosDays),
    fhs: new Uint8Array(oosDays),
    fhs_evt: new Uint8Array(oosDays),
  };
  const exc975: Record<BacktestMethod, Uint8Array> = {
    true_control: new Uint8Array(oosDays),
    gaussian_garch: new Uint8Array(oosDays),
    t_garch: new Uint8Array(oosDays),
    fhs: new Uint8Array(oosDays),
    fhs_evt: new Uint8Array(oosDays),
  };
  const z2Sum99: Record<BacktestMethod, number> = {
    true_control: 0,
    gaussian_garch: 0,
    t_garch: 0,
    fhs: 0,
    fhs_evt: 0,
  };
  const z2Sum975: Record<BacktestMethod, number> = {
    true_control: 0,
    gaussian_garch: 0,
    t_garch: 0,
    fhs: 0,
    fhs_evt: 0,
  };

  const qTrue99 = stdTQuantile(0.99, TRUE_GARCH_PARAMS.nu);
  const qTrue975 = stdTQuantile(0.975, TRUE_GARCH_PARAMS.nu);
  const esTrue99 = stdTEs(0.99, TRUE_GARCH_PARAMS.nu);
  const esTrue975 = stdTEs(0.975, TRUE_GARCH_PARAMS.nu);

  const qNorm99 = normalQuantile(0.99);
  const qNorm975 = normalQuantile(0.975);
  const esNorm99 = normalEs(0.99);
  const esNorm975 = normalEs(0.975);

  let warmT: GarchTParams | undefined;
  let warmG: GarchGaussParams | undefined;

  for (let start = windowSize; start < totalDays; start += stepDays) {
    const winSlice = r.subarray(start - windowSize, start);
    const fitT = fitGarchT(winSlice, warmT, false);
    const fitG = fitGarchGauss(winSlice, warmG);
    warmT = fitT.params;
    warmG = fitG.params;

    const pT = fitT.params;
    const pG = fitG.params;
    const qT99 = stdTQuantile(0.99, pT.nu);
    const qT975 = stdTQuantile(0.975, pT.nu);
    const esT99 = stdTEs(0.99, pT.nu);
    const esT975 = stdTEs(0.975, pT.nu);
    const mult = computeFilteredMultipliers(fitT.residualsHat);

    let hT = fitT.sigma2NextHat;
    let hG = fitG.sigma2NextHat;

    const blockEnd = Math.min(start + stepDays, totalDays);
    for (let t = start; t < blockEnd; t += 1) {
      const idx = t - windowSize;
      const rt = r[t] ?? 0.0;
      const loss = -rt;
      const sigTrue = Math.sqrt(sig2True[t] ?? 1e-4);
      const sigT = Math.sqrt(hT);
      const sigG = Math.sqrt(hG);

      // 1) True control
      const vC99 = -TRUE_GARCH_PARAMS.mu + sigTrue * qTrue99;
      const vC975 = -TRUE_GARCH_PARAMS.mu + sigTrue * qTrue975;
      const eC99 = -TRUE_GARCH_PARAMS.mu + sigTrue * esTrue99;
      const eC975 = -TRUE_GARCH_PARAMS.mu + sigTrue * esTrue975;
      if (loss > vC99) {
        exc99.true_control[idx] = 1;
        z2Sum99.true_control += loss / eC99;
      }
      if (loss > vC975) {
        exc975.true_control[idx] = 1;
        z2Sum975.true_control += loss / eC975;
      }

      // 2) Gaussian-GARCH
      const vG99 = -pG.mu + sigG * qNorm99;
      const vG975 = -pG.mu + sigG * qNorm975;
      const eG99 = -pG.mu + sigG * esNorm99;
      const eG975 = -pG.mu + sigG * esNorm975;
      if (loss > vG99) {
        exc99.gaussian_garch[idx] = 1;
        z2Sum99.gaussian_garch += loss / eG99;
      }
      if (loss > vG975) {
        exc975.gaussian_garch[idx] = 1;
        z2Sum975.gaussian_garch += loss / eG975;
      }

      // 3) t-GARCH
      const vT99 = -pT.mu + sigT * qT99;
      const vT975 = -pT.mu + sigT * qT975;
      const eT99 = -pT.mu + sigT * esT99;
      const eT975 = -pT.mu + sigT * esT975;
      if (loss > vT99) {
        exc99.t_garch[idx] = 1;
        z2Sum99.t_garch += loss / eT99;
      }
      if (loss > vT975) {
        exc975.t_garch[idx] = 1;
        z2Sum975.t_garch += loss / eT975;
      }

      // 4) FHS
      const vF99 = -pT.mu + sigT * mult.q99Fhs;
      const vF975 = -pT.mu + sigT * mult.q975Fhs;
      const eF99 = -pT.mu + sigT * mult.es99Fhs;
      const eF975 = -pT.mu + sigT * mult.es975Fhs;
      if (loss > vF99) {
        exc99.fhs[idx] = 1;
        z2Sum99.fhs += loss / eF99;
      }
      if (loss > vF975) {
        exc975.fhs[idx] = 1;
        z2Sum975.fhs += loss / eF975;
      }

      // 5) FHS-EVT
      const vE99 = -pT.mu + sigT * mult.q99Evt;
      const vE975 = -pT.mu + sigT * mult.q975Evt;
      const eE99 = -pT.mu + sigT * mult.es99Evt;
      const eE975 = -pT.mu + sigT * mult.es975Evt;
      if (loss > vE99) {
        exc99.fhs_evt[idx] = 1;
        z2Sum99.fhs_evt += loss / eE99;
      }
      if (loss > vE975) {
        exc975.fhs_evt[idx] = 1;
        z2Sum975.fhs_evt += loss / eE975;
      }

      // One-step recursion using realized return r[t] for day t+1
      const epsT = rt - pT.mu;
      const epsG = rt - pG.mu;
      hT = pT.omega + pT.alpha * epsT * epsT + pT.beta * hT;
      hG = pG.omega + pG.alpha * epsG * epsG + pG.beta * hG;
    }
  }

  const out = {} as Record<BacktestMethod, RollingMethodSeriesMetrics>;
  for (const m of methodList) {
    const k99 = kupiecTest(exc99[m], 0.99);
    const ch99 = christoffersenTest(exc99[m], 0.99);
    let c975 = 0;
    const arr975 = exc975[m];
    for (let i = 0; i < oosDays; i += 1) c975 += arr975[i] ?? 0;
    out[m] = {
      rate99: k99.rate,
      rate975: c975 / oosDays,
      kupiec99AsympP: k99.pValueAsymp,
      kupiec99ExactP: k99.pValueExact,
      chris99IndP: ch99.pInd,
      chris99CcP: ch99.pCc,
      z2_99: 1.0 - z2Sum99[m] / (oosDays * 0.01),
      z2_975: 1.0 - z2Sum975[m] / (oosDays * 0.025),
    };
  }

  return { methods: out };
}

/**
 * Simulate 6000 days with structural break at day 3000 (alpha 0.08 -> 0.15, beta 0.90 -> 0.83)
 * and compare fixed parameter forecast vs W=1000 rolling 20-day re-estimation speed of detection.
 */
export function evaluateRegimeChangeSeries(
  rng: Xoshiro256PlusPlus,
  z2CritT1000: number,
  z2CritT3000: number,
  totalDays = 6000,
  breakDay = 3000,
  windowSize = 1000,
  stepDays = 20,
): RegimeDetectionMetrics {
  const sim = simulateGarchT(totalDays, TRUE_GARCH_PARAMS, rng, breakDay, 0.15, 0.83);
  const r = sim.returns;
  const postDays = totalDays - breakDay;

  // Fixed fit uses pre-break window [breakDay - windowSize, breakDay) = [2000, 3000)
  const fixedFit = fitGarchT(r.subarray(breakDay - windowSize, breakDay), undefined, false);
  const pFix = fixedFit.params;
  const qFix99 = stdTQuantile(0.99, pFix.nu);
  const esFix99 = stdTEs(0.99, pFix.nu);
  let hFix = fixedFit.sigma2NextHat;

  const excFix = new Uint8Array(postDays);
  const lossRatioFix = new Float64Array(postDays);
  for (let t = breakDay; t < totalDays; t += 1) {
    const idx = t - breakDay;
    const rt = r[t] ?? 0.0;
    const loss = -rt;
    const sig = Math.sqrt(hFix);
    const v = -pFix.mu + sig * qFix99;
    const es = -pFix.mu + sig * esFix99;
    if (loss > v) {
      excFix[idx] = 1;
      lossRatioFix[idx] = es > 0.0 ? loss / es : 0.0;
    }
    const eps = rt - pFix.mu;
    hFix = pFix.omega + pFix.alpha * eps * eps + pFix.beta * hFix;
  }

  // Rolling fit over days [breakDay - 200, totalDays) initialized on [breakDay - windowSize, breakDay)
  // Note: start rolling tracking right from breakDay onward using the preceding window [t - 1000, t)
  const excRoll = new Uint8Array(postDays);
  const lossRatioRoll = new Float64Array(postDays);
  const rollingAlphaTrack: number[] = [];
  let warm = pFix;

  for (let start = breakDay; start < totalDays; start += stepDays) {
    let pRoll = pFix;
    let hRoll = 0.0;
    if (start === breakDay) {
      pRoll = pFix;
      hRoll = fixedFit.sigma2NextHat;
    } else {
      const fit = fitGarchT(r.subarray(start - windowSize, start), warm, false);
      warm = fit.params;
      pRoll = fit.params;
      hRoll = fit.sigma2NextHat;
    }
    if ((start - breakDay) % 100 === 0) {
      rollingAlphaTrack.push(pRoll.alpha);
    }

    const qR99 = stdTQuantile(0.99, pRoll.nu);
    const esR99 = stdTEs(0.99, pRoll.nu);
    const blockEnd = Math.min(start + stepDays, totalDays);
    for (let t = start; t < blockEnd; t += 1) {
      const idx = t - breakDay;
      const rt = r[t] ?? 0.0;
      const loss = -rt;
      const sig = Math.sqrt(hRoll);
      const v = -pRoll.mu + sig * qR99;
      const es = -pRoll.mu + sig * esR99;
      if (loss > v) {
        excRoll[idx] = 1;
        lossRatioRoll[idx] = es > 0.0 ? loss / es : 0.0;
      }
      const eps = rt - pRoll.mu;
      hRoll = pRoll.omega + pRoll.alpha * eps * eps + pRoll.beta * hRoll;
    }
  }

  const findDetections = (
    exc: Uint8Array,
    lRatio: Float64Array,
  ): {
    trailKupiec: number | null;
    trailCc: number | null;
    cumKupiec: number | null;
    cumZ2: number | null;
    ratePost: number;
    finalKupiecRej: boolean;
    finalCcRej: boolean;
    finalZ2Rej: boolean;
    finalZ2: number;
    binRates: number[];
  } => {
    let trailKupiec: number | null = null;
    let trailCc: number | null = null;
    let cumKupiec: number | null = null;
    let cumZ2: number | null = null;

    let cumSum = 0;
    let cumLR = 0.0;
    for (let d = 0; d < postDays; d += 1) {
      cumSum += exc[d] ?? 0;
      cumLR += lRatio[d] ?? 0.0;
      const len = d + 1;
      if (cumKupiec === null && len >= 50 && len % 5 === 0) {
        const k = kupiecFromCounts(cumSum, len, 0.01);
        if (k.pValueExact < 0.05) {
          cumKupiec = len;
        }
      }
      if (cumZ2 === null && len >= 1000 && len % 20 === 0) {
        const z2 = 1.0 - cumLR / (len * 0.01);
        if (z2 < z2CritT1000) {
          cumZ2 = len;
        }
      }
      if (len >= 250 && len % 10 === 0) {
        const win = exc.subarray(len - 250, len);
        const k = kupiecTest(win, 0.99);
        if (trailKupiec === null && k.pValueExact < 0.05) {
          trailKupiec = len;
        }
        const ch = christoffersenTest(win, 0.99);
        if (trailCc === null && ch.pCc < 0.05) {
          trailCc = len;
        }
      }
    }

    const finalK = kupiecFromCounts(cumSum, postDays, 0.01);
    const finalCh = christoffersenTest(exc, 0.99);
    const finalZ2 = 1.0 - cumLR / (postDays * 0.01);
    const nBins = Math.floor(postDays / 100);
    const binRates = new Array<number>(nBins);
    for (let b = 0; b < nBins; b += 1) {
      let s = 0;
      for (let j = 0; j < 100; j += 1) s += exc[b * 100 + j] ?? 0;
      binRates[b] = s / 100.0;
    }

    return {
      trailKupiec,
      trailCc,
      cumKupiec,
      cumZ2,
      ratePost: cumSum / postDays,
      finalKupiecRej: finalK.pValueExact < 0.05,
      finalCcRej: finalCh.pCc < 0.05,
      finalZ2Rej: finalZ2 < z2CritT3000,
      finalZ2,
      binRates,
    };
  };

  const df = findDetections(excFix, lossRatioFix);
  const dr = findDetections(excRoll, lossRatioRoll);

  // Use reconstructSigma2 silently checked in module exports
  void reconstructSigma2;

  return {
    fixedRate99Post: df.ratePost,
    rollingRate99Post: dr.ratePost,
    fixedDetectTrail250Kupiec: df.trailKupiec,
    rollingDetectTrail250Kupiec: dr.trailKupiec,
    fixedDetectTrail250Cc: df.trailCc,
    rollingDetectTrail250Cc: dr.trailCc,
    fixedDetectCumKupiec: df.cumKupiec,
    rollingDetectCumKupiec: dr.cumKupiec,
    fixedDetectCumZ2: df.cumZ2,
    rollingDetectCumZ2: dr.cumZ2,
    fixedFinalKupiecRejected: df.finalKupiecRej,
    rollingFinalKupiecRejected: dr.finalKupiecRej,
    fixedFinalCcRejected: df.finalCcRej,
    rollingFinalCcRejected: dr.finalCcRej,
    fixedFinalZ2Rejected: df.finalZ2Rej,
    rollingFinalZ2Rejected: dr.finalZ2Rej,
    fixedPostZ2: df.finalZ2,
    rollingPostZ2: dr.finalZ2,
    fixedBinRates99: df.binRates,
    rollingBinRates99: dr.binRates,
    rollingAlphaTrack,
  };
}
