/**
 * Extreme Value Theory (EVT) estimators:
 * - Generalized Pareto Distribution (GPD) MLE with analytical gradient + feasibility guard
 * - Hill tail-index estimator on upper k order statistics
 * - Peaks-Over-Threshold (POT) pipeline returning threshold u, N_u, GPD params, VaR, and ES
 */

import { gpdPotEs, gpdPotVar } from "./distributions.js";
import { bfgs, type ValueGrad } from "./optim.js";
import { quantile } from "./stats.js";

export interface GpdFitResult {
  xi: number;
  beta: number;
  logLik: number;
  nExceed: number;
  converged: boolean;
  seXi: number;
  seBeta: number;
}

export interface PotSummary {
  qThreshold: number;
  u: number;
  nTotal: number;
  nExceed: number;
  pExceed: number;
  gpd: GpdFitResult;
  vars: Record<string, number>;
  ess: Record<string, number>;
}

/**
 * Evaluate GPD negative log-likelihood and gradient w.r.t. (xi, log(beta)) on exceedances y >= 0.
 * Guarded smoothly around xi = 0 using exact Taylor expansion to 3rd order.
 */
function gpdNegLogLikAndGrad(
  exceedances: ArrayLike<number>,
  xi: number,
  logBeta: number,
): ValueGrad {
  const n = exceedances.length;
  const beta = Math.exp(logBeta);
  if (!(beta > 1e-15 && Number.isFinite(beta))) {
    return { value: 1e15, grad: [0, 0] };
  }

  // Exact exponential limit neighborhood |xi| < 1e-6
  if (Math.abs(xi) < 1e-6) {
    let sumY = 0.0;
    let sumY2 = 0.0;
    let sumY3 = 0.0;
    const invBeta = 1.0 / beta;
    for (let i = 0; i < n; i += 1) {
      const y = (exceedances[i] ?? 0.0) * invBeta;
      const y2 = y * y;
      sumY += y;
      sumY2 += y2;
      sumY3 += y2 * y;
    }
    // log L(xi, beta) = -n*log(beta) - sum(y) + xi * (-sum(y) + 0.5*sum(y^2)) + xi^2 * (0.5*sum(y^2) - 1/3*sum(y^3)) + O(xi^3)
    const ll =
      -n * logBeta - sumY + xi * (-sumY + 0.5 * sumY2) + xi * xi * (0.5 * sumY2 - sumY3 / 3.0);
    const dLogLDxi = -sumY + 0.5 * sumY2 + xi * (sumY2 - (2.0 * sumY3) / 3.0);
    const dLogLDLogBeta = -n + (1.0 + xi) * (sumY - xi * sumY2);
    return { value: -ll / n, grad: [-dLogLDxi / n, -dLogLDLogBeta / n] };
  }

  let sumLogU = 0.0;
  let sumW = 0.0;
  const invBeta = 1.0 / beta;

  for (let i = 0; i < n; i += 1) {
    const r = (exceedances[i] ?? 0.0) * invBeta;
    const u = 1.0 + xi * r;
    if (!(u > 1e-14)) {
      return { value: 1e15, grad: [0, 0] };
    }
    sumLogU += Math.log(u);
    sumW += r / u;
  }

  const invXi = 1.0 / xi;
  const ll = -n * logBeta - (invXi + 1.0) * sumLogU;
  const dLogLDxi = invXi * invXi * sumLogU - (invXi + 1.0) * sumW;
  const dLogLDLogBeta = -n + (1.0 + xi) * sumW;

  return { value: -ll / n, grad: [-dLogLDxi / n, -dLogLDLogBeta / n] };
}

/** Fit GPD (xi, beta) by Maximum Likelihood on positive exceedances y_i > 0. */
export function fitGpdMle(exceedances: ArrayLike<number>): GpdFitResult {
  const n = exceedances.length;
  if (n < 5) {
    throw new Error(`fitGpdMle requires at least 5 exceedances, got ${n}`);
  }

  let sum = 0.0;
  let yMax = 0.0;
  for (let i = 0; i < n; i += 1) {
    const y = exceedances[i] ?? 0.0;
    sum += y;
    if (y > yMax) yMax = y;
  }
  const yBar = Math.max(sum / n, 1e-12);
  let ss = 0.0;
  for (let i = 0; i < n; i += 1) {
    const d = (exceedances[i] ?? 0.0) - yBar;
    ss += d * d;
  }
  const s2 = Math.max(ss / (n - 1), yBar * yBar * 1e-4);
  const ratio = (yBar * yBar) / s2;

  // Method-of-Moments start clamped cleanly inside strictly feasible region 1 + xi * yMax / beta > 0
  const xi0 = Math.max(-0.25, Math.min(0.6, 0.5 * (1.0 - ratio)));
  let beta0 = Math.max(0.5 * yBar * (1.0 + ratio), 1e-8);
  if (xi0 < 0.0 && 1.0 + (xi0 * yMax) / beta0 <= 0.05) {
    beta0 = Math.max(beta0, (-xi0 * yMax) / 0.9);
  }

  const candidates: Array<[number, number]> = [
    [xi0, Math.log(beta0)],
    [0.15, Math.log(Math.max(yBar * 0.85, (-0.15 * yMax) / 0.9, 1e-8))],
    [0.0, Math.log(yBar)],
  ];

  let bestRes = bfgs(
    (theta) => {
      const xi = theta[0] ?? 0.0;
      const lb = theta[1] ?? 0.0;
      if (xi < -0.49 || xi > 5.0) return { value: 1e15, grad: [0, 0] };
      return gpdNegLogLikAndGrad(exceedances, xi, lb);
    },
    candidates[0] ?? [0.1, Math.log(yBar)],
    { maxIter: 150, gradTol: 1e-9, fTol: 1e-13 },
  );

  for (let c = 1; c < candidates.length && !bestRes.converged; c += 1) {
    const cand = bfgs(
      (theta) => {
        const xi = theta[0] ?? 0.0;
        const lb = theta[1] ?? 0.0;
        if (xi < -0.49 || xi > 5.0) return { value: 1e15, grad: [0, 0] };
        return gpdNegLogLikAndGrad(exceedances, xi, lb);
      },
      candidates[c] ?? [0.0, Math.log(yBar)],
      { maxIter: 150, gradTol: 1e-9, fTol: 1e-13 },
    );
    if (cand.value < bestRes.value) {
      bestRes = cand;
    }
  }

  const xiHat = bestRes.x[0] ?? 0.0;
  const betaHat = Math.exp(bestRes.x[1] ?? 0.0);

  // Asymptotic Fisher information standard errors (Smith 1985 for xi > -0.5):
  // Var([xi, beta]) = (1 + xi)/n * [[1 + xi, -beta], [-beta, 2*beta^2]]
  const seXi = xiHat > -0.5 ? Math.sqrt(((1.0 + xiHat) * (1.0 + xiHat)) / n) : NaN;
  const seBeta = xiHat > -0.5 ? betaHat * Math.sqrt((2.0 * (1.0 + xiHat)) / n) : NaN;

  return {
    xi: xiHat,
    beta: betaHat,
    logLik: -bestRes.value * n,
    nExceed: n,
    converged: bestRes.converged,
    seXi,
    seBeta,
  };
}

/**
 * Hill estimator on the upper k order statistics X_{(1)} >= ... >= X_{(k)} > X_{(k+1)}:
 *   xi_Hill = (1/k) * sum_{i=1}^k log(X_{(i)} / X_{(k+1)})
 * Caller may pass pre-sorted descending values via alreadySortedDesc=true.
 */
export function hillEstimator(
  values: ArrayLike<number>,
  k: number,
  alreadySortedDesc = false,
): { xiHat: number; se: number; threshold: number; k: number } {
  const n = values.length;
  if (k < 1 || k >= n) {
    throw new Error(`hillEstimator requires 1 <= k < n (${n}), got ${k}`);
  }
  const sortedDesc = alreadySortedDesc
    ? values
    : Array.from({ length: n }, (_, i) => values[i] ?? 0.0).sort((a, b) => b - a);

  const anchor = sortedDesc[k] ?? 0.0;
  if (!(anchor > 0.0)) {
    throw new Error(`Hill anchor order statistic X_(k+1) must be positive, got ${anchor}`);
  }
  const logAnchor = Math.log(anchor);
  let sumLog = 0.0;
  for (let i = 0; i < k; i += 1) {
    const v = sortedDesc[i] ?? 0.0;
    sumLog += Math.log(v) - logAnchor;
  }
  const xiHat = sumLog / k;
  return {
    xiHat,
    se: xiHat / Math.sqrt(k),
    threshold: anchor,
    k,
  };
}

/** Run Peaks-Over-Threshold GPD fit at quantile threshold qThreshold and compute VaR/ES levels. */
export function potFromSample(
  values: ArrayLike<number>,
  qThreshold: number,
  varLevels: number[],
  esLevels: number[],
  alreadySortedAsc = false,
): PotSummary {
  const n = values.length;
  const sortedAsc = alreadySortedAsc
    ? values
    : Array.from({ length: n }, (_, i) => values[i] ?? 0.0).sort((a, b) => a - b);

  const u = quantile(sortedAsc, qThreshold, true);
  const exceedances: number[] = [];
  for (let i = 0; i < n; i += 1) {
    const v = sortedAsc[i] ?? 0.0;
    if (v > u) {
      exceedances.push(v - u);
    }
  }
  const nExceed = exceedances.length;
  const pExceed = nExceed / n;
  const gpd = fitGpdMle(exceedances);

  const vars: Record<string, number> = {};
  for (const p of varLevels) {
    vars[p.toString()] = gpdPotVar(u, pExceed, gpd.xi, gpd.beta, p);
  }
  const ess: Record<string, number> = {};
  for (const p of esLevels) {
    const vp = gpdPotVar(u, pExceed, gpd.xi, gpd.beta, p);
    ess[p.toString()] = gpdPotEs(u, vp, gpd.xi, gpd.beta);
  }

  return {
    qThreshold,
    u,
    nTotal: n,
    nExceed,
    pExceed,
    gpd,
    vars,
    ess,
  };
}
