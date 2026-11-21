/**
 * Probability distributions implemented from scratch:
 * - Raw Student-t, standardized Student-t (variance 1), location-scale Student-t
 * - Exact Expected Shortfall (ES) for standardized Student-t & Gaussian
 * - Generalized Pareto Distribution (GPD) for Peaks Over Threshold (POT)
 */

import {
  adaptiveSimpson,
  logBeta,
  logGamma,
  normalPdf,
  normalQuantile,
  regularizedIncompleteBeta,
} from "./special.js";

const SQRT_PI = Math.sqrt(Math.PI);

/** Student-t PDF f(x; nu) for x in R and nu > 0. */
export function studentTPdf(x: number, nu: number): number {
  if (!(nu > 0.0)) {
    throw new Error(`studentTPdf requires nu > 0, got ${nu}`);
  }
  const logNorm = -0.5 * Math.log(nu) - logBeta(0.5 * nu, 0.5);
  return Math.exp(logNorm - 0.5 * (nu + 1.0) * Math.log1p((x * x) / nu));
}

/** Student-t CDF P(T <= x) via regularized incomplete beta (overflow-safe for large |x|). */
export function studentTCdf(x: number, nu: number): number {
  if (!(nu > 0.0)) {
    throw new Error(`studentTCdf requires nu > 0, got ${nu}`);
  }
  if (x === 0.0) return 0.5;
  const ax = Math.abs(x);
  const u = ax > 1e100 ? nu / ax / ax : nu / (nu + ax * ax);
  const tail = 0.5 * regularizedIncompleteBeta(u, 0.5 * nu, 0.5);
  return x > 0.0 ? 1.0 - tail : tail;
}

/** Student-t survival function P(T > x) without upper-tail cancellation (overflow-safe for large |x|). */
export function studentTSf(x: number, nu: number): number {
  if (!(nu > 0.0)) {
    throw new Error(`studentTSf requires nu > 0, got ${nu}`);
  }
  if (x === 0.0) return 0.5;
  const ax = Math.abs(x);
  const u = ax > 1e100 ? nu / ax / ax : nu / (nu + ax * ax);
  const halfBeta = 0.5 * regularizedIncompleteBeta(u, 0.5 * nu, 0.5);
  return x > 0.0 ? halfBeta : 1.0 - halfBeta;
}

/** Upper-tail Student-t quantile solving P(T > x) = qTail for qTail in (0, 0.5] without 1-p cancellation. */
export function studentTSfQuantile(qTail: number, nu: number): number {
  if (!(qTail > 0.0 && qTail <= 0.5)) {
    if (qTail === 0.0) return Infinity;
    throw new Error(`studentTSfQuantile requires qTail in (0, 0.5], got ${qTail}`);
  }
  if (!(nu > 0.0)) {
    throw new Error(`studentTSfQuantile requires nu > 0, got ${nu}`);
  }
  if (qTail === 0.5) return 0.0;
  if (nu === 1.0) {
    return 1.0 / Math.tan(Math.PI * qTail);
  }
  if (nu === 2.0) {
    const a = 0.5 - qTail;
    return (2.0 * a) / Math.sqrt(1.0 - 4.0 * a * a);
  }

  // Initial guess: Cornish-Fisher for central tail (qTail > 1e-4, avoiding 1-qTail underflow), power-law for extreme tail
  let x = 0.0;
  if (qTail > 1e-4) {
    const z = normalQuantile(1.0 - qTail);
    const z3 = z * z * z;
    const z5 = z3 * z * z;
    const g1 = (z3 + z) / (4.0 * nu);
    const g2 = (5.0 * z5 + 16.0 * z3 + 3.0 * z) / (96.0 * nu * nu);
    const g3 = (3.0 * z5 * z * z + 19.0 * z5 + 17.0 * z3 - 15.0 * z) / (384.0 * nu * nu * nu);
    x = Math.max(1e-12, z + g1 + g2 + g3);
  } else {
    const logC =
      logGamma(0.5 * (nu + 1.0)) -
      logGamma(0.5 * nu) -
      0.5 * Math.log(Math.PI) +
      0.5 * nu * Math.log(nu) -
      Math.log(nu);
    x = Math.exp((logC - Math.log(qTail)) / nu);
    if (!Number.isFinite(x) || x <= 0.0) {
      x = 1.0;
    }
  }

  let lo = 0.0;
  let hi = Math.max(x * 1.5, 2.0);
  while (studentTSf(hi, nu) > qTail) {
    lo = hi;
    hi *= 2.0;
    if (hi > 1e140) break;
  }
  if (!(x > lo && x < hi)) {
    x = 0.5 * (lo + hi);
  }

  for (let iter = 0; iter < 60; iter += 1) {
    const sf = studentTSf(x, nu);
    const diff = sf - qTail;
    if (Math.abs(diff) <= 1e-15 * qTail) {
      return x;
    }
    if (diff > 0.0) {
      lo = x;
    } else {
      hi = x;
    }
    const pdf = studentTPdf(x, nu);
    const xNext = pdf > 0.0 ? x + diff / pdf : 0.5 * (lo + hi);
    if (!(Number.isFinite(xNext) && xNext > lo && xNext < hi)) {
      x = 0.5 * (lo + hi);
    } else {
      x = xNext;
    }
    if (hi - lo <= 1e-14 * (1.0 + Math.abs(x))) {
      break;
    }
  }
  return x;
}

/**
 * Student-t quantile T^{-1}(p; nu) solved via Newton-Raphson with strict bisection fallback
 * on the survival function in the upper tail (exact symmetry around p = 0.5).
 */
export function studentTQuantile(p: number, nu: number): number {
  if (!(p > 0.0 && p < 1.0)) {
    if (p === 0.0) return -Infinity;
    if (p === 1.0) return Infinity;
    throw new Error(`studentTQuantile probability out of (0,1): ${p}`);
  }
  if (p === 0.5) return 0.0;
  return p < 0.5 ? -studentTSfQuantile(p, nu) : studentTSfQuantile(1.0 - p, nu);
}

/** Variance scale factor s(nu) = sqrt((nu - 2) / nu) so X = s(nu) * T has variance 1. */
export function stdTScale(nu: number): number {
  if (!(nu > 2.0)) {
    throw new Error(`Standardized Student-t requires nu > 2, got ${nu}`);
  }
  return Math.sqrt((nu - 2.0) / nu);
}

/** Standardized Student-t(nu) PDF with mean 0 and variance 1 (nu > 2). */
export function stdTPdf(z: number, nu: number): number {
  const s = stdTScale(nu);
  return studentTPdf(z / s, nu) / s;
}

/** Standardized Student-t(nu) CDF with mean 0 and variance 1 (nu > 2). */
export function stdTCdf(z: number, nu: number): number {
  const s = stdTScale(nu);
  return studentTCdf(z / s, nu);
}

/** Standardized Student-t(nu) survival function P(Z > z) with variance 1 (nu > 2). */
export function stdTSf(z: number, nu: number): number {
  const s = stdTScale(nu);
  return studentTSf(z / s, nu);
}

/** Standardized Student-t(nu) quantile q_std(p; nu) with mean 0 and variance 1 (nu > 2). */
export function stdTQuantile(p: number, nu: number): number {
  return stdTScale(nu) * studentTQuantile(p, nu);
}

/** Standardized Student-t(nu) upper-tail quantile for survival probability qTail in (0, 0.5]. */
export function stdTQuantileSf(qTail: number, nu: number): number {
  return stdTScale(nu) * studentTSfQuantile(qTail, nu);
}

/**
 * Exact closed-form Expected Shortfall ES_p(Z) = E[Z | Z > q_p] for standardized Student-t(nu):
 *   ES_p(Z) = [f_std(q_p; nu) / (1 - p)] * [(nu - 2 + q_p^2) / (nu - 1)]
 * Derived directly from integral_{q_p}^infty z f_std(z) dz / (1 - p).
 */
export function stdTEs(p: number, nu: number): number {
  if (!(p > 0.0 && p < 1.0)) {
    throw new Error(`stdTEs requires p in (0,1), got ${p}`);
  }
  if (!(nu > 2.0)) {
    throw new Error(`stdTEs requires nu > 2, got ${nu}`);
  }
  const qp = stdTQuantile(p, nu);
  const fqp = stdTPdf(qp, nu);
  return (fqp / (1.0 - p)) * ((nu - 2.0 + qp * qp) / (nu - 1.0));
}

/** Independent numerical ES via adaptive Simpson quadrature of upper-tail quantile integral (v in [0,1], qTail = (1-p)*v^6). */
export function stdTEsNumerical(p: number, nu: number, tol = 1e-11): number {
  return adaptiveSimpson(
    (v) => {
      if (v <= 0.0) return 0.0;
      const v2 = v * v;
      const v5 = v2 * v2 * v;
      const qTail = (1.0 - p) * v5 * v;
      if (qTail <= 1e-300) return 0.0;
      return 6.0 * v5 * stdTQuantileSf(qTail, nu);
    },
    0.0,
    1.0,
    tol,
    16,
  );
}

/** Location-scale Student-t PDF (X = mu + sigma * T_raw or unit-variance when std=true). */
export function locationScaleTPdf(
  x: number,
  mu: number,
  sigma: number,
  nu: number,
  unitVariance = false,
): number {
  const z = (x - mu) / sigma;
  return (unitVariance ? stdTPdf(z, nu) : studentTPdf(z, nu)) / sigma;
}

/** Location-scale Student-t CDF. */
export function locationScaleTCdf(
  x: number,
  mu: number,
  sigma: number,
  nu: number,
  unitVariance = false,
): number {
  const z = (x - mu) / sigma;
  return unitVariance ? stdTCdf(z, nu) : studentTCdf(z, nu);
}

/** Location-scale Student-t quantile. */
export function locationScaleTQuantile(
  p: number,
  mu: number,
  sigma: number,
  nu: number,
  unitVariance = false,
): number {
  return mu + sigma * (unitVariance ? stdTQuantile(p, nu) : studentTQuantile(p, nu));
}

/** Independent numerical Student-t CDF by adaptive Simpson integration from 0 to x. */
export function studentTCdfNumerical(x: number, nu: number, tol = 1e-13): number {
  if (x === 0.0) return 0.5;
  const sign = x >= 0.0 ? 1.0 : -1.0;
  const integral = adaptiveSimpson((t) => studentTPdf(t, nu), 0.0, Math.abs(x), tol, 24);
  return 0.5 + sign * integral;
}

/** Exact standard normal Expected Shortfall ES_p(Z) = phi(Phi^{-1}(p)) / (1 - p). */
export function normalEs(p: number): number {
  const z = normalQuantile(p);
  return normalPdf(z) / (1.0 - p);
}

/** Generalized Pareto Distribution (GPD) PDF for exceedance y >= 0, shape xi, scale beta > 0. */
export function gpdPdf(y: number, xi: number, beta: number): number {
  if (!(beta > 0.0)) throw new Error(`gpdPdf requires beta > 0, got ${beta}`);
  if (y < 0.0) return 0.0;
  if (Math.abs(xi) < 1e-8) {
    return Math.exp(-y / beta) / beta;
  }
  const arg = 1.0 + (xi * y) / beta;
  if (arg <= 0.0) return 0.0;
  return Math.exp(-Math.log(beta) - (1.0 / xi + 1.0) * Math.log(arg));
}

/** GPD CDF G_{xi, beta}(y) for exceedance y >= 0. */
export function gpdCdf(y: number, xi: number, beta: number): number {
  if (!(beta > 0.0)) throw new Error(`gpdCdf requires beta > 0, got ${beta}`);
  if (y <= 0.0) return 0.0;
  if (Math.abs(xi) < 1e-8) {
    return 1.0 - Math.exp(-y / beta);
  }
  const arg = 1.0 + (xi * y) / beta;
  if (arg <= 0.0) return 1.0;
  return 1.0 - Math.pow(arg, -1.0 / xi);
}

/** GPD conditional quantile y_q (exceedance over threshold u) for q in [0, 1). */
export function gpdQuantile(q: number, xi: number, beta: number): number {
  if (!(beta > 0.0)) throw new Error(`gpdQuantile requires beta > 0, got ${beta}`);
  if (!(q >= 0.0 && q < 1.0)) throw new Error(`gpdQuantile requires q in [0, 1), got ${q}`);
  if (q === 0.0) return 0.0;
  if (Math.abs(xi) < 1e-8) {
    return -beta * Math.log1p(-q);
  }
  return (beta / xi) * (Math.pow(1.0 - q, -xi) - 1.0);
}

/**
 * POT VaR at confidence level p > 1 - p_u where p_u = P(X > u) (empirically N_u / n):
 *   VaR_p = u + (beta / xi) * ( (p_u / (1 - p))^xi - 1 )
 */
export function gpdPotVar(u: number, pExceed: number, xi: number, beta: number, p: number): number {
  const ratio = pExceed / (1.0 - p);
  if (Math.abs(xi) < 1e-8) {
    return u + beta * Math.log(ratio);
  }
  return u + (beta / xi) * (Math.pow(ratio, xi) - 1.0);
}

/**
 * POT Expected Shortfall for tail above threshold u (requires xi < 1):
 *   ES_p = (VaR_p + beta - xi * u) / (1 - xi)
 */
export function gpdPotEs(u: number, varP: number, xi: number, beta: number): number {
  if (xi >= 1.0) {
    return Infinity;
  }
  return (varP + beta - xi * u) / (1.0 - xi);
}

export const DIST_CONSTANTS = { SQRT_PI };
