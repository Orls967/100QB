/**
 * GARCH(1,1)-t and Gaussian GARCH(1,1) (QMLE) simulation, likelihood & exact analytical score,
 * variance path reconstruction, BFGS estimation with smooth unconstrained reparameterization,
 * and OPG / Hessian-based Sandwich & Inverse-Information standard errors.
 */

import { bfgs, centralDiffGradient, invertMatrix, type ValueGrad } from "./optim.js";
import { Xoshiro256PlusPlus } from "./rng.js";
import { digamma, logGamma } from "./special.js";

export interface GarchTParams {
  mu: number;
  omega: number;
  alpha: number;
  beta: number;
  nu: number;
}

export interface GarchGaussParams {
  mu: number;
  omega: number;
  alpha: number;
  beta: number;
}

export const TRUE_GARCH_PARAMS: GarchTParams = {
  mu: 0.0003,
  omega: 2e-6,
  alpha: 0.08,
  beta: 0.9,
  nu: 5.0,
};

export interface SimulatedGarchSeries {
  returns: Float64Array;
  sigma2: Float64Array;
  z: Float64Array;
  sigma2NextTrue: number;
}

export interface GarchFitResult {
  params: GarchTParams;
  logLik: number;
  seInfo: GarchTParams;
  seSandwich: GarchTParams;
  converged: boolean;
  iters: number;
  gradNormInf: number;
  hitBoundary: boolean;
  boundaryFlags: string[];
  sigma2Hat: Float64Array;
  residualsHat: Float64Array;
  sigma2NextHat: number;
}

export interface GarchGaussFitResult {
  params: GarchGaussParams;
  logLik: number;
  converged: boolean;
  iters: number;
  hitBoundary: boolean;
  sigma2Hat: Float64Array;
  residualsHat: Float64Array;
  sigma2NextHat: number;
}

export interface GarchTheoreticalMoments {
  unconditionalVariance: number;
  innovationKurtosisZ: number;
  innovationExcessKurtosisZ: number;
  fourthMomentDenomD: number;
  fourthMomentFinite: boolean;
  formalAcfLag1: number;
  formalAcfLag2: number;
}

const LOG_PI = Math.log(Math.PI);
const LOG_SQRT_2PI = 0.5 * Math.log(2.0 * Math.PI);
const NU_CAP = 300.0;

/**
 * Compute theoretical moments for GARCH(1,1)-t:
 *   E[z_t^4] = 3 + 6/(nu - 4) (nu > 4)
 *   E[eps_t^4] < infty iff D = 1 - E[z_t^4]*alpha^2 - 2*alpha*beta - beta^2 > 0
 * When D > 0, Bollerslev/He-Terasvirta ACF of eps_t^2 is:
 *   rho_1 = alpha * (1 - beta^2 - alpha*beta) / (1 - beta^2 - 2*alpha*beta)
 *   rho_2 = (alpha + beta) * rho_1
 */
export function theoreticalGarchMoments(p: GarchTParams): GarchTheoreticalMoments {
  const uncondVar = p.omega / (1.0 - p.alpha - p.beta);
  const excKurt = p.nu > 4.0 ? 6.0 / (p.nu - 4.0) : Infinity;
  const kurtZ = p.nu > 4.0 ? 3.0 + excKurt : Infinity;
  const d =
    p.nu > 4.0
      ? 1.0 - kurtZ * p.alpha * p.alpha - 2.0 * p.alpha * p.beta - p.beta * p.beta
      : -Infinity;
  const a2 = p.alpha * p.alpha;
  const b2 = p.beta * p.beta;
  const ab = p.alpha * p.beta;
  const rho1 = (p.alpha * (1.0 - b2 - ab)) / (1.0 - b2 - 2.0 * ab);
  const rho2 = (p.alpha + p.beta) * rho1;
  // Note: He & Terasvirta (1999) eq. for GARCH(1,1) squared-error ACF depends on E[z^4] ONLY via finiteness condition D>0!
  // Wait: verify whether (1 - beta^2 - 2 alpha beta) vs (1 - E[z^4] alpha^2 - 2 alpha beta - beta^2) cancels:
  // Var(eps^2) = (E[z^4] - 1)*E[h^2], Cov(eps_t^2, eps_{t-1}^2) = E[eps_{t-1}^2 h_t] - E[h]^2 = (w E[h] + alpha E[eps^4] + beta E[h^2]) - E[h]^2.
  // Actually let's compute the exact He-Terasvirta formula:
  //   E[h^2] = (w^2 (1 + a + b) / (1 - a - b)) / D
  //   Var(eps^2) = k_z E[h^2] - E[h]^2
  // Let's keep both forms clear: when E[eps^4] exists (D>0), standard Bollerslev (1988) formula rho_1 = alpha*(1 - b^2 - a*b)/(1 - b^2 - 2*a*b) for Gaussian GARCH;
  // wait: with non-Gaussian k_z != 3, does k_z enter rho_1(eps^2)?
  // Check: Cov(eps_t^2, eps_{t-1}^2) = w E[h] + alpha k_z E[h^2] + beta E[h^2] - E[h]^2.
  // Using E[h]^2 * (1 - k_z a^2 - 2ab - b^2) = ... let's evaluate numerically or exact!
  void a2;
  return {
    unconditionalVariance: uncondVar,
    innovationKurtosisZ: kurtZ,
    innovationExcessKurtosisZ: excKurt,
    fourthMomentDenomD: d,
    fourthMomentFinite: d > 0.0,
    formalAcfLag1: rho1,
    formalAcfLag2: rho2,
  };
}

/**
 * Exact theoretical rho_1(eps^2) and rho_2(eps^2) for GARCH(1,1) with general kurtosis k_z = E[z^4]
 * when D = 1 - k_z*alpha^2 - 2*alpha*beta - beta^2 > 0:
 *   Note that Cov(eps_t^2, eps_{t-1}^2) = w E[eps^2] + (alpha*k_z + beta) E[h^2] - E[h]^2,
 *   while Var(eps_t^2) = k_z E[h^2] - E[h]^2.
 */
export function exactGarchSqAcf(
  alpha: number,
  beta: number,
  kurtZ: number,
): { rho1: number; rho2: number; D: number } {
  const sum = alpha + beta;
  const D = 1.0 - kurtZ * alpha * alpha - 2.0 * alpha * beta - beta * beta;
  if (!(D > 0.0) || !(sum < 1.0)) {
    return { rho1: NaN, rho2: NaN, D };
  }
  // Normalize w = 1 - alpha - beta so E[h] = 1:
  // E[h^2] = ( (1 - a - b)^2 + 2*(1 - a - b)*(a + b) ) / D = (1 - (a+b)^2) / D
  const eh2 = (1.0 - sum * sum) / D;
  const varEps2 = kurtZ * eh2 - 1.0;
  // Since eps_t^2 = h_t * z_t^2 and E[z_t^2]=1,
  // E[eps_t^2 eps_{t-1}^2] = E[(w + alpha eps_{t-1}^2 + beta h_{t-1}) eps_{t-1}^2]
  //                       = (1 - sum)*1 + alpha*E[eps_{t-1}^4] + beta*E[h_{t-1}^2]
  //                       = (1 - sum) + (alpha*kurtZ + beta)*eh2
  const cov1 = 1.0 - sum + (alpha * kurtZ + beta) * eh2 - 1.0;
  const rho1 = cov1 / varEps2;
  const rho2 = sum * rho1;
  return { rho1, rho2, D };
}

/**
 * Simulate GARCH(1,1)-t series with sigma_1^2 = omega / (1 - alpha - beta).
 * Optionally switches to (alphaPost, betaPost) starting at index `breakIdx`.
 */
export function simulateGarchT(
  n: number,
  params: GarchTParams,
  rng: Xoshiro256PlusPlus,
  breakIdx?: number,
  alphaPost?: number,
  betaPost?: number,
): SimulatedGarchSeries {
  const returns = new Float64Array(n);
  const sigma2 = new Float64Array(n);
  const z = new Float64Array(n);

  let h = params.omega / (1.0 - params.alpha - params.beta);
  let epsPrev = 0.0;

  for (let t = 0; t < n; t += 1) {
    const a =
      breakIdx !== undefined && t >= breakIdx && alphaPost !== undefined ? alphaPost : params.alpha;
    const b =
      breakIdx !== undefined && t >= breakIdx && betaPost !== undefined ? betaPost : params.beta;
    if (t > 0) {
      h = params.omega + a * epsPrev * epsPrev + b * h;
    }
    const zt = rng.nextStdT(params.nu);
    const eps = Math.sqrt(h) * zt;
    sigma2[t] = h;
    z[t] = zt;
    returns[t] = params.mu + eps;
    epsPrev = eps;
  }

  const aLast =
    breakIdx !== undefined && n >= breakIdx && alphaPost !== undefined ? alphaPost : params.alpha;
  const bLast =
    breakIdx !== undefined && n >= breakIdx && betaPost !== undefined ? betaPost : params.beta;
  const sigma2NextTrue = params.omega + aLast * epsPrev * epsPrev + bLast * h;

  return { returns, sigma2, z, sigma2NextTrue };
}

/** Deterministic conditional variance reconstruction h_t and standardized residuals z_t. */
export function reconstructSigma2(
  returns: ArrayLike<number>,
  params: { mu: number; omega: number; alpha: number; beta: number },
): { sigma2: Float64Array; residuals: Float64Array; sigma2Next: number } {
  const n = returns.length;
  const sigma2 = new Float64Array(n);
  const residuals = new Float64Array(n);

  const persistence = params.alpha + params.beta;
  let h = params.omega / Math.max(1.0 - persistence, 1e-12);
  let eps = (returns[0] ?? 0.0) - params.mu;
  sigma2[0] = h;
  residuals[0] = eps / Math.sqrt(h);

  for (let t = 1; t < n; t += 1) {
    h = params.omega + params.alpha * eps * eps + params.beta * h;
    eps = (returns[t] ?? 0.0) - params.mu;
    sigma2[t] = h;
    residuals[t] = eps / Math.sqrt(h);
  }

  const sigma2Next = params.omega + params.alpha * eps * eps + params.beta * h;
  return { sigma2, residuals, sigma2Next };
}

/**
 * Total log-likelihood and analytical score [d/dmu, d/domega, d/dalpha, d/dbeta, d/dnu]
 * for GARCH(1,1) with standardized Student-t(nu) innovations (nu > 2).
 */
export function garchTLogLikAndScore(
  returns: ArrayLike<number>,
  params: GarchTParams,
  perObsOut?: number[][],
): { logLik: number; score: [number, number, number, number, number] } {
  const n = returns.length;
  const { mu, omega, alpha, beta, nu } = params;
  const delta = 1.0 - alpha - beta;
  if (!(omega > 0.0 && alpha >= 0.0 && beta >= 0.0 && delta > 1e-12 && nu > 2.000001)) {
    return { logLik: -Infinity, score: [0, 0, 0, 0, 0] };
  }

  const nm2 = nu - 2.0;
  const halfNp1 = 0.5 * (nu + 1.0);
  const constTerm = logGamma(halfNp1) - logGamma(0.5 * nu) - 0.5 * LOG_PI - 0.5 * Math.log(nm2);
  const dConstDnu = 0.5 * digamma(halfNp1) - 0.5 * digamma(0.5 * nu) - 0.5 / nm2;

  let h = omega / delta;
  let dhDmu = 0.0;
  let dhDomega = 1.0 / delta;
  let dhDalpha = h / delta;
  let dhDbeta = h / delta;

  let logLik = 0.0;
  let gMu = 0.0;
  let gOmega = 0.0;
  let gAlpha = 0.0;
  let gBeta = 0.0;
  let gNu = 0.0;

  let epsPrev = 0.0;
  let hPrev = h;

  for (let t = 0; t < n; t += 1) {
    if (t > 0) {
      hPrev = h;
      h = omega + alpha * epsPrev * epsPrev + beta * hPrev;
      dhDmu = -2.0 * alpha * epsPrev + beta * dhDmu;
      dhDomega = 1.0 + beta * dhDomega;
      dhDalpha = epsPrev * epsPrev + beta * dhDalpha;
      dhDbeta = hPrev + beta * dhDbeta;
    }

    const eps = (returns[t] ?? 0.0) - mu;
    const eps2 = eps * eps;
    const z2OverNm2 = eps2 / (nm2 * h);
    const Dt = 1.0 + z2OverNm2;
    const logDt = Math.log1p(z2OverNm2);

    logLik += constTerm - 0.5 * Math.log(h) - halfNp1 * logDt;

    const wt = (nu + 1.0) / (nm2 * Dt);
    const z2 = eps2 / h;
    const dlDmuDirect = (wt * eps) / h;
    const dlDh = (0.5 / h) * (wt * z2 - 1.0);

    const gtMu = dlDmuDirect + dlDh * dhDmu;
    const gtOmega = dlDh * dhDomega;
    const gtAlpha = dlDh * dhDalpha;
    const gtBeta = dlDh * dhDbeta;
    const gtNu = dConstDnu - 0.5 * logDt + (wt * z2) / (2.0 * nm2);

    gMu += gtMu;
    gOmega += gtOmega;
    gAlpha += gtAlpha;
    gBeta += gtBeta;
    gNu += gtNu;

    if (perObsOut) {
      perObsOut.push([gtMu, gtOmega, gtAlpha, gtBeta, gtNu]);
    }

    epsPrev = eps;
  }

  return { logLik, score: [gMu, gOmega, gAlpha, gBeta, gNu] };
}

/**
 * Total log-likelihood and analytical score [d/dmu, d/domega, d/dalpha, d/dbeta]
 * for Gaussian GARCH(1,1) (also used as QMLE).
 */
export function garchGaussLogLikAndScore(
  returns: ArrayLike<number>,
  params: GarchGaussParams,
): { logLik: number; score: [number, number, number, number] } {
  const n = returns.length;
  const { mu, omega, alpha, beta } = params;
  const delta = 1.0 - alpha - beta;
  if (!(omega > 0.0 && alpha >= 0.0 && beta >= 0.0 && delta > 1e-12)) {
    return { logLik: -Infinity, score: [0, 0, 0, 0] };
  }

  let h = omega / delta;
  let dhDmu = 0.0;
  let dhDomega = 1.0 / delta;
  let dhDalpha = h / delta;
  let dhDbeta = h / delta;

  let logLik = 0.0;
  let gMu = 0.0;
  let gOmega = 0.0;
  let gAlpha = 0.0;
  let gBeta = 0.0;

  let epsPrev = 0.0;
  let hPrev = h;

  for (let t = 0; t < n; t += 1) {
    if (t > 0) {
      hPrev = h;
      h = omega + alpha * epsPrev * epsPrev + beta * hPrev;
      dhDmu = -2.0 * alpha * epsPrev + beta * dhDmu;
      dhDomega = 1.0 + beta * dhDomega;
      dhDalpha = epsPrev * epsPrev + beta * dhDalpha;
      dhDbeta = hPrev + beta * dhDbeta;
    }

    const eps = (returns[t] ?? 0.0) - mu;
    const z2 = (eps * eps) / h;

    logLik += -LOG_SQRT_2PI - 0.5 * Math.log(h) - 0.5 * z2;

    const dlDmuDirect = eps / h;
    const dlDh = (0.5 / h) * (z2 - 1.0);

    gMu += dlDmuDirect + dlDh * dhDmu;
    gOmega += dlDh * dhDomega;
    gAlpha += dlDh * dhDalpha;
    gBeta += dlDh * dhDbeta;

    epsPrev = eps;
  }

  return { logLik, score: [gMu, gOmega, gAlpha, gBeta] };
}

function sampleMeanVar(returns: ArrayLike<number>): { mu0: number; v0: number; sMu: number } {
  const n = returns.length;
  let sum = 0.0;
  for (let i = 0; i < n; i += 1) sum += returns[i] ?? 0.0;
  const mu0 = sum / n;
  let ss = 0.0;
  for (let i = 0; i < n; i += 1) {
    const d = (returns[i] ?? 0.0) - mu0;
    ss += d * d;
  }
  const v0 = Math.max(ss / n, 1e-10);
  return { mu0, v0, sMu: Math.sqrt(v0) };
}

function decodeGarchTEta(eta: number[], sMu: number): { params: GarchTParams; dNuDeta: number } {
  const e0 = eta[0] ?? 0.0;
  const e1 = Math.max(-45.0, Math.min(10.0, eta[1] ?? -12.0));
  const e2 = Math.max(-35.0, Math.min(25.0, eta[2] ?? -1.0));
  const e3 = Math.max(-35.0, Math.min(25.0, eta[3] ?? 2.0));
  const e4 = Math.max(-18.0, Math.min(Math.log(NU_CAP - 2.0), eta[4] ?? 1.0));

  const mu = sMu * e0;
  const omega = Math.exp(e1);
  // Numerically stable softmax for (alpha, beta, delta=1-alpha-beta)
  const m = Math.max(0.0, e2, e3);
  const aScaled = Math.exp(e2 - m);
  const bScaled = Math.exp(e3 - m);
  const oneScaled = Math.exp(-m);
  const sum = oneScaled + aScaled + bScaled;
  const alpha = aScaled / sum;
  const beta = bScaled / sum;
  const nu = 2.0 + Math.exp(e4);

  return {
    params: { mu, omega, alpha, beta, nu },
    dNuDeta: nu - 2.0,
  };
}

function encodeGarchTParams(p: GarchTParams, sMu: number): number[] {
  const delta = Math.max(1.0 - p.alpha - p.beta, 1e-7);
  const a = Math.max(p.alpha, 1e-7);
  const b = Math.max(p.beta, 1e-7);
  return [
    p.mu / sMu,
    Math.log(Math.max(p.omega, 1e-14)),
    Math.log(a / delta),
    Math.log(b / delta),
    Math.log(Math.max(Math.min(p.nu, NU_CAP - 1.0) - 2.0, 1e-5)),
  ];
}

function solveGarchTWithStart(returns: ArrayLike<number>, start: GarchTParams, sMu: number) {
  const n = returns.length;
  const eta0 = encodeGarchTParams(start, sMu);

  const obj = (eta: number[]): ValueGrad => {
    const { params, dNuDeta } = decodeGarchTEta(eta, sMu);
    const delta = 1.0 - params.alpha - params.beta;
    if (!(delta > 1e-10 && params.omega > 1e-15)) {
      return { value: 1e12, grad: [0, 0, 0, 0, 0] };
    }
    const { logLik, score } = garchTLogLikAndScore(returns, params);
    if (!Number.isFinite(logLik)) {
      return { value: 1e12, grad: [0, 0, 0, 0, 0] };
    }
    const [gMu, gOmega, gAlpha, gBeta, gNu] = score;
    const invN = -1.0 / n;
    const gEta0 = invN * gMu * sMu;
    const gEta1 = invN * gOmega * params.omega;
    const gEta2 =
      invN * (gAlpha * params.alpha * (1.0 - params.alpha) - gBeta * params.alpha * params.beta);
    const gEta3 =
      invN * (-gAlpha * params.alpha * params.beta + gBeta * params.beta * (1.0 - params.beta));
    const gEta4 = invN * gNu * dNuDeta;
    return {
      value: -logLik / n,
      grad: [gEta0, gEta1, gEta2, gEta3, gEta4],
    };
  };

  const res = bfgs(obj, eta0, { maxIter: 200, gradTol: 1e-7, fTol: 1e-12, xTol: 1e-11 });
  const { params } = decodeGarchTEta(res.x, sMu);
  const finalLL = -res.value * n;
  return { res, params, finalLL };
}

/**
 * Fit GARCH(1,1)-t via BFGS with optional warm start and automatic fallback start.
 * Computes both inverse-Hessian SE and OPG/Hessian Sandwich SE in original parameter space.
 */
export function fitGarchT(
  returns: ArrayLike<number>,
  warmStart?: GarchTParams,
  computeSes = true,
): GarchFitResult {
  const { mu0, v0, sMu } = sampleMeanVar(returns);
  const defaultStart: GarchTParams = {
    mu: mu0,
    omega: 0.03 * v0,
    alpha: 0.07,
    beta: 0.9,
    nu: 7.0,
  };

  let best = warmStart
    ? solveGarchTWithStart(returns, warmStart, sMu)
    : solveGarchTWithStart(returns, defaultStart, sMu);

  if (warmStart && (!best.res.converged || best.res.gradNorm > 1e-4)) {
    const fresh = solveGarchTWithStart(returns, defaultStart, sMu);
    if (fresh.finalLL > best.finalLL) {
      best = fresh;
    }
  } else if (!warmStart && !best.res.converged) {
    const altStart: GarchTParams = {
      mu: mu0,
      omega: 0.08 * v0,
      alpha: 0.12,
      beta: 0.8,
      nu: 4.5,
    };
    const alt = solveGarchTWithStart(returns, altStart, sMu);
    if (alt.finalLL > best.finalLL) {
      best = alt;
    }
  }

  const pHat = best.params;
  const flags: string[] = [];
  if (pHat.omega <= 1e-12) flags.push("omega_lower");
  if (pHat.alpha <= 1e-6) flags.push("alpha_lower");
  if (pHat.beta <= 1e-6) flags.push("beta_lower");
  if (pHat.alpha + pHat.beta >= 1.0 - 1e-5) flags.push("persistence_upper");
  if (pHat.nu <= 2.02) flags.push("nu_lower");
  if (pHat.nu >= NU_CAP - 10.0) flags.push("nu_upper");

  let seInfo: GarchTParams = { mu: NaN, omega: NaN, alpha: NaN, beta: NaN, nu: NaN };
  let seSandwich: GarchTParams = { mu: NaN, omega: NaN, alpha: NaN, beta: NaN, nu: NaN };

  if (computeSes) {
    const perObs: number[][] = [];
    garchTLogLikAndScore(returns, pHat, perObs);
    // OPG matrix B = sum g_t g_t^T
    const B: number[][] = Array.from({ length: 5 }, () => new Array<number>(5).fill(0.0));
    for (let t = 0; t < perObs.length; t += 1) {
      const gt = perObs[t];
      if (!gt) continue;
      for (let i = 0; i < 5; i += 1) {
        const bi = B[i];
        const gi = gt[i] ?? 0.0;
        if (!bi) continue;
        for (let j = 0; j < 5; j += 1) {
          bi[j] = (bi[j] ?? 0.0) + gi * (gt[j] ?? 0.0);
        }
      }
    }

    // Negative Hessian A = -d^2 L / dtheta^2 via central differences of analytical total score
    const thetaHat = [pHat.mu, pHat.omega, pHat.alpha, pHat.beta, pHat.nu];
    const hSteps = [
      1e-6 * Math.max(Math.abs(pHat.mu), 1e-3),
      1e-5 * Math.max(pHat.omega, 1e-7),
      1e-5 * Math.max(pHat.alpha, 1e-3),
      1e-5 * Math.max(pHat.beta, 1e-3),
      1e-5 * Math.max(pHat.nu, 0.5),
    ];
    // Clip step so nu stays > 2 and persistence stays < 1
    hSteps[1] = Math.min(hSteps[1] ?? 1e-7, 0.25 * pHat.omega);
    hSteps[2] = Math.min(
      hSteps[2] ?? 1e-5,
      0.25 * (1.0 - pHat.alpha - pHat.beta),
      0.25 * pHat.alpha,
    );
    hSteps[3] = Math.min(
      hSteps[3] ?? 1e-5,
      0.25 * (1.0 - pHat.alpha - pHat.beta),
      0.25 * pHat.beta,
    );
    hSteps[4] = Math.min(hSteps[4] ?? 1e-4, 0.25 * (pHat.nu - 2.0));

    const A: number[][] = Array.from({ length: 5 }, () => new Array<number>(5).fill(0.0));
    for (let j = 0; j < 5; j += 1) {
      const hj = hSteps[j] ?? 1e-6;
      const tp = [...thetaHat];
      const tm = [...thetaHat];
      tp[j] = (tp[j] ?? 0.0) + hj;
      tm[j] = (tm[j] ?? 0.0) - hj;
      const gp = garchTLogLikAndScore(returns, {
        mu: tp[0] ?? 0,
        omega: tp[1] ?? 0,
        alpha: tp[2] ?? 0,
        beta: tp[3] ?? 0,
        nu: tp[4] ?? 5,
      }).score;
      const gm = garchTLogLikAndScore(returns, {
        mu: tm[0] ?? 0,
        omega: tm[1] ?? 0,
        alpha: tm[2] ?? 0,
        beta: tm[3] ?? 0,
        nu: tm[4] ?? 5,
      }).score;
      for (let i = 0; i < 5; i += 1) {
        const ai = A[i];
        if (ai) ai[j] = -((gp[i] ?? 0.0) - (gm[i] ?? 0.0)) / (2.0 * hj);
      }
    }
    // Symmetrize A
    for (let i = 0; i < 5; i += 1) {
      for (let j = i + 1; j < 5; j += 1) {
        const avg = 0.5 * ((A[i]?.[j] ?? 0.0) + (A[j]?.[i] ?? 0.0));
        const ri = A[i];
        const rj = A[j];
        if (ri) ri[j] = avg;
        if (rj) rj[i] = avg;
      }
    }

    try {
      const invA = invertMatrix(A);
      // Sandwich = invA * B * invA
      const invAB: number[][] = Array.from({ length: 5 }, (_, i) =>
        Array.from({ length: 5 }, (_, j) => {
          let s = 0.0;
          for (let k = 0; k < 5; k += 1) {
            s += (invA[i]?.[k] ?? 0.0) * (B[k]?.[j] ?? 0.0);
          }
          return s;
        }),
      );
      const sand: number[][] = Array.from({ length: 5 }, (_, i) =>
        Array.from({ length: 5 }, (_, j) => {
          let s = 0.0;
          for (let k = 0; k < 5; k += 1) {
            s += (invAB[i]?.[k] ?? 0.0) * (invA[k]?.[j] ?? 0.0);
          }
          return s;
        }),
      );
      seInfo = {
        mu: Math.sqrt(Math.max(invA[0]?.[0] ?? NaN, 0.0)),
        omega: Math.sqrt(Math.max(invA[1]?.[1] ?? NaN, 0.0)),
        alpha: Math.sqrt(Math.max(invA[2]?.[2] ?? NaN, 0.0)),
        beta: Math.sqrt(Math.max(invA[3]?.[3] ?? NaN, 0.0)),
        nu: Math.sqrt(Math.max(invA[4]?.[4] ?? NaN, 0.0)),
      };
      seSandwich = {
        mu: Math.sqrt(Math.max(sand[0]?.[0] ?? NaN, 0.0)),
        omega: Math.sqrt(Math.max(sand[1]?.[1] ?? NaN, 0.0)),
        alpha: Math.sqrt(Math.max(sand[2]?.[2] ?? NaN, 0.0)),
        beta: Math.sqrt(Math.max(sand[3]?.[3] ?? NaN, 0.0)),
        nu: Math.sqrt(Math.max(sand[4]?.[4] ?? NaN, 0.0)),
      };
    } catch {
      // Singular Hessian left as NaN and recorded honestly
    }
  }

  const { sigma2, residuals, sigma2Next } = reconstructSigma2(returns, pHat);

  return {
    params: pHat,
    logLik: best.finalLL,
    seInfo,
    seSandwich,
    converged: best.res.converged,
    iters: best.res.iters,
    gradNormInf: best.res.gradNorm,
    hitBoundary: flags.length > 0,
    boundaryFlags: flags,
    sigma2Hat: sigma2,
    residualsHat: residuals,
    sigma2NextHat: sigma2Next,
  };
}

/** Fit Gaussian GARCH(1,1) (QMLE) via BFGS with optional warm start. */
export function fitGarchGauss(
  returns: ArrayLike<number>,
  warmStart?: GarchGaussParams,
): GarchGaussFitResult {
  const n = returns.length;
  const { mu0, v0, sMu } = sampleMeanVar(returns);
  const start = warmStart ?? {
    mu: mu0,
    omega: 0.03 * v0,
    alpha: 0.08,
    beta: 0.9,
  };

  const decode = (eta: number[]): GarchGaussParams => {
    const e1 = Math.max(-45.0, Math.min(10.0, eta[1] ?? -12.0));
    const e2 = Math.max(-35.0, Math.min(25.0, eta[2] ?? -1.0));
    const e3 = Math.max(-35.0, Math.min(25.0, eta[3] ?? 2.0));
    const m = Math.max(0.0, e2, e3);
    const aS = Math.exp(e2 - m);
    const bS = Math.exp(e3 - m);
    const oS = Math.exp(-m);
    const sum = oS + aS + bS;
    return {
      mu: sMu * (eta[0] ?? 0.0),
      omega: Math.exp(e1),
      alpha: aS / sum,
      beta: bS / sum,
    };
  };

  const encode = (p: GarchGaussParams): number[] => {
    const delta = Math.max(1.0 - p.alpha - p.beta, 1e-7);
    return [
      p.mu / sMu,
      Math.log(Math.max(p.omega, 1e-14)),
      Math.log(Math.max(p.alpha, 1e-7) / delta),
      Math.log(Math.max(p.beta, 1e-7) / delta),
    ];
  };

  const obj = (eta: number[]): ValueGrad => {
    const p = decode(eta);
    if (!(1.0 - p.alpha - p.beta > 1e-10 && p.omega > 1e-15)) {
      return { value: 1e12, grad: [0, 0, 0, 0] };
    }
    const { logLik, score } = garchGaussLogLikAndScore(returns, p);
    if (!Number.isFinite(logLik)) {
      return { value: 1e12, grad: [0, 0, 0, 0] };
    }
    const [gMu, gOmega, gAlpha, gBeta] = score;
    const invN = -1.0 / n;
    return {
      value: -logLik / n,
      grad: [
        invN * gMu * sMu,
        invN * gOmega * p.omega,
        invN * (gAlpha * p.alpha * (1.0 - p.alpha) - gBeta * p.alpha * p.beta),
        invN * (-gAlpha * p.alpha * paramsOrZero(p.beta) + gBeta * p.beta * (1.0 - p.beta)),
      ],
    };
  };

  let res = bfgs(obj, encode(start), { maxIter: 200, gradTol: 1e-7, fTol: 1e-12, xTol: 1e-11 });
  let pHat = decode(res.x);
  let finalLL = -res.value * n;

  if (warmStart && (!res.converged || res.gradNorm > 1e-4)) {
    const fresh = bfgs(obj, encode({ mu: mu0, omega: 0.03 * v0, alpha: 0.08, beta: 0.9 }), {
      maxIter: 200,
      gradTol: 1e-7,
    });
    const freshLL = -fresh.value * n;
    if (freshLL > finalLL) {
      res = fresh;
      pHat = decode(res.x);
      finalLL = freshLL;
    }
  }

  const hitBoundary =
    pHat.omega <= 1e-12 ||
    pHat.alpha <= 1e-6 ||
    pHat.beta <= 1e-6 ||
    pHat.alpha + pHat.beta >= 1.0 - 1e-5;
  const { sigma2, residuals, sigma2Next } = reconstructSigma2(returns, pHat);

  return {
    params: pHat,
    logLik: finalLL,
    converged: res.converged,
    iters: res.iters,
    hitBoundary,
    sigma2Hat: sigma2,
    residualsHat: residuals,
    sigma2NextHat: sigma2Next,
  };
}

function paramsOrZero(v: number): number {
  return Number.isFinite(v) ? v : 0.0;
}

/** Verify analytical total score against central finite differences. */
export function scoreVsCentralFd(
  returns: ArrayLike<number>,
  params: GarchTParams,
): { maxRelErr: number; relErrs: number[]; analytical: number[]; fd: number[] } {
  const { score } = garchTLogLikAndScore(returns, params);
  const theta = [params.mu, params.omega, params.alpha, params.beta, params.nu];
  const steps = [1e-7, 1e-10, 1e-7, 1e-7, 1e-5];
  const fd = centralDiffGradient(
    (x) =>
      garchTLogLikAndScore(returns, {
        mu: x[0] ?? 0,
        omega: x[1] ?? 1e-6,
        alpha: x[2] ?? 0.05,
        beta: x[3] ?? 0.9,
        nu: x[4] ?? 5,
      }).logLik,
    theta,
    steps,
  );
  const relErrs = score.map((s, i) => {
    const f = fd[i] ?? 0.0;
    return Math.abs(s - f) / Math.max(Math.abs(f), 1.0);
  });
  return {
    maxRelErr: Math.max(...relErrs),
    relErrs,
    analytical: [...score],
    fd,
  };
}
