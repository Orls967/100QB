/**
 * Conditional one-day-ahead loss VaR and Expected Shortfall estimators:
 * - True GARCH-t(nu=5) loss VaR / ES (positive loss convention L_t = -r_t)
 * - (a) Gaussian
 * - (b) t-parametric (estimated nu_hat)
 * - (c) Filtered Historical Simulation (FHS empirical quantiles & tail mean)
 * - (d) FHS-EVT (GPD POT at threshold quantile q_u = 0.90)
 */

import { gpdPotEs, gpdPotVar, normalEs, stdTEs, stdTQuantile } from "./distributions.js";
import { fitGpdMle } from "./evt.js";
import { TRUE_GARCH_PARAMS, type GarchTParams } from "./garch.js";
import { normalQuantile } from "./special.js";
import { quantile } from "./stats.js";

export type RiskMethodName = "gaussian" | "t_parametric" | "fhs" | "fhs_evt";

export const RISK_VAR_LEVELS = [0.95, 0.99, 0.999] as const;
export const RISK_ES_LEVELS = [0.975, 0.99] as const;

export interface MethodVarEs {
  var95: number;
  var99: number;
  var999: number;
  es975: number;
  es99: number;
}

export interface ThresholdEvtCheck {
  q: number;
  u: number;
  nExceed: number;
  xiHat: number;
  betaHat: number;
  var99: number;
  var999: number;
  es975: number;
  es99: number;
}

export interface SampleRiskEstimates {
  trueValues: MethodVarEs;
  methods: Record<RiskMethodName, MethodVarEs>;
  evtXiHat: number;
  invNuHat: number;
  thresholdSensitivity: ThresholdEvtCheck[];
}

/**
 * Compute one-day-ahead positive loss VaR and ES for standardized t(nu) conditional volatility:
 *   VaR_p(t) = -mu + sigma_t * q_std(p; nu)
 *   ES_p(t)  = -mu + sigma_t * ES_std(p; nu)
 */
export function lossVarEsFromTParams(p: GarchTParams, sigmaNext: number): MethodVarEs {
  return {
    var95: -p.mu + sigmaNext * stdTQuantile(0.95, p.nu),
    var99: -p.mu + sigmaNext * stdTQuantile(0.99, p.nu),
    var999: -p.mu + sigmaNext * stdTQuantile(0.999, p.nu),
    es975: -p.mu + sigmaNext * stdTEs(0.975, p.nu),
    es99: -p.mu + sigmaNext * stdTEs(0.99, p.nu),
  };
}

/** Conditional Gaussian VaR/ES given (mu, sigmaNext). */
export function lossVarEsFromGaussian(mu: number, sigmaNext: number): MethodVarEs {
  return {
    var95: -mu + sigmaNext * normalQuantile(0.95),
    var99: -mu + sigmaNext * normalQuantile(0.99),
    var999: -mu + sigmaNext * normalQuantile(0.999),
    es975: -mu + sigmaNext * normalEs(0.975),
    es99: -mu + sigmaNext * normalEs(0.99),
  };
}

/** Empirical tail mean E[X | X > q_p] using top exceedances (falls back to >= q_p if needed). */
function empiricalTailMean(sortedAsc: number[], thresh: number): number {
  let sum = 0.0;
  let cnt = 0;
  for (let i = sortedAsc.length - 1; i >= 0; i -= 1) {
    const v = sortedAsc[i] ?? 0.0;
    if (v > thresh) {
      sum += v;
      cnt += 1;
    } else {
      break;
    }
  }
  if (cnt === 0) {
    return sortedAsc[sortedAsc.length - 1] ?? thresh;
  }
  return sum / cnt;
}

/** Compute FHS-EVT VaR & ES on pre-sorted ascending standardized loss residuals ell_t = -z_t. */
export function evtLossQuantilesFromSortedLosses(
  sortedLossResAsc: number[],
  qThreshold: number,
  muHat: number,
  sigmaNextHat: number,
): ThresholdEvtCheck {
  const n = sortedLossResAsc.length;
  const u = quantile(sortedLossResAsc, qThreshold, true);
  const exc: number[] = [];
  for (let i = 0; i < n; i += 1) {
    const v = sortedLossResAsc[i] ?? 0.0;
    if (v > u) exc.push(v - u);
  }
  const nExceed = exc.length;
  const pExceed = nExceed / n;
  const gpd = fitGpdMle(exc);
  const q99 = gpdPotVar(u, pExceed, gpd.xi, gpd.beta, 0.99);
  const q999 = gpdPotVar(u, pExceed, gpd.xi, gpd.beta, 0.999);
  const q975 = gpdPotVar(u, pExceed, gpd.xi, gpd.beta, 0.975);
  const e975 = gpdPotEs(u, q975, gpd.xi, gpd.beta);
  const e99 = gpdPotEs(u, q99, gpd.xi, gpd.beta);

  return {
    q: qThreshold,
    u,
    nExceed,
    xiHat: gpd.xi,
    betaHat: gpd.beta,
    var99: -muHat + sigmaNextHat * q99,
    var999: -muHat + sigmaNextHat * q999,
    es975: -muHat + sigmaNextHat * e975,
    es99: -muHat + sigmaNextHat * e99,
  };
}

/**
 * Evaluate all four VaR/ES methods (Gaussian, t-parametric, FHS, FHS-EVT) on GARCH-t filtered
 * standardized residuals z_hat, plus GPD threshold sensitivity at q in {0.85, 0.90, 0.95}.
 */
export function estimateRiskMethodsOnSample(
  fitParams: GarchTParams,
  residualsHat: ArrayLike<number>,
  sigmaNextHat: number,
  sigmaNextTrue: number,
  trueParams: GarchTParams = TRUE_GARCH_PARAMS,
): SampleRiskEstimates {
  const trueValues = lossVarEsFromTParams(trueParams, sigmaNextTrue);
  const gaussian = lossVarEsFromGaussian(fitParams.mu, sigmaNextHat);
  const tParametric = lossVarEsFromTParams(fitParams, sigmaNextHat);

  const n = residualsHat.length;
  // Loss standardized residuals ell_t = -(r_t - mu_hat) / sigma_hat = -z_hat_t
  const sortedLossAsc = new Array<number>(n);
  for (let i = 0; i < n; i += 1) {
    sortedLossAsc[i] = -(residualsHat[i] ?? 0.0);
  }
  sortedLossAsc.sort((a, b) => a - b);

  const q95 = quantile(sortedLossAsc, 0.95, true);
  const q975 = quantile(sortedLossAsc, 0.975, true);
  const q99 = quantile(sortedLossAsc, 0.99, true);
  const q999 = quantile(sortedLossAsc, 0.999, true);

  const fhs: MethodVarEs = {
    var95: -fitParams.mu + sigmaNextHat * q95,
    var99: -fitParams.mu + sigmaNextHat * q99,
    var999: -fitParams.mu + sigmaNextHat * q999,
    es975: -fitParams.mu + sigmaNextHat * empiricalTailMean(sortedLossAsc, q975),
    es99: -fitParams.mu + sigmaNextHat * empiricalTailMean(sortedLossAsc, q99),
  };

  // Threshold sensitivity at q in {0.85, 0.90, 0.95}; baseline FHS-EVT uses q = 0.90
  const evt85 = evtLossQuantilesFromSortedLosses(sortedLossAsc, 0.85, fitParams.mu, sigmaNextHat);
  const evt90 = evtLossQuantilesFromSortedLosses(sortedLossAsc, 0.9, fitParams.mu, sigmaNextHat);
  const evt95 = evtLossQuantilesFromSortedLosses(sortedLossAsc, 0.95, fitParams.mu, sigmaNextHat);

  // Also compute p=0.95 VaR for FHS-EVT at q=0.90 via GPD POT
  const u90 = evt90.u;
  const pExc90 = evt90.nExceed / n;
  const q95Evt = gpdPotVar(u90, pExc90, evt90.xiHat, evt90.betaHat, 0.95);

  const fhsEvt: MethodVarEs = {
    var95: -fitParams.mu + sigmaNextHat * q95Evt,
    var99: evt90.var99,
    var999: evt90.var999,
    es975: evt90.es975,
    es99: evt90.es99,
  };

  return {
    trueValues,
    methods: {
      gaussian,
      t_parametric: tParametric,
      fhs,
      fhs_evt: fhsEvt,
    },
    evtXiHat: evt90.xiHat,
    invNuHat: 1.0 / fitParams.nu,
    thresholdSensitivity: [evt85, evt90, evt95],
  };
}
