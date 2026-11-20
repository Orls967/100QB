/** High-precision statistical helpers: Cody erf, normal CDF, moments, and log-log OLS slope. */

const SQRT1_2 = Math.SQRT1_2;
const CI_Z_95 = 1.96;

export interface EstimatorResult {
  price: number;
  se: number;
  ci_low: number;
  ci_high: number;
}

const A0 = 3.1611237438705655;
const A1 = 113.86415415105016;
const A2 = 377.485237685302;
const A3 = 3209.3775891384694;
const A4 = 0.18577770618460315;

const B0 = 23.60129095234412;
const B1 = 244.02463793444417;
const B2 = 1282.6165260773723;
const B3 = 2844.236833439171;

const CODY_C = [
  0.5641884969886701, 8.88314979438838, 66.11919063714163, 298.6351381974001, 881.952221241769,
  1712.047612634071, 2051.078377826071, 1230.339354797997, 2.153115354744039e-8,
] as const;
const CODY_D = [
  15.744926110709835, 117.6939508913125, 537.1811018620099, 1621.38957456669, 3290.79923573346,
  4362.619090143247, 3439.367674143722, 1230.339354803749,
] as const;

const CODY_P = [
  0.3053266349612323, 0.3603448999498044, 0.1257817261112292, 0.01608378514874228,
  6.58749161529838e-4, 0.0163153871373021,
] as const;
const CODY_Q = [
  2.568520192289822, 1.87295284992346, 0.527905102951428, 0.06051834131244132, 0.002335204976268692,
] as const;

/** W. J. Cody rational Chebyshev approximation to erf(x) (relative error < 1e-15). */
export function erf(x: number): number {
  if (Number.isNaN(x)) return Number.NaN;
  const absX = Math.abs(x);
  let res: number;

  if (absX < 0.46875) {
    const x2 = absX * absX;
    const num = (((A4 * x2 + A0) * x2 + A1) * x2 + A2) * x2 + A3;
    const den = (((x2 + B0) * x2 + B1) * x2 + B2) * x2 + B3;
    res = absX * (num / den);
  } else if (absX < 4.0) {
    let num = CODY_C[8];
    let den = 1.0;
    for (let i = 0; i < 8; i += 1) {
      num = num * absX + (CODY_C[i] ?? 0);
      den = den * absX + (CODY_D[i] ?? 0);
    }
    const erfcVal = Math.exp(-absX * absX) * (num / den);
    res = 1.0 - erfcVal;
  } else if (absX < 26.5) {
    const invX2 = 1.0 / (absX * absX);
    let num = CODY_P[5];
    let den = 1.0;
    for (let i = 0; i < 5; i += 1) {
      num = num * invX2 + (CODY_P[i] ?? 0);
      den = den * invX2 + (CODY_Q[i] ?? 0);
    }
    const erfcVal = (Math.exp(-absX * absX) / absX) * (0.5641895835477563 - invX2 * (num / den));
    res = 1.0 - erfcVal;
  } else {
    res = 1.0;
  }

  return x < 0.0 ? -res : res;
}

export function normCdf(x: number): number {
  return 0.5 * (1.0 + erf(x * SQRT1_2));
}

export function summarizeMoments(count: number, sum: number, sumSq: number): EstimatorResult {
  const mean = sum / count;
  const variance = Math.max((sumSq - count * mean * mean) / (count - 1), 0.0);
  const se = Math.sqrt(variance / count);
  return {
    price: mean,
    se,
    ci_low: mean - CI_Z_95 * se,
    ci_high: mean + CI_Z_95 * se,
  };
}

export function fitLogLogSlope(ns: readonly number[], ses: readonly number[]): number {
  if (ns.length !== ses.length || ns.length < 2) {
    throw new Error("ns and ses must have equal length >= 2");
  }
  const n = ns.length;
  const xs = new Float64Array(n);
  const ys = new Float64Array(n);
  let sumX = 0.0;
  let sumY = 0.0;
  for (let i = 0; i < n; i += 1) {
    const x = Math.log10(ns[i] ?? 1);
    const y = Math.log10(ses[i] ?? 1);
    xs[i] = x;
    ys[i] = y;
    sumX += x;
    sumY += y;
  }
  const mx = sumX / n;
  const my = sumY / n;
  let num = 0.0;
  let den = 0.0;
  for (let i = 0; i < n; i += 1) {
    const dx = (xs[i] ?? 0) - mx;
    num += dx * ((ys[i] ?? 0) - my);
    den += dx * dx;
  }
  if (den <= 0.0) {
    throw new Error("zero variance in x");
  }
  return num / den;
}
