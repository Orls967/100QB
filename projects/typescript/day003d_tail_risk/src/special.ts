/**
 * Special mathematical functions and numerical quadrature implemented from scratch:
 * log-gamma, digamma, trigamma, regularized incomplete beta, error functions,
 * normal PDF/CDF/SF/quantile, chi-square survival, and adaptive Simpson quadrature.
 */

const LOG_SQRT_2PI = 0.5 * Math.log(2 * Math.PI);

// Lanczos g=7, n=9 coefficients (GNU Scientific Library / Boost standard double precision)
const LANCZOS_G = 7;
const LANCZOS_COEFF = [
  0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406,
  12.507343278686905, -0.13857109526572012, 0.000009984369578019572, 1.5056327351493116e-7,
];

/** Natural logarithm of Gamma(z) for real z > 0 (or non-integer via reflection). */
export function logGamma(z: number): number {
  if (!(z > 0)) {
    throw new Error(`logGamma requires positive argument, got ${z}`);
  }
  if (z < 0.5) {
    return Math.log(Math.PI) - Math.log(Math.abs(Math.sin(Math.PI * z))) - logGamma(1 - z);
  }
  const zm1 = z - 1;
  let x = LANCZOS_COEFF[0] ?? 1;
  for (let i = 1; i < LANCZOS_COEFF.length; i += 1) {
    x += (LANCZOS_COEFF[i] ?? 0) / (zm1 + i);
  }
  const t = zm1 + LANCZOS_G + 0.5;
  return LOG_SQRT_2PI + (zm1 + 0.5) * Math.log(t) - t + Math.log(x);
}

/** Digamma function psi(x) = d/dx log(Gamma(x)) for x > 0. */
export function digamma(xInput: number): number {
  if (!(xInput > 0)) {
    throw new Error(`digamma requires x > 0, got ${xInput}`);
  }
  let x = xInput;
  let result = 0;
  while (x < 10) {
    result -= 1 / x;
    x += 1;
  }
  const invX = 1 / x;
  const invX2 = invX * invX;
  result += Math.log(x) - 0.5 * invX;
  // Asymptotic Bernoulli series: -sum B_{2k} / (2k * x^{2k})
  result -=
    invX2 *
    (1 / 12 -
      invX2 *
        (1 / 120 -
          invX2 * (1 / 252 - invX2 * (1 / 240 - invX2 * (1 / 132 - invX2 * (691 / 32760))))));
  return result;
}

/** Trigamma function psi'(x) = d^2/dx^2 log(Gamma(x)) for x > 0. */
export function trigamma(xInput: number): number {
  if (!(xInput > 0)) {
    throw new Error(`trigamma requires x > 0, got ${xInput}`);
  }
  let x = xInput;
  let result = 0;
  while (x < 10) {
    result += 1 / (x * x);
    x += 1;
  }
  const invX = 1 / x;
  const invX2 = invX * invX;
  result +=
    invX +
    0.5 * invX2 +
    invX *
      invX2 *
      (1 / 6 - invX2 * (1 / 30 - invX2 * (1 / 42 - invX2 * (1 / 30 - invX2 * (5 / 66)))));
  return result;
}

/** Log Beta function log B(a, b) for a, b > 0. */
export function logBeta(a: number, b: number): number {
  return logGamma(a) + logGamma(b) - logGamma(a + b);
}

/** Evaluate continued fraction for I_x(a, b) via modified Lentz method (Numerical Recipes). */
function betacf(a: number, b: number, x: number): number {
  const MAXIT = 300;
  const EPS = 3e-16;
  const FPMIN = 1e-300;
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - (qab * x) / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m += 1) {
    const m2 = 2 * m;
    let aa = (m * (b - m) * x) / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    h *= d * c;
    aa = (-(a + m) * (qab + m) * x) / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) <= EPS) {
      break;
    }
  }
  return h;
}

/** Regularized incomplete beta function I_x(a, b) for x in [0, 1], a > 0, b > 0. */
export function regularizedIncompleteBeta(x: number, a: number, b: number): number {
  if (!(a > 0 && b > 0)) {
    throw new Error(`Incomplete beta requires a, b > 0, got a=${a}, b=${b}`);
  }
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log(1 - x),
  );
  if (x < (a + 1) / (a + b + 2)) {
    return (bt * betacf(a, b, x)) / a;
  }
  return 1 - (bt * betacf(b, a, 1 - x)) / b;
}

// Cody (1969) rational approximations for erf/erfc (16-digit IEEE-754 double literals)
const ERF_A = [
  3.1611237438705655, 113.86415415105016, 377.485237685302, 3209.3775891384694, 0.18577770618460315,
];
const ERF_B = [23.60129095234412, 244.02463793444417, 1282.6165260773723, 2844.236833439171];

const ERFC_C = [
  0.5641884969886701, 8.883149794388377, 66.11919063714163, 298.6351381974001, 881.952221241769,
  1712.0476126340707, 2051.0783778260716, 1230.3393547979972, 2.1531153547440383e-8,
];
const ERFC_D = [
  15.744926110709835, 117.6939508913125, 537.1811018620099, 1621.3895745666903, 3290.7992357334597,
  4362.619090143247, 3439.3676741437216, 1230.3393548037493,
];

const ERFC_P = [
  0.30532663496123236, 0.36034489994980445, 0.12578172611122926, 0.01608378514874228,
  0.0006587491615298378, 0.016315387137302097,
];
const ERFC_Q = [
  2.568520192289822, 1.8729528499234604, 0.5279051029514285, 0.06051834131244132,
  0.002335204976268692,
];

const SQRPI = 0.5641895835477563; // 1 / sqrt(pi)

function codyCalerf(arg: number, jint: 0 | 1): number {
  const x = Math.abs(arg);
  let result = 0;

  if (x <= 0.46875) {
    const xsq = x > 1.11e-16 ? x * x : 0;
    let xnum = (ERF_A[4] ?? 0) * xsq;
    let xden = xsq;
    for (let i = 0; i < 3; i += 1) {
      xnum = (xnum + (ERF_A[i] ?? 0)) * xsq;
      xden = (xden + (ERF_B[i] ?? 0)) * xsq;
    }
    result = (x * (xnum + (ERF_A[3] ?? 0))) / (xden + (ERF_B[3] ?? 0));
    if (jint !== 0) result = 1 - result;
    if (jint === 0 && arg < 0) result = -result;
    return result;
  }

  if (x <= 4) {
    let xnum = (ERFC_C[8] ?? 0) * x;
    let xden = x;
    for (let i = 0; i < 7; i += 1) {
      xnum = (xnum + (ERFC_C[i] ?? 0)) * x;
      xden = (xden + (ERFC_D[i] ?? 0)) * x;
    }
    result = (xnum + (ERFC_C[7] ?? 0)) / (xden + (ERFC_D[7] ?? 0));
    const xsq = Math.trunc(x * 16) / 16;
    const del = (x - xsq) * (x + xsq);
    result *= Math.exp(-xsq * xsq) * Math.exp(-del);
  } else {
    if (x >= 26.6) {
      result = 0;
    } else {
      const xsq = 1 / (x * x);
      let xnum = (ERFC_P[5] ?? 0) * xsq;
      let xden = xsq;
      for (let i = 0; i < 4; i += 1) {
        xnum = (xnum + (ERFC_P[i] ?? 0)) * xsq;
        xden = (xden + (ERFC_Q[i] ?? 0)) * xsq;
      }
      result = (xsq * (xnum + (ERFC_P[4] ?? 0))) / (xden + (ERFC_Q[4] ?? 0));
      result = (SQRPI - result) / x;
      const xsq2 = Math.trunc(x * 16) / 16;
      const del = (x - xsq2) * (x + xsq2);
      result *= Math.exp(-xsq2 * xsq2) * Math.exp(-del);
    }
  }

  if (jint === 0) {
    result = 1 - result;
    if (arg < 0) result = -result;
  } else if (arg < 0) {
    result = 2 - result;
  }
  return result;
}

/** Error function erf(x) to machine precision. */
export function erf(x: number): number {
  return codyCalerf(x, 0);
}

/** Complementary error function erfc(x) = 1 - erf(x) without tail cancellation. */
export function erfc(x: number): number {
  return codyCalerf(x, 1);
}

/** Standard normal PDF phi(x). */
export function normalPdf(x: number): number {
  return Math.exp(-0.5 * x * x - LOG_SQRT_2PI);
}

/** Standard normal CDF Phi(x). */
export function normalCdf(x: number): number {
  return 0.5 * erfc(-x * Math.SQRT1_2);
}

/** Standard normal survival function 1 - Phi(x) without cancellation in upper tail. */
export function normalSf(x: number): number {
  return 0.5 * erfc(x * Math.SQRT1_2);
}

/** Standard normal quantile Phi^{-1}(p) using Acklam's approximation + Halley refinement. */
export function normalQuantile(p: number): number {
  if (!(p > 0 && p < 1)) {
    if (p === 0) return -Infinity;
    if (p === 1) return Infinity;
    throw new Error(`normalQuantile probability out of (0,1): ${p}`);
  }
  const a1 = -39.69683028665376;
  const a2 = 220.9460984245205;
  const a3 = -275.9285104469687;
  const a4 = 138.357751867269;
  const a5 = -30.66479806614716;
  const a6 = 2.506628277459239;

  const b1 = -54.47609879822406;
  const b2 = 161.5858368580409;
  const b3 = -155.6989798598866;
  const b4 = 66.80131188771972;
  const b5 = -13.28068155288572;

  const c1 = -0.007784894002430293;
  const c2 = -0.3223964580411365;
  const c3 = -2.400758277161838;
  const c4 = -2.549732539343734;
  const c5 = 4.374664141464968;
  const c6 = 2.938163982698783;

  const d1 = 0.007784695709041462;
  const d2 = 0.3224671290700398;
  const d3 = 2.445134137142996;
  const d4 = 3.754408661907416;

  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let x = 0;
  if (p < pLow) {
    const q = Math.sqrt(-2 * Math.log(p));
    x =
      (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
      ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  } else if (p <= pHigh) {
    const q = p - 0.5;
    const r = q * q;
    x =
      ((((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q) /
      (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
  } else {
    const q = Math.sqrt(-2 * Math.log(1 - p));
    x =
      -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) /
      ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
  }

  // Two Halley refinements to machine precision
  for (let i = 0; i < 2; i += 1) {
    const err = p <= 0.5 ? normalCdf(x) - p : p - normalCdf(x) === 0 ? 0 : -(normalSf(x) - (1 - p));
    const u = err / normalPdf(x);
    x -= u / (1 + 0.5 * x * u);
  }
  return x;
}

/** Upper regularized incomplete gamma function Q(a, x) = Gamma(a, x) / Gamma(a). */
export function regularizedGammaQ(a: number, x: number): number {
  if (!(a > 0) || x < 0) {
    throw new Error(`regularizedGammaQ invalid args a=${a}, x=${x}`);
  }
  if (x === 0) return 1;
  if (x < a + 1) {
    // Series representation for P(a, x), return 1 - P
    let ap = a;
    let del = 1 / a;
    let sum = del;
    for (let n = 1; n <= 400; n += 1) {
      ap += 1;
      del *= x / ap;
      sum += del;
      if (Math.abs(del) <= Math.abs(sum) * 3e-16) {
        break;
      }
    }
    const pVal = sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
    return Math.max(0, 1 - pVal);
  }
  // Modified Lentz continued fraction for Q(a, x)
  const FPMIN = 1e-300;
  let b = x + 1 - a;
  let c = 1 / FPMIN;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i <= 400; i += 1) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) <= 3e-16) {
      break;
    }
  }
  return Math.exp(-x + a * Math.log(x) - logGamma(a)) * h;
}

/** Chi-square upper-tail survival p-value P(X >= x) for df > 0. */
export function chi2Sf(x: number, df: number): number {
  if (x <= 0) return 1;
  if (df === 1) {
    return 2 * normalSf(Math.sqrt(x));
  }
  if (df === 2) {
    return Math.exp(-0.5 * x);
  }
  return regularizedGammaQ(0.5 * df, 0.5 * x);
}

/** Recursive adaptive Simpson quadrature on [a, b]. */
export function adaptiveSimpson(
  f: (x: number) => number,
  a: number,
  b: number,
  tol = 1e-12,
  maxDepth = 20,
): number {
  const c = 0.5 * (a + b);
  const fa = f(a);
  const fb = f(b);
  const fc = f(c);
  const whole = ((b - a) / 6) * (fa + 4 * fc + fb);

  function recurse(
    lo: number,
    hi: number,
    fLo: number,
    fMid: number,
    fHi: number,
    sWhole: number,
    eps: number,
    depth: number,
  ): number {
    const mid = 0.5 * (lo + hi);
    const lMid = 0.5 * (lo + mid);
    const rMid = 0.5 * (mid + hi);
    const fL = f(lMid);
    const fR = f(rMid);
    const h = 0.5 * (hi - lo);
    const sLeft = (h / 6) * (fLo + 4 * fL + fMid);
    const sRight = (h / 6) * (fMid + 4 * fR + fHi);
    const delta = sLeft + sRight - sWhole;
    if (depth <= 0 || Math.abs(delta) <= 15 * eps) {
      return sLeft + sRight + delta / 15;
    }
    return (
      recurse(lo, mid, fLo, fL, fMid, sLeft, 0.5 * eps, depth - 1) +
      recurse(mid, hi, fMid, fR, fHi, sRight, 0.5 * eps, depth - 1)
    );
  }

  return recurse(a, b, fa, fc, fb, whole, tol, maxDepth);
}
