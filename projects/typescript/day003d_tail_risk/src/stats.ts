/**
 * Summary statistics, empirical quantiles (Type 7), sample ACF, excess kurtosis,
 * batch-means standard errors, and correlation utilities.
 */

export interface MomentsSummary {
  mean: number;
  sd: number;
  se: number;
  ci_low: number;
  ci_high: number;
  median: number;
  q025: number;
  q975: number;
}

export function mean(xs: ArrayLike<number>): number {
  const n = xs.length;
  if (n === 0) return NaN;
  let sum = 0.0;
  for (let i = 0; i < n; i += 1) {
    sum += xs[i] ?? 0.0;
  }
  return sum / n;
}

/** Sample variance (ddof = 1 by default). */
export function variance(xs: ArrayLike<number>, ddof = 1): number {
  const n = xs.length;
  if (n <= ddof) return NaN;
  const m = mean(xs);
  let ss = 0.0;
  for (let i = 0; i < n; i += 1) {
    const d = (xs[i] ?? 0.0) - m;
    ss += d * d;
  }
  return ss / (n - ddof);
}

export function sd(xs: ArrayLike<number>, ddof = 1): number {
  return Math.sqrt(variance(xs, ddof));
}

export function se(xs: ArrayLike<number>): number {
  return sd(xs, 1) / Math.sqrt(xs.length);
}

/** Empirical quantile using Type-7 linear interpolation on sorted copy (unless alreadySorted=true). */
export function quantile(xs: ArrayLike<number>, p: number, alreadySorted = false): number {
  const n = xs.length;
  if (n === 0) return NaN;
  if (p <= 0.0) {
    let mn = Infinity;
    for (let i = 0; i < n; i += 1) {
      const v = xs[i] ?? Infinity;
      if (v < mn) mn = v;
    }
    return mn;
  }
  if (p >= 1.0) {
    let mx = -Infinity;
    for (let i = 0; i < n; i += 1) {
      const v = xs[i] ?? -Infinity;
      if (v > mx) mx = v;
    }
    return mx;
  }
  const arr = alreadySorted
    ? xs
    : Array.from({ length: n }, (_, i) => xs[i] ?? 0.0).sort((a, b) => a - b);
  const idx = p * (n - 1);
  const lo = Math.floor(idx);
  const hi = Math.min(lo + 1, n - 1);
  const frac = idx - lo;
  return (1.0 - frac) * (arr[lo] ?? 0.0) + frac * (arr[hi] ?? 0.0);
}

export function median(xs: ArrayLike<number>, alreadySorted = false): number {
  return quantile(xs, 0.5, alreadySorted);
}

export function summarizeSeries(xs: number[]): MomentsSummary {
  const sorted = [...xs].sort((a, b) => a - b);
  const m = mean(xs);
  const s = sd(xs, 1);
  const stdErr = s / Math.sqrt(xs.length);
  return {
    mean: m,
    sd: s,
    se: stdErr,
    ci_low: m - 1.959963984540054 * stdErr,
    ci_high: m + 1.959963984540054 * stdErr,
    median: median(sorted, true),
    q025: quantile(sorted, 0.025, true),
    q975: quantile(sorted, 0.975, true),
  };
}

/** Sample autocorrelation at integer lag >= 1. */
export function sampleAcf(xs: ArrayLike<number>, lag: number): number {
  const n = xs.length;
  if (lag < 1 || lag >= n) {
    throw new Error(`Invalid ACF lag ${lag} for length ${n}`);
  }
  const m = mean(xs);
  let num = 0.0;
  let den = 0.0;
  for (let i = 0; i < n; i += 1) {
    const d = (xs[i] ?? 0.0) - m;
    den += d * d;
    if (i >= lag) {
      num += ((xs[i - lag] ?? 0.0) - m) * d;
    }
  }
  return den > 0.0 ? num / den : NaN;
}

/** Excess kurtosis E[(X - mu)^4] / sigma^4 - 3 (population plug-in). */
export function sampleExcessKurtosis(xs: ArrayLike<number>): number {
  const n = xs.length;
  const m = mean(xs);
  let m2 = 0.0;
  let m4 = 0.0;
  for (let i = 0; i < n; i += 1) {
    const d = (xs[i] ?? 0.0) - m;
    const d2 = d * d;
    m2 += d2;
    m4 += d2 * d2;
  }
  m2 /= n;
  m4 /= n;
  return m4 / (m2 * m2) - 3.0;
}

/** Batch-means standard error for dependent stationary series. */
export function batchMeansSe(xs: ArrayLike<number>, nBatches: number): number {
  const n = xs.length;
  const bSize = Math.floor(n / nBatches);
  if (bSize < 2 || nBatches < 2) {
    return NaN;
  }
  const bMeans = new Array<number>(nBatches);
  for (let b = 0; b < nBatches; b += 1) {
    let s = 0.0;
    const start = b * bSize;
    for (let i = 0; i < bSize; i += 1) {
      s += xs[start + i] ?? 0.0;
    }
    bMeans[b] = s / bSize;
  }
  return se(bMeans);
}

/** Sample Pearson correlation between two equal-length arrays. */
export function correlation(xs: ArrayLike<number>, ys: ArrayLike<number>): number {
  const n = xs.length;
  if (n !== ys.length || n < 2) return NaN;
  const mx = mean(xs);
  const my = mean(ys);
  let sxy = 0.0;
  let sxx = 0.0;
  let syy = 0.0;
  for (let i = 0; i < n; i += 1) {
    const dx = (xs[i] ?? 0.0) - mx;
    const dy = (ys[i] ?? 0.0) - my;
    sxy += dx * dy;
    sxx += dx * dx;
    syy += dy * dy;
  }
  const den = Math.sqrt(sxx * syy);
  return den > 0.0 ? sxy / den : NaN;
}

export interface HistogramData {
  edges: number[];
  counts: number[];
  density: number[];
}

export function computeHistogram(
  xs: ArrayLike<number>,
  nBins: number,
  xMin?: number,
  xMax?: number,
): HistogramData {
  const n = xs.length;
  let lo = xMin ?? Infinity;
  let hi = xMax ?? -Infinity;
  if (xMin === undefined || xMax === undefined) {
    for (let i = 0; i < n; i += 1) {
      const v = xs[i] ?? 0.0;
      if (xMin === undefined && v < lo) lo = v;
      if (xMax === undefined && v > hi) hi = v;
    }
  }
  if (!(hi > lo)) {
    lo -= 0.5;
    hi += 0.5;
  }
  const width = (hi - lo) / nBins;
  const edges = Array.from({ length: nBins + 1 }, (_, i) => lo + i * width);
  const counts = new Array<number>(nBins).fill(0);
  for (let i = 0; i < n; i += 1) {
    const v = xs[i] ?? 0.0;
    let bin = Math.floor((v - lo) / width);
    if (bin < 0) bin = 0;
    if (bin >= nBins) bin = nBins - 1;
    counts[bin] = (counts[bin] ?? 0) + 1;
  }
  const density = counts.map((c) => c / (Math.max(n, 1) * width));
  return { edges, counts, density };
}
