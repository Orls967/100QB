/**
 * Unconstrained multivariate optimization & linear algebra utilities from scratch:
 * - BFGS with Nocedal & Wright Strong Wolfe line search
 * - Nelder-Mead simplex optimizer (for independent derivative-free cross-checks)
 * - Dense matrix inversion (Gauss-Jordan with partial pivoting) and central finite differences
 */

export interface ValueGrad {
  value: number;
  grad: number[];
}

export interface BfgsOptions {
  maxIter?: number;
  gradTol?: number;
  fTol?: number;
  xTol?: number;
}

export interface BfgsResult {
  x: number[];
  value: number;
  grad: number[];
  gradNorm: number;
  iters: number;
  evals: number;
  converged: boolean;
}

function dot(a: number[], b: number[]): number {
  let s = 0.0;
  for (let i = 0; i < a.length; i += 1) {
    s += (a[i] ?? 0.0) * (b[i] ?? 0.0);
  }
  return s;
}

function infNorm(a: number[]): number {
  let m = 0.0;
  for (let i = 0; i < a.length; i += 1) {
    const v = Math.abs(a[i] ?? 0.0);
    if (v > m) m = v;
  }
  return m;
}

function addScaled(x: number[], d: number[], alpha: number): number[] {
  const n = x.length;
  const out = new Array<number>(n);
  for (let i = 0; i < n; i += 1) {
    out[i] = (x[i] ?? 0.0) + alpha * (d[i] ?? 0.0);
  }
  return out;
}

/**
 * Strong Wolfe line search (Nocedal & Wright Algorithms 3.5 & 3.6) with finite-value guards.
 */
function strongWolfeLineSearch(
  fnGrad: (x: number[]) => ValueGrad,
  x: number[],
  f0: number,
  g0: number[],
  dir: number[],
  alphaInit: number,
  c1 = 1e-4,
  c2 = 0.9,
): { alpha: number; xNew: number[]; fg: ValueGrad; evals: number } {
  const slope0 = dot(g0, dir);
  let evals = 0;

  const evalAt = (alpha: number): { xCand: number[]; fg: ValueGrad; slope: number } => {
    const xCand = addScaled(x, dir, alpha);
    const fg = fnGrad(xCand);
    evals += 1;
    return { xCand, fg, slope: dot(fg.grad, dir) };
  };

  let alphaLo = 0.0;
  let fLo = f0;
  let slopeLo = slope0;

  const zoom = (
    aLoInit: number,
    fLoInit: number,
    dLoInit: number,
    aHiInit: number,
    fHiInit: number,
  ): { alpha: number; xNew: number[]; fg: ValueGrad } => {
    let aLo = aLoInit;
    let fL = fLoInit;
    let dL = dLoInit;
    let aHi = aHiInit;
    let fH = fHiInit;

    for (let j = 0; j < 25; j += 1) {
      // Quadratic interpolation inside [aLo, aHi] with safeguard toward midpoint
      let aJ = 0.5 * (aLo + aHi);
      const denom = 2.0 * (fH - fL - dL * (aHi - aLo));
      if (Number.isFinite(fH) && Number.isFinite(denom) && Math.abs(denom) > 1e-18) {
        const quad = aLo - (dL * (aHi - aLo) * (aHi - aLo)) / denom;
        const loB = Math.min(aLo, aHi);
        const hiB = Math.max(aLo, aHi);
        const span = hiB - loB;
        if (quad > loB + 0.1 * span && quad < hiB - 0.1 * span) {
          aJ = quad;
        }
      }
      const ev = evalAt(aJ);
      if (
        !Number.isFinite(ev.fg.value) ||
        ev.fg.value > f0 + c1 * aJ * slope0 ||
        ev.fg.value >= fL
      ) {
        aHi = aJ;
        fH = ev.fg.value;
      } else {
        if (Math.abs(ev.slope) <= -c2 * slope0) {
          return { alpha: aJ, xNew: ev.xCand, fg: ev.fg };
        }
        if (ev.slope * (aHi - aLo) >= 0.0) {
          aHi = aLo;
          fH = fL;
        }
        aLo = aJ;
        fL = ev.fg.value;
        dL = ev.slope;
      }
      if (Math.abs(aHi - aLo) <= 1e-14 * Math.max(1.0, Math.abs(aLo))) {
        return { alpha: aJ, xNew: ev.xCand, fg: ev.fg };
      }
    }
    const mid = 0.5 * (aLo + aHi);
    const evMid = evalAt(mid);
    if (Number.isFinite(evMid.fg.value) && evMid.fg.value < f0) {
      return { alpha: mid, xNew: evMid.xCand, fg: evMid.fg };
    }
    return { alpha: 0.0, xNew: [...x], fg: { value: f0, grad: [...g0] } };
  };

  let alphaPrev = 0.0;
  let fPrev = f0;
  let alphaCur = alphaInit;
  const alphaMax = 64.0;

  for (let i = 1; i <= 25; i += 1) {
    const ev = evalAt(alphaCur);
    if (!Number.isFinite(ev.fg.value)) {
      alphaCur *= 0.5;
      continue;
    }
    if (ev.fg.value > f0 + c1 * alphaCur * slope0 || (i > 1 && ev.fg.value >= fPrev)) {
      const z = zoom(alphaPrev, fPrev, slopeLo, alphaCur, ev.fg.value);
      return { ...z, evals };
    }
    if (Math.abs(ev.slope) <= -c2 * slope0) {
      return { alpha: alphaCur, xNew: ev.xCand, fg: ev.fg, evals };
    }
    if (ev.slope >= 0.0) {
      const z = zoom(alphaCur, ev.fg.value, ev.slope, alphaPrev, fPrev);
      return { ...z, evals };
    }
    alphaLo = alphaCur;
    fLo = ev.fg.value;
    slopeLo = ev.slope;
    alphaPrev = alphaLo;
    fPrev = fLo;
    alphaCur = Math.min(2.0 * alphaCur, alphaMax);
    if (alphaCur >= alphaMax) {
      return { alpha: alphaLo, xNew: ev.xCand, fg: ev.fg, evals };
    }
  }

  // Fallback Armijo backtracking
  let a = 0.25;
  for (let k = 0; k < 40; k += 1) {
    const ev = evalAt(a);
    if (Number.isFinite(ev.fg.value) && ev.fg.value <= f0 + c1 * a * slope0) {
      return { alpha: a, xNew: ev.xCand, fg: ev.fg, evals };
    }
    a *= 0.5;
  }
  return { alpha: 0.0, xNew: [...x], fg: { value: f0, grad: [...g0] }, evals };
}

/**
 * Minimize smooth scalar function using BFGS with inverse Hessian updates and Strong Wolfe line search.
 */
export function bfgs(
  fnGrad: (x: number[]) => ValueGrad,
  x0: number[],
  options: BfgsOptions = {},
): BfgsResult {
  const maxIter = options.maxIter ?? 250;
  const gradTol = options.gradTol ?? 1e-8;
  const fTol = options.fTol ?? 1e-13;
  const xTol = options.xTol ?? 1e-12;

  const n = x0.length;
  let x = [...x0];
  let fg = fnGrad(x);
  let evals = 1;

  // Identity inverse Hessian
  const H: number[][] = Array.from({ length: n }, (_, i) =>
    Array.from({ length: n }, (_, j) => (i === j ? 1.0 : 0.0)),
  );

  let iters = 0;
  let converged = false;

  for (let iter = 0; iter < maxIter; iter += 1) {
    const gNorm = infNorm(fg.grad);
    if (gNorm <= gradTol) {
      converged = true;
      break;
    }

    // Search direction d = -H g
    const dir = new Array<number>(n).fill(0.0);
    for (let i = 0; i < n; i += 1) {
      const row = H[i];
      if (!row) continue;
      let sum = 0.0;
      for (let j = 0; j < n; j += 1) {
        sum += (row[j] ?? 0.0) * (fg.grad[j] ?? 0.0);
      }
      dir[i] = -sum;
    }

    // Reset H if not a descent direction
    if (dot(fg.grad, dir) >= 0.0) {
      for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
          const r = H[i];
          if (r) r[j] = i === j ? 1.0 : 0.0;
        }
        dir[i] = -(fg.grad[i] ?? 0.0);
      }
    }

    const ls = strongWolfeLineSearch(fnGrad, x, fg.value, fg.grad, dir, 1.0);
    evals += ls.evals;

    if (ls.alpha === 0.0) {
      // Try steepest descent restart once
      const sdDir = fg.grad.map((v) => -v);
      const sdLs = strongWolfeLineSearch(fnGrad, x, fg.value, fg.grad, sdDir, 0.1);
      evals += sdLs.evals;
      if (sdLs.alpha === 0.0) {
        break;
      }
      ls.alpha = sdLs.alpha;
      ls.xNew = sdLs.xNew;
      ls.fg = sdLs.fg;
    }

    iters = iter + 1;
    const s = ls.xNew.map((v, i) => v - (x[i] ?? 0.0));
    const y = ls.fg.grad.map((v, i) => v - (fg.grad[i] ?? 0.0));
    const fDrop = Math.abs(fg.value - ls.fg.value);
    const stepNorm = infNorm(s);

    x = ls.xNew;
    fg = ls.fg;

    if (
      infNorm(fg.grad) <= gradTol ||
      (stepNorm <= xTol && fDrop <= fTol * (1.0 + Math.abs(fg.value)))
    ) {
      converged = infNorm(fg.grad) <= Math.max(gradTol * 100.0, 1e-5);
      if (converged) break;
    }

    const sy = dot(s, y);
    const sn = Math.sqrt(dot(s, s));
    const yn = Math.sqrt(dot(y, y));
    if (sy > 1e-12 * sn * yn) {
      const rho = 1.0 / sy;
      // Scale H on first iteration per Nocedal & Wright eq. (6.20)
      if (iter === 0) {
        const scale = sy / Math.max(dot(y, y), 1e-30);
        for (let i = 0; i < n; i += 1) {
          const r = H[i];
          if (r) {
            for (let j = 0; j < n; j += 1) {
              r[j] = i === j ? scale : 0.0;
            }
          }
        }
      }

      const Hy = new Array<number>(n).fill(0.0);
      for (let i = 0; i < n; i += 1) {
        const r = H[i];
        if (!r) continue;
        let sum = 0.0;
        for (let j = 0; j < n; j += 1) {
          sum += (r[j] ?? 0.0) * (y[j] ?? 0.0);
        }
        Hy[i] = sum;
      }
      const yHy = dot(y, Hy);
      const factor = (1.0 + rho * yHy) * rho;
      for (let i = 0; i < n; i += 1) {
        const r = H[i];
        if (!r) continue;
        const si = s[i] ?? 0.0;
        const hyi = Hy[i] ?? 0.0;
        for (let j = 0; j < n; j += 1) {
          const sj = s[j] ?? 0.0;
          const hyj = Hy[j] ?? 0.0;
          r[j] = (r[j] ?? 0.0) + factor * si * sj - rho * (hyi * sj + si * hyj);
        }
      }
    }
  }

  return {
    x,
    value: fg.value,
    grad: fg.grad,
    gradNorm: infNorm(fg.grad),
    iters,
    evals,
    converged,
  };
}

/**
 * Derivative-free Nelder-Mead simplex minimizer for independent optimizer cross-checking.
 */
export function nelderMead(
  fn: (x: number[]) => number,
  x0: number[],
  maxIter = 1500,
  tol = 1e-9,
): { x: number[]; value: number; iters: number; converged: boolean } {
  const n = x0.length;
  const simplex: Array<{ x: number[]; fx: number }> = [];
  simplex.push({ x: [...x0], fx: fn(x0) });
  for (let i = 0; i < n; i += 1) {
    const pt = [...x0];
    const xi = pt[i] ?? 0.0;
    const step = Math.abs(xi) > 1e-8 ? 0.05 * xi : 0.025;
    pt[i] = xi + step;
    simplex.push({ x: pt, fx: fn(pt) });
  }

  let iters = 0;
  let converged = false;

  for (let iter = 0; iter < maxIter; iter += 1) {
    simplex.sort((a, b) => a.fx - b.fx);
    const best = simplex[0];
    const worst = simplex[n];
    const secondWorst = simplex[n - 1];
    if (!best || !worst || !secondWorst) break;

    let maxSpan = 0.0;
    for (let i = 1; i <= n; i += 1) {
      const si = simplex[i];
      if (!si) continue;
      for (let j = 0; j < n; j += 1) {
        maxSpan = Math.max(maxSpan, Math.abs((si.x[j] ?? 0.0) - (best.x[j] ?? 0.0)));
      }
    }
    if (Math.abs(worst.fx - best.fx) <= tol && maxSpan <= tol) {
      converged = true;
      break;
    }

    const centroid = new Array<number>(n).fill(0.0);
    for (let i = 0; i < n; i += 1) {
      const si = simplex[i];
      if (!si) continue;
      for (let j = 0; j < n; j += 1) {
        centroid[j] = (centroid[j] ?? 0.0) + (si.x[j] ?? 0.0) / n;
      }
    }

    const reflect = centroid.map((c, j) => c + (c - (worst.x[j] ?? 0.0)));
    const fReflect = fn(reflect);

    if (fReflect >= best.fx && fReflect < secondWorst.fx) {
      simplex[n] = { x: reflect, fx: fReflect };
    } else if (fReflect < best.fx) {
      const expand = centroid.map((c, j) => c + 2.0 * ((reflect[j] ?? 0.0) - c));
      const fExpand = fn(expand);
      simplex[n] = fExpand < fReflect ? { x: expand, fx: fExpand } : { x: reflect, fx: fReflect };
    } else {
      const contract =
        fReflect < worst.fx
          ? centroid.map((c, j) => c + 0.5 * ((reflect[j] ?? 0.0) - c))
          : centroid.map((c, j) => c + 0.5 * ((worst.x[j] ?? 0.0) - c));
      const fContract = fn(contract);
      if (fContract < Math.min(fReflect, worst.fx)) {
        simplex[n] = { x: contract, fx: fContract };
      } else {
        for (let i = 1; i <= n; i += 1) {
          const si = simplex[i];
          if (!si) continue;
          const shr = best.x.map((b, j) => b + 0.5 * ((si.x[j] ?? 0.0) - b));
          simplex[i] = { x: shr, fx: fn(shr) };
        }
      }
    }
    iters = iter + 1;
  }

  simplex.sort((a, b) => a.fx - b.fx);
  const top = simplex[0] ?? { x: x0, fx: fn(x0) };
  return { x: top.x, value: top.fx, iters, converged };
}

/** Gauss-Jordan matrix inversion with partial pivoting for small dense matrices. */
export function invertMatrix(mat: number[][]): number[][] {
  const n = mat.length;
  const aug: number[][] = Array.from({ length: n }, (_, i) => {
    const row = new Array<number>(2 * n).fill(0.0);
    const src = mat[i] ?? [];
    for (let j = 0; j < n; j += 1) {
      row[j] = src[j] ?? 0.0;
    }
    row[n + i] = 1.0;
    return row;
  });

  for (let col = 0; col < n; col += 1) {
    let pivotRow = col;
    let pivotAbs = Math.abs(aug[col]?.[col] ?? 0.0);
    for (let r = col + 1; r < n; r += 1) {
      const v = Math.abs(aug[r]?.[col] ?? 0.0);
      if (v > pivotAbs) {
        pivotAbs = v;
        pivotRow = r;
      }
    }
    if (pivotAbs < 1e-18) {
      throw new Error(`Singular matrix during inversion at column ${col}`);
    }
    if (pivotRow !== col) {
      const tmp = aug[col];
      const sw = aug[pivotRow];
      if (tmp && sw) {
        aug[col] = sw;
        aug[pivotRow] = tmp;
      }
    }
    const pRow = aug[col];
    if (!pRow) continue;
    const invPivot = 1.0 / (pRow[col] ?? 1.0);
    for (let j = 0; j < 2 * n; j += 1) {
      pRow[j] = (pRow[j] ?? 0.0) * invPivot;
    }
    for (let r = 0; r < n; r += 1) {
      if (r === col) continue;
      const rRow = aug[r];
      if (!rRow) continue;
      const factor = rRow[col] ?? 0.0;
      if (factor !== 0.0) {
        for (let j = 0; j < 2 * n; j += 1) {
          rRow[j] = (rRow[j] ?? 0.0) - factor * (pRow[j] ?? 0.0);
        }
      }
    }
  }

  return aug.map((row) => row.slice(n));
}

/** Central finite-difference gradient with per-coordinate step sizes. */
export function centralDiffGradient(
  fn: (x: number[]) => number,
  x: number[],
  steps: number[],
): number[] {
  const n = x.length;
  const g = new Array<number>(n).fill(0.0);
  for (let i = 0; i < n; i += 1) {
    const h = steps[i] ?? 1e-6;
    const xp = [...x];
    const xm = [...x];
    xp[i] = (xp[i] ?? 0.0) + h;
    xm[i] = (xm[i] ?? 0.0) - h;
    g[i] = (fn(xp) - fn(xm)) / (2.0 * h);
  }
  return g;
}
