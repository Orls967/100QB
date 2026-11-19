"""Risiko kredit Vasicek: distribusi kerugian eksak (DP + Gauss-Hermite), LHP & importance sampling."""

from __future__ import annotations

import time

import numpy as np
from scipy import integrate
from scipy.interpolate import CubicSpline
from scipy.stats import binom, multivariate_normal, norm

DEFAULT_SEED = 42
DEFAULT_PD = 0.02
DEFAULT_RHO = 0.20
DEFAULT_N = 100


def make_exposures(case: str, n_obligors: int = DEFAULT_N) -> np.ndarray:
    """Bangkitkan vektor eksposur integer: ``'A'`` (homogen 1) atau ``'B'`` (heterogen 1..5)."""
    if case.upper() == "A":
        return np.ones(n_obligors, dtype=int)
    if case.upper() == "B":
        return 1 + (np.arange(n_obligors, dtype=int) % 5)
    raise ValueError(f"Kasus tidak dikenal: {case}")


def conditional_pd(z: np.ndarray | float, pd: float = DEFAULT_PD, rho: float = DEFAULT_RHO) -> np.ndarray | float:
    """Peluang gagal bayar bersyarat faktor sistemik ``z``."""
    thresh = norm.ppf(pd)
    return norm.cdf((thresh - np.sqrt(rho) * z) / np.sqrt(1.0 - rho))


def conditional_loss_pmf_dp(exposures: np.ndarray, p_z: np.ndarray) -> np.ndarray:
    """Hitung PMF kerugian bersyarat via dynamic programming untuk vektor ``p_z`` berukuran ``(M,)``."""
    p_arr = np.atleast_1d(np.asarray(p_z, dtype=float))
    total = int(np.sum(exposures))

    if np.all(exposures == 1):
        n_ob = len(exposures)
        k_grid = np.arange(total + 1)
        return binom.pmf(k_grid[None, :], n_ob, p_arr[:, None])

    pmf = np.zeros((p_arr.size, total + 1), dtype=float)
    pmf[:, 0] = 1.0
    for exp_i in exposures:
        e = int(exp_i)
        p_col = p_arr[:, None]
        pmf[:, e:] = (1.0 - p_col) * pmf[:, e:] + p_col * pmf[:, :-e]
        pmf[:, :e] = (1.0 - p_col) * pmf[:, :e]
    return pmf


def exact_loss_pmf(
    exposures: np.ndarray,
    pd: float = DEFAULT_PD,
    rho: float = DEFAULT_RHO,
    n_gh: int = 200,
) -> np.ndarray:
    """PMF kerugian eksak melalui kuadratur Gauss-Hermite probabilis ``n_gh`` node."""
    z_nodes, w_nodes = np.polynomial.hermite_e.hermegauss(n_gh)
    w_norm = w_nodes / np.sum(w_nodes)
    p_z = conditional_pd(z_nodes, pd=pd, rho=rho)
    pmf_cond = conditional_loss_pmf_dp(exposures, p_z)
    return w_norm @ pmf_cond


def exact_risk_metrics(pmf: np.ndarray, alphas: tuple[float, ...] = (0.99, 0.999)) -> dict[str, object]:
    """Hitung mean, varians, VaR diskret, dan Expected Shortfall dengan koreksi massa diskret."""
    losses = np.arange(len(pmf), dtype=float)
    mean = float(np.sum(losses * pmf))
    var = float(np.sum(((losses - mean) ** 2) * pmf))
    cdf = np.cumsum(pmf)

    vars_dict: dict[float, int] = {}
    es_dict: dict[float, float] = {}
    for alpha in alphas:
        idx = int(np.searchsorted(cdf, alpha, side="left"))
        idx = min(idx, len(pmf) - 1)
        vars_dict[alpha] = idx
        tail_mass = float(np.sum(losses[idx + 1 :] * pmf[idx + 1 :]))
        es_dict[alpha] = ((cdf[idx] - alpha) * idx + tail_mass) / (1.0 - alpha)

    return {"mean": mean, "var": var, "std": float(np.sqrt(var)), "cdf": cdf, "var_q": vars_dict, "es_q": es_dict}


def bivariate_loss_variance(exposures: np.ndarray, pd: float = DEFAULT_PD, rho: float = DEFAULT_RHO) -> float:
    """Verifikasi varians kerugian eksak memakai ``Phi_2(a, a; rho)``."""
    a = float(norm.ppf(pd))
    cov_mat = np.array([[1.0, rho], [rho, 1.0]], dtype=float)
    p_joint = float(multivariate_normal.cdf([a, a], mean=[0.0, 0.0], cov=cov_mat, abseps=1e-12, maxpts=1_000_000))
    cov_ii = pd * (1.0 - pd)
    cov_ij = p_joint - (pd * pd)
    sum_e2 = float(np.sum(exposures.astype(float) ** 2))
    sum_e_sq = float(np.sum(exposures.astype(float)) ** 2)
    return sum_e2 * cov_ii + (sum_e_sq - sum_e2) * cov_ij


def exact_pmf_quad_homogeneous(
    k_values: list[int] | tuple[int, ...],
    n_obligors: int = DEFAULT_N,
    pd: float = DEFAULT_PD,
    rho: float = DEFAULT_RHO,
) -> np.ndarray:
    """Hitung ``P(L = k)`` kasus homogen via ``scipy.integrate.quad`` sebagai pembanding independen."""
    out = np.empty(len(k_values), dtype=float)
    for idx, k in enumerate(k_values):
        val, _ = integrate.quad(
            lambda z, kk=k: float(binom.pmf(kk, n_obligors, conditional_pd(z, pd=pd, rho=rho)) * norm.pdf(z)),
            -10.0,
            10.0,
            epsabs=1e-13,
            epsrel=1e-13,
            limit=200,
        )
        out[idx] = val
    return out


def lhp_tail_prob(x_frac: np.ndarray | float, pd: float = DEFAULT_PD, rho: float = DEFAULT_RHO) -> np.ndarray | float:
    """Peluang ekor ``P(L/N >= x)`` menurut pendekatan Large Homogeneous Portfolio (Vasicek)."""
    x_arr = np.asarray(x_frac, dtype=float)
    z_x = norm.ppf(x_arr)
    z_p = norm.ppf(pd)
    cdf_val = norm.cdf((np.sqrt(1.0 - rho) * z_x - z_p) / np.sqrt(rho))
    return 1.0 - cdf_val


def lhp_var(alpha: float, pd: float = DEFAULT_PD, rho: float = DEFAULT_RHO) -> float:
    """VaR kuantil ``alpha`` (fraksi) menurut pendekatan LHP."""
    return float(norm.cdf((norm.ppf(pd) + np.sqrt(rho) * norm.ppf(alpha)) / np.sqrt(1.0 - rho)))


def build_conditional_tail_interpolator(
    exposures: np.ndarray,
    thresholds: list[int] | tuple[int, ...],
    pd: float = DEFAULT_PD,
    rho: float = DEFAULT_RHO,
    z_min: float = -8.0,
    z_max: float = 8.0,
    n_grid: int = 2001,
) -> dict[int, CubicSpline]:
    """Prekomputasi ``P(L >= x | z)`` pada grid z halus dengan DP lalu bangun ``CubicSpline``."""
    z_grid = np.linspace(z_min, z_max, n_grid)
    p_z = conditional_pd(z_grid, pd=pd, rho=rho)
    pmf_cond = conditional_loss_pmf_dp(exposures, p_z)
    # Tail prob P(L >= k) = reverse cumsum dari kanan ke kiri
    tail_cond = np.cumsum(pmf_cond[:, ::-1], axis=1)[:, ::-1]
    return {int(thr): CubicSpline(z_grid, np.clip(tail_cond[:, int(thr)], 0.0, 1.0), extrapolate=True) for thr in thresholds}


def _eval_tail_spline(spline: CubicSpline, z: np.ndarray) -> np.ndarray:
    """Evaluasi spline ``P(L >= x | z)`` dengan penjagaan batas ekor ``[-8, 8]``."""
    z_clamped = np.clip(z, -8.0, 8.0)
    return np.clip(spline(z_clamped), 0.0, 1.0)


def _summarize_mc(samples: np.ndarray, exact_val: float, elapsed_s: float) -> dict[str, float]:
    n = len(samples)
    mean = float(np.mean(samples))
    var = float(np.var(samples, ddof=1))
    se = float(np.sqrt(var / n))
    z_score = abs(mean - exact_val) / se if se > 0.0 else (0.0 if mean == exact_val else np.inf)
    return {
        "estimate": mean,
        "se": se,
        "var": var,
        "ci_low": mean - 1.96 * se,
        "ci_high": mean + 1.96 * se,
        "z_score": z_score,
        "time_sec": elapsed_s,
    }


def estimate_naive_mc(
    exposures: np.ndarray,
    threshold: int,
    n_samples: int = 200_000,
    chunk_size: int = 20_000,
    pd: float = DEFAULT_PD,
    rho: float = DEFAULT_RHO,
    seed: int = DEFAULT_SEED,
    exact_val: float = 0.0,
) -> dict[str, float]:
    """Estimator Monte Carlo naif (simulasi ``Z`` dan ``eps_i`` penuh dengan chunking)."""
    t0 = time.perf_counter()
    n_ob = len(exposures)
    a = float(norm.ppf(pd))
    sq_rho = np.sqrt(rho)
    sq_1m = np.sqrt(1.0 - rho)
    exp_f = exposures.astype(float)

    n_chunks = (n_samples + chunk_size - 1) // chunk_size
    child_seeds = np.random.SeedSequence(seed).spawn(n_chunks)
    indicators = np.empty(n_samples, dtype=float)

    offset = 0
    for c_seed in child_seeds:
        cur = min(chunk_size, n_samples - offset)
        rng = np.random.default_rng(c_seed)
        z = rng.standard_normal((cur, 1))
        eps = rng.standard_normal((cur, n_ob))
        assets = sq_rho * z + sq_1m * eps
        losses = (assets < a) @ exp_f
        indicators[offset : offset + cur] = (losses >= threshold).astype(float)
        offset += cur

    return _summarize_mc(indicators, exact_val, time.perf_counter() - t0)


def estimate_conditional_mc(
    spline: CubicSpline,
    n_samples: int = 200_000,
    seed: int = DEFAULT_SEED,
    exact_val: float = 0.0,
) -> dict[str, float]:
    """Estimator Conditional Monte Carlo (hanya simulasi ``Z``, pakai ``P(L >= x | Z)`` eksak)."""
    t0 = time.perf_counter()
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(n_samples)
    vals = _eval_tail_spline(spline, z)
    return _summarize_mc(vals, exact_val, time.perf_counter() - t0)


def find_optimal_is_theta(
    spline: CubicSpline,
    theta_grid: np.ndarray | None = None,
    n_pilot: int = 20_000,
    seed: int = 1042,
) -> tuple[float, np.ndarray, np.ndarray]:
    """Cari pergeseran ``theta`` optimal (``Z ~ N(-theta, 1)``) pada sampel pilot independen."""
    if theta_grid is None:
        theta_grid = np.linspace(0.0, 4.0, 41)
    rng = np.random.default_rng(seed)
    z_base = rng.standard_normal(n_pilot)

    variances = np.empty(len(theta_grid), dtype=float)
    for idx, th in enumerate(theta_grid):
        z_shifted = z_base - th
        lr = np.exp(th * z_shifted + 0.5 * th * th)
        vals = _eval_tail_spline(spline, z_shifted) * lr
        variances[idx] = float(np.var(vals, ddof=1))

    best_idx = int(np.argmin(variances))
    return float(theta_grid[best_idx]), theta_grid, variances


def estimate_importance_sampling_mc(
    spline: CubicSpline,
    theta: float,
    n_samples: int = 200_000,
    seed: int = DEFAULT_SEED,
    exact_val: float = 0.0,
) -> dict[str, float]:
    """Estimator Importance Sampling pada ``Z ~ N(-theta, 1)`` + Conditional MC."""
    t0 = time.perf_counter()
    rng = np.random.default_rng(seed)
    z = rng.normal(loc=-theta, scale=1.0, size=n_samples)
    lr = np.exp(theta * z + 0.5 * theta * theta)
    vals = _eval_tail_spline(spline, z) * lr
    return _summarize_mc(vals, exact_val, time.perf_counter() - t0)


def fit_loglog_slope(n_values: list[int] | np.ndarray, se_values: list[float] | np.ndarray) -> float:
    """Fit slope regresi linear ``log10(SE)`` terhadap ``log10(n)``."""
    xs = np.log10(np.asarray(n_values, dtype=float))
    ys = np.log10(np.maximum(np.asarray(se_values, dtype=float), 1e-16))
    slope, _ = np.polyfit(xs, ys, 1)
    return float(slope)
