"""Kalibrasi model Heston ke permukaan implied volatility dan analisis identifiabilitas."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .cf import HestonParams
from .pricing import heston_cos, implied_vol

LOWER_BOUNDS = np.array([1e-4, 0.1, 1e-4, 0.05, -0.99], dtype=float)
UPPER_BOUNDS = np.array([1.0, 10.0, 1.0, 2.0, 0.99], dtype=float)
DEFAULT_STARTS = (
    np.array([0.09, 1.0, 0.09, 0.3, -0.3], dtype=float),
    np.array([0.02, 3.0, 0.06, 0.7, -0.5], dtype=float),
    np.array([0.06, 1.5, 0.03, 0.4, -0.8], dtype=float),
)


def calibrate_heston(
    market_iv: np.ndarray,
    s0: float,
    strikes: np.ndarray,
    r: float,
    maturities: np.ndarray,
    n_cos: int = 256,
    x0: np.ndarray | None = None,
    multi_start: bool = True,
) -> dict[str, object]:
    """Kalibrasi parameter Heston ``(v0, kappa, theta, xi, rho)`` via ``least_squares`` ('trf')."""
    target_iv = market_iv.ravel()

    def residuals(x: np.ndarray) -> np.ndarray:
        p = HestonParams(v0=float(x[0]), kappa=float(x[1]), theta=float(x[2]), xi=float(x[3]), rho=float(x[4]))
        prices = heston_cos(s0, strikes, r, maturities, p, n_cos=n_cos, l_trunc=10.0)
        iv_mod = implied_vol(prices, s0, strikes, r, maturities)
        return np.asarray(iv_mod, dtype=float).ravel() - target_iv

    starts = [np.asarray(x0, dtype=float)] if x0 is not None else [DEFAULT_STARTS[0]]
    if multi_start and x0 is None:
        first_res = least_squares(residuals, starts[0], bounds=(LOWER_BOUNDS, UPPER_BOUNDS), method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12)
        best = first_res
        rmse_first = float(np.sqrt(np.mean(first_res.fun**2)))
        runs = [(starts[0], rmse_first)]
        if rmse_first > 1e-6:
            for alt in DEFAULT_STARTS[1:]:
                res_alt = least_squares(residuals, alt, bounds=(LOWER_BOUNDS, UPPER_BOUNDS), method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12)
                rmse_alt = float(np.sqrt(np.mean(res_alt.fun**2)))
                runs.append((alt, rmse_alt))
                if res_alt.cost < best.cost:
                    best = res_alt
    else:
        best = least_squares(residuals, starts[0], bounds=(LOWER_BOUNDS, UPPER_BOUNDS), method="trf", ftol=1e-12, xtol=1e-12, gtol=1e-12)
        runs = [(starts[0], float(np.sqrt(np.mean(best.fun**2))))]

    opt_p = HestonParams(*[float(v) for v in best.x])
    rmse = float(np.sqrt(np.mean(best.fun**2)))
    return {"params": opt_p, "x": best.x, "rmse": rmse, "jac": best.jac, "residuals": best.fun, "runs": runs}


def identifiability_from_jacobian(jac: np.ndarray, noise_std: float = 0.002) -> dict[str, np.ndarray | float]:
    """Hitung kovarians asimtotik ``sigma^2 (J^T J)^-1``, SE asimtotik, dan bilangan kondisi."""
    jtj = jac.T @ jac
    cond_num = float(np.linalg.cond(jtj))
    cov = (noise_std**2) * np.linalg.inv(jtj)
    se = np.sqrt(np.maximum(np.diag(cov), 0.0))
    return {"cov": cov, "se": se, "cond": cond_num}


def run_noisy_calibration_mc(
    true_iv: np.ndarray,
    s0: float,
    strikes: np.ndarray,
    r: float,
    maturities: np.ndarray,
    n_reps: int = 100,
    noise_std: float = 0.002,
    n_cos: int = 192,
    seed: int = 42,
) -> np.ndarray:
    """Jalankan ``n_reps`` replikasi kalibrasi dengan noise ``N(0, noise_std^2)`` pada implied vol."""
    child_seeds = np.random.SeedSequence(seed).spawn(n_reps)
    est_all = np.empty((n_reps, 5), dtype=float)
    x0 = DEFAULT_STARTS[2]
    for idx, c_seed in enumerate(child_seeds):
        rng = np.random.default_rng(c_seed)
        noisy_iv = np.maximum(true_iv + rng.normal(0.0, noise_std, size=true_iv.shape), 1e-4)
        res = calibrate_heston(noisy_iv, s0, strikes, r, maturities, n_cos=n_cos, x0=x0, multi_start=False)
        est_all[idx] = res["x"]
    return est_all
