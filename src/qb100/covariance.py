"""Estimasi kovarians & portofolio varians minimum: Wishart, Marchenko-Pastur, Ledoit-Wolf."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

DEFAULT_SEED = 42


def generate_factor_model(
    n_assets: int = 50,
    n_factors: int = 3,
    seed: int = DEFAULT_SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Bangkitkan ``(Sigma_true, mu)`` sintetis dari model K-faktor."""
    if n_assets < 2 or n_factors < 1:
        raise ValueError("n_assets >= 2 dan n_factors >= 1")
    rng = np.random.default_rng(seed)
    b_loadings = rng.normal(0.0, 0.15, size=(n_assets, n_factors))
    d_idio = rng.uniform(0.02, 0.05, size=n_assets)
    sigma_true = b_loadings @ b_loadings.T + np.diag(d_idio)
    mu = rng.uniform(0.0002, 0.0010, size=n_assets)
    return sigma_true, mu


def simulate_returns(
    mu: np.ndarray,
    sigma: np.ndarray,
    n_obs: int,
    n_trials: int = 1,
    seed: int | np.random.Generator = DEFAULT_SEED,
) -> np.ndarray:
    """Simulasikan return Gaussian ``r_t ~ N(mu, sigma)`` via Cholesky.

    Mengembalikan array bentuk ``(n_trials, n_obs, n_assets)`` jika ``n_trials > 1``,
    atau ``(n_obs, n_assets)`` jika ``n_trials == 1`` dan dipanggil langsung dengan
    ``squeeze=False`` tidak diperlukan — di sini selalu 3D untuk konsistensi batch.
    """
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    chol = np.linalg.cholesky(sigma)
    n_assets = sigma.shape[0]
    z = rng.standard_normal(size=(n_trials, n_obs, n_assets))
    return mu[None, None, :] + z @ chol.T


def sample_covariance(returns: np.ndarray, ddof: int = 1) -> np.ndarray:
    """Kovarians sampel ternormalisasi ``1 / (T - ddof)`` untuk input 2D atau 3D batched."""
    x = returns - returns.mean(axis=-2, keepdims=True)
    n_obs = returns.shape[-2]
    return np.swapaxes(x, -1, -2) @ x / float(n_obs - ddof)


def ledoit_wolf(returns: np.ndarray) -> tuple[np.ndarray, np.ndarray | float]:
    """Estimator shrinkage Ledoit-Wolf (target identitas) untuk 2D ``(T, N)`` atau 3D ``(M, T, N)``."""
    x = returns - returns.mean(axis=-2, keepdims=True)
    n_obs, n_assets = x.shape[-2], x.shape[-1]
    s_t = np.swapaxes(x, -1, -2) @ x / float(n_obs)

    eye = np.eye(n_assets, dtype=float)
    trace_st = np.trace(s_t, axis1=-2, axis2=-1)
    m_target = trace_st / float(n_assets)

    if returns.ndim == 2:
        diff = s_t - m_target * eye
        d_sq = float(np.sum(diff * diff) / n_assets)
        st_frob_sq = float(np.sum(s_t * s_t))
        x_norm_sq = np.sum(x * x, axis=-1)
        sum_b = float(np.sum(x_norm_sq * x_norm_sq)) - n_obs * st_frob_sq
        b_bar_sq = max(0.0, sum_b / (float(n_obs * n_obs) * n_assets))
        b_sq = min(b_bar_sq, d_sq)
        intensity = (b_sq / d_sq) if d_sq > 0.0 else 0.0
        sigma_lw = intensity * m_target * eye + (1.0 - intensity) * s_t
        return sigma_lw, float(intensity)

    diff = s_t - m_target[:, None, None] * eye[None, :, :]
    d_sq = np.sum(diff * diff, axis=(-2, -1)) / float(n_assets)
    st_frob_sq = np.sum(s_t * s_t, axis=(-2, -1))
    x_norm_sq = np.sum(x * x, axis=-1)
    sum_b = np.sum(x_norm_sq * x_norm_sq, axis=-1) - n_obs * st_frob_sq
    b_bar_sq = np.maximum(0.0, sum_b / (float(n_obs * n_obs) * n_assets))
    b_sq = np.minimum(b_bar_sq, d_sq)
    intensity = np.where(d_sq > 0.0, b_sq / np.where(d_sq > 0.0, d_sq, 1.0), 0.0)
    sigma_lw = (
        intensity[:, None, None] * m_target[:, None, None] * eye[None, :, :]
        + (1.0 - intensity)[:, None, None] * s_t
    )
    return sigma_lw, intensity


def gmv_closed_form(sigma: np.ndarray) -> tuple[np.ndarray, np.ndarray | float]:
    """Bobot & varians GMV tak berkendala dalam bentuk tertutup (mendukung 2D dan 3D)."""
    n_assets = sigma.shape[-1]
    ones = np.ones(n_assets, dtype=float)
    if sigma.ndim == 2:
        inv_ones = np.linalg.solve(sigma, ones)
        denom = float(ones @ inv_ones)
        return inv_ones / denom, 1.0 / denom

    inv_ones = np.linalg.solve(sigma, ones)
    denom = inv_ones.sum(axis=-1)
    weights = inv_ones / denom[:, None]
    return weights, 1.0 / denom


def gmv_numerical(sigma: np.ndarray) -> tuple[np.ndarray, float]:
    """Bobot & varians GMV tak berkendala via ``scipy.optimize.minimize`` (SLSQP)."""
    n_assets = sigma.shape[0]
    w0 = np.full(n_assets, 1.0 / n_assets)
    scale = 1.0 / max(float(np.mean(np.diag(sigma))), 1e-12)
    res = minimize(
        lambda w: scale * float(w @ sigma @ w),
        w0,
        jac=lambda w: 2.0 * scale * (sigma @ w),
        method="SLSQP",
        constraints={"type": "eq", "fun": lambda w: float(np.sum(w) - 1.0), "jac": lambda _w: np.ones(n_assets)},
        options={"ftol": 1e-15, "maxiter": 2000},
    )
    if not res.success:
        raise RuntimeError(f"SLSQP gagal konvergen: {res.message}")
    w_opt = np.asarray(res.x, dtype=float)
    return w_opt, float(w_opt @ sigma @ w_opt)


def long_only_gmv(sigma: np.ndarray) -> tuple[np.ndarray, float]:
    """Bobot & varians GMV berkendala long-only (``w >= 0``, ``sum(w) = 1``)."""
    n_assets = sigma.shape[0]
    w0 = np.full(n_assets, 1.0 / n_assets)
    bounds = [(0.0, 1.0)] * n_assets
    res = minimize(
        lambda w: float(w @ sigma @ w),
        w0,
        jac=lambda w: 2.0 * (sigma @ w),
        method="SLSQP",
        bounds=bounds,
        constraints={"type": "eq", "fun": lambda w: float(np.sum(w) - 1.0), "jac": lambda _w: np.ones(n_assets)},
        options={"ftol": 1e-13, "maxiter": 1000},
    )
    if not res.success:
        raise RuntimeError(f"SLSQP long-only GMV gagal: {res.message}")
    w_opt = np.asarray(res.x, dtype=float)
    return w_opt, float(w_opt @ sigma @ w_opt)


def long_only_frontier(
    sigma_est: np.ndarray,
    mu_est: np.ndarray,
    sigma_true: np.ndarray,
    mu_true: np.ndarray,
    n_points: int = 10,
) -> dict[str, np.ndarray]:
    """Hitung frontier long-only untuk ``n_points`` target return dan evaluasi di ``(Sigma_true, mu_true)``."""
    n_assets = sigma_est.shape[0]
    w_gmv, _ = long_only_gmv(sigma_est)
    ret_min = float(w_gmv @ mu_est)
    ret_max = float(np.max(mu_est)) * 0.999
    targets = np.linspace(ret_min, ret_max, n_points)

    bounds = [(0.0, 1.0)] * n_assets
    w_init = w_gmv.copy()
    weights = np.zeros((n_points, n_assets), dtype=float)
    realized_var = np.zeros(n_points, dtype=float)
    realized_ret = np.zeros(n_points, dtype=float)

    for idx, target in enumerate(targets):
        constraints = [
            {"type": "eq", "fun": lambda w: float(np.sum(w) - 1.0), "jac": lambda _w: np.ones(n_assets)},
            {
                "type": "eq",
                "fun": lambda w, t=target: float(w @ mu_est - t),
                "jac": lambda _w: mu_est,
            },
        ]
        res = minimize(
            lambda w: float(w @ sigma_est @ w),
            w_init,
            jac=lambda w: 2.0 * (sigma_est @ w),
            method="SLSQP",
            bounds=bounds,
            constraints=constraints,
            options={"ftol": 1e-12, "maxiter": 1000},
        )
        w_opt = np.asarray(res.x, dtype=float)
        w_init = w_opt
        weights[idx] = w_opt
        realized_var[idx] = float(w_opt @ sigma_true @ w_opt)
        realized_ret[idx] = float(w_opt @ mu_true)

    return {
        "targets": targets,
        "weights": weights,
        "realized_std": np.sqrt(realized_var),
        "realized_var": realized_var,
        "realized_ret": realized_ret,
    }


def marchenko_pastur_bounds(q: float) -> tuple[float, float]:
    """Batas bawah dan atas spektrum Marchenko-Pastur ``[(1 - sqrt(q))^2, (1 + sqrt(q))^2]``."""
    sq = np.sqrt(q)
    return float((1.0 - sq) ** 2), float((1.0 + sq) ** 2)


def marchenko_pastur_pdf(lam: np.ndarray, q: float) -> np.ndarray:
    """Kerapatan probabilitas teoretis Marchenko-Pastur untuk kasus ``Sigma = I``."""
    lam_min, lam_max = marchenko_pastur_bounds(q)
    arr = np.asarray(lam, dtype=float)
    inside = (arr >= lam_min) & (arr <= lam_max) & (arr > 0.0)
    pdf = np.zeros_like(arr)
    vals = arr[inside]
    pdf[inside] = np.sqrt((lam_max - vals) * (vals - lam_min)) / (2.0 * np.pi * q * vals)
    return pdf


def sample_correlation_eigenvalues(
    n_assets: int,
    n_obs: int,
    seed: int = DEFAULT_SEED,
) -> np.ndarray:
    """Eigenvalue matriks korelasi sampel dari data Gaussian ``N(0, I)``."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(size=(n_obs, n_assets))
    corr = np.corrcoef(z, rowvar=False)
    return np.linalg.eigvalsh(corr)


def run_wishart_lw_experiment(
    sigma_true: np.ndarray,
    mu: np.ndarray,
    n_obs: int = 200,
    n_trials: int = 5_000,
    batch_size: int = 500,
    seed: int = DEFAULT_SEED,
) -> dict[str, np.ndarray | float]:
    """Eksperimen Monte Carlo batched untuk rasio in-sample, realized sampel, dan Ledoit-Wolf."""
    _, sigma_min = gmv_closed_form(sigma_true)
    sigma_min_val = float(sigma_min)

    n_batches = (n_trials + batch_size - 1) // batch_size
    child_seeds = np.random.SeedSequence(seed).spawn(n_batches)

    in_sample_ratios = np.empty(n_trials, dtype=float)
    realized_sample_ratios = np.empty(n_trials, dtype=float)
    realized_lw_ratios = np.empty(n_trials, dtype=float)
    lw_intensities = np.empty(n_trials, dtype=float)

    offset = 0
    for b_idx, c_seed in enumerate(child_seeds):
        cur_m = min(batch_size, n_trials - offset)
        rng = np.random.default_rng(c_seed)
        rets = simulate_returns(mu, sigma_true, n_obs=n_obs, n_trials=cur_m, seed=rng)

        s_mat = sample_covariance(rets, ddof=1)
        w_s, var_in = gmv_closed_form(s_mat)
        var_real_s = np.einsum("bi,ij,bj->b", w_s, sigma_true, w_s)

        sigma_lw, intensity = ledoit_wolf(rets)
        w_lw, _ = gmv_closed_form(sigma_lw)
        var_real_lw = np.einsum("bi,ij,bj->b", w_lw, sigma_true, w_lw)

        sl = slice(offset, offset + cur_m)
        in_sample_ratios[sl] = var_in / sigma_min_val
        realized_sample_ratios[sl] = var_real_s / sigma_min_val
        realized_lw_ratios[sl] = var_real_lw / sigma_min_val
        lw_intensities[sl] = intensity
        offset += cur_m

    return {
        "sigma_min": sigma_min_val,
        "in_sample_ratio": in_sample_ratios,
        "realized_sample_ratio": realized_sample_ratios,
        "realized_lw_ratio": realized_lw_ratios,
        "lw_intensity": lw_intensities,
    }


def run_q_sweep(
    n_grid: list[int] | tuple[int, ...] = (10, 25, 50, 100, 150),
    n_obs: int = 200,
    n_trials: int = 2_000,
    seed: int = DEFAULT_SEED,
) -> list[dict[str, float]]:
    """Sapuan rasio ``q = N / T`` untuk membandingkan rasio empiris vs teoretis."""
    results: list[dict[str, float]] = []
    child_seeds = np.random.SeedSequence(seed).spawn(len(n_grid))
    for n_assets, c_seed in zip(n_grid, child_seeds):
        sigma_true, mu = generate_factor_model(n_assets=n_assets, n_factors=3, seed=DEFAULT_SEED)
        exp = run_wishart_lw_experiment(
            sigma_true, mu, n_obs=n_obs, n_trials=n_trials, batch_size=500, seed=int(c_seed.generate_state(1)[0])
        )
        in_r = exp["in_sample_ratio"]
        re_r = exp["realized_sample_ratio"]
        lw_r = exp["realized_lw_ratio"]
        theo_in = float(n_obs - n_assets) / float(n_obs - 1)
        theo_re = float(n_obs - 2) / float(n_obs - n_assets - 1)

        one_ret = simulate_returns(mu, sigma_true, n_obs=n_obs, n_trials=1, seed=int(c_seed.generate_state(1)[0]))
        s_one = sample_covariance(one_ret)[0]

        results.append(
            {
                "n": int(n_assets),
                "t": int(n_obs),
                "q": float(n_assets) / float(n_obs),
                "in_sample_mean": float(np.mean(in_r)),
                "in_sample_se": float(np.std(in_r, ddof=1) / np.sqrt(n_trials)),
                "in_sample_theo": theo_in,
                "realized_mean": float(np.mean(re_r)),
                "realized_se": float(np.std(re_r, ddof=1) / np.sqrt(n_trials)),
                "realized_theo": theo_re,
                "lw_realized_mean": float(np.mean(lw_r)),
                "lw_realized_se": float(np.std(lw_r, ddof=1) / np.sqrt(n_trials)),
                "cond_true": float(np.linalg.cond(sigma_true)),
                "cond_sample": float(np.linalg.cond(s_one)),
            }
        )
    return results
