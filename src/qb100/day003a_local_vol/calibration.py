"""Kalibrasi balik SSVI dengan Levenberg-Marquardt buatan sendiri dan studi Monte Carlo R=50."""

from __future__ import annotations

import numpy as np

from .ssvi import SSVIParams, ssvi_theta, ssvi_total_variance


def generate_noisy_ssvi_surface(
    params: SSVIParams = SSVIParams(),
    m_grid: np.ndarray | None = None,
    T_grid: np.ndarray | None = None,
    noise_std: float = 5e-4,
    seed: int = 20261121,
) -> dict[str, np.ndarray]:
    """Bangkitkan grid implied vol SSVI (7 strike x 5 maturitas) + noise N(0, noise_std^2)."""
    if m_grid is None:
        m_grid = np.array([0.8, 0.85, 0.9, 1.0, 1.1, 1.15, 1.2], dtype=np.float64)
    if T_grid is None:
        T_grid = np.array([0.25, 0.5, 1.0, 1.5, 2.0], dtype=np.float64)

    y_vals = np.log(m_grid)
    Y, TT = np.meshgrid(y_vals, T_grid, indexing="ij")
    w_true = ssvi_total_variance(Y, TT, params)
    iv_true = np.sqrt(w_true / TT)

    rng = np.random.Generator(np.random.PCG64(seed))
    iv_noisy = iv_true + rng.normal(0.0, noise_std, size=iv_true.shape)
    return {
        "m_grid": m_grid,
        "y_vals": y_vals,
        "T_grid": T_grid,
        "iv_true": iv_true,
        "iv_noisy": iv_noisy,
        "theta_true": np.asarray(ssvi_theta(T_grid, params), dtype=np.float64),
    }


def _ssvi_iv_and_jacobian(
    beta: np.ndarray,
    y_vals: np.ndarray,
    T_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Hitung vektor implied vol model (n_k * n_t) dan Jacobian analitik (n_obs, 3 + n_t).

    Parameter ``beta = [rho, eta, gamma, theta_0, ..., theta_{n_t-1}]``.
    """
    rho, eta, gamma = float(beta[0]), float(beta[1]), float(beta[2])
    thetas = beta[3:]
    n_k = len(y_vals)
    n_t = len(T_grid)
    n_obs = n_k * n_t
    n_par = 3 + n_t

    iv_vec = np.empty(n_obs, dtype=np.float64)
    jac = np.zeros((n_obs, n_par), dtype=np.float64)

    for j in range(n_t):
        T_j = float(T_grid[j])
        th = float(thetas[j])
        phi = eta / ((th ** gamma) * ((1.0 + th) ** (1.0 - gamma)))
        dphi_deta = phi / eta
        dphi_dgamma = phi * np.log((1.0 + th) / th)
        dphi_dth = -phi * (gamma + th) / (th * (1.0 + th))

        z = phi * y_vals + rho
        R = np.sqrt(z * z + (1.0 - rho * rho))
        A = 1.0 + rho * phi * y_vals + R
        w = 0.5 * th * A
        iv = np.sqrt(w / T_j)

        # Turunan parsial iv terhadap w: d(iv)/dw = 1 / (2 * T_j * iv)
        div_dw = 1.0 / (2.0 * T_j * iv)

        # Turunan A terhadap phi: (rho + z / R) * y
        dA_dphi = (rho + z / R) * y_vals
        # Turunan A terhadap rho: phi * y + (z - rho) / R = phi * y + phi * y / R = phi * y * (1 + 1/R)
        dA_drho = phi * y_vals + (z - rho) / R

        dw_drho = 0.5 * th * dA_drho
        dw_deta = 0.5 * th * dA_dphi * dphi_deta
        dw_dgamma = 0.5 * th * dA_dphi * dphi_dgamma
        dw_dth = 0.5 * A + 0.5 * th * dA_dphi * dphi_dth

        sl = slice(j * n_k, (j + 1) * n_k)
        iv_vec[sl] = iv
        jac[sl, 0] = div_dw * dw_drho
        jac[sl, 1] = div_dw * dw_deta
        jac[sl, 2] = div_dw * dw_dgamma
        jac[sl, 3 + j] = div_dw * dw_dth

    return iv_vec, jac


def _project_params(beta: np.ndarray) -> np.ndarray:
    """Proyeksi parameter ke domain valid SSVI."""
    out = beta.copy()
    out[0] = float(np.clip(out[0], -0.995, 0.995))
    out[1] = float(max(out[1], 1e-5))
    out[2] = float(np.clip(out[2], 0.01, 0.99))
    out[3:] = np.maximum(out[3:], 1e-6)
    return out


def calibrate_ssvi_lm(
    iv_obs: np.ndarray,
    y_vals: np.ndarray,
    T_grid: np.ndarray,
    initial_guess: np.ndarray | None = None,
    max_iter: int = 200,
    tol_step: float = 1e-10,
    tol_cost: float = 1e-14,
) -> dict[str, object]:
    """Kalibrasi (rho, eta, gamma, theta_1..theta_m) dengan Levenberg-Marquardt buatan sendiri."""
    n_k, n_t = iv_obs.shape
    # Urutkan observasi sesuai (j, i) agar konsisten dengan _ssvi_iv_and_jacobian
    obs_vec = iv_obs.T.reshape(-1)

    if initial_guess is None:
        atm_idx = int(np.argmin(np.abs(y_vals)))
        theta_init = np.square(iv_obs[atm_idx, :]) * T_grid * 0.9
        beta = np.concatenate([[-0.3, 0.6, 0.3], theta_init])
    else:
        beta = np.asarray(initial_guess, dtype=np.float64).copy()
    beta = _project_params(beta)

    lam = 1e-3
    iv_mod, jac = _ssvi_iv_and_jacobian(beta, y_vals, T_grid)
    res = iv_mod - obs_vec
    cost = 0.5 * float(np.dot(res, res))

    for it in range(max_iter):
        jtj = jac.T @ jac
        jtr = jac.T @ res
        diag = np.maximum(np.diag(jtj), 1e-12)

        accepted = False
        for _ in range(25):
            lhs = jtj + lam * np.diag(diag)
            try:
                step = -np.linalg.solve(lhs, jtr)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue

            cand = _project_params(beta + step)
            iv_cand, jac_cand = _ssvi_iv_and_jacobian(cand, y_vals, T_grid)
            res_cand = iv_cand - obs_vec
            cost_cand = 0.5 * float(np.dot(res_cand, res_cand))

            if cost_cand < cost:
                dcost = cost - cost_cand
                dstep = float(np.max(np.abs(cand - beta)))
                beta, iv_mod, jac, res, cost = cand, iv_cand, jac_cand, res_cand, cost_cand
                lam = max(lam * 0.3, 1e-12)
                accepted = True
                if dstep < tol_step or dcost < tol_cost:
                    break
                break
            lam *= 5.0

        if not accepted or float(np.max(np.abs(jtr))) < 1e-12:
            break
        if dstep < tol_step or dcost < tol_cost:
            break

    n_obs = len(obs_vec)
    n_par = len(beta)
    rms_err = float(np.sqrt(np.mean(np.square(res))))
    sigma2_hat = float(np.dot(res, res)) / max(n_obs - n_par, 1)
    cov = np.linalg.inv(jac.T @ jac) * sigma2_hat
    sd = np.sqrt(np.diag(cov))
    corr = cov / np.outer(sd, sd)

    return {
        "beta": beta,
        "rho": float(beta[0]),
        "eta": float(beta[1]),
        "gamma": float(beta[2]),
        "thetas": beta[3:],
        "rms_vol_error": rms_err,
        "jacobian": jac,
        "cov": cov,
        "corr": corr,
        "se": sd,
        "iterations": it + 1,
    }


def run_calibration_mc_study(
    params: SSVIParams = SSVIParams(),
    n_rep: int = 50,
    noise_std: float = 5e-4,
    base_seed: int = 20261121,
) -> dict[str, object]:
    """Ulangi kalibrasi LM pada ``n_rep`` noise seed untuk distribusi estimasi parameter."""
    estimates = []
    rms_list = []
    for rep in range(n_rep):
        data = generate_noisy_ssvi_surface(
            params=params, noise_std=noise_std, seed=base_seed + rep
        )
        fit = calibrate_ssvi_lm(data["iv_noisy"], data["y_vals"], data["T_grid"])
        estimates.append(fit["beta"])
        rms_list.append(fit["rms_vol_error"])

    est_mat = np.vstack(estimates)
    true_vec = np.concatenate([[params.rho, params.eta, params.gamma], data["theta_true"]])
    mean_est = np.mean(est_mat, axis=0)
    sd_est = np.std(est_mat, axis=0, ddof=1)
    bias = mean_est - true_vec

    return {
        "n_rep": n_rep,
        "true_params": true_vec,
        "estimates": est_mat,
        "mean": mean_est,
        "sd": sd_est,
        "bias": bias,
        "bias_over_sd": bias / np.maximum(sd_est, 1e-15),
        "mean_rms": float(np.mean(rms_list)),
    }
