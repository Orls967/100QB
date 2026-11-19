"""Hidden Markov Model (HMM) dua/multi-rezim Gaussian dengan NumPy murni."""

from __future__ import annotations

import itertools

import numpy as np

DEFAULT_SEED = 42
TRUE_MU = np.array([0.0005, -0.0010], dtype=float)
TRUE_SIGMA = np.array([0.008, 0.025], dtype=float)
TRUE_A = np.array([[0.98, 0.02], [0.05, 0.95]], dtype=float)
TRUE_PI = np.array([0.05 / 0.07, 0.02 / 0.07], dtype=float)


def stationary_distribution(trans_a: np.ndarray) -> np.ndarray:
    """Hitung distribusi stasioner ``pi`` dari eigenvector kiri matriks transisi ``A``."""
    eigvals, eigvecs = np.linalg.eig(trans_a.T)
    idx = int(np.argmin(np.abs(eigvals - 1.0)))
    vec = np.real(eigvecs[:, idx])
    vec = np.maximum(vec, 0.0)
    return vec / np.sum(vec)


def simulate_hmm(
    n_obs: int = 2000,
    mu: np.ndarray = TRUE_MU,
    sigma: np.ndarray = TRUE_SIGMA,
    trans_a: np.ndarray = TRUE_A,
    pi: np.ndarray = TRUE_PI,
    seed: int | np.random.Generator = DEFAULT_SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulasikan satu deret return dan state tersembunyi HMM Gaussian."""
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    n_states = len(pi)
    states = np.empty(n_obs, dtype=int)
    states[0] = int(rng.choice(n_states, p=pi))
    cdf_a = np.cumsum(trans_a, axis=1)
    u_trans = rng.random(n_obs - 1)
    for t in range(1, n_obs):
        states[t] = int(np.searchsorted(cdf_a[states[t - 1]], u_trans[t - 1], side="right"))
        if states[t] >= n_states:
            states[t] = n_states - 1
    obs = rng.normal(loc=mu[states], scale=sigma[states])
    return obs, states


def simulate_hmm_batch(
    n_series: int = 200,
    n_obs: int = 2000,
    mu: np.ndarray = TRUE_MU,
    sigma: np.ndarray = TRUE_SIGMA,
    trans_a: np.ndarray = TRUE_A,
    pi: np.ndarray = TRUE_PI,
    seed: int = DEFAULT_SEED,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulasikan ``n_series`` deret independen HMM Gaussian dengan ``SeedSequence.spawn``."""
    child_seeds = np.random.SeedSequence(seed).spawn(n_series)
    obs_all = np.empty((n_series, n_obs), dtype=float)
    states_all = np.empty((n_series, n_obs), dtype=int)
    for idx, c_seed in enumerate(child_seeds):
        obs_all[idx], states_all[idx] = simulate_hmm(
            n_obs=n_obs, mu=mu, sigma=sigma, trans_a=trans_a, pi=pi, seed=np.random.default_rng(c_seed)
        )
    return obs_all, states_all


def gaussian_emission(obs: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Densitas emisi Gaussian untuk 1D ``(T,)`` -> ``(T, K)`` atau batched ``(S, T)`` -> ``(S, T, K)``."""
    diff = obs[..., None] - mu[..., None, :]
    var = sigma[..., None, :] ** 2
    norm_const = 1.0 / (np.sqrt(2.0 * np.pi) * sigma[..., None, :])
    return norm_const * np.exp(-0.5 * (diff * diff) / var)


def forward_naive(obs: np.ndarray, pi: np.ndarray, trans_a: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Forward algorithm naif tanpa skala (rawan overflow/underflow untuk T besar)."""
    b = gaussian_emission(obs, mu, sigma)
    n_obs, n_states = b.shape
    alpha = np.empty((n_obs, n_states), dtype=float)
    alpha[0] = pi * b[0]
    for t in range(1, n_obs):
        alpha[t] = (alpha[t - 1] @ trans_a) * b[t]
    return alpha


def forward_backward_scaled(
    obs: np.ndarray,
    pi: np.ndarray,
    trans_a: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> dict[str, np.ndarray | float]:
    """Forward-backward ber-skala (Rabiner) untuk 1D ``(T,)`` atau batched ``(S, T)``."""
    is_1d = obs.ndim == 1
    obs_b = obs[None, :] if is_1d else obs
    pi_b = pi[None, :] if pi.ndim == 1 else pi
    a_b = trans_a[None, :, :] if trans_a.ndim == 2 else trans_a
    mu_b = mu[None, :] if mu.ndim == 1 else mu
    sig_b = sigma[None, :] if sigma.ndim == 1 else sigma

    b = gaussian_emission(obs_b, mu_b, sig_b)
    n_series, n_obs, n_states = b.shape

    alpha_hat = np.empty((n_series, n_obs, n_states), dtype=float)
    c_scale = np.empty((n_series, n_obs), dtype=float)

    alpha_0 = pi_b * b[:, 0, :]
    c_0 = np.maximum(alpha_0.sum(axis=-1, keepdims=True), 1e-300)
    c_scale[:, 0] = c_0[:, 0]
    alpha_hat[:, 0, :] = alpha_0 / c_0

    for t in range(1, n_obs):
        prev = np.einsum("si,sij->sj", alpha_hat[:, t - 1, :], a_b)
        alpha_t = prev * b[:, t, :]
        c_t = np.maximum(alpha_t.sum(axis=-1, keepdims=True), 1e-300)
        c_scale[:, t] = c_t[:, 0]
        alpha_hat[:, t, :] = alpha_t / c_t

    beta_hat = np.ones((n_series, n_obs, n_states), dtype=float)
    for t in range(n_obs - 2, -1, -1):
        weighted = b[:, t + 1, :] * beta_hat[:, t + 1, :]
        back = np.einsum("sij,sj->si", a_b, weighted)
        beta_hat[:, t, :] = back / c_scale[:, t + 1 : t + 2]

    gamma = alpha_hat * beta_hat
    gamma /= np.maximum(gamma.sum(axis=-1, keepdims=True), 1e-300)

    # Ekspektasi jumlah transisi xi_sum[s, i, j] = sum_{t=0}^{T-2} xi[s, t, i, j]
    rhs = (b[:, 1:, :] * beta_hat[:, 1:, :]) / c_scale[:, 1:, None]
    xi_sum = np.einsum("sti,sij,stj->sij", alpha_hat[:, :-1, :], a_b, rhs)
    log_lik = np.sum(np.log(c_scale), axis=-1)

    if is_1d:
        return {
            "log_likelihood": float(log_lik[0]),
            "gamma": gamma[0],
            "xi_sum": xi_sum[0],
            "alpha_hat": alpha_hat[0],
            "beta_hat": beta_hat[0],
            "c_scale": c_scale[0],
        }
    return {"log_likelihood": log_lik, "gamma": gamma, "xi_sum": xi_sum}


def viterbi(
    obs: np.ndarray,
    pi: np.ndarray,
    trans_a: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Algoritma Viterbi di ruang log untuk mendekode jalur state MAP."""
    b = gaussian_emission(obs, mu, sigma)
    log_b = np.log(np.maximum(b, 1e-300))
    log_pi = np.log(np.maximum(pi, 1e-300))
    log_a = np.log(np.maximum(trans_a, 1e-300))

    n_obs, n_states = log_b.shape
    delta = np.empty((n_obs, n_states), dtype=float)
    psi = np.zeros((n_obs, n_states), dtype=int)

    delta[0] = log_pi + log_b[0]
    for t in range(1, n_obs):
        scores = delta[t - 1][:, None] + log_a
        psi[t] = np.argmax(scores, axis=0)
        delta[t] = np.max(scores, axis=0) + log_b[t]

    path = np.empty(n_obs, dtype=int)
    path[-1] = int(np.argmax(delta[-1]))
    best_log_prob = float(delta[-1, path[-1]])
    for t in range(n_obs - 2, -1, -1):
        path[t] = psi[t + 1, path[t + 1]]
    return path, best_log_prob


def brute_force_hmm(
    obs: np.ndarray,
    pi: np.ndarray,
    trans_a: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
) -> dict[str, np.ndarray | float]:
    """Enumerasi eksak semua ``K^T`` jalur state untuk verifikasi (khusus ``T`` kecil)."""
    n_obs = len(obs)
    n_states = len(pi)
    if n_states**n_obs > 100_000:
        raise ValueError("K^T terlalu besar untuk brute force")

    log_b = np.log(np.maximum(gaussian_emission(obs, mu, sigma), 1e-300))
    log_pi = np.log(np.maximum(pi, 1e-300))
    log_a = np.log(np.maximum(trans_a, 1e-300))

    paths = np.array(list(itertools.product(range(n_states), repeat=n_obs)), dtype=int)
    log_joints = log_pi[paths[:, 0]] + log_b[0, paths[:, 0]]
    for t in range(1, n_obs):
        log_joints += log_a[paths[:, t - 1], paths[:, t]] + log_b[t, paths[:, t]]

    max_lp = float(np.max(log_joints))
    log_lik = max_lp + float(np.log(np.sum(np.exp(log_joints - max_lp))))
    posteriors = np.exp(log_joints - log_lik)

    gamma = np.zeros((n_obs, n_states), dtype=float)
    for t in range(n_obs):
        for k in range(n_states):
            gamma[t, k] = float(np.sum(posteriors[paths[:, t] == k]))

    xi_sum = np.zeros((n_states, n_states), dtype=float)
    for t in range(n_obs - 1):
        for i in range(n_states):
            for j in range(n_states):
                mask = (paths[:, t] == i) & (paths[:, t + 1] == j)
                xi_sum[i, j] += float(np.sum(posteriors[mask]))

    best_idx = int(np.argmax(log_joints))
    return {
        "log_likelihood": log_lik,
        "gamma": gamma,
        "xi_sum": xi_sum,
        "viterbi_path": paths[best_idx],
        "viterbi_log_prob": float(log_joints[best_idx]),
    }


def _init_em_params(obs_b: np.ndarray, n_states: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Inisialisasi deterministik berbasis kuantil ``|r|`` untuk 1 atau banyak deret."""
    n_series = obs_b.shape[0]
    pi = np.full((n_series, n_states), 1.0 / n_states, dtype=float)

    if n_states == 1:
        mu = np.mean(obs_b, axis=1, keepdims=True)
        sigma = np.maximum(np.std(obs_b, axis=1, keepdims=True), 1e-6)
        trans_a = np.ones((n_series, 1, 1), dtype=float)
        return pi, trans_a, mu, sigma

    a_init = np.full((n_states, n_states), 0.1 / (n_states - 1), dtype=float)
    np.fill_diagonal(a_init, 0.9)
    trans_a = np.broadcast_to(a_init, (n_series, n_states, n_states)).copy()

    abs_obs = np.abs(obs_b)
    quantiles = np.linspace(0.0, 1.0, n_states + 1)[1:-1]
    edges = np.quantile(abs_obs, quantiles, axis=1).T  # (S, K-1)

    mu = np.empty((n_series, n_states), dtype=float)
    sigma = np.empty((n_series, n_states), dtype=float)
    for k in range(n_states):
        if k == 0:
            mask = abs_obs <= edges[:, 0:1]
        elif k == n_states - 1:
            mask = abs_obs > edges[:, -1:]
        else:
            mask = (abs_obs > edges[:, k - 1 : k]) & (abs_obs <= edges[:, k : k + 1])
        counts = np.maximum(mask.sum(axis=1), 1)
        mu[:, k] = (obs_b * mask).sum(axis=1) / counts
        var_k = ((obs_b - mu[:, k : k + 1]) ** 2 * mask).sum(axis=1) / counts
        sigma[:, k] = np.maximum(np.sqrt(var_k), 1e-6)

    return pi, trans_a, mu, sigma


def fit_em(
    obs: np.ndarray,
    n_states: int = 2,
    max_iter: int = 500,
    tol: float = 1e-8,
) -> dict[str, np.ndarray | float | bool | int | list[float]]:
    """Fit Baum-Welch (EM) pada satu deret ``(T,)`` dan urutkan state menurut ``sigma`` menaik."""
    res = fit_em_batch(obs[None, :], n_states=n_states, max_iter=max_iter, tol=tol, track_history=True)
    return {
        "pi": res["pi"][0],
        "trans_a": res["trans_a"][0],
        "mu": res["mu"][0],
        "sigma": res["sigma"][0],
        "log_likelihood": float(res["log_likelihood"][0]),
        "n_iter": int(res["n_iter"][0]),
        "converged": bool(res["converged"][0]),
        "ll_history": res["ll_history"][0],
    }


def fit_em_batch(
    obs_batch: np.ndarray,
    n_states: int = 2,
    max_iter: int = 500,
    tol: float = 1e-8,
    track_history: bool = False,
) -> dict[str, np.ndarray | list[list[float]]]:
    """Fit Baum-Welch (EM) tervektorisasi pada batch ``(S, T)`` dan urutkan state per deret."""
    n_series = obs_batch.shape[0]
    pi, trans_a, mu, sigma = _init_em_params(obs_batch, n_states)

    active = np.ones(n_series, dtype=bool)
    prev_ll = np.full(n_series, -np.inf, dtype=float)
    converged = np.zeros(n_series, dtype=bool)
    n_iter = np.zeros(n_series, dtype=int)
    histories: list[list[float]] = [[] for _ in range(n_series)] if track_history else []

    for it in range(1, max_iter + 1):
        idx = np.where(active)[0]
        if idx.size == 0:
            break

        fb = forward_backward_scaled(obs_batch[idx], pi[idx], trans_a[idx], mu[idx], sigma[idx])
        ll_cur = fb["log_likelihood"]
        gamma = fb["gamma"]
        xi_sum = fb["xi_sum"]

        if track_history:
            for pos, s_idx in enumerate(idx):
                histories[s_idx].append(float(ll_cur[pos]))

        # M-step pada deret yang masih aktif
        pi[idx] = gamma[:, 0, :]
        denom_a = np.maximum(gamma[:, :-1, :].sum(axis=1)[:, :, None], 1e-300)
        trans_a[idx] = xi_sum / denom_a

        gamma_sum = np.maximum(gamma.sum(axis=1), 1e-300)
        mu_new = np.einsum("stk,st->sk", gamma, obs_batch[idx]) / gamma_sum
        diff_sq = (obs_batch[idx, :, None] - mu_new[:, None, :]) ** 2
        var_new = np.einsum("stk,stk->sk", gamma, diff_sq) / gamma_sum
        mu[idx] = mu_new
        sigma[idx] = np.maximum(np.sqrt(var_new), 1e-6)

        n_iter[idx] = it
        done = np.abs(ll_cur - prev_ll[idx]) < tol
        converged[idx[done]] = True
        active[idx[done]] = False
        prev_ll[idx] = ll_cur

    # Hitung log-likelihood akhir setelah M-step terakhir
    final_ll = forward_backward_scaled(obs_batch, pi, trans_a, mu, sigma)["log_likelihood"]

    # Urutkan label state menurut sigma menaik agar bebas label switching
    order = np.argsort(sigma, axis=1)
    rows = np.arange(n_series)[:, None]
    pi = pi[rows, order]
    mu = mu[rows, order]
    sigma = sigma[rows, order]
    trans_a = trans_a[rows[:, :, None], order[:, :, None], order[:, None, :]]

    out: dict[str, np.ndarray | list[list[float]]] = {
        "pi": pi,
        "trans_a": trans_a,
        "mu": mu,
        "sigma": sigma,
        "log_likelihood": final_ll,
        "n_iter": n_iter,
        "converged": converged,
    }
    if track_history:
        out["ll_history"] = histories
    return out


def bic_score(log_likelihood: float, n_obs: int, n_states: int) -> tuple[float, int]:
    """Hitung BIC dan jumlah parameter bebas: ``2`` (K=1) atau ``(K-1) + K(K-1) + 2K`` (K>1)."""
    n_params = 2 if n_states == 1 else (n_states - 1) + n_states * (n_states - 1) + 2 * n_states
    bic = -2.0 * float(log_likelihood) + n_params * float(np.log(n_obs))
    return bic, n_params


def viterbi_diagnostics(
    true_states: np.ndarray,
    pred_states: np.ndarray,
    n_states: int = 2,
) -> dict[str, np.ndarray | float | list[float]]:
    """Akurasi, matriks konfusi, durasi run rata-rata, dan fraksi waktu empiris."""
    acc = float(np.mean(true_states == pred_states))
    conf = np.zeros((n_states, n_states), dtype=int)
    for t_val, p_val in zip(true_states, pred_states):
        conf[int(t_val), int(p_val)] += 1

    durations: list[float] = []
    for k in range(n_states):
        runs: list[int] = []
        cur = 0
        for s in pred_states:
            if s == k:
                cur += 1
            elif cur > 0:
                runs.append(cur)
                cur = 0
        if cur > 0:
            runs.append(cur)
        durations.append(float(np.mean(runs)) if runs else 0.0)

    fracs = np.array([float(np.mean(pred_states == k)) for k in range(n_states)])
    return {"accuracy": acc, "confusion": conf, "mean_durations": durations, "state_fractions": fracs}


def reliability_diagram(
    gamma_state1: np.ndarray,
    true_states: np.ndarray,
    n_bins: int = 10,
) -> dict[str, np.ndarray]:
    """Kelompokkan ``gamma(state 1)`` ke ``n_bins`` desil vs frekuensi empiris state 1."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_pred = np.full(n_bins, np.nan)
    bin_true = np.full(n_bins, np.nan)
    bin_counts = np.zeros(n_bins, dtype=int)

    for b in range(n_bins):
        if b == n_bins - 1:
            mask = (gamma_state1 >= edges[b]) & (gamma_state1 <= edges[b + 1])
        else:
            mask = (gamma_state1 >= edges[b]) & (gamma_state1 < edges[b + 1])
        cnt = int(np.sum(mask))
        bin_counts[b] = cnt
        if cnt > 0:
            bin_pred[b] = float(np.mean(gamma_state1[mask]))
            bin_true[b] = float(np.mean(true_states[mask] == 1))
    return {"pred_mean": bin_pred, "true_freq": bin_true, "counts": bin_counts}
