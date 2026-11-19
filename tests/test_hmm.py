"""Unit test untuk modul qb100.hmm."""

from __future__ import annotations

import numpy as np

from qb100 import hmm


def test_brute_force_matches_forward_backward_and_viterbi():
    for seed in (7, 42, 99):
        obs, _ = hmm.simulate_hmm(n_obs=10, seed=seed)
        bf = hmm.brute_force_hmm(obs, hmm.TRUE_PI, hmm.TRUE_A, hmm.TRUE_MU, hmm.TRUE_SIGMA)
        fb = hmm.forward_backward_scaled(obs, hmm.TRUE_PI, hmm.TRUE_A, hmm.TRUE_MU, hmm.TRUE_SIGMA)
        vit_path, vit_lp = hmm.viterbi(obs, hmm.TRUE_PI, hmm.TRUE_A, hmm.TRUE_MU, hmm.TRUE_SIGMA)

        assert abs(fb["log_likelihood"] - bf["log_likelihood"]) < 1e-12
        assert np.allclose(fb["gamma"], bf["gamma"], atol=1e-12)
        assert np.allclose(fb["xi_sum"], bf["xi_sum"], atol=1e-12)
        assert np.array_equal(vit_path, bf["viterbi_path"])
        assert abs(vit_lp - bf["viterbi_log_prob"]) < 1e-12


def test_naive_forward_overflows_while_scaled_is_stable():
    obs, _ = hmm.simulate_hmm(n_obs=2000, seed=42)
    alpha_naive = hmm.forward_naive(obs, hmm.TRUE_PI, hmm.TRUE_A, hmm.TRUE_MU, hmm.TRUE_SIGMA)
    assert not np.all(np.isfinite(alpha_naive[-1]))

    fb = hmm.forward_backward_scaled(obs, hmm.TRUE_PI, hmm.TRUE_A, hmm.TRUE_MU, hmm.TRUE_SIGMA)
    assert np.isfinite(fb["log_likelihood"])
    assert np.allclose(fb["gamma"].sum(axis=1), 1.0, atol=1e-12)


def test_em_monotonic_sorted_and_reproducible():
    obs, _ = hmm.simulate_hmm(n_obs=1000, seed=42)
    fit1 = hmm.fit_em(obs, n_states=2, max_iter=200, tol=1e-8)
    fit2 = hmm.fit_em(obs, n_states=2, max_iter=200, tol=1e-8)

    assert fit1["converged"]
    assert np.allclose(fit1["trans_a"].sum(axis=1), 1.0, atol=1e-12)
    assert fit1["sigma"][0] < fit1["sigma"][1]

    diffs = np.diff(fit1["ll_history"])
    assert np.all(diffs >= -1e-9)

    assert np.allclose(fit1["mu"], fit2["mu"], atol=1e-15)
    assert np.allclose(fit1["sigma"], fit2["sigma"], atol=1e-15)
    assert np.allclose(fit1["trans_a"], fit2["trans_a"], atol=1e-15)
