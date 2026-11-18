import numpy as np
import pytest

from qb100 import simulation as sim


def test_get_rng_same_seed_is_reproducible():
    a = sim.get_rng(123).normal(size=5)
    b = sim.get_rng(123).normal(size=5)
    np.testing.assert_array_equal(a, b)


def test_get_rng_returns_given_generator():
    g = np.random.default_rng(1)
    assert sim.get_rng(g) is g


def test_get_rng_none_uses_env_seed(monkeypatch):
    monkeypatch.setenv("QB100_SEED", "7")
    np.testing.assert_array_equal(
        sim.get_rng(None).random(3), np.random.default_rng(7).random(3)
    )


def test_get_rng_none_default_seed(monkeypatch):
    monkeypatch.delenv("QB100_SEED", raising=False)
    np.testing.assert_array_equal(
        sim.get_rng(None).random(3), np.random.default_rng(sim.DEFAULT_SEED).random(3)
    )


def test_random_walk_shape_and_start():
    paths = sim.random_walk(100, n_paths=5, start=10.0, seed=0)
    assert paths.shape == (101, 5)
    np.testing.assert_array_equal(paths[0], 10.0)


def test_random_walk_reproducible():
    np.testing.assert_array_equal(sim.random_walk(50, seed=3), sim.random_walk(50, seed=3))


def test_random_walk_zero_std_is_pure_drift():
    paths = sim.random_walk(4, step_std=0.0, drift=2.0, start=1.0, seed=0)
    np.testing.assert_allclose(paths[:, 0], [1.0, 3.0, 5.0, 7.0, 9.0])


def test_random_walk_step_statistics():
    paths = sim.random_walk(20_000, step_std=2.0, drift=0.5, seed=1)
    steps = np.diff(paths[:, 0])
    assert steps.mean() == pytest.approx(0.5, abs=0.05)
    assert steps.std() == pytest.approx(2.0, rel=0.03)


def test_random_walk_invalid_args():
    with pytest.raises(ValueError):
        sim.random_walk(-1)
    with pytest.raises(ValueError):
        sim.random_walk(10, n_paths=0)


def _pi_trial(rng):
    x, y = rng.random(2)
    return 4.0 * (x * x + y * y <= 1.0)


def test_monte_carlo_estimates_pi():
    samples = sim.monte_carlo(_pi_trial, 20_000, seed=0)
    assert samples.shape == (20_000,)
    assert samples.mean() == pytest.approx(np.pi, abs=0.05)


def test_monte_carlo_reproducible():
    a = sim.monte_carlo(lambda r: r.normal(), 100, seed=5)
    b = sim.monte_carlo(lambda r: r.normal(), 100, seed=5)
    np.testing.assert_array_equal(a, b)


def test_monte_carlo_invalid_n():
    with pytest.raises(ValueError):
        sim.monte_carlo(lambda r: 0.0, 0)


def test_mc_summary_keys_and_values():
    x = np.arange(1, 101, dtype=float)
    s = sim.mc_summary(x, ci=0.9)
    assert set(s) == {"mean", "std", "se", "ci_low", "ci_high"}
    assert s["mean"] == pytest.approx(50.5)
    assert s["std"] == pytest.approx(x.std(ddof=1))
    assert s["se"] == pytest.approx(x.std(ddof=1) / 10)
    assert s["ci_low"] == pytest.approx(np.quantile(x, 0.05))
    assert s["ci_high"] == pytest.approx(np.quantile(x, 0.95))


def test_mc_summary_needs_two_samples():
    with pytest.raises(ValueError):
        sim.mc_summary([1.0])


def test_running_mean_values():
    rm = sim.running_mean([2.0, 4.0, 6.0, 8.0])
    np.testing.assert_allclose(rm, [2.0, 3.0, 4.0, 5.0])


def test_running_mean_empty_raises():
    with pytest.raises(ValueError):
        sim.running_mean([])
