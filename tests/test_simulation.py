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
