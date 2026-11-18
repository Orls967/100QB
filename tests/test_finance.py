import numpy as np
import pandas as pd
import pytest

from qb100 import finance as fin


@pytest.fixture
def prices():
    return pd.Series([100.0, 110.0, 99.0, 108.9, 120.0])


def test_simple_returns_values(prices):
    r = fin.simple_returns(prices)
    assert len(r) == len(prices) - 1
    np.testing.assert_allclose(r.values, [0.10, -0.10, 0.10, 120.0 / 108.9 - 1])


def test_simple_returns_accepts_numpy_array():
    r = fin.simple_returns(np.array([1.0, 2.0, 1.0]))
    np.testing.assert_allclose(r.values, [1.0, -0.5])


def test_simple_returns_dataframe():
    df = pd.DataFrame({"a": [1.0, 2.0], "b": [10.0, 5.0]})
    r = fin.simple_returns(df)
    assert list(r.columns) == ["a", "b"]
    np.testing.assert_allclose(r.iloc[0].values, [1.0, -0.5])


def test_log_returns_values(prices):
    r = fin.log_returns(prices)
    np.testing.assert_allclose(r.values, np.log(prices.values[1:] / prices.values[:-1]))


def test_log_returns_sum_equals_total_log_growth(prices):
    r = fin.log_returns(prices)
    assert r.sum() == pytest.approx(np.log(prices.iloc[-1] / prices.iloc[0]))


def test_log_and_simple_returns_relation(prices):
    np.testing.assert_allclose(
        np.exp(fin.log_returns(prices).values) - 1, fin.simple_returns(prices).values
    )
