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


def test_volatility_annualized():
    r = pd.Series([0.01, -0.02, 0.015, 0.0, -0.005])
    expected = r.std(ddof=1) * np.sqrt(252)
    assert fin.volatility(r) == pytest.approx(expected)


def test_volatility_not_annualized_and_custom_periods():
    r = pd.Series([0.01, -0.02, 0.015, 0.0, -0.005])
    assert fin.volatility(r, annualize=False) == pytest.approx(r.std(ddof=1))
    assert fin.volatility(r, periods_per_year=12) == pytest.approx(r.std(ddof=1) * np.sqrt(12))


def test_volatility_constant_returns_is_zero():
    assert fin.volatility([0.01, 0.01, 0.01]) == pytest.approx(0.0)


def test_sharpe_ratio_matches_formula():
    r = pd.Series([0.01, -0.005, 0.02, 0.003, -0.01, 0.007])
    rf = 0.05
    ex = r - rf / 252
    expected = ex.mean() / ex.std(ddof=1) * np.sqrt(252)
    assert fin.sharpe_ratio(r, risk_free=rf) == pytest.approx(expected)


def test_sharpe_ratio_zero_std_is_nan():
    assert np.isnan(fin.sharpe_ratio([0.01, 0.01, 0.01]))


def test_sharpe_ratio_dataframe_per_column():
    df = pd.DataFrame({"a": [0.01, 0.02, -0.01, 0.0], "b": [0.01, 0.01, 0.01, 0.01]})
    sr = fin.sharpe_ratio(df)
    assert sr["a"] == pytest.approx(df["a"].mean() / df["a"].std() * np.sqrt(252))
    assert np.isnan(sr["b"])


def test_drawdown_series(prices):
    dd = fin.drawdown(prices)
    np.testing.assert_allclose(dd.values, [0.0, 0.0, -0.1, 108.9 / 110 - 1, 0.0])
    assert (dd <= 0).all()


def test_max_drawdown_value(prices):
    assert fin.max_drawdown(prices) == pytest.approx(-0.1)


def test_max_drawdown_monotonic_increase_is_zero():
    assert fin.max_drawdown([1.0, 2.0, 3.0, 4.0]) == 0.0


def test_max_drawdown_dataframe():
    df = pd.DataFrame({"a": [100.0, 50.0, 75.0], "b": [10.0, 20.0, 15.0]})
    mdd = fin.max_drawdown(df)
    assert mdd["a"] == pytest.approx(-0.5)
    assert mdd["b"] == pytest.approx(-0.25)
