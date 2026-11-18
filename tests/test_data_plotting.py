import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from qb100 import data, plotting  # noqa: E402


def test_cache_path_is_deterministic_and_safe(tmp_path):
    p = data.cache_path("bbca.jk", "2024-01-01", None, cache_dir=tmp_path)
    assert p == tmp_path / "BBCA.JK_2024-01-01_max_1d.csv"
    assert data.cache_path("^GSPC", cache_dir=tmp_path).name == "_GSPC_min_max_1d.csv"


def test_get_prices_reads_cache_without_network(tmp_path):
    idx = pd.date_range("2024-01-01", periods=3, freq="D")
    df = pd.DataFrame({"Close": [1.0, 2.0, 3.0]}, index=idx)
    df.to_csv(data.cache_path("FAKE", cache_dir=tmp_path))
    out = data.get_prices("FAKE", cache_dir=tmp_path)
    assert out["Close"].tolist() == [1.0, 2.0, 3.0]


def test_save_figure_writes_file(tmp_path):
    plotting.set_style()
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    path = plotting.save_figure(fig, "test_fig", folder=tmp_path)
    plt.close(fig)
    assert path.exists() and path.suffix == ".png"
