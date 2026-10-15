import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest

import pydada2 as dada2
from pydada2 import align
from pydada2.datasets import err_balanced_f, example_path
from pydada2.io import quality_profile
from pydada2.plotting import (
    error_plot_data,
    plot_complexity,
    plot_errors,
    plot_quality_profile,
    quality_profile_stats,
)

from rparity import load


@pytest.mark.parity
@pytest.mark.parametrize("fname", ["sam1F.fastq.gz", "samPB.fastq.gz"])
def test_quality_profile_stats_parity(fname):
    st = quality_profile_stats(quality_profile(example_path(fname)))
    r = load("quality_profile")[fname]["stats"]
    for py, rk in [("cycle", "Cycle"), ("mean", "Mean"), ("q25", "Q25"), ("q50", "Q50"),
                   ("q75", "Q75"), ("cum", "Cum")]:  # fmt: skip
        np.testing.assert_allclose(st[py].to_numpy(), np.array(r[rk], float), rtol=1e-12)


def test_plots_smoke():
    f1, f2 = example_path("sam1F.fastq.gz"), example_path("samPB.fastq.gz")
    fig = plot_quality_profile([f1, f2])
    assert len(fig.axes) == 2
    assert len(plot_quality_profile([f1, f2], aggregate=True).axes) == 1
    le = dada2.learn_errors(f1, n_bases=1)
    fig = plot_errors(le, nominal_q=True, err_in=True)
    assert len(fig.axes) == 16
    assert len(plot_errors(le, nti=["A"], ntj=["C", "G"]).axes) == 2
    with pytest.raises(ValueError, match="nti and ntj"):
        plot_errors(le, nti=["X"])
    assert len(plot_complexity([f1, f1], aggregate=True).axes) == 1
    assert len(plot_complexity(f1, n=50).axes) == 1


def test_error_plot_data():
    df = error_plot_data(err_balanced_f())
    assert df["observed"].isna().all() and df["estimated"].notna().any()
    with pytest.raises(ValueError, match="quality scores"):
        error_plot_data(np.full((16, 1), 0.01))
    with pytest.raises(ValueError, match="must be available"):
        error_plot_data({"err_out": None, "err_in": None, "trans": None})


def test_align_api():
    a, b = align.nwalign("ACGTACGT", "ACGACGT")
    assert (a, b) == ("ACGTACGT", "ACG-ACGT")
    assert align.nwalign("ACGTACGT", "ACGACGT", vec=True) == (a, b)
    assert align.eval_pair(a, b) == (7, 0, 1)
    assert align.pair_consensus(a, b) == "ACGTACGT"
    assert align.nwhamming("ACGTACGTTT", "ACGAACGTTT") == 1
    assert align.nwhamming(["ACGTACGTTT", "ACGT"], ["ACGAACGTTT", "ACGT"]).tolist() == [1, 0]
    assert align.nwhamming(["ACGTACGTTT", "ACGAACGTTT"], "ACGAACGTTT").tolist() == [1, 0]
    with pytest.raises(ValueError, match="same length"):
        align.nwhamming(["A", "C", "G"], ["A", "C"])
    with pytest.raises(ValueError, match="A/C/G/T"):
        align.nwalign("ACGN", "ACGT")
    with pytest.raises(ValueError, match="Homopolymer"):
        align.nwalign("ACGT", "ACGT", vec=True, homo_gap=-1)
    with pytest.raises(TypeError):
        align.nwalign(1, "ACGT")  # type: ignore[arg-type]
