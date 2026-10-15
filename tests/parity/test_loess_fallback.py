"""Loess error model on few observed quality scores (REVIEW.md 1.3): parity with R."""

import numpy as np
import pytest

import pydada2 as dada2
from pydada2.datasets import tperr1
from pydada2.error_models import LoessErrorModel, _skmisc_loess

from rparity import INPUTS, load, matrix

pytestmark = pytest.mark.parity


@pytest.fixture(scope="module")
def rv():
    return load("review_python")


@pytest.mark.parametrize("k", range(len(load("loess_degenerate"))))
def test_loess_few_quality_scores(k):
    case = load("loess_degenerate")[k]
    trans = matrix(case["trans"])
    if not case["ok"][0]:  # R's loess fails too (fewer than 3 quality scores for a nucleotide)
        with pytest.raises(ValueError):
            LoessErrorModel()(trans)
        return
    np.testing.assert_allclose(LoessErrorModel()(trans), matrix(case["err"]), rtol=1e-10)


def test_loess_fallback_is_exercised():
    """The five binned-quality configurations make skmisc abort (R only warns)."""
    for case in load("loess_degenerate")[:5]:
        trans = matrix(case["trans"])
        tot = trans[0:4].sum(axis=0)
        ok = tot > 0
        q = np.arange(trans.shape[1], dtype=float)
        with pytest.raises(ValueError):
            _skmisc_loess(q[ok], np.log10((trans[1, ok] + 1) / tot[ok]), tot[ok], q[ok])


def test_learn_errors_binned_singletons(rv):
    fx = rv["binned"]
    path = INPUTS / "binned_singletons.fastq.gz"
    le = dada2.learn_errors(path)
    np.testing.assert_allclose(le.err_, matrix(fx["learn"]["err_out"]), rtol=1e-10)
    np.testing.assert_array_equal(le.trans_, matrix(fx["learn"]["trans"]))
    assert len(le.err_in_) == len(fx["learn"]["err_in"])
    for a, b in zip(le.err_in_, fx["learn"]["err_in"], strict=True):
        np.testing.assert_allclose(a, matrix(b), rtol=1e-10)
    dd = dada2.DADA(err=tperr1()).fit(dada2.derep_fastq(path))
    np.testing.assert_allclose(dd.err_out_, matrix(fx["dada_err_out"]), rtol=1e-10)
