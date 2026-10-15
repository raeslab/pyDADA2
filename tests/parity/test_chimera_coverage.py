"""Chimera methods and options on a table where they disagree (chimera_coverage.json).

On the original chimera fixture every method and option gave the same calls. The generator
builds a 5-sample table (6 parents, 5 chimeras: parents in all / 3 / 4 / 1 samples, an abundant
chimera whose parents occur in one sample only, and a one-off chimera) and checks that the
methods and options disagree on it. It also runs DADA on reads with injected chimeras (so
``is_bimera_denovo`` on a DADA result has positives) and isShiftDenovo on shifted / contained
sequences with short and long overlaps.
"""

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
import pydada2._core as core
from pydada2.chimeras import (
    BimeraRemover,
    is_bimera_denovo,
    is_bimera_denovo_table,
    is_shift_denovo,
    remove_bimera_denovo,
)
from pydada2.datasets import err_balanced_f
from pydada2.denoise import DADA

from rparity import INPUTS, load, matrix

pytestmark = pytest.mark.parity

CH = load("chimera_coverage")
R2PY = {
    "minSampleFraction": "min_sample_fraction",
    "ignoreNNegatives": "ignore_n_negatives",
    "minFoldParentOverAbundance": "min_fold_parent_over_abundance",
    "minParentAbundance": "min_parent_abundance",
    "allowOneOff": "allow_one_off",
    "maxShift": "max_shift",
}


def _kwargs(opts):
    if not isinstance(opts, dict):  # empty R list -> []
        return {}
    return {R2PY[k]: (v[0] if isinstance(v, list) else v) for k, v in opts.items()}


def _table(obj):
    return pd.DataFrame(
        matrix(obj, dtype=np.int64), index=obj.get("rownames"), columns=obj["colnames"]
    )


def _assert_table(got, obj):
    exp = _table(obj)
    assert got.columns.tolist() == exp.columns.tolist()
    assert [str(i) for i in got.index] == list(exp.index)
    np.testing.assert_array_equal(got.to_numpy(), exp.to_numpy())


ST = _table(CH["table"])


def test_table_bimera_counts():
    nflag, nsam = core.table_bimera(
        np.ascontiguousarray(ST.to_numpy(dtype=np.int32)), ST.columns.tolist(),
        1.5, 2, False, 4, 5, -4, -8, 16,
    )  # fmt: skip
    assert nflag.tolist() == CH["table_bimera"]["nflag"]
    assert nsam.tolist() == CH["table_bimera"]["nsam"]


@pytest.mark.parametrize("opt", sorted(CH["table_opts"]))
def test_consensus_calls(opt):
    got = is_bimera_denovo_table(ST, **_kwargs(CH["table_opts"][opt]))
    assert got.tolist() == CH["table_calls"][opt]
    assert got.index.tolist() == ST.columns.tolist()


@pytest.mark.parametrize("opt", sorted(CH["pooled_opts"]))
def test_pooled_calls(opt):
    got = is_bimera_denovo(ST.sum(axis=0), **_kwargs(CH["pooled_opts"][opt]))
    assert got.tolist() == CH["pooled_calls"][opt]


@pytest.mark.parametrize(
    "method,kwargs,key",
    [
        ("consensus", {}, "remove_consensus"),
        ("pooled", {}, "remove_pooled"),
        ("per-sample", {}, "remove_per_sample"),
        ("per-sample", {"allow_one_off": True}, "remove_per_sample_oneoff"),
    ],
)
def test_remove_methods(method, kwargs, key):
    _assert_table(remove_bimera_denovo(ST, method=method, **kwargs), CH[key])
    _assert_table(BimeraRemover(method=method, **kwargs).fit_transform(ST), CH[key])


def test_methods_disagree_on_this_table():
    kept = {
        k: CH[k]["colnames"] for k in ("remove_consensus", "remove_pooled", "remove_per_sample")
    }
    assert len({tuple(v) for v in kept.values()}) == 3


@pytest.fixture(scope="module")
def dada_chim():
    drp = dada2.derep_fastq(INPUTS / "chimeric_sam1F.fastq.gz")
    return DADA(err=err_balanced_f()).fit(drp)[0]


def test_dada_with_chimeras(dada_chim):
    fx = CH["dada_chim"]
    assert dada_chim.denoised.index.tolist() == fx["denoised"]["sequence"]
    assert dada_chim.denoised.tolist() == fx["denoised"]["abundance"]
    assert any(fx["is_bimera"])
    assert is_bimera_denovo(dada_chim).tolist() == fx["is_bimera"]
    assert is_bimera_denovo(dada_chim, allow_one_off=True).tolist() == fx["is_bimera_oneoff"]
    removed = remove_bimera_denovo(dada_chim)
    assert removed.index.tolist() == fx["removed"]["sequence"]
    assert removed.tolist() == fx["removed"]["abundance"]


# The fixture also holds R's flagSubseqs=TRUE calls ("subseqs"): R's isShiftDenovo never passes
# flagSubseqs on (chimeras.R), so they equal the default calls. pydada2 currently honours the
# flag (REVIEW.md §1.6b); that case is added once the behaviour is settled.
SHIFT = {"default": {}, "ov100": {"min_overlap": 100}}


@pytest.mark.parametrize("case", sorted(SHIFT))
def test_is_shift_denovo(case):
    d = CH["shift"]
    unq = pd.Series(d["uniques"]["abundance"], index=d["uniques"]["sequence"])
    assert is_shift_denovo(unq, **SHIFT[case]).tolist() == d["calls"][case]
