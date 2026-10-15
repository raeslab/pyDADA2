"""Phase 6: mergePairs, makeSequenceTable, mergeSequenceTables, collapseNoMismatch vs R."""

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import err_balanced_f, err_balanced_r
from pydada2.denoise import DADA
from pydada2.merge import merge_pairs
from pydada2.seqtab import collapse_no_mismatch, make_sequence_table, merge_sequence_tables

from rparity import assert_close_r, load, matrix

pytestmark = pytest.mark.parity


@pytest.fixture(scope="module")
def denoised(dereps):
    ddF = DADA(err=err_balanced_f()).fit(dereps["dF_P"])
    ddR = DADA(err=err_balanced_r()).fit(dereps["dR"])
    return ddF, ddR, list(dereps["dF_P"].values()), list(dereps["dR"].values())


MERGE_CASES = {
    "default": {},
    "rejects": {"return_rejects": True},
    "mismatch1": {"max_mismatch": 1, "return_rejects": True},
    "min_overlap20": {"min_overlap": 20},
    "concat": {"just_concatenate": True},
    "trim_overhang": {"trim_overhang": True},
    "propagate": {"propagate_col": ["n0", "birth_fold"]},
}


def _seq_col(values):
    return [v[0] if isinstance(v, list) else v for v in values]


@pytest.mark.parametrize("case", sorted(MERGE_CASES))
def test_merge_pairs_parity(case, denoised):
    ddF, ddR, dF, dR = denoised
    got = merge_pairs(ddF, dF, ddR, dR, **MERGE_CASES[case])
    exp = load("merge_pairs")[case]
    assert len(got) == len(exp)
    for g, e in zip(got, exp.values(), strict=True):
        assert g["sequence"].tolist() == _seq_col(e["sequence"])
        assert g["abundance"].tolist() == e["abundance"]
        assert (g["forward"] + 1).tolist() == e["forward"]
        assert (g["reverse"] + 1).tolist() == e["reverse"]
        for col in ("nmatch", "nmismatch", "nindel", "accept"):
            assert g[col].tolist() == e[col], col
        assert [None if pd.isna(v) else int(v) for v in g["prefer"]] == e["prefer"]
        for col in ("F.n0", "R.n0"):
            if col in e:
                assert g[col].tolist() == e[col]
        for col in ("F.birth_fold", "R.birth_fold"):
            if col in e:
                assert_close_r(g[col].to_numpy(dtype=float), e[col], rtol=1e-10)


def _table(obj):
    return pd.DataFrame(
        matrix(obj, dtype=np.int64), index=obj.get("rownames"), columns=obj["colnames"]
    )


def _assert_table(got, obj):
    exp = _table(obj)
    assert got.columns.tolist() == exp.columns.tolist()
    if obj.get("rownames") is not None:
        assert [str(i) for i in got.index] == list(exp.index)
    np.testing.assert_array_equal(got.to_numpy(), exp.to_numpy())


def test_make_sequence_table_parity(denoised):
    ddF, ddR, dF, dR = denoised
    st = load("seqtab")
    mergers = merge_pairs(ddF, dF, ddR, dR)
    named = dict(zip(["sam1", "sam2"], mergers, strict=True))
    _assert_table(make_sequence_table(named), st["seqtab"])
    _assert_table(make_sequence_table(named, order_by="nsamples"), st["seqtab_nsamples"])
    _assert_table(make_sequence_table(named, order_by=None), st["seqtab_unordered"])
    seqtab_f = make_sequence_table(ddF)
    _assert_table(seqtab_f, st["seqtabF"])
    assert list(seqtab_f.index) == ["sam1", "sam2"]  # names from the fitted DADA


def test_collapse_parity():
    st = load("seqtab")
    st_var = _table(st["st_var"])
    _assert_table(collapse_no_mismatch(st_var), st["collapse"])
    _assert_table(collapse_no_mismatch(st_var, min_overlap=100), st["collapse_ov100"])
    _assert_table(collapse_no_mismatch(st_var, vec=False), st["collapse_novec"])
    _assert_table(
        collapse_no_mismatch(_table(st["st_dup"]), identical_only=True), st["collapse_identical"]
    )
    for n in (2, 4):
        _assert_table(collapse_no_mismatch(st_var, n_jobs=n), st["collapse"])


def test_merge_sequence_tables_parity():
    st = load("seqtab")
    a = _table(st["seqtab"])
    _assert_table(merge_sequence_tables(a, _table(st["stB"])), st["merge_AB"])
    _assert_table(merge_sequence_tables(a, _table(st["stC"]), repeats="sum"), st["merge_AC_sum"])
    _assert_table(
        merge_sequence_tables(a, _table(st["stB"]), order_by=None), st["merge_AB_unordered"]
    )
    _assert_table(merge_sequence_tables(a, _table(st["stRC"]), try_rc=True), st["merge_ARC_tryrc"])
    with pytest.raises(ValueError, match="Duplicated sample names"):
        merge_sequence_tables(a, _table(st["stC"]))


def test_get_uniques_of_mergers(denoised):
    ddF, ddR, dF, dR = denoised
    m = merge_pairs(ddF, dF, ddR, dR)[0]
    u = dada2.get_uniques(m)
    assert u.sum() == m["abundance"].sum()
