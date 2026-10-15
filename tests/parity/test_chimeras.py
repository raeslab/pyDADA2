"""Phase 7: chimera detection/removal and shift detection vs R."""

import numpy as np
import pandas as pd
import pytest

from pydada2.chimeras import (
    BimeraRemover,
    is_bimera,
    is_bimera_denovo,
    is_bimera_denovo_table,
    is_shift_denovo,
    remove_bimera_denovo,
)
from pydada2.datasets import err_balanced_f
from pydada2.denoise import DADA

from rparity import load, matrix

pytestmark = pytest.mark.parity


@pytest.fixture(scope="module")
def ch():
    return load("chimeras")


def _table(obj):
    return pd.DataFrame(
        matrix(obj, dtype=np.int64), index=obj.get("rownames"), columns=obj["colnames"]
    )


def test_is_bimera(ch):
    d = ch["is_bimera"]
    for k, s in enumerate(d["seqs"]):
        pars = [p for p in d["parents"] if p != s]
        assert is_bimera(s, pars) == d["default"][k]
        assert is_bimera(s, pars, allow_one_off=True) == d["one_off"][k]
        assert is_bimera(s, pars, max_shift=2) == d["shift2"][k]


def test_is_bimera_denovo_table(ch):
    st = _table(ch["st_chim"])
    assert st.columns.tolist()[-len(ch["chims"]) :] == ch["chims"]
    assert is_bimera_denovo_table(st).tolist() == ch["is_bimera_denovo_table"]
    got = is_bimera_denovo_table(st, allow_one_off=True, min_sample_fraction=0.5)
    assert got.tolist() == ch["is_bimera_denovo_table_oneoff"]
    assert is_bimera_denovo_table(st, n_jobs=4).tolist() == ch["is_bimera_denovo_table"]


def test_is_bimera_denovo(ch):
    st = _table(ch["st_chim"])
    unq = st.sum(axis=0)
    assert is_bimera_denovo(unq).tolist() == ch["is_bimera_denovo"]
    got = is_bimera_denovo(
        unq, allow_one_off=True, min_fold_parent_over_abundance=1, min_parent_abundance=2
    )
    assert got.tolist() == ch["is_bimera_denovo_oneoff"]


def test_is_bimera_denovo_dada(ch, dereps):
    dd = DADA(err=err_balanced_f()).fit(dereps["dF_P"])
    assert is_bimera_denovo(dd[0]).tolist() == ch["is_bimera_denovo_dada"]


def test_remove_bimera_denovo(ch):
    st = _table(ch["st_chim"])
    assert remove_bimera_denovo(st).columns.tolist() == ch["remove_consensus"]
    assert remove_bimera_denovo(st, method="pooled").columns.tolist() == ch["remove_pooled"]
    ps = remove_bimera_denovo(st, method="per-sample")
    exp = _table(ch["remove_per_sample"])
    assert ps.columns.tolist() == exp.columns.tolist()
    np.testing.assert_array_equal(ps.to_numpy(), exp.to_numpy())
    br = BimeraRemover(method="per-sample").fit(st)
    assert br.is_bimera_.shape == st.shape


def test_is_shift_denovo(ch):
    d = ch["shift"]
    unq = pd.Series(d["uniques"]["abundance"], index=d["uniques"]["sequence"])
    assert is_shift_denovo(unq).tolist() == d["default"]
    assert is_shift_denovo(unq, min_overlap=100).tolist() == d["ov100"]
    assert is_shift_denovo(unq, flag_subseqs=True).tolist() == d["subseqs"]
