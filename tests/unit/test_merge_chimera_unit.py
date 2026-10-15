import logging

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2 import DADA, BimeraRemover
from pydada2.datasets import err_balanced_f, err_balanced_r, example_path

A = "ACGTACGTAC" * 3
B = "TTGGCCAATT" * 3
CHIM = A[:15] + B[15:]


@pytest.fixture(scope="module")
def pair():
    drp_f = dada2.derep_fastq(example_path("sam1F.fastq.gz"))
    drp_r = dada2.derep_fastq(example_path("sam1R.fastq.gz"))
    ddf = DADA(err=err_balanced_f()).fit(drp_f)
    ddr = DADA(err=err_balanced_r()).fit(drp_r)
    return ddf, drp_f, ddr, drp_r


def test_merge_pairs_single_and_validation(pair, caplog):
    ddf, drp_f, ddr, drp_r = pair
    with caplog.at_level(logging.INFO, logger="pydada2"):
        m = dada2.merge_pairs(ddf, drp_f, ddr, drp_r, verbose=1, n_jobs=2)
    assert isinstance(m, pd.DataFrame) and m["accept"].all()
    assert any("successfully merged" in msg for msg in caplog.messages)
    m2 = dada2.merge_pairs(ddf[0], drp_f, ddr[0], drp_r)
    assert m2.equals(m)
    with pytest.raises(ValueError, match="same length"):
        dada2.merge_pairs([ddf[0]], [drp_f, drp_f], [ddr[0]], [drp_r])
    with pytest.raises(ValueError, match="Non-corresponding"):
        dada2.merge_pairs(ddf, drp_r, ddr, drp_r)
    with pytest.raises(TypeError, match="dada_f"):
        dada2.merge_pairs([3], [drp_f], [ddr[0]], [drp_r])
    with pytest.raises(TypeError, match="derep_f"):
        dada2.merge_pairs(ddf, [3], ddr, [drp_r])


def test_merge_pairs_no_pairs(pair, caplog):
    ddf, drp_f, ddr, drp_r = pair
    res = ddf[0]
    import dataclasses

    empty = dataclasses.replace(res, map=np.full_like(res.map, -1))
    with caplog.at_level(logging.INFO, logger="pydada2"):
        out = dada2.merge_pairs(empty, drp_f, ddr[0], drp_r, verbose=1)
    assert out.empty and "accept" in out.columns


def test_make_sequence_table_inputs():
    t = dada2.make_sequence_table([{"AC": 1, "GG": 5}, {"AC": 7}])
    assert t.columns.tolist() == ["AC", "GG"] and t.index.tolist() == [0, 1]
    t1 = dada2.make_sequence_table(pd.Series({"AC": 2}))
    assert t1.shape == (1, 1)
    with pytest.raises(ValueError, match="order_by"):
        dada2.make_sequence_table([{"AC": 1}], order_by="bogus")  # type: ignore[arg-type]


def test_merge_sequence_tables_errors():
    t = pd.DataFrame({"AC": [1]}, index=["s1"])
    with pytest.raises(ValueError, match="At least two"):
        dada2.merge_sequence_tables(t)
    with pytest.raises(TypeError, match="not a DataFrame"):
        dada2.merge_sequence_tables(t, 3)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="negative"):
        dada2.merge_sequence_tables(t, pd.DataFrame({"AC": [-1]}, index=["s2"]))
    with pytest.raises(ValueError, match="blank"):
        dada2.merge_sequence_tables(t, pd.DataFrame({"": [1]}, index=["s2"]))
    m = dada2.merge_sequence_tables([t, pd.DataFrame({"GG": [2]}, index=["s2"])])
    assert m.shape == (2, 2)


def test_collapse_duplicates_and_verbose(caplog):
    tab = pd.DataFrame([[1, 2, 3]], columns=["ACGT" * 6, "ACGT" * 6, "ACGT" * 6 + "A"])
    with caplog.at_level(logging.INFO, logger="pydada2"):
        out = dada2.collapse_no_mismatch(tab, verbose=1)
    assert out.shape == (1, 1) and out.iloc[0, 0] == 6


def test_bimera_remover_api():
    tab = pd.DataFrame({A: [100, 80], B: [90, 70], CHIM: [5, 4]}, index=["s1", "s2"])
    br = BimeraRemover()
    with pytest.raises(dada2.NotFittedError):
        br.transform(tab)
    assert br.fit_transform(tab).columns.tolist() == [A, B]
    assert BimeraRemover(method="pooled", min_parent_abundance=1).fit_transform(tab).shape[1] == 2
    with pytest.raises(ValueError, match="Valid values"):
        BimeraRemover(method="nope").fit(tab)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="DataFrame"):
        BimeraRemover().fit({"a": 1})  # type: ignore[arg-type]


def test_remove_bimera_denovo_types(pair):
    ddf, drp_f, ddr, drp_r = pair
    u = pd.Series({A: 100, B: 90, CHIM: 5})
    assert remove_out(dada2.remove_bimera_denovo(u, min_parent_abundance=1)) == [A, B]
    assert isinstance(dada2.remove_bimera_denovo(ddf[0]), pd.Series)
    assert isinstance(dada2.remove_bimera_denovo(drp_f), pd.Series)
    m = dada2.merge_pairs(ddf, drp_f, ddr, drp_r)
    assert isinstance(dada2.remove_bimera_denovo(m), pd.DataFrame)
    assert len(dada2.remove_bimera_denovo([u, u], min_parent_abundance=1)) == 2


def remove_out(s):
    return list(s.index)


def test_is_bimera_denovo_table_validation(caplog):
    with pytest.raises(ValueError, match="valid sequence table"):
        dada2.is_bimera_denovo_table({"a": 1})  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Duplicate"):
        dada2.is_bimera_denovo_table(pd.DataFrame([[1, 2]], columns=["AC", "AC"]))
    with pytest.raises(ValueError, match="integer"):
        dada2.is_bimera_denovo_table(pd.DataFrame([[1.5]], columns=["AC"]))
    tab = pd.DataFrame({A: [100], B: [90], CHIM: [5]})
    with caplog.at_level(logging.INFO, logger="pydada2"):
        assert dada2.is_bimera_denovo_table(tab, verbose=1).tolist() == [False, False, True]
        dada2.is_bimera_denovo(pd.Series({A: 100, CHIM: 1}), verbose=1)
    assert any("bimeras out of" in m for m in caplog.messages)


def test_is_bimera_denovo_duplicates(caplog):
    df = pd.DataFrame({"sequence": [A, A, B], "abundance": [5, 6, 7]})
    with caplog.at_level(logging.INFO, logger="pydada2"):
        out = dada2.is_bimera_denovo(df)
    assert len(out) == 3
