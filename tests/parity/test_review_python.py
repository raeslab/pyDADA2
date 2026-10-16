"""Parity of the review fixes with R (fixtures: review_python, loess_degenerate)."""

import hashlib
import warnings

import dnaio
import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import example_path, tperr1

from rparity import INPUTS, colnames, load, matrix, rownames, scalar

pytestmark = pytest.mark.parity


@pytest.fixture(scope="module")
def rv():
    return load("review_python")


# ------------------------------------------- 1.4 duplicate sequences / tables ---
def _table(obj):
    return pd.DataFrame(matrix(obj, dtype=np.int64), index=rownames(obj), columns=colnames(obj))


def test_make_sequence_table_duplicate_mergers(rv):
    fx = rv["dup_mergers"]
    mergers = {
        name: pd.DataFrame({"sequence": m["sequence"], "abundance": m["abundance"]})
        for name, m in fx["mergers"].items()
    }
    assert any(m["sequence"].duplicated().any() for m in mergers.values())
    for name, m in mergers.items():
        u = dada2.get_uniques(m)
        assert u.index.tolist() == fx["uniques"][name]["sequence"]
        assert u.tolist() == fx["uniques"][name]["abundance"]
    for key, order_by in (
        ("seqtab", "abundance"),
        ("seqtab_nsamples", "nsamples"),
        ("seqtab_unordered", None),
    ):
        got = dada2.make_sequence_table(mergers, order_by=order_by)
        exp = _table(fx[key])
        assert got.columns.tolist() == exp.columns.tolist(), key
        assert got.index.tolist() == exp.index.tolist()
        np.testing.assert_array_equal(got.to_numpy(), exp.to_numpy())


def test_make_sequence_table_tied_duplicates(rv):
    """Duplicates merged in sorted order (R's tapply) decide the order of tied columns."""
    fx = rv["dup_synthetic"]
    mergers = {
        name: pd.DataFrame({"sequence": m["sequence"], "abundance": m["abundance"]})
        for name, m in fx["mergers"].items()
    }
    got = dada2.make_sequence_table(mergers)
    exp = _table(fx["seqtab"])
    assert got.columns.tolist() == exp.columns.tolist()
    np.testing.assert_array_equal(got.to_numpy(), exp.to_numpy())


def test_get_sequences_collapse_and_taxonomy_table(rv):
    assert (
        dada2.get_sequences(["acgt", "ACGT", "GG", "acgt"], collapse=True)
        == (rv["get_sequences_collapse"])
    )
    tax = pd.DataFrame(
        {"Kingdom": ["Bacteria", "Bacteria"], "Phylum": ["Firmicutes", None]},
        index=["ACGTAC", "GGTTAA"],
    )
    assert dada2.get_sequences(tax) == rv["get_sequences_taxtab"]


# ---------------------------------------------------------------- 1.5 priors ---
def _as_dict(denoised):
    return dict(zip(denoised.index.tolist(), denoised.tolist(), strict=True))


@pytest.mark.parametrize("form", ["string", "fasta", "derep_sam2"])
def test_dada_priors_forms(rv, form):
    fx = rv["priors"]
    d1 = dada2.derep_fastq(example_path("sam1F.fastq.gz"))
    priors = {
        "string": fx["prior"][0],
        "fasta": str(INPUTS / "priors.fa"),
        "derep_sam2": dada2.derep_fastq(example_path("sam2F.fastq.gz")),
    }[form]
    params = {"omega_a": 1e-60}
    got = dada2.DADA(err=tperr1(), priors=priors, params=params).fit(d1).results_[0].denoised
    exp = dict(zip(fx[form]["sequence"], fx[form]["abundance"], strict=True))
    assert _as_dict(got) == exp
    none = dada2.DADA(err=tperr1(), params=params).fit(d1).results_[0].denoised
    assert _as_dict(none) == dict(zip(fx["none"]["sequence"], fx["none"]["abundance"], strict=True))
    assert len(exp) != len(fx["none"]["sequence"])  # the prior matters


# ------------------------------------------------------------- 1.6 b, c, f ---
def test_is_shift_denovo_ignores_flag_subseqs(rv):
    d = rv["shift"]
    u = pd.Series(d["uniques"]["abundance"], index=d["uniques"]["sequence"])
    assert dada2.is_shift_denovo(u).tolist() == d["default"]
    with pytest.warns(UserWarning, match="no effect"):
        assert dada2.is_shift_denovo(u, flag_subseqs=True).tolist() == d["subseqs"]


@pytest.mark.parametrize(("key", "allow"), [("single", False), ("multiple", True)])
def test_assign_species_tokenisation(rv, key, allow):
    fx = rv["species_tokens"]
    got = dada2.assign_species(fx["queries"], INPUTS / "species_tokens.fa", allow_multiple=allow)
    assert got["Genus"].tolist() == fx[key]["V1"]
    assert got["Species"].tolist() == fx[key]["V2"]


@pytest.mark.parametrize(
    ("key", "kwargs"),
    [
        ("maxN_inf", {"max_n": np.inf, "min_len": -np.inf}),
        ("maxN_1", {"max_n": 1}),
        ("minQ_inf", {"min_q": np.inf, "min_len": np.inf}),
    ],
)
def test_filter_infinite_thresholds(rv, tmp_path, key, kwargs):
    fx = rv["filter_inf"][key]
    out = tmp_path / "out.fastq"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stats = dada2.filter_and_trim(
            INPUTS / "with_n.fastq.gz", out, trunc_q=0, rm_phix=False, compress=False, **kwargs
        )
    assert stats["reads_in"].tolist() == fx["reads_in"]
    assert stats["reads_out"].tolist() == fx["reads_out"]
    exp = fx["outputs"][0]
    assert out.exists() == scalar(exp["exists"])
    if out.exists():
        with dnaio.open(out) as fh:
            recs = list(fh)
        assert [r.name for r in recs] == exp["ids"]
        txt = "".join(f"{r.name}\n{r.sequence}\n{r.qualities}\n" for r in recs)
        assert hashlib.md5(txt.encode()).hexdigest() == scalar(exp["md5"])
