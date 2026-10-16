"""Unit tests for the review fixes of the Python layer."""

import gzip

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2._parallel import effective_n_jobs
from pydada2.datasets import example_path
from pydada2.io import read_fasta


# ----------------------------------------------------------------- n_jobs ---
def test_n_jobs_numpy_integers():
    assert effective_n_jobs(np.int64(2)) == 2
    assert effective_n_jobs(np.int32(-1)) >= 1
    with pytest.raises(ValueError):
        effective_n_jobs(np.True_)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        effective_n_jobs(2.0)  # type: ignore[arg-type]


# ------------------------------------------------------------------- derep ---
def test_derep_empty_file_raises_r_error(tmp_path):
    empty = tmp_path / "empty.fastq.gz"
    with gzip.open(empty, "wb"):
        pass
    with pytest.raises(ValueError, match="Only zero-length sequences detected"):
        dada2.derep_fastq(empty)


# ------------------------------------------------------------ ErrorLearner ---
def test_error_learner_mapping_and_param_conflicts():
    d1 = dada2.derep_fastq(example_path("sam1F.fastq.gz"))
    d2 = dada2.derep_fastq(example_path("sam2F.fastq.gz"))
    a = dada2.learn_errors({"s1": d1, "s2": d2}, max_consist=2)
    b = dada2.learn_errors([d1, d2], max_consist=2)
    np.testing.assert_array_equal(a.err_, b.err_)
    assert a.n_samples_ == 2
    with pytest.raises(ValueError, match="omega_c"):
        dada2.learn_errors([d1], params={"omega_c": 1e-40})
    with pytest.raises(ValueError, match="max_consist"):
        dada2.learn_errors([d1], params=dada2.DadaParams(max_consist=3))
    # the same value in both places is not a conflict
    c = dada2.learn_errors([d1, d2], max_consist=2, params={"omega_c": 0, "max_consist": 2})
    np.testing.assert_array_equal(c.err_, b.err_)


# ------------------------------------------------------------ get_sequences ---
def test_get_sequences_taxonomy_table_and_dedup_order():
    tax = dada2.assign_taxonomy(
        ["ACGT" * 20, "TTGA" * 20], example_path("example_train_set.fa.gz"), random_state=1
    )
    assert dada2.get_sequences(tax) == ["ACGT" * 20, "TTGA" * 20]
    assert dada2.get_sequences(["acgt", "ACGT"], collapse=True) == ["ACGT", "ACGT"]
    assert dada2.get_uniques({"TT": 1, "AA": 2}).index.tolist() == ["TT", "AA"]
    dup = pd.DataFrame({"sequence": ["TT", "AA", "TT"], "abundance": [1, 2, 1]})
    assert dada2.get_uniques(dup).to_dict() == {"AA": 2, "TT": 2}


def test_read_fasta_raw_ids(tmp_path):
    fa = tmp_path / "x.fa"
    fa.write_text(">a b  \nACGT\n>c\td\t\nGG\nTT\n")
    assert read_fasta(fa) == [("a b", "ACGT"), ("c\td", "GGTT")]
    assert read_fasta(fa, raw_ids=True) == [("a b  ", "ACGT"), ("c\td\t", "GGTT")]
    fq = tmp_path / "x.fq"
    fq.write_text("@r1 x \nACGT\n+\n>>>>\n")
    assert read_fasta(fq, raw_ids=True) == read_fasta(fq)


# --------------------------------------------------------------- filtering ---
def test_filter_infinite_only_for_thresholds(tmp_path):
    out = tmp_path / "o.fastq"
    with pytest.raises(ValueError, match="trim_left"):
        dada2.filter_and_trim(example_path("sam1F.fastq.gz"), out, trim_left=np.inf)
    with pytest.raises(ValueError, match="trunc_q"):
        dada2.filter_and_trim(example_path("sam1F.fastq.gz"), out, trunc_q=-np.inf)
