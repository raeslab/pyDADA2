"""The behaviour changes listed in DEVIATIONS.md ("no effect on valid inputs"), each pinned.

R would warn, error or read out of bounds on these inputs; pydada2 behaves as documented.
"""

import warnings

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2 import align
from pydada2.datasets import err_balanced_f, example_path
from pydada2.taxonomy import TaxonomyClassifier


def _fq(path, records):
    path.write_text("".join(f"@{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in records))
    return path


def test_unequal_alignment_strings_raise():
    # R: C_eval_pair / C_pair_consensus warn and return NULL
    with pytest.raises(ValueError, match="same length"):
        align.eval_pair("ACG-T", "ACGT")
    with pytest.raises(ValueError, match="same length"):
        align.pair_consensus("ACG-T", "ACGT")


def test_empty_alignment_strings_give_zero_counts():
    # R: reads out of bounds
    assert tuple(align.eval_pair("", "")) == (0, 0, 0)


def test_single_end_orient_short_reads_do_not_match(tmp_path):
    # R: narrow() errors and aborts the whole file; pydada2 treats short reads as not matching
    f = _fq(tmp_path / "in.fq", [("short", "TAC"), ("fwd", "TACGTACGTACGTACGTACGTACG")])
    stats = dada2.fastq_filter(f, tmp_path / "o.fq", orient_fwd="TACG", min_len=1,
                               rm_phix=False, compress=False)  # fmt: skip
    assert (stats.reads_in, stats.reads_out) == (2, 1)
    assert "fwd" in (tmp_path / "o.fq").read_text()


def test_rm_lowcomplex_removes_reads_without_acgt_kmers(tmp_path):
    # R: seqComplexity is NA for these reads and filterAndTrim indexes with NA
    f = _fq(tmp_path / "in.fq", [("n", "N" * 30), ("ok", "TACGGAGGATCCGAGCGTTATCCGGATTTA")])
    stats = dada2.fastq_filter(f, tmp_path / "o.fq", rm_lowcomplex=1, min_len=1, max_n=100,
                               rm_phix=False, compress=False)  # fmt: skip
    assert (stats.reads_in, stats.reads_out) == (2, 1)
    assert "@ok" in (tmp_path / "o.fq").read_text()


def test_remove_primers_empty_primer_raises(tmp_path):
    with pytest.raises(ValueError, match="Primer must not be empty"):
        dada2.remove_primers(example_path("samPBprimers.fastq.gz"), tmp_path / "o.fq", "")


def test_dada_rejects_sequences_at_or_above_seqlen():
    # C++ fixed buffers of SEQLEN (9999); upstream checked this too
    long = "ACGT" * 2500
    drp = dada2.Dereplication(
        uniques=pd.Series([5, 3], index=[long, long[:-1] + "A"]),
        quals=np.full((2, len(long)), 30.0),
        map=np.zeros(8, dtype=np.int64),
    )
    with pytest.raises(ValueError, match="maximum allowed string length"):
        dada2.DADA(err=err_balanced_f()).fit(drp)


@pytest.fixture(scope="module")
def classifier():
    return TaxonomyClassifier(min_boot=0, random_state=1).fit(
        example_path("example_train_set.fa.gz")
    )


def test_taxonomy_long_and_short_queries(classifier):
    # upstream: fixed-size stack arrays overflow above 9999 nt; queries shorter than the k-mer
    # size blew up the uniform-buffer size
    long = "".join(classifier_ref_seq() for _ in range(8))
    assert len(long) > 9999
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tax, boot = classifier.predict([long, "ACGT", long[:300]], return_bootstrap=True)
    assert tax.shape[0] == 3 and boot.shape[0] == 3
    assert tax.iloc[1].isna().all()  # too short to classify
    assert tax.iloc[0, 0] is not None and tax.iloc[2, 0] is not None


def classifier_ref_seq():
    import dnaio

    with dnaio.open(example_path("example_train_set.fa.gz")) as fh:
        return next(iter(fh)).sequence
