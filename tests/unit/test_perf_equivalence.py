"""Faster code paths give the same results (REVIEW.md 5.3)."""

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import example_path


def test_derep_threads_same_result():
    files = [example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")] * 2
    a = dada2.derep_fastq(files)
    b = dada2.derep_fastq(files, n_jobs=3)
    for x, y in zip(a, b, strict=True):
        assert x.uniques.equals(y.uniques)
        np.testing.assert_array_equal(x.quals, y.quals)
        np.testing.assert_array_equal(x.map, y.map)


def test_per_sample_bimeras_skip_absent_sequences():
    rng = np.random.default_rng(1)
    base = "".join(rng.choice(list("ACGT"), 120))
    seqs = set()
    while len(seqs) < 40:
        s = list(base)
        for p in rng.choice(120, 6, replace=False):
            s[p] = rng.choice(list("ACGT"))
        seqs.add("".join(s))
    seqs = sorted(seqs)
    a, b = seqs[0], seqs[1]
    chim = a[:60] + b[60:]
    cols = [*seqs, chim]
    mat = np.zeros((6, len(cols)), dtype=np.int64)
    for i in range(6):
        idx = rng.choice(len(seqs), 15, replace=False)
        mat[i, idx] = rng.integers(1, 500, 15)
        mat[i, :2] = [900, 800]
        mat[i, -1] = 0 if i % 2 else 7
    tab = pd.DataFrame(mat, columns=cols)
    br = dada2.BimeraRemover(method="per-sample").fit(tab)
    flags = br.is_bimera_.to_numpy()
    assert not flags[mat == 0].any()
    for i in range(6):  # the full-row test (zeros included, as R does) agrees on present cells
        full = dada2.is_bimera_denovo(
            pd.Series(mat[i], index=cols), min_fold_parent_over_abundance=2, min_parent_abundance=8
        ).to_numpy()
        np.testing.assert_array_equal(flags[i, mat[i] > 0], full[mat[i] > 0])
    assert flags[0, -1]
    out = br.transform(tab)
    masked = tab.mask(br.is_bimera_.to_numpy(), 0)
    pd.testing.assert_frame_equal(out, masked.loc[:, masked.sum(axis=0).to_numpy() != 0])


def test_filter_writer_thread_errors_and_chunks(tmp_path):
    with pytest.raises(FileNotFoundError):
        dada2.fastq_filter(
            example_path("sam1F.fastq.gz"), tmp_path / "missing_dir" / "o.fastq", compress=False
        )
    one = tmp_path / "one.fastq.gz"
    many = tmp_path / "many.fastq.gz"
    a = dada2.fastq_filter(example_path("sam1F.fastq.gz"), one, trunc_len=240, max_ee=2)
    b = dada2.fastq_filter(example_path("sam1F.fastq.gz"), many, trunc_len=240, max_ee=2, n=7)
    assert a == b
    assert dada2.derep_fastq(one).uniques.equals(dada2.derep_fastq(many).uniques)


def test_paired_filter_forward_error_wins(tmp_path):
    fn = (example_path("sam1F.fastq.gz"), example_path("sam1R.fastq.gz"))
    fo = (tmp_path / "f.fastq", tmp_path / "r.fastq")
    with pytest.raises(ValueError, match="trunc_q"):
        dada2.fastq_paired_filter(fn, fo, trunc_q=(500, 600), compress=False)
    out = dada2.fastq_paired_filter(fn, fo, trunc_len=(240, 160), compress=False, n_threads=2)
    assert out.reads_out > 0
