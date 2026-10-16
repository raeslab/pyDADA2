"""n_jobs determinism: every kernel that uses parallel_for gives identical results
for 1 and 4 threads. (dada_uniques, DADA pseudo-pooling, table_bimera, collapse_no_mismatch,
filter_and_trim and taxonomy are covered in the parity tests.)"""

import warnings

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2.datasets import err_balanced_f, err_balanced_r, example_path

from rparity import INPUTS, load, matrix

F27 = "AGRGTTYGATYMTGGCTCAG"
R1492_RC = dada2.reverse_complement("RGYTACCTTGTTACGACTT")


def _st_chim():
    for name in ("chimera_coverage", "chimeras"):
        key = "table" if name == "chimera_coverage" else "st_chim"
        st = load(name)[key]
        yield pd.DataFrame(matrix(st, dtype=np.int64), index=st["rownames"], columns=st["colnames"])


def _cov_seqs():
    import dnaio

    with dnaio.open(INPUTS / "cov_F.fastq.gz") as fh:
        return [r.sequence.upper() for r in fh]


def test_seq_complexity_n_jobs():
    seqs = _cov_seqs() + load("sequences")["complexity"]["seqs"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # non-ACGT letters
        for kw in ({}, {"kmer_size": 3}, {"window": 50, "by": 5}):
            a = dada2.seq_complexity(seqs, **kw)
            b = dada2.seq_complexity(seqs, n_jobs=4, **kw)
            np.testing.assert_array_equal(a, b)


@pytest.mark.parametrize("allow_indels", [False, True])
def test_remove_primers_n_jobs(tmp_path, allow_indels):
    src = example_path("samPBprimers.fastq.gz")
    kw = {"primer_fwd": F27, "primer_rev": R1492_RC, "allow_indels": allow_indels}
    a = dada2.remove_primers(src, tmp_path / "a.fq", compress=False, **kw)
    b = dada2.remove_primers(src, tmp_path / "b.fq", compress=False, n_jobs=4, **kw)
    assert a.equals(b)
    assert (tmp_path / "a.fq").read_bytes() == (tmp_path / "b.fq").read_bytes()


def test_filter_kernel_n_jobs(tmp_path):
    # one file -> all threads go to the K1 kernel; inputs where every filter bites
    kw = {"trunc_q": 2, "max_n": 1, "rm_phix": True, "rm_lowcomplex": 8, "max_len": 255}
    pairs = (INPUTS / "cov_F.fastq.gz", INPUTS / "cov_R.fastq.gz")
    a = dada2.filter_and_trim(pairs[0], tmp_path / "a.fq", pairs[1], tmp_path / "ar.fq",
                              compress=False, **kw)  # fmt: skip
    b = dada2.filter_and_trim(pairs[0], tmp_path / "b.fq", pairs[1], tmp_path / "br.fq",
                              compress=False, n_jobs=4, **kw)  # fmt: skip
    assert a.equals(b)
    for x, y in (("a.fq", "b.fq"), ("ar.fq", "br.fq")):
        assert (tmp_path / x).read_bytes() == (tmp_path / y).read_bytes()


@pytest.mark.parametrize("allow_multiple", [False, True])
def test_assign_species_n_jobs(allow_multiple):
    tx = load("taxonomy")
    seqs = tx["queries_species"]
    ref = example_path("example_species_assignment.fa.gz")
    a = dada2.assign_species(seqs, ref, allow_multiple=allow_multiple, try_rc=True)
    b = dada2.assign_species(seqs, ref, allow_multiple=allow_multiple, try_rc=True, n_jobs=4)
    assert a.equals(b)


@pytest.mark.parametrize("st", list(_st_chim()), ids=["coverage", "chimeras"])
def test_chimera_functions_n_jobs(st):
    unq = st.sum(axis=0)
    for kw in ({}, {"allow_one_off": True}):
        a = dada2.is_bimera_denovo(unq, **kw)
        assert a.equals(dada2.is_bimera_denovo(unq, n_jobs=4, **kw))
    for method in ("consensus", "pooled", "per-sample"):
        a = dada2.remove_bimera_denovo(st, method=method)
        assert a.equals(dada2.remove_bimera_denovo(st, method=method, n_jobs=4))
    for kw in ({}, {"min_overlap": 100}, {"flag_subseqs": True}):
        a = dada2.is_shift_denovo(unq, **kw)
        assert a.equals(dada2.is_shift_denovo(unq, n_jobs=4, **kw))


def test_shift_and_dada_chimera_n_jobs():
    d = load("chimera_coverage")["shift"]
    unq = pd.Series(d["uniques"]["abundance"], index=d["uniques"]["sequence"])
    assert dada2.is_shift_denovo(unq).equals(dada2.is_shift_denovo(unq, n_jobs=4))
    drp = dada2.derep_fastq(INPUTS / "chimeric_sam1F.fastq.gz")
    dd = dada2.DADA(err=err_balanced_f()).fit(drp)[0]
    assert dada2.is_bimera_denovo(dd).equals(dada2.is_bimera_denovo(dd, n_jobs=4))


def test_merge_pairs_and_pooled_dada_n_jobs():
    drp_f = dada2.derep_fastq([example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")])
    drp_r = dada2.derep_fastq([example_path("sam1R.fastq.gz"), example_path("sam2R.fastq.gz")])
    for pool in (False, True):
        a = dada2.DADA(err=err_balanced_f(), pool=pool).fit(drp_f)
        b = dada2.DADA(err=err_balanced_f(), pool=pool, n_jobs=4).fit(drp_f)
        for ra, rb in zip(a.results_, b.results_, strict=True):
            pd.testing.assert_frame_equal(ra.clusters, rb.clusters)
            np.testing.assert_array_equal(ra.map, rb.map)
    dd_f = dada2.DADA(err=err_balanced_f()).fit(drp_f)
    dd_r = dada2.DADA(err=err_balanced_r()).fit(drp_r)
    for kw in ({}, {"max_mismatch": 1, "return_rejects": True}, {"just_concatenate": True}):
        a = dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r, **kw)
        b = dada2.merge_pairs(dd_f, drp_f, dd_r, drp_r, n_jobs=4, **kw)
        for x, y in zip(a, b, strict=True):
            pd.testing.assert_frame_equal(x, y)
