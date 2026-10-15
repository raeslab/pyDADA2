"""Phase 3: dereplication, sequence utilities and quality profiles vs R."""

import numpy as np
import pytest

import pydada2 as dada2
from pydada2.io import quality_profile

from rparity import EXTDATA, assert_close_r, derep_fixture, load, scalar

pytestmark = pytest.mark.parity

DEREP_CASES = {
    "raw_sam1F": ("sam1F.fastq.gz", {}),
    "raw_samPB": ("samPB.fastq.gz", {}),
    "raw_samPB_n100": ("samPB.fastq.gz", {"n": 100}),
}


@pytest.mark.parametrize("case", sorted(DEREP_CASES))
def test_derep_raw_files(case):
    fname, kwargs = DEREP_CASES[case]
    drp = dada2.derep_fastq(EXTDATA / fname, **kwargs)
    seqs, abund, quals, m = derep_fixture(case)
    assert drp.sequences == seqs
    np.testing.assert_array_equal(drp.uniques.to_numpy(), abund)
    np.testing.assert_array_equal(drp.map, m)
    assert drp.quals.shape == quals.shape
    np.testing.assert_allclose(drp.quals, quals, rtol=1e-12, equal_nan=True)


def test_derep_filtered_and_chunked(tmp_path):
    # Reproduce the filtered file of fixture F1 (sam1F), then dereplicate with n=default and 200.
    out = tmp_path / "sam1F.fastq"
    dada2.filter_and_trim(EXTDATA / "sam1F.fastq.gz", out, trunc_len=240, max_n=0, max_ee=2,
                          trunc_q=2, rm_phix=True, compress=False)  # fmt: skip
    for case, kwargs in (("filtF_sam1", {}), ("filtF_sam1_n200", {"n": 200})):
        drp = dada2.derep_fastq(out, **kwargs)
        seqs, abund, quals, m = derep_fixture(case)
        assert drp.sequences == seqs, case
        np.testing.assert_array_equal(drp.uniques.to_numpy(), abund)
        np.testing.assert_array_equal(drp.map, m)
        np.testing.assert_allclose(drp.quals, quals, rtol=1e-12, equal_nan=True)


def test_derep_list_and_repr():
    drps = dada2.derep_fastq([EXTDATA / "sam1F.fastq.gz", EXTDATA / "sam2F.fastq.gz"])
    assert len(drps) == 2 and drps[0].n_reads == 1500
    assert "1500 reads" in repr(drps[0])


def test_seq_complexity_parity():
    d = load("sequences")["complexity"]
    with pytest.warns(UserWarning, match="A/C/G/T"):
        assert_close_r(dada2.seq_complexity(d["seqs"]), d["k2"], rtol=1e-12)
    with pytest.warns(UserWarning):
        assert_close_r(dada2.seq_complexity(d["seqs"], kmer_size=3), d["k3"], rtol=1e-12)
    with pytest.warns(UserWarning):
        assert_close_r(dada2.seq_complexity(d["seqs"], window=50, by=5), d["win50"], rtol=1e-12)


def test_is_phix_parity():
    d = load("sequences")["is_phix"]
    assert dada2.is_phix(d["seqs"]).tolist() == d["default"]
    assert dada2.is_phix(d["seqs"], word_size=12, min_matches=1).tolist() == d["ws12_mm1"]


def test_reverse_complement_parity():
    d = load("sequences")["rc"]
    assert dada2.reverse_complement(d["seqs"]) == d["result"]


@pytest.mark.parametrize("fname", ["sam1F.fastq.gz", "samPB.fastq.gz"])
def test_quality_profile_parity(fname):
    d = load("quality_profile")[fname]
    df = quality_profile(EXTDATA / fname)
    pc = d["per_cycle"]
    exp = sorted(zip(pc["Cycle"], pc["Score"], pc["Count"], strict=True))
    got = sorted(zip(df["cycle"], df["score"], df["count"], strict=True))
    assert got == exp
    assert df.attrs["read_count"] == scalar(d["read_count"])


def test_uniques_to_fasta_parity(tmp_path):
    misc = load("misc")
    fx = load("dada")["default_F"]["sam1"]
    unqs = dict(zip(fx["denoised"]["sequence"], fx["denoised"]["abundance"], strict=True))
    out = tmp_path / "u.fasta"
    dada2.uniques_to_fasta(unqs, out)
    assert out.read_text().splitlines() == misc["uniques_to_fasta"]
