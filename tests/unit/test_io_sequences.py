import gzip

import numpy as np
import pandas as pd
import pytest

import pydada2 as dada2
from pydada2 import datasets
from pydada2.io import (
    format_fastq,
    list_fastq,
    quality_offset,
    quality_profile,
    read_fasta,
    read_fastq_chunks,
    write_fasta,
)


def _write_fq(path, recs):
    path.write_text("".join(f"@{i}\n{s}\n+\n{q}\n" for i, s, q in recs))
    return path


def test_datasets():
    assert "samPB.fastq.gz" in datasets.list_examples()
    with pytest.raises(FileNotFoundError, match="No example file"):
        datasets.example_path("nope.fq")
    for f in (datasets.err_balanced_f, datasets.err_balanced_r, datasets.tperr1):
        a = f()
        assert a.shape == (16, 41) and not a.flags.writeable


def test_quality_offset():
    assert quality_offset(["IIII"]) == 33
    assert quality_offset(["hhhh", "]]]]"]) == 64  # min > ':' and max > 'J'
    assert quality_offset(["IIII"], "SFastqQuality") == 64
    assert quality_offset(["hhhh"], "FastqQuality") == 33
    with pytest.raises(ValueError, match="quality_type"):
        quality_offset(["II"], "bogus")  # type: ignore[arg-type]


def test_read_fastq_chunks(tmp_path):
    f = _write_fq(tmp_path / "a.fq", [(f"r{i} x", "acgt", "IIII") for i in range(5)])
    chunks = list(read_fastq_chunks(f, n=2))
    assert [len(c) for c in chunks] == [2, 2, 1]
    assert chunks[0].seqs == ["ACGT", "ACGT"]  # upper-cased like ShortRead
    assert chunks[0].ids[0] == "r0 x"
    with pytest.raises(ValueError, match="positive"):
        list(read_fastq_chunks(f, n=0))


def test_format_fastq():
    assert format_fastq(["a"], ["AC"], ["II"]) == b"@a\nAC\n+\nII\n"


def test_list_fastq(tmp_path):
    for name in ("b.fastq.gz", "a.fastq", "c.txt"):
        (tmp_path / name).write_text("")
    assert [p.name for p in list_fastq(tmp_path)] == ["a.fastq", "b.fastq.gz"]
    assert [p.name for p in list_fastq(tmp_path, r"\.txt$")] == ["c.txt"]
    with pytest.raises(NotADirectoryError):
        list_fastq(tmp_path / "a.fastq")


def test_fasta_roundtrip(tmp_path):
    out = tmp_path / "x.fa.gz"
    write_fasta([("s1", "ACGT" * 5), ("s2", "GG")], out, width=8)
    assert read_fasta(out) == [("s1", "ACGT" * 5), ("s2", "GG")]
    with gzip.open(out, "rt") as fh:
        assert fh.read().splitlines()[1] == "ACGTACGT"
    write_fasta([("s3", "TT")], out, mode="a")
    assert [r[0] for r in read_fasta(out)] == ["s1", "s2", "s3"]
    with pytest.raises(ValueError, match="mode"):
        write_fasta([], out, mode="x")


def test_quality_profile_empty(tmp_path):
    f = tmp_path / "e.fq"
    f.write_text("")
    df = quality_profile(f)
    assert df.empty and df.attrs["read_count"] == 0


def test_reverse_complement():
    assert dada2.reverse_complement("") == ""
    assert dada2.reverse_complement(["ACGT-"]) == ["-ACGT"]


def test_get_uniques_and_sequences(tmp_path):
    s = dada2.get_uniques(pd.Series({"AC": 2, "GG": 1}))
    assert s.to_dict() == {"AC": 2, "GG": 1} and s.dtype == np.int64
    df = pd.DataFrame({"sequence": ["AC", "AC", "GT"], "abundance": [1, 2, 5]})
    assert dada2.get_uniques(df).to_dict() == {"AC": 3, "GT": 5}
    assert dada2.get_uniques(df, collapse=False).tolist() == [1, 2, 5]
    tab = pd.DataFrame([[1, 2], [3, 4]], columns=["AA", "CC"], index=["s1", "s2"])
    assert dada2.get_uniques(tab).to_dict() == {"AA": 4, "CC": 6}
    with pytest.raises(TypeError, match="Unrecognized"):
        dada2.get_uniques(3.5)
    drp = dada2.derep_fastq(datasets.example_path("sam1F.fastq.gz"))
    assert dada2.get_uniques(drp).equals(drp.uniques)
    assert dada2.get_sequences(drp)[:2] == drp.sequences[:2]
    assert dada2.get_uniques(datasets.example_path("sam1F.fastq.gz")).sum() == 1500
    fa = tmp_path / "x.fa"
    write_fasta([("a", "acg"), ("b", "TTA")], fa)
    assert dada2.get_sequences(fa) == ["ACG", "TTA"]
    assert dada2.get_sequences("acg") == ["ACG"]
    assert dada2.get_sequences(np.array(["AC", "AC"]), collapse=True) == ["AC"]


def test_uniques_to_fasta_ids(tmp_path):
    out = tmp_path / "u.fa"
    dada2.uniques_to_fasta({"AC": 2}, out, ids=["x"])
    assert out.read_text() == ">x\nAC\n"
    with pytest.raises(ValueError, match="one entry"):
        dada2.uniques_to_fasta({"AC": 2}, out, ids=["x", "y"])


def test_seq_complexity_errors():
    with pytest.raises(ValueError, match="larger than the kmer_size"):
        dada2.seq_complexity(["ACGTACGT"], kmer_size=4, window=3)
    assert np.isnan(dada2.seq_complexity(["A"])[0])
