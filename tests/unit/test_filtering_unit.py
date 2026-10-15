import gzip

import pytest

import pydada2 as dada2
from pydada2.datasets import example_path


def _fq(path, recs):
    path.write_text("".join(f"@{i}\n{s}\n+\n{q}\n" for i, s, q in recs))
    return path


def test_fastq_filter_same_file(tmp_path):
    f = _fq(tmp_path / "a.fq", [("r1", "ACGT" * 10, "I" * 40)])
    with pytest.raises(ValueError, match="must be different"):
        dada2.fastq_filter(f, f)


def test_fastq_filter_none_pass_removes_output(tmp_path):
    f = _fq(tmp_path / "a.fq", [("r1", "ACGT" * 10, "I" * 40)])
    out = tmp_path / "o.fq.gz"
    stats = dada2.fastq_filter(f, out, min_len=100, rm_phix=False)
    assert stats == (1, 0) and not out.exists()


def test_fastq_filter_uncompressed_and_trimming(tmp_path):
    f = _fq(tmp_path / "a.fq", [("r1", "ACGTAC" * 10, "I" * 30 + "#" + "I" * 29)])
    out = tmp_path / "o.fq"
    stats = dada2.fastq_filter(f, out, rm_phix=False, trim_left=2, compress=False, verbose=1)
    assert stats == (1, 1)
    lines = out.read_text().splitlines()
    assert lines[1] == ("ACGTAC" * 10)[2:30]  # truncated before the '#' (Q2) base


def test_trunc_q_encoding_error(tmp_path):
    f = _fq(tmp_path / "a.fq", [("r1", "ACGT" * 10, "I" * 40)])
    with pytest.raises(ValueError, match="Encoding for this trunc_q"):
        dada2.fastq_filter(f, tmp_path / "o.fq", trunc_q=200)


def test_orient_requires_acgt(tmp_path):
    f = _fq(tmp_path / "a.fq", [("r1", "ACGT" * 10, "I" * 40)])
    with pytest.raises(ValueError, match="Non-ACGT"):
        dada2.fastq_filter(f, tmp_path / "o.fq", orient_fwd="ACGN")


def test_filter_and_trim_validation(tmp_path):
    src = example_path("sam1F.fastq.gz")
    with pytest.raises(ValueError, match="corresponding output"):
        dada2.filter_and_trim([src, src], [tmp_path / "a"])
    with pytest.raises(FileNotFoundError):
        dada2.filter_and_trim(tmp_path / "missing.fq", tmp_path / "a")
    with pytest.raises(ValueError, match="distinct"):
        dada2.filter_and_trim(
            [src, example_path("sam2F.fastq.gz")], [tmp_path / "a", tmp_path / "a"]
        )
    with pytest.raises(ValueError, match="scalars for single-end"):
        dada2.filter_and_trim(src, tmp_path / "a", trunc_len=(240, 200))
    with pytest.raises(ValueError, match="reverse reads are required"):
        dada2.filter_and_trim(src, tmp_path / "a", rev=src)
    with pytest.raises(ValueError, match="must correspond"):
        dada2.filter_and_trim(src, tmp_path / "a", [src, src], [tmp_path / "b", tmp_path / "c"])
    with pytest.raises(ValueError, match="integer"):
        dada2.filter_and_trim(src, tmp_path / "a", trunc_len=2.5)


def test_filter_and_trim_creates_dirs_and_warns(tmp_path):
    out = tmp_path / "new" / "dir" / "o.fq.gz"
    with pytest.warns(UserWarning, match="No reads passed"):
        df = dada2.filter_and_trim(example_path("sam1F.fastq.gz"), out, trunc_len=400)
    assert df["reads_out"].tolist() == [0] and out.parent.is_dir()


def test_paired_mismatched_counts(tmp_path):
    f = _fq(tmp_path / "f.fq", [("r1", "ACGT" * 10, "I" * 40)] * 2)
    r = _fq(tmp_path / "r.fq", [("r1", "ACGT" * 10, "I" * 40)])
    with pytest.raises(ValueError, match="Mismatched forward and reverse"):
        dada2.fastq_paired_filter((f, r), (tmp_path / "a.fq", tmp_path / "b.fq"), rm_phix=False)
    with pytest.raises(ValueError, match="must be different"):
        dada2.fastq_paired_filter((f, r), (f, tmp_path / "b.fq"))
    with pytest.raises(ValueError, match="pair"):
        dada2.fastq_paired_filter((f, r), (tmp_path / "a", tmp_path / "b"), trunc_len=(1, 2, 3))


def test_paired_match_ids_undetectable(tmp_path):
    f = _fq(tmp_path / "f.fq", [("read1", "ACGT" * 10, "I" * 40)])
    r = _fq(tmp_path / "r.fq", [("read1", "ACGT" * 10, "I" * 40)])
    with pytest.raises(ValueError, match="identifier field"):
        dada2.fastq_paired_filter((f, r), (tmp_path / "a", tmp_path / "b"), match_ids=True)
    stats = dada2.fastq_paired_filter(
        (f, r), (tmp_path / "a", tmp_path / "b"), match_ids=True, id_field=0, rm_phix=False,
        compress=False, verbose=1,
    )  # fmt: skip
    assert stats == (1, 1)


def test_paired_old_casava_ids(tmp_path):
    ids_f = ["HWUSI-EAS100R:6:73:941:1973#0/1", "HWUSI-EAS100R:6:73:941:1974#0/1"]
    f = _fq(tmp_path / "f.fq", [(i, "ACGT" * 10, "I" * 40) for i in ids_f])
    r = _fq(tmp_path / "r.fq", [(i.split("#")[0], "ACGT" * 10, "I" * 40) for i in ids_f[1:]])
    stats = dada2.fastq_paired_filter(
        (f, r), (tmp_path / "a.fq", tmp_path / "b.fq"), match_ids=True, rm_phix=False
    )
    assert stats == (2, 1)


def test_remove_primers_validation_and_no_hits(tmp_path):
    src = example_path("samPBprimers.fastq.gz")
    with pytest.raises(ValueError, match="distinct from the input"):
        dada2.remove_primers(src, src, "ACGT")
    with pytest.raises(TypeError, match="strings"):
        dada2.remove_primers(src, tmp_path / "o.fq", 3)  # type: ignore[arg-type]
    with pytest.warns(UserWarning, match="No reads passed the primer"):
        df = dada2.remove_primers(src, tmp_path / "o.fq", "GGGGGGGGGGGGGGGGGGGGGGG", max_mismatch=0)
    assert df["reads_out"].tolist() == [0] and not (tmp_path / "o.fq").exists()


def test_remove_primers_verbose_uncompressed(tmp_path):
    out = tmp_path / "o.fq"
    df = dada2.remove_primers(
        example_path("samPBprimers.fastq.gz"), out, "AGRGTTYGATYMTGGCTCAG",
        dada2.reverse_complement("RGYTACCTTGTTACGACTT"), compress=False, verbose=1,
    )  # fmt: skip
    assert df["reads_out"].iloc[0] == 45
    with pytest.raises(OSError), gzip.open(out) as fh:
        fh.read()
