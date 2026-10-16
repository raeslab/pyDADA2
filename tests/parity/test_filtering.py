"""filter_and_trim / remove_primers parity with R (output records via canonical md5)."""

import pytest

import pydada2 as dada2

from rparity import EXTDATA, check_fastq_output, load, run_filter_case, scalar

pytestmark = pytest.mark.parity


CASES = sorted(load("filter_and_trim").keys())


@pytest.mark.parametrize("case", CASES)
def test_filter_and_trim_parity(case, tmp_path):
    run_filter_case(load("filter_and_trim")[case], tmp_path)


def test_filter_n_jobs_identical(tmp_path):
    fwd = [EXTDATA / "sam1F.fastq.gz", EXTDATA / "sam2F.fastq.gz"]
    a = [tmp_path / "a1.fq", tmp_path / "a2.fq"]
    b = [tmp_path / "b1.fq", tmp_path / "b2.fq"]
    d1 = dada2.filter_and_trim(fwd, a, trunc_len=240, max_ee=2, compress=False)
    d2 = dada2.filter_and_trim(fwd, b, trunc_len=240, max_ee=2, compress=False, n_jobs=4)
    assert d1.equals(d2)
    for x, y in zip(a, b, strict=True):
        assert x.read_bytes() == y.read_bytes()


PRIMER_CASES = ["default", "indels", "fwd_only", "no_trim", "mm0"]
PRIMER_ARGS = {
    "primer.fwd": "primer_fwd", "primer.rev": "primer_rev", "orient": "orient",
    "allow.indels": "allow_indels", "trim.fwd": "trim_fwd", "trim.rev": "trim_rev",
    "max.mismatch": "max_mismatch",
}  # fmt: skip


@pytest.mark.parametrize("case", PRIMER_CASES)
def test_remove_primers_parity(case, tmp_path):
    fx = load("remove_primers")[case]
    kwargs = {PRIMER_ARGS[k]: scalar(v) for k, v in fx["params"].items()}
    assert kwargs["primer_fwd"] == scalar(load("remove_primers")["F27"])
    if "primer_rev" in kwargs:
        assert kwargs["primer_rev"] == scalar(load("remove_primers")["R1492_rc"])
    fout = tmp_path / "out.fastq.gz"
    df = dada2.remove_primers(EXTDATA / "samPBprimers.fastq.gz", fout, **kwargs)
    assert df["reads_in"].tolist() == fx["reads_in"]
    assert df["reads_out"].tolist() == fx["reads_out"]
    check_fastq_output(fout, fx["output"])
