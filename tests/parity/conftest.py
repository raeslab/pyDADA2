"""Shared inputs for parity tests: filtered + dereplicated example data (as in the R generator)."""

import pytest

import pydada2 as dada2

from rparity import EXTDATA, load, r_warnings


@pytest.fixture(scope="session")
def filtered(tmp_path_factory):
    d = tmp_path_factory.mktemp("filtered")
    fx = load("filter_and_trim")
    # Same configurations as fixtures F1, P1 and FPB; pydada2 must warn exactly when R does.
    with r_warnings(fx["F1"]["warnings"] + fx["P1"]["warnings"] + fx["FPB"]["warnings"]):
        f1 = [d / "F1_sam1F.fastq", d / "F1_sam2F.fastq"]
        dada2.filter_and_trim(
            [EXTDATA / "sam1F.fastq.gz", EXTDATA / "sam2F.fastq.gz"], f1, trunc_len=240,
            max_n=0, max_ee=2, trunc_q=2, rm_phix=True, compress=False,
        )  # fmt: skip
        pf = [d / "P1_sam1F.fastq", d / "P1_sam2F.fastq"]
        pr = [d / "P1_sam1R.fastq", d / "P1_sam2R.fastq"]
        dada2.filter_and_trim(
            [EXTDATA / "sam1F.fastq.gz", EXTDATA / "sam2F.fastq.gz"], pf,
            [EXTDATA / "sam1R.fastq.gz", EXTDATA / "sam2R.fastq.gz"], pr, trunc_len=(240, 160),
            max_n=0, max_ee=(2, 2), trunc_q=2, rm_phix=True, compress=False,
        )  # fmt: skip
        pb = d / "FPB_samPB.fastq"
        dada2.filter_and_trim(
            EXTDATA / "samPB.fastq.gz", pb, min_q=3, min_len=1000, max_len=1600, max_n=0,
            rm_phix=False, max_ee=2, compress=False,
        )  # fmt: skip
    return {"F1": f1, "PF": pf, "PR": pr, "PB": pb}


@pytest.fixture(scope="session")
def dereps(filtered):
    return {
        "dF": dict(zip(["sam1", "sam2"], dada2.derep_fastq(filtered["F1"]), strict=True)),
        "dF_P": dict(zip(["sam1", "sam2"], dada2.derep_fastq(filtered["PF"]), strict=True)),
        "dR": dict(zip(["sam1", "sam2"], dada2.derep_fastq(filtered["PR"]), strict=True)),
        "dPB": dada2.derep_fastq(filtered["PB"]),
    }
