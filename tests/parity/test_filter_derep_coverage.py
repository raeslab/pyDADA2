"""Filters and dereplication on synthetic inputs where every filter bites (filter_coverage.json).

The example FASTQs contain no N/IUPAC base, no phiX and no low-complexity read, so maxN, maxLen,
rm.phix and rm.lowcomplex never removed a read in filter_and_trim.json. The generator builds
inputs/cov_F|R.fastq.gz (N's, IUPAC codes, lowercase stretches, phiX fragments, low-complexity,
shorter and longer reads) and checks, for each case, that every filter in ``bites`` changes R's
output when it is disabled. derep_coverage.json pins the dereplication tie order of reads that
differ only by IUPAC codes / N / case (Biostrings letter order, chunk-dependent).
"""

import numpy as np
import pytest

import pydada2 as dada2

from rparity import INPUTS, derep_fixture, load, run_filter_case

pytestmark = pytest.mark.parity

FCOV = load("filter_coverage")


@pytest.mark.parametrize("case", sorted(FCOV))
def test_filter_coverage_parity(case, tmp_path):
    fx = FCOV[case]
    assert fx["bites"], "the generator records which filters remove reads in this case"
    run_filter_case(fx, tmp_path)


def test_filter_coverage_inputs_have_ambiguous_and_lowercase_bases():
    import dnaio

    with dnaio.open(INPUTS / "cov_F.fastq.gz") as fh:
        raw = [r.sequence for r in fh]
    assert any("N" in s for s in raw)
    assert any(set(s) & set("RYSWKMBDHV") for s in raw)
    assert any(s != s.upper() for s in raw)


DEREP_CASES = {
    "default": ("derep_iupac.fastq.gz", {}),
    "n5": ("derep_iupac.fastq.gz", {"n": 5}),
    "cov_F": ("cov_F.fastq.gz", {}),
    "cov_F_n37": ("cov_F.fastq.gz", {"n": 37}),
}


@pytest.mark.parametrize("case", sorted(DEREP_CASES))
def test_derep_iupac_order(case):
    fname, kwargs = DEREP_CASES[case]
    drp = dada2.derep_fastq(INPUTS / fname, **kwargs)
    seqs, abund, quals, m = derep_fixture(case, "derep_coverage")
    assert drp.sequences == seqs
    np.testing.assert_array_equal(drp.uniques.to_numpy(), abund)
    np.testing.assert_array_equal(drp.map, m)
    assert drp.quals.shape == quals.shape
    np.testing.assert_allclose(drp.quals, quals, rtol=1e-12, equal_nan=True)


def test_derep_iupac_fixture_has_ties_between_letters():
    seqs, abund, _, _ = derep_fixture("default", "derep_coverage")
    by_abund: dict[int, set[str]] = {}
    for s, a in zip(seqs, abund, strict=True):
        if len(s) == 60:
            by_abund.setdefault(int(a), set()).add(s[29])
    assert {"C", "G", "R", "S", "W", "B", "V"} <= by_abund[1]  # ties among letter variants
    assert {"N", "Y"} <= by_abund[3] | by_abund[2]  # duplicated / lowercase copies
