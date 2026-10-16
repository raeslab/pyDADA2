"""Hypothesis property tests."""

import tempfile
from pathlib import Path

import numpy as np
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

import pydada2 as dada2
from pydada2 import align

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])

IUPAC = "ACGTRYSWKMBDHVN"
acgt = st.text(alphabet="ACGT", min_size=1, max_size=60)
iupac = st.text(alphabet=IUPAC + IUPAC.lower(), max_size=80)


# ------------------------------------------------------------- reverse complement ---
@SETTINGS
@given(iupac)
def test_reverse_complement_is_an_involution(s):
    rc = dada2.reverse_complement(s)
    assert len(rc) == len(s)
    assert dada2.reverse_complement(rc) == s.upper()  # output is upper case, like R's rc()


@SETTINGS
@given(st.lists(iupac, max_size=5))
def test_reverse_complement_vectorised(seqs):
    assert dada2.reverse_complement(seqs) == [dada2.reverse_complement(s) for s in seqs]


# ---------------------------------------------------------------------- alignment ---
def _endsfree_score(al1, al2, match=5, mismatch=-4, gap=-8):
    """Score of an alignment under nwalign's model: leading/trailing gaps are free."""
    n = len(al1)
    lead = 0
    if al1.startswith("-"):
        lead = len(al1) - len(al1.lstrip("-"))
    elif al2.startswith("-"):
        lead = len(al2) - len(al2.lstrip("-"))
    trail = 0
    if al1.endswith("-"):
        trail = len(al1) - len(al1.rstrip("-"))
    elif al2.endswith("-"):
        trail = len(al2) - len(al2.rstrip("-"))
    score = 0
    for a, b in zip(al1[lead : n - trail], al2[lead : n - trail], strict=True):
        if a == "-" or b == "-":
            score += gap
        else:
            score += match if a == b else mismatch
    return score


@st.composite
def related_pair(draw):
    """A random sequence and a mutated / shifted copy (realistic alignment problems)."""
    s1 = draw(st.text(alphabet="ACGT", min_size=5, max_size=80))
    s2 = list(s1)
    for _ in range(draw(st.integers(0, 6))):
        op = draw(st.sampled_from(["sub", "ins", "del"]))
        i = draw(st.integers(0, max(0, len(s2) - 1)))
        if op == "sub" and s2:
            s2[i] = draw(st.sampled_from("ACGT"))
        elif op == "ins":
            s2.insert(i, draw(st.sampled_from("ACGT")))
        elif op == "del" and len(s2) > 1:
            del s2[i]
    lo = draw(st.integers(0, min(10, len(s2) - 1)))
    hi = len(s2) - draw(st.integers(0, min(10, len(s2) - lo - 1)))
    return s1, "".join(s2[lo:hi])


@SETTINGS
@given(st.one_of(related_pair(), st.tuples(acgt, acgt)))
def test_nwalign_properties(pair):
    s1, s2 = pair
    al1, al2 = align.nwalign(s1, s2)
    assert len(al1) == len(al2)
    assert al1.replace("-", "") == s1 and al2.replace("-", "") == s2
    assert not any(a == "-" and b == "-" for a, b in zip(al1, al2, strict=True))
    score = _endsfree_score(al1, al2)
    assert 0 <= score <= 5 * min(len(s1), len(s2))
    # the optimal score does not depend on the argument order
    rev1, rev2 = align.nwalign(s2, s1)
    assert _endsfree_score(rev1, rev2) == score
    m, mm, indel = align.eval_pair(al1, al2)
    assert m + mm <= min(len(s1), len(s2)) and min(m, mm, indel) >= 0
    assert align.nwhamming(s1, s2) == mm + indel


@SETTINGS
@given(acgt)
def test_nwalign_identity(s):
    assert align.nwalign(s, s) == (s, s)
    assert align.nwhamming(s, s) == 0
    assert tuple(align.eval_pair(s, s)) == (len(s), 0, 0)
    assert align.pair_consensus(s, s) == s


# ------------------------------------------------------------------ dereplication ---
read = st.tuples(
    st.text(alphabet="ACGTN", min_size=0, max_size=12),
    st.integers(min_value=2, max_value=40),
)


def _write_fastq(path: Path, reads):
    with open(path, "w") as fh:
        for i, (seq, q) in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{chr(33 + q) * len(seq)}\n")


@SETTINGS
@given(st.lists(read, min_size=1, max_size=60), st.integers(min_value=1, max_value=25))
def test_derep_counts_and_chunk_independence(reads, n):
    # As R's derepFastq, a chunk of only zero-length reads is an error (pinned in
    # test_derep_zero_length_chunk_is_an_error); the property covers the other inputs.
    assume(all(any(seq for seq, _ in reads[i : i + n]) for i in range(0, len(reads), n)))
    with tempfile.TemporaryDirectory() as d:
        fq = Path(d) / "r.fastq"
        _write_fastq(fq, reads)
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # zero-length reads are reported with a warning
            a = dada2.derep_fastq(fq)
            b = dada2.derep_fastq(fq, n=n)
    nonempty = [seq for seq, _ in reads if seq]
    # every non-empty read is counted once; zero-length reads map to -1
    assert int(a.uniques.sum()) == len(nonempty) == a.n_reads
    assert len(a.map) == len(reads)
    for (seq, _), m in zip(reads, a.map, strict=True):
        assert (m == -1) if not seq else (a.sequences[m] == seq)
    # abundances are non-increasing; chunking only changes the order of equal abundances
    for drp in (a, b):
        assert (np.diff(drp.uniques.to_numpy()) <= 0).all()
    assert dict(a.uniques) == dict(b.uniques)
    pos_b = {s: i for i, s in enumerate(b.sequences)}
    for i, s in enumerate(a.sequences):
        np.testing.assert_array_equal(a.quals[i], b.quals[pos_b[s]])
    by_ab = {}
    for drp_name, drp in (("a", a), ("b", b)):
        for s, ab in drp.uniques.items():
            by_ab.setdefault((drp_name, ab), set()).add(s)
    for (name, ab), seqs in by_ab.items():
        if name == "a":
            assert by_ab[("b", ab)] == seqs


def test_derep_zero_length_chunk_is_an_error():
    # R (sequenceIO.R): derepFastq stops with this message when any chunk of `n` reads has
    # only zero-length reads, and warns (then ignores them) otherwise.
    import warnings

    import pytest

    with tempfile.TemporaryDirectory() as d:
        fq = Path(d) / "r.fastq"
        _write_fastq(fq, [("A", 20), ("", 20)])
        with pytest.raises(ValueError, match="Only zero-length sequences detected"):
            dada2.derep_fastq(fq, n=1)
        with pytest.warns(UserWarning, match="Zero-length sequences detected"):
            drp = dada2.derep_fastq(fq, n=10)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert dict(drp.uniques) == {"A": 1} and drp.map.tolist() == [0, -1]


# ---------------------------------------------------------------------- filtering ---
filter_read = st.tuples(
    st.text(alphabet="ACGTN", min_size=1, max_size=60),
    st.lists(st.integers(min_value=0, max_value=40), min_size=60, max_size=60),
)


@SETTINGS
@given(
    st.lists(filter_read, min_size=1, max_size=40),
    st.fixed_dictionaries(
        {
            "trunc_q": st.integers(0, 15),
            "trim_left": st.integers(0, 5),
            "trim_right": st.integers(0, 5),
            "trunc_len": st.sampled_from([0, 10, 25, 40]),
            "min_len": st.integers(0, 20),
            "max_n": st.integers(0, 3),
            "max_ee": st.sampled_from([float("inf"), 1.0, 3.0]),
        }
    ),
)
def test_filtering_never_lengthens_reads(reads, params):
    with tempfile.TemporaryDirectory() as d:
        fin, fout = Path(d) / "in.fastq", Path(d) / "out.fastq"
        with open(fin, "w") as fh:
            for i, (seq, quals) in enumerate(reads):
                q = "".join(chr(33 + x) for x in quals[: len(seq)])
                fh.write(f"@r{i}\n{seq}\n+\n{q}\n")
        stats = dada2.fastq_filter(fin, fout, rm_phix=False, compress=False, **params)
        assert stats.reads_in == len(reads)
        out = []
        if fout.exists():
            import dnaio

            with dnaio.open(fout) as fh:
                out = [(r.name, r.sequence, r.qualities) for r in fh]
    assert stats.reads_out == len(out)
    src = {f"r{i}": (seq, q) for i, (seq, q) in enumerate(reads)}
    order = [int(name[1:]) for name, _, _ in out]
    assert order == sorted(order)  # reads keep their input order
    for name, seq, qual in out:
        orig, _ = src[name]
        assert len(seq) <= len(orig) and len(qual) == len(seq)
        assert seq in orig  # a contiguous piece of the input read
        assert len(seq) >= params["min_len"]
