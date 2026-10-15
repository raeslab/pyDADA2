"""Phase 2b: low-level parity of the C++ core (pydada2._core) against R dada2 1.36."""

import dnaio
import numpy as np
import pytest

import pydada2._core as core

from rparity import EXTDATA, load, matrix, scalar

pytestmark = pytest.mark.parity

MATCH, MISMATCH, GAP = 5, -4, -8


@pytest.fixture(scope="module")
def align():
    return load("align")


def _nw(s1, s2, v):
    if v.get("vec"):
        band = scalar(v["band"]) if "band" in v else -1
        return core.nwvec([s1], [s2], MATCH, MISMATCH, GAP, band, True)[0]
    match = scalar(v.get("match", [MATCH]))
    mismatch = scalar(v.get("mismatch", [MISMATCH]))
    gap = scalar(v.get("gap", [GAP]))
    homo = scalar(v.get("homo_gap", [gap]))
    band = scalar(v.get("band", [-1]))
    endsfree = scalar(v.get("endsfree", [True]))
    return core.nwalign(s1, s2, match, mismatch, gap, homo, band, endsfree)


@pytest.mark.parametrize(
    "variant",
    ["default", "band16", "noends", "homo", "vec", "vec_band16", "merge0", "merge8"],
)
def test_nwalign_and_eval_pair(align, variant):
    panel = align["panel"]
    pairs = np.array(align["pairs"]) - 1
    v = align["variants"][variant]
    if isinstance(v, list):  # empty R list -> []
        v = {}
    res = align["results"][variant]
    for k, (i, j) in enumerate(pairs):
        a1, a2 = _nw(panel[i], panel[j], v)
        assert a1 == res["al1"][k], (variant, k)
        assert a2 == res["al2"][k], (variant, k)
        assert list(core.eval_pair(a1, a2)) == [
            res["match"][k],
            res["mismatch"][k],
            res["indel"][k],
        ]


def test_nwhamming(align):
    panel = align["panel"]
    pairs = np.array(align["pairs"]) - 1
    for key, band in (("nwhamming", -1), ("nwhamming_band16", 16)):
        got = []
        for i, j in pairs:
            a1, a2 = core.nwalign(panel[i], panel[j], MATCH, MISMATCH, GAP, GAP, band, True)
            _, mm, indel = core.eval_pair(a1, a2)
            got.append(mm + indel)
        assert got == align[key]


def test_kmer_functions(align):
    panel = align["panel"]
    pairs = np.array(align["pairs"]) - 1
    s1 = [panel[i] for i in pairs[:, 0]]
    s2 = [panel[j] for j in pairs[:, 1]]
    km = align["kmer"]
    np.testing.assert_allclose(core.kmer_dist(s1, s2, 5), km["kmer_dist5"], rtol=1e-12)
    for sse in (0, 1, 2):
        np.testing.assert_allclose(
            core.kord_dist(s1, s2, 5, sse), km[f"kord_dist5_sse{sse}"], rtol=1e-12
        )
    assert core.kmer_matches(s1, s2, 5) == km["kmer_matches5"]
    assert core.kdist_matches(s1, s2, 5) == km["kdist_matches5"]


def test_pair_consensus(align):
    pc = align["pair_consensus"]
    for r in pc["results"]:
        a1, a2 = scalar(r["al1"]), scalar(r["al2"])
        assert list(core.eval_pair(a1, a2)) == r["eval"]
        assert core.pair_consensus(a1, a2, 1, False) == scalar(r["p1"])
        assert core.pair_consensus(a1, a2, 2, False) == scalar(r["p2"])
        assert core.pair_consensus(a1, a2, 1, True) == scalar(r["p1_trim"])
        assert core.pair_consensus(a1, a2, 2, True) == scalar(r["p2_trim"])


def test_is_acgt(align):
    d = align["is_acgt"]
    assert core.is_acgt(d["seqs"]).tolist() == d["result"]


def _read_quals(path, n):
    rows = []
    with dnaio.open(path) as fh:
        for i, rec in enumerate(fh):
            if i == n:
                break
            rows.append(np.frombuffer(rec.qualities.encode(), dtype=np.uint8).astype(np.int32) - 33)
    width = max(len(r) for r in rows)
    mat = np.full((len(rows), width), core.NA_QUAL, dtype=np.int32)
    for i, r in enumerate(rows):
        mat[i, : len(r)] = r
    return mat


def test_matrix_ee():
    d = load("matrix_ee")
    np.testing.assert_allclose(
        core.matrix_ee(_read_quals(EXTDATA / "sam1F.fastq.gz", 300)), d["sam1F_300"], rtol=1e-12
    )
    q = _read_quals(EXTDATA / "samPB.fastq.gz", 100)
    assert q.shape[1] == scalar(d["samPB_100_ncol"])
    np.testing.assert_allclose(core.matrix_ee(q), d["samPB_100"], rtol=1e-12)


def test_match_ref():
    d = load("sequences")["is_phix"]
    with dnaio.open(EXTDATA / "phix_genome.fa") as fh:
        phix = next(iter(fh)).sequence
    assert core.match_ref(d["seqs"], phix, 16, True) == d["match_ref"]
    assert core.match_ref(d["seqs"], phix, 16, False) == d["match_ref_overlap"]


def test_ppois_upper():
    d = load("ppois")
    for n, lam, p in zip(d["n"], d["lambda"], d["p_upper"], strict=True):
        got = core.ppois_upper(int(n), lam)
        # rel 1e-12; subnormal results (< ~1e-305) only carry a few significant digits
        assert got == pytest.approx(p, rel=1e-12, abs=1e-305), (n, lam)


def test_ppois_upper_edges():
    assert core.ppois_upper(-1, 3.0) == 1.0
    assert core.ppois_upper(3, 0.0) == 0.0
    assert core.ppois_upper(3, float("inf")) == 1.0
    assert np.isnan(core.ppois_upper(3, -1.0))


def test_table_bimera():
    ch = load("chimeras")
    st = matrix(ch["st_chim"], dtype=np.int32)
    seqs = ch["st_chim"]["colnames"]
    nflag, nsam = core.table_bimera(st, seqs, 1.5, 2, False, 4, MATCH, MISMATCH, GAP, 16)
    assert nflag.tolist() == ch["table_bimera"]["nflag"]
    assert nsam.tolist() == ch["table_bimera"]["nsam"]
    for n in (2, 4):
        f2, s2 = core.table_bimera(st, seqs, 1.5, 2, False, 4, MATCH, MISMATCH, GAP, 16, n)
        assert f2.tolist() == nflag.tolist() and s2.tolist() == nsam.tolist()


def test_is_bimera():
    d = load("chimeras")["is_bimera"]
    parents = d["parents"]
    for k, s in enumerate(d["seqs"]):
        pars = [p for p in parents if p != s]
        assert core.is_bimera(s, pars, False, 4, MATCH, MISMATCH, GAP, 16) == d["default"][k]
        assert core.is_bimera(s, pars, True, 4, MATCH, MISMATCH, GAP, 16) == d["one_off"][k]
        assert core.is_bimera(s, pars, False, 4, MATCH, MISMATCH, GAP, 2) == d["shift2"][k]


def test_input_validation():
    with pytest.raises(ValueError, match="same length"):
        core.eval_pair("AC-", "AC")
    with pytest.raises(ValueError, match="same length"):
        core.pair_consensus("AC-", "AC", 1, False)
    assert list(core.eval_pair("", "")) == [0, 0, 0]
    with pytest.raises(ValueError, match="Mismatched"):
        core.kmer_dist(["ACGTACGT"], [], 5)
    with pytest.raises(ValueError, match="equal length"):
        core.nwvec(["ACGT"], [], 5, -4, -8, -1, True)
