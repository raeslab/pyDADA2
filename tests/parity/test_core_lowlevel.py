"""Low-level parity of the C++ core (pydada2._core) against R dada2 1.36."""

import sys

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


# One unit in the last place of a subnormal double (2**-1074).
SUBNORMAL_ULP = 5e-324


def _check_ppois(ns, lams, ps):
    """ppois_upper vs R's ppois(lower.tail=FALSE).

    * zero and non-zero results must agree exactly (R keeps subnormal tails);
    * p >= 1e-292: relative difference <= 1e-12;
    * p < 1e-292: R recomputes these as exp(log p) and is itself off by up to 1.9e-12 from a
      45-digit reference there (pydada2: <= 6e-13), so <= 2.5e-12 relative; subnormal results
      carry fewer bits, so at least one subnormal ulp is allowed.
    """
    for n, lam, p in zip(ns, lams, ps, strict=True):
        got = core.ppois_upper(int(n), lam)
        assert (got == 0.0) == (p == 0.0), (n, lam, got, p)
        tol = 1e-12 * p if p >= 1e-292 else max(2.5e-12 * p, SUBNORMAL_ULP)
        assert abs(got - p) <= tol, (n, lam, got, p)


def test_ppois_upper():
    d = load("ppois")
    _check_ppois(d["n"], d["lambda"], d["p_upper"])


def test_ppois_upper_extended():
    """Large counts (n up to 2e6, lambda/n in [0.5, 1.5]), small counts, lambda down to 1e-320
    and a dense grid around the underflow boundary (Poisson tail underflow and large-n accuracy)."""
    d = load("ppois_extended")
    assert sum(p == 0 for p in d["p_upper"]) > 500  # exact zeros near the boundary
    assert sum(0 < p < 2.2250738585072014e-308 for p in d["p_upper"]) > 1000  # subnormals
    _check_ppois(d["n"], d["lambda"], d["p_upper"])


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


def test_seq_complexity_exact():
    """seqComplexity sums like R's long-double sum(): uniform k-mer spectra give exact integers
    (NGAGTACCNN -> 6, not 5.999999999999998), so rm_lowcomplex thresholds match R exactly.

    The general sequences are exact against R (glibc); on Windows and macOS the platform's log/exp
    can differ from glibc by an ulp (deviation #9 in docs/porting-notes.md)."""
    import pydada2 as dada2

    d = load("complexity_thresholds")
    assert d["k2"]["seqs"][0] == "NGAGTACCNN" and d["k2"]["value"][0] == 6.0
    with pytest.warns(UserWarning):
        got = dada2.seq_complexity(d["k2"]["seqs"])
    np.testing.assert_array_equal(got, d["k2"]["value"])
    with pytest.warns(UserWarning):
        got = dada2.seq_complexity(d["k3"]["seqs"], kmer_size=3)
    np.testing.assert_array_equal(got, d["k3"]["value"])
    got = dada2.seq_complexity(d["window"]["seqs"], window=12, by=3)
    np.testing.assert_array_equal(got, d["window"]["value"])
    cx = load("sequences")["complexity"]
    for key, kw in (("k2", {}), ("k3", {"kmer_size": 3}), ("win50", {"window": 50, "by": 5})):
        with pytest.warns(UserWarning):
            got = dada2.seq_complexity(cx["seqs"], **kw)
        exp = np.array([np.nan if v is None else v for v in cx[key]], dtype=float)
        if sys.platform in ("win32", "darwin"):  # not glibc
            np.testing.assert_array_equal(np.isnan(got), np.isnan(exp), err_msg=key)
            ok = ~np.isnan(exp)
            np.testing.assert_array_max_ulp(got[ok], exp[ok], maxulp=2)
        else:
            np.testing.assert_array_equal(got, exp, err_msg=key)


def test_match_pattern_max_mismatch_ge_length():
    """Biostrings reports every overlapping window when max.mismatch >= pattern length."""
    d = load("match_pattern_degenerate")
    for case in d["cases"]:
        subject = scalar(case["subject"]) if case["subject"] else ""
        pattern, k, fixed = (
            scalar(case["pattern"]),
            scalar(case["max_mismatch"]),
            scalar(case["fixed"]),
        )
        for mode, indels in (("mismatch", False), ("indels", True)):
            exp = [(s - 1, e) for s, e in zip(case[mode]["start"], case[mode]["end"], strict=True)]
            got = [tuple(x) for x in core.match_pattern(subject, pattern, k, indels, fixed)]
            assert got == exp, (subject, pattern, k, mode)
            hits = core.match_primer([subject], pattern, k, indels, fixed)
            assert int(hits["nhits"][0]) == len(exp)
            if exp:
                assert (int(hits["first_start"][0]), int(hits["first_end"][0])) == exp[0]
                assert (int(hits["last_start"][0]), int(hits["last_end"][0])) == exp[-1]


def test_remove_primers_max_mismatch_ge_length(tmp_path):
    import hashlib

    import pydada2 as dada2

    fx = load("match_pattern_degenerate")["remove_primers"]
    fout = tmp_path / "out.fastq"
    df = dada2.remove_primers(
        EXTDATA / "samPBprimers.fastq.gz", fout, primer_fwd="AGRG", primer_rev="GTCG",
        max_mismatch=4, orient=False, compress=False,
    )  # fmt: skip
    assert df["reads_in"].tolist() == fx["reads_in"]
    assert df["reads_out"].tolist() == fx["reads_out"]
    with dnaio.open(fout) as fh:
        recs = list(fh)
    assert [r.name for r in recs] == fx["output"]["ids"]
    assert [len(r.sequence) for r in recs] == fx["output"]["widths"]
    txt = "".join(f"{r.name}\n{r.sequence}\n{r.qualities}\n" for r in recs)
    assert hashlib.md5(txt.encode()).hexdigest() == scalar(fx["output"]["md5"])


def test_kmer_size_and_acgt_validation():
    s1, s2 = ["ACGTACGTAC"], ["ACGTACGTAA"]
    for k in (2, 9, 15, 16, -1):
        with pytest.raises(ValueError, match="kmer-size"):
            core.kmer_dist(s1, s2, k)
        with pytest.raises(ValueError, match="kmer-size"):
            core.kdist_matches(s1, s2, k)
    for k in (0, 9):
        with pytest.raises(ValueError, match="kmer-size"):
            core.kmer_matches(s1, s2, k)
        with pytest.raises(ValueError, match="kmer-size"):
            core.kord_dist(s1, s2, k, 0)
    with pytest.raises(ValueError, match="A/C/G/T"):
        core.nwalign("ACNT", "ACGT", MATCH, MISMATCH, GAP, GAP, -1, True)
    # the vectorized aligner only compares letters; R accepts any letters there
    assert core.nwvec(["ACNNGT"], ["ACGGT"], MATCH, MISMATCH, GAP, -1, True) == [
        ("ACNNGT", "ACG-GT")
    ]


def test_kord_dist_sse_paths():
    # SSE=1 uses the minimum length (kord_dist_SSEi, emulated exactly off x86); SSE=0 returns -1
    # for unequal lengths. The result must not depend on the platform.
    assert core.kord_dist(["ACGTACGTAC"], ["ACGTACGTAAGG"], 5, 0) == [-1.0]
    assert core.kord_dist(["ACGTACGTAC"], ["ACGTACGTAAGG"], 5, 1) == [pytest.approx(1 / 6)]


def test_int32_matrices_reject_out_of_range():
    seqs = ["ACGTACGTACGTAAAACCCCGGGGTTTT", "ACGTACGTACGTAAAACCCCGGGGTTTA"]
    args = (seqs, 1.5, 2, False, 4, MATCH, MISMATCH, GAP, 16)
    ref = core.table_bimera(np.array([[100, 3]], dtype=np.int32), *args)
    for mat in (
        np.array([[100, 3]], dtype=np.int64),
        np.asfortranarray([[100, 3]], dtype=np.int32),
    ):
        got = core.table_bimera(mat, *args)
        assert got[0].tolist() == ref[0].tolist() and got[1].tolist() == ref[1].tolist()
    for bad in (np.array([[2**32 + 5, 3]], dtype=np.int64), np.array([[1e12, 3.0]])):
        with pytest.raises(ValueError, match="int32 range"):
            core.table_bimera(bad, *args)
    q = np.array([[30, 20, core.NA_QUAL]], dtype=np.int64)
    np.testing.assert_array_equal(core.matrix_ee(q), core.matrix_ee(q.astype(np.int32)))
    with pytest.raises(ValueError, match="int32 range"):
        core.matrix_ee(np.array([[30, 2**33]], dtype=np.int64))


def test_dereplicator_concurrent_add_chunk():
    import threading

    d = core.Dereplicator()
    seqs, quals = ["ACGT" * 30, "ACGA" * 30] * 50, ["I" * 120] * 100

    def work():
        for _ in range(100):
            d.add_chunk(seqs, quals, 33)

    threads = [threading.Thread(target=work) for _ in range(6)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    res = d.finish()
    assert res["abundance"].tolist() == [30000, 30000]
    assert len(res["map"]) == 60000
