"""dada_uniques parity (single-call fixtures) and thread-count determinism."""

import numpy as np
import pytest

import pydada2._core as core

from rparity import assert_close_r, derep_fixture, load, matrix, r_index, r_int, scalar

pytestmark = pytest.mark.parity

# dada fixture -> (derep fixtures per sample, err source)
CASES = {
    "single_filtF": (["filtF_sam1"], "errBalancedF"),
    "default_F": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "default_R": (["filtR_P_sam1", "filtR_P_sam2"], "errBalancedR"),
    "no_quals": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "not_gapless": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "band32": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "priors": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "detect_singletons": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "omega_a_1e10": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "no_kmers": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "not_vectorized": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "homo_gap": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "sse0": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "sse1": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "not_greedy": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "min_fold_hamming": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "max_clust5": (["filtF_P_sam1", "filtF_P_sam2"], "errBalancedF"),
    "pacbio": (["filt_samPB"], "pacbio"),
}


def _err(source):
    if source == "pacbio":
        return np.ascontiguousarray(matrix(load("learn_errors")["pacbio"]["err_out"]))
    return np.ascontiguousarray(matrix(load("bundled_err")[source]))


def _opts_kwargs(opts):
    o = {k: scalar(v) if isinstance(v, list) and len(v) == 1 else v for k, v in opts.items()}
    gap = o["GAP_PENALTY"]
    homo = o["HOMOPOLYMER_GAP_PENALTY"]
    return dict(
        match=o["MATCH"], mismatch=o["MISMATCH"], gap=gap, use_kmers=o["USE_KMERS"],
        kdist_cutoff=o["KDIST_CUTOFF"], band_size=int(o["BAND_SIZE"]), omega_a=o["OMEGA_A"],
        omega_p=o["OMEGA_P"], omega_c=o["OMEGA_C"], detect_singletons=o["DETECT_SINGLETONS"],
        max_clust=int(o["MAX_CLUST"]), min_fold=o["MIN_FOLD"], min_hamming=int(o["MIN_HAMMING"]),
        min_abund=int(o["MIN_ABUNDANCE"]), use_quals=True, final_consensus=False,
        vectorized_alignment=o["VECTORIZED_ALIGNMENT"], homo_gap=homo, sse=int(o["SSE"]),
        gapless=o["GAPLESS"], greedy=o["GREEDY"],
    )  # fmt: skip


def _run(case, sample, n_threads=1):
    derep_names, err_src = CASES[case]
    entry = load("dada")[case]
    fx = list(entry.values())[sample] if isinstance(entry, dict) else entry[sample]
    seqs, abund, quals, _ = derep_fixture(derep_names[sample])
    priors = set(load("dada")["prior_seqs"]) if case == "priors" else set()
    err = _err(err_src)
    qmax = int(np.ceil(np.nanmax(quals)))
    if err.shape[1] < qmax + 1:  # R extends err by repeating the last column
        err = np.hstack([err] + [err[:, -1:]] * (qmax + 1 - err.shape[1]))
    out = core.dada_uniques(
        seqs, abund.tolist(), [s in priors for s in seqs], np.ascontiguousarray(err),
        np.ascontiguousarray(quals), n_threads=n_threads, **_opts_kwargs(fx["opts"]),
    )  # fmt: skip
    return out, fx


def _compare(out, fx):
    cl = fx["clustering"]
    assert out["sequence"] == cl["sequence"]
    for col in ("abundance", "n0", "n1", "nunq"):
        assert out[col].tolist() == cl[col], col
    np.testing.assert_allclose(out["pval"], np.array(cl["pval"], float), rtol=1e-10)
    np.testing.assert_array_equal(out["birth_from"], r_index(cl["birth_from"]))
    np.testing.assert_array_equal(out["birth_ham"], r_int(cl["birth_ham"]))
    for col in ("birth_pval", "birth_fold", "birth_qave"):
        assert_close_r(out[col], cl[col], rtol=1e-10, err_msg=col)
    np.testing.assert_array_equal(out["map"], r_index(fx["map"]))
    np.testing.assert_allclose(out["pval_raw"], np.array(fx["pval"], float), rtol=1e-10)
    bs = fx["birth_subs"]
    np.testing.assert_array_equal(out["bs_pos"], np.array(bs["pos"], int) - 1)
    assert out["bs_ref"] == bs["ref"] and out["bs_sub"] == bs["sub"]
    np.testing.assert_allclose(out["bs_qual"], np.array(bs["qual"], float), rtol=1e-12)
    np.testing.assert_array_equal(out["bs_clust"], np.array(bs["clust"], int) - 1)
    np.testing.assert_array_equal(out["trans"], matrix(fx["trans"], dtype=np.int64))
    np.testing.assert_allclose(out["cluster_quals"], matrix(fx["quality"]), rtol=1e-12)


@pytest.mark.parametrize("case", sorted(CASES))
def test_dada_uniques_parity(case):
    for sample in range(len(CASES[case][0])):
        out, fx = _run(case, sample)
        _compare(out, fx)


@pytest.mark.parametrize("case", ["default_F", "not_greedy", "pacbio"])
def test_dada_uniques_thread_determinism(case):
    ref, _ = _run(case, 0, n_threads=1)
    for n in (2, 4):
        out, _ = _run(case, 0, n_threads=n)
        for key, val in ref.items():
            if isinstance(val, np.ndarray):
                np.testing.assert_array_equal(out[key], val, err_msg=key)
            else:
                assert out[key] == val, key


def test_dada_uniques_validation():
    err = np.ones((16, 41))
    with pytest.raises(ValueError, match="Zero input"):
        core.dada_uniques([], [], [], err, None)
    with pytest.raises(ValueError, match="A/C/G/T"):
        core.dada_uniques(["ACGTNACGTACG"], [1], [False], err, None)
    with pytest.raises(ValueError, match="16 rows"):
        core.dada_uniques(["ACGTACGTACG"], [1], [False], np.ones((4, 41)), None)
    with pytest.raises(ValueError, match="kmer-size"):
        core.dada_uniques(["ACGT"], [1], [False], err, None)
    with pytest.raises(ValueError, match="different lengths"):
        core.dada_uniques(["ACGTACGTACG"], [1, 2], [False], err, None)
    with pytest.raises(ValueError, match="qualities for each"):
        core.dada_uniques(["ACGTACGTACG"], [1], [False], err, np.ones((1, 5)))


def test_dada_uniques_log_callback():
    seqs, abund, quals, _ = derep_fixture("filtF_sam1")
    msgs = []
    err = np.ascontiguousarray(matrix(load("bundled_err")["errBalancedF"]))
    core.dada_uniques(
        seqs, abund.tolist(), [False] * len(seqs), err, np.ascontiguousarray(quals),
        verbose=True, log=msgs.append,
    )  # fmt: skip
    assert any("ALIGN" in m for m in msgs)
