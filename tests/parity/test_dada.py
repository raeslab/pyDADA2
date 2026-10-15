"""Phase 5: DADA estimator parity with R's dada() for all fixture variants; learn_errors."""

import numpy as np
import pytest

from pydada2.datasets import err_balanced_f, err_balanced_r
from pydada2.denoise import DADA, DadaParams
from pydada2.error_models import NoQualErrorModel, PacBioErrorModel, learn_errors

from rparity import assert_close_r, load, matrix, r_index, r_int

pytestmark = pytest.mark.parity

PB_ERR = None


def _pb_err():
    return matrix(load("learn_errors")["pacbio"]["err_out"])


# fixture name -> (dereps key, err, DADA kwargs)
CASES = {
    "default_F": ("dF_P", "F", {}),
    "default_R": ("dR", "R", {}),
    "single_filtF": ("dF.sam1", "F", {}),
    "pool_true": ("dF_P", "F", {"pool": True}),
    "pool_pseudo": ("dF_P", "F", {"pool": "pseudo"}),
    "no_quals": ("dF_P", "F", {"params": {"use_quals": False}}),
    "not_gapless": ("dF_P", "F", {"params": {"gapless": False}}),
    "band32": ("dF_P", "F", {"params": {"band_size": 32}}),
    "priors": ("dF_P", "F", {"priors": "PRIORS"}),
    "detect_singletons": ("dF_P", "F", {"params": {"detect_singletons": True}}),
    "omega_a_1e10": ("dF_P", "F", {"params": {"omega_a": 1e-10}}),
    "no_kmers": ("dF_P", "F", {"params": {"use_kmers": False}}),
    "not_vectorized": ("dF_P", "F", {"params": {"vectorized_alignment": False}}),
    "homo_gap": (
        "dF_P", "F", {"params": {"homopolymer_gap_penalty": -1, "vectorized_alignment": False}}
    ),
    "sse0": ("dF_P", "F", {"params": {"sse": 0}}),
    "sse1": ("dF_P", "F", {"params": {"sse": 1}}),
    "not_greedy": ("dF_P", "F", {"params": {"greedy": False}}),
    "min_fold_hamming": (
        "dF_P", "F", {"params": {"min_fold": 2, "min_hamming": 2, "min_abundance": 2}}
    ),
    "max_clust5": ("dF_P", "F", {"params": {"max_clust": 5}}),
    "self_consist": ("dF_P", "F", {"self_consistent": True}),
    "pacbio": ("dPB", "PB", {"params": {"band_size": 32}}),
}  # fmt: skip


def _inputs(dereps, key):
    if key == "dF.sam1":
        return dereps["dF"]["sam1"]
    return dereps[key]


def _err(code):
    return {"F": err_balanced_f, "R": err_balanced_r, "PB": _pb_err}[code]()


def _fx_samples(case):
    e = load("dada")[case]
    return list(e.values()) if isinstance(e, dict) else e


def _compare(res, fx, case):
    cl = fx["clustering"]
    c = res.clusters
    assert c["sequence"].tolist() == cl["sequence"]
    for col in ("abundance", "n0", "n1", "nunq"):
        assert c[col].tolist() == cl[col], col
    assert_close_r(c["pval"], cl["pval"], err_msg="pval")
    np.testing.assert_array_equal(c["birth_from"], r_index(cl["birth_from"]))
    np.testing.assert_array_equal(c["birth_ham"], r_int(cl["birth_ham"]))
    for col in ("birth_pval", "birth_fold", "birth_qave"):
        assert_close_r(c[col], cl[col], err_msg=col)
    assert res.denoised.index.tolist() == fx["denoised"]["sequence"]
    assert res.denoised.tolist() == fx["denoised"]["abundance"]
    np.testing.assert_array_equal(res.map, r_index(fx["map"]))
    assert_close_r(res.pval, fx["pval"], err_msg="pval_raw")
    bs = fx["birth_subs"]
    np.testing.assert_array_equal(res.birth_subs["pos"], np.array(bs["pos"], int) - 1)
    assert res.birth_subs["ref"].tolist() == bs["ref"]
    np.testing.assert_array_equal(res.birth_subs["clust"], np.array(bs["clust"], int) - 1)
    np.testing.assert_array_equal(res.trans, matrix(fx["trans"], dtype=np.int64))
    np.testing.assert_allclose(res.quality, matrix(fx["quality"]), rtol=1e-12)
    np.testing.assert_allclose(res.err_out, matrix(fx["err_out"]), rtol=1e-10)
    if isinstance(fx["err_in"], list):
        assert len(res.err_in) == len(fx["err_in"])
        for a, b in zip(res.err_in, fx["err_in"], strict=True):
            np.testing.assert_allclose(a, matrix(b), rtol=1e-10)
    else:
        np.testing.assert_allclose(res.err_in, matrix(fx["err_in"]), rtol=1e-10)


def _compare_unordered(res, fx):
    """Cluster-order-independent comparison (DEVIATIONS.md #4: self-consistent mode)."""
    cl = fx["clustering"]
    r_seqs = cl["sequence"]
    c = res.clusters.set_index("sequence")
    assert sorted(c.index) == sorted(r_seqs)
    for k, sq in enumerate(r_seqs):
        for col in ("abundance", "n0", "n1", "nunq", "birth_ham"):
            exp = -1 if cl[col][k] is None else cl[col][k]
            assert c.loc[sq, col] == exp, (col, k)
    # same cluster membership of every unique sequence
    r_map = r_index(fx["map"])
    got = [res.clusters["sequence"].iloc[m] if m >= 0 else None for m in res.map]
    exp = [r_seqs[m] if m >= 0 else None for m in r_map]
    assert got == exp
    np.testing.assert_array_equal(res.trans, matrix(fx["trans"], dtype=np.int64))
    np.testing.assert_allclose(res.err_out, matrix(fx["err_out"]), rtol=1e-10)
    for a, b in zip(res.err_in, fx["err_in"], strict=True):
        np.testing.assert_allclose(a, matrix(b), rtol=1e-10)


UNORDERED = {"self_consist"}


@pytest.mark.parametrize("case", sorted(CASES))
def test_dada_parity(case, dereps):
    key, err, kwargs = CASES[case]
    kwargs = dict(kwargs)
    if kwargs.get("priors") == "PRIORS":
        kwargs["priors"] = load("dada")["prior_seqs"]
    dd = DADA(err=_err(err), **kwargs).fit(_inputs(dereps, key))
    fxs = _fx_samples(case)
    assert len(dd.results_) == len(fxs)
    for res, fx in zip(dd.results_, fxs, strict=True):
        if case in UNORDERED:
            _compare_unordered(res, fx)
        else:
            _compare(res, fx, case)


def _check_learned(le, fx):
    np.testing.assert_allclose(le.err_, matrix(fx["err_out"]), rtol=1e-10)
    np.testing.assert_array_equal(le.trans_, matrix(fx["trans"]))
    exp_in = fx["err_in"]
    assert len(le.err_in_) == len(exp_in)
    for a, b in zip(le.err_in_, exp_in, strict=True):
        np.testing.assert_allclose(a, matrix(b), rtol=1e-10)


def test_learn_errors_fwd_rev(filtered):
    fx = load("learn_errors")
    _check_learned(learn_errors(filtered["PF"]), fx["fwd"])
    _check_learned(learn_errors(filtered["PR"]), fx["rev"])


def test_learn_errors_noqual_and_nbases(filtered):
    fx = load("learn_errors")
    _check_learned(learn_errors(filtered["PF"], error_model=NoQualErrorModel()), fx["noqual"])
    le = learn_errors(filtered["PF"], n_bases=1e5)
    assert le.n_samples_ == 1
    _check_learned(le, fx["nbases_1e5"])


def test_learn_errors_pacbio(filtered):
    fx = load("learn_errors")["pacbio"]
    le = learn_errors(filtered["PB"], error_model=PacBioErrorModel(), params={"band_size": 32})
    _check_learned(le, fx)


def test_dada_n_jobs_determinism(dereps):
    a = DADA(err=err_balanced_f(), pool="pseudo").fit(dereps["dF_P"])
    b = DADA(err=err_balanced_f(), pool="pseudo", n_jobs=4).fit(dereps["dF_P"])
    for ra, rb in zip(a.results_, b.results_, strict=True):
        assert ra.clusters.equals(rb.clusters)
        np.testing.assert_array_equal(ra.map, rb.map)


def test_dada_params_object(dereps):
    a = DADA(err=err_balanced_f(), params=DadaParams(band_size=32)).fit(dereps["dF_P"])
    b = DADA(err=err_balanced_f(), params={"band_size": 32}).fit(dereps["dF_P"])
    assert a.results_[0].clusters.equals(b.results_[0].clusters)
