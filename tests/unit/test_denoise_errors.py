import logging

import numpy as np
import pandas as pd
import pytest
from sklearn.base import clone

import pydada2 as dada2
from pydada2 import DADA, DadaParams, NotFittedError
from pydada2.datasets import err_balanced_f, example_path
from pydada2.denoise import combine_dereps, make_table
from pydada2.error_models import (
    TRANSITIONS,
    BinnedQualErrorModel,
    ErrorLearner,
    LoessErrorModel,
    NoQualErrorModel,
    PacBioErrorModel,
    accumulate_trans,
    get_errors,
)


@pytest.fixture(scope="module")
def drp():
    return dada2.derep_fastq(example_path("sam1F.fastq.gz"))


def _trans():
    rng = np.random.default_rng(1)
    t = np.zeros((16, 41))
    t[[0, 5, 10, 15]] = rng.integers(1000, 5000, size=(4, 41))
    for r in range(16):
        if r not in (0, 5, 10, 15):
            t[r] = rng.integers(0, 20, size=41)
    return t


def test_dada_params_replace_and_resolve():
    p = DadaParams().replace(gap_penalty=8, homopolymer_gap_penalty=1)
    r = p.resolved()
    assert r.gap_penalty == -8 and r.homopolymer_gap_penalty == -1
    assert not r.vectorized_alignment
    assert not DadaParams(band_size=0).resolved().vectorized_alignment
    with pytest.raises(TypeError, match="Unknown DADA option"):
        DadaParams().replace(nope=1)


def test_dada_estimator_api(drp):
    est = DADA(err=err_balanced_f(), params={"band_size": 32})
    assert "band_size" in repr(est)
    assert clone(est).get_params()["params"] == {"band_size": 32}
    est.fit(drp)
    assert est.n_iter_ == 1 and est.converged_
    r = est[0]
    assert "sequence variants" in repr(r)
    assert r.sequence == r.denoised.index.tolist()
    assert dada2.get_uniques(r).equals(r.denoised)
    assert np.array_equal(get_errors(est), est.err_out_)


def test_dada_inputs(drp, tmp_path):
    a = DADA(err=err_balanced_f()).fit({"s1": drp})
    assert a.sample_names_ == ["s1"]
    b = DADA(err=err_balanced_f()).fit(example_path("sam1F.fastq.gz"))
    assert b.sample_names_ == ["sam1F.fastq.gz"]
    c = DADA(err=err_balanced_f()).fit([example_path("sam1F.fastq.gz")])
    assert c.results_[0].clusters.equals(b.results_[0].clusters)
    with pytest.raises(TypeError, match="Dereplication objects"):
        DADA(err=err_balanced_f()).fit([3])
    with pytest.raises(ValueError, match="No samples"):
        DADA(err=err_balanced_f()).fit([])


def test_dada_validation(drp):
    with pytest.raises(ValueError, match="pool"):
        DADA(err=err_balanced_f(), pool="yes").fit([drp, drp])
    with pytest.raises(ValueError, match="omega_a"):
        DADA(err=err_balanced_f(), params={"omega_a": 2}).fit(drp)
    with pytest.raises(ValueError, match="omega_p"):
        DADA(err=err_balanced_f(), params={"omega_p": -1}).fit(drp)
    with (
        pytest.raises(ValueError, match="Some error correction"),
        pytest.warns(UserWarning, match="Strict error correction"),
    ):
        DADA(self_consistent=True, params={"omega_c": 1}).fit(drp)
    with pytest.warns(UserWarning, match="omega_p should"):
        DADA(err=err_balanced_f(), priors=[drp.sequences[0]], params={"omega_p": 1e-50}).fit(drp)
    with pytest.raises(TypeError, match="params must be"):
        DADA(err=err_balanced_f(), params=3).fit(drp)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="Error matrix"):
        DADA(err=None).fit(drp)
    bad = dada2.Dereplication(
        uniques=pd.Series([2], index=["ACGTNACGTACGT"]), quals=np.full((1, 13), 30.0),
        map=np.zeros(2, dtype=np.int64),
    )  # fmt: skip
    with pytest.raises(ValueError, match="A/C/G/T"):
        DADA(err=err_balanced_f()).fit(bad)


def test_dada_verbose_logging(drp, caplog):
    with caplog.at_level(logging.DEBUG, logger="pydada2"):
        DADA(err=err_balanced_f()[:, :30], verbose=2).fit(drp)
    assert any("repeating its last column" in m for m in caplog.messages)


def test_combine_dereps_and_table(drp):
    pooled = combine_dereps([drp, drp])
    assert pooled.n_reads == 2 * drp.n_reads
    np.testing.assert_allclose(pooled.quals, drp.quals, rtol=1e-12, equal_nan=True)
    with pytest.raises(ValueError):
        combine_dereps([])
    t = make_table([pd.Series({"A": 1, "C": 2}), pd.Series({"C": 3})])
    assert t.columns.tolist() == ["A", "C"] and t.to_numpy().tolist() == [[1, 2], [0, 3]]


def test_error_models_api():
    t = _trans()
    for m in (LoessErrorModel(), NoQualErrorModel(), PacBioErrorModel()):
        e = m.fit(t).err_
        assert e.shape == (16, 41) and ((e >= 0) & (e <= 1)).all()
        assert m.to_frame().index.tolist() == list(TRANSITIONS)
    with pytest.raises(NotFittedError):
        LoessErrorModel().to_frame()
    with pytest.raises(ValueError, match="16 rows"):
        LoessErrorModel()(np.zeros((4, 5)))
    df = pd.DataFrame(t, index=list(TRANSITIONS), columns=range(41))
    np.testing.assert_array_equal(NoQualErrorModel()(df), NoQualErrorModel()(t))
    with pytest.raises(NotImplementedError):
        dada2.error_models.ErrorModel()(t)


def test_binned_errors():
    t = _trans()
    with pytest.raises(ValueError, match="must be provided"):
        BinnedQualErrorModel()(t)
    with pytest.raises(ValueError, match="higher quality"):
        BinnedQualErrorModel([0, 20])(t)
    with pytest.raises(ValueError, match="lower quality"):
        BinnedQualErrorModel([5, 40])(t)
    with pytest.warns(UserWarning, match="Maximum observed"):
        BinnedQualErrorModel([0, 20, 41])(np.hstack([t, np.zeros((16, 1))]))


def test_pacbio_q93_position():
    t = np.hstack([_trans()] * 3)[:, :94]
    t2 = np.hstack([t[:, :93], t[:, 93:]])
    assert PacBioErrorModel()(t2).shape == (16, 94)


def test_loess_failure_and_accumulate():
    with pytest.raises(ValueError, match="could not be estimated"):
        LoessErrorModel()(np.zeros((16, 3)))
    with pytest.raises(ValueError):
        accumulate_trans([])
    assert accumulate_trans([np.ones((16, 2)), np.ones((16, 3))])[:, 2].tolist() == [1] * 16


def test_get_errors_variants():
    e = err_balanced_f()
    assert np.array_equal(get_errors(pd.DataFrame(e)), e)
    d = get_errors({"err_out": e, "err_in": None, "trans": None}, detailed=True)
    assert d["err_in"] is None
    with pytest.raises(ValueError, match="None"):
        get_errors(None)
    with pytest.raises(ValueError, match="16 rows"):
        get_errors(np.ones((4, 4)))
    with pytest.raises(ValueError, match=">= 0"):
        get_errors(-np.ones((16, 4)))
    with pytest.raises(ValueError, match="<= 1"):
        get_errors(2 * np.ones((16, 4)))
    with pytest.warns(UserWarning, match="Zero"):
        get_errors(np.zeros((16, 4)))
    with pytest.raises(TypeError):
        get_errors(3)

    class R:
        def __init__(self, err):
            self.err_out, self.err_in, self.trans = err, err, np.ones((16, 41))

    assert get_errors([R(e), R(e)], detailed=True)["trans"][0, 0] == 2
    with pytest.raises(ValueError, match="same output error"):
        get_errors([R(e), R(e * 0.5)])
    learner = type("L", (), {"err_": e})()
    assert np.array_equal(get_errors(learner), e)


def test_error_learner_randomize_and_derep_input(caplog):
    drps = dada2.derep_fastq([example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")])
    with caplog.at_level(logging.INFO, logger="pydada2"):
        le = ErrorLearner(randomize=True, random_state=0, n_bases=1, verbose=1).fit(drps)
    assert le.n_samples_ == 1 and le.err_.shape[0] == 16
    assert any("will be used for learning" in m for m in caplog.messages)
    with pytest.raises(ValueError, match="No samples"):
        ErrorLearner().fit([])
    le2 = ErrorLearner(params=DadaParams(), n_bases=1).fit(drps[0])
    assert le2.n_iter_ >= 2
