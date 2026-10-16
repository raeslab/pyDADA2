"""Error models vs R (loessErrfun, noqualErrfun, makeBinnedQualErrfun, inflateErr...)."""

import numpy as np
import pytest

from pydada2.datasets import err_balanced_f
from pydada2.denoise import DADA
from pydada2.error_models import (
    BinnedQualErrorModel,
    LoessErrorModel,
    NoQualErrorModel,
    PacBioErrorModel,
    get_errors,
    inflate_err,
)

from rparity import load, matrix

pytestmark = pytest.mark.parity


@pytest.fixture(scope="module")
def em():
    return load("error_models")


def test_loess(em):
    np.testing.assert_allclose(
        LoessErrorModel()(matrix(em["trans"])), matrix(em["loess"]), rtol=1e-6
    )
    # tighter in practice: skmisc ships the same netlib dloess code as R
    np.testing.assert_allclose(
        LoessErrorModel()(matrix(em["trans"])), matrix(em["loess"]), rtol=1e-12
    )


def test_noqual(em):
    trans = matrix(em["trans"])
    np.testing.assert_allclose(NoQualErrorModel()(trans), matrix(em["noqual"]), rtol=1e-12)
    np.testing.assert_allclose(NoQualErrorModel(5)(trans), matrix(em["noqual_pc5"]), rtol=1e-12)


def test_binned(em):
    trans = matrix(em["trans"])
    got = BinnedQualErrorModel(em["binned_q"])(trans)
    np.testing.assert_allclose(got, matrix(em["binned"]), rtol=1e-12)


def test_inflate(em):
    err = matrix(load("bundled_err")["errBalancedF"])
    np.testing.assert_allclose(inflate_err(err, 3), matrix(em["inflate3"]), rtol=1e-12)
    np.testing.assert_allclose(
        inflate_err(err, 3, inflate_self_transitions=True), matrix(em["inflate3_self"]), rtol=1e-12
    )


def test_pacbio_err_from_trans():
    le = load("learn_errors")["pacbio"]
    trans = matrix(le["trans"])
    assert trans.shape[1] == 94  # Q0..Q93
    np.testing.assert_allclose(PacBioErrorModel()(trans), matrix(le["err_out"]), rtol=1e-6)


@pytest.mark.parametrize("key,model", [("fwd", LoessErrorModel()), ("noqual", NoQualErrorModel())])
def test_error_model_on_learned_trans(key, model):
    le = load("learn_errors")[key]
    np.testing.assert_allclose(model(matrix(le["trans"])), matrix(le["err_out"]), rtol=1e-6)


def test_get_errors_from_dada(em, dereps):
    # R: getErrors(d0) and getErrors(d0, detailed=TRUE), d0 <- dada(dF$sam1, err=errBalancedF)
    dd = DADA(err=err_balanced_f()).fit(dereps["dF"]["sam1"])
    exp = matrix(em["get_errors"])
    np.testing.assert_allclose(get_errors(dd), exp, rtol=1e-10)
    np.testing.assert_allclose(get_errors(dd[0]), exp, rtol=1e-10)
    det = get_errors(dd[0], detailed=True)
    np.testing.assert_allclose(
        det["err_out"], matrix(em["get_errors_detailed"]["err_out"]), rtol=1e-10
    )
    np.testing.assert_array_equal(
        det["trans"], matrix(em["get_errors_detailed"]["trans"], dtype=np.int64)
    )
    assert np.array_equal(get_errors(exp), exp)  # a plain matrix passes through
