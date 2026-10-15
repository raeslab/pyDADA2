import pytest
from sklearn.base import clone

from pydada2._base import BaseEstimator


class Inner(BaseEstimator):
    def __init__(self, alpha=1.0, beta="x"):
        self.alpha = alpha
        self.beta = beta


class Outer(BaseEstimator):
    def __init__(self, inner=None, n=3, values=(1, 2)):
        self.inner = inner
        self.n = n
        self.values = values


class Bad(BaseEstimator):
    def __init__(self, **kwargs):
        self.kwargs = kwargs


class NoInit(BaseEstimator):
    pass


def test_get_params_shallow_and_deep():
    est = Outer(inner=Inner(alpha=2.0), n=5)
    assert est.get_params(deep=False) == {"inner": est.inner, "n": 5, "values": (1, 2)}
    deep = est.get_params()
    assert deep["inner__alpha"] == 2.0
    assert deep["inner__beta"] == "x"


def test_set_params_nested_and_plain():
    est = Outer(inner=Inner())
    assert est.set_params(n=7, inner__alpha=9.0) is est
    assert est.n == 7
    assert est.inner.alpha == 9.0
    assert est.set_params() is est


def test_set_params_invalid():
    with pytest.raises(ValueError, match="Invalid parameter 'nope'"):
        Outer().set_params(nope=1)


def test_repr_shows_non_defaults_only():
    assert repr(Inner()) == "Inner()"
    assert repr(Inner(alpha=3)) == "Inner(alpha=3)"
    assert repr(Outer(inner=Inner(beta="y"))) == "Outer(inner=Inner(beta='y'))"


def test_varkw_rejected():
    with pytest.raises(TypeError, match=r"\*args/\*\*kwargs"):
        Bad().get_params()


def test_no_init():
    assert NoInit().get_params() == {}
    assert repr(NoInit()) == "NoInit()"


def test_sklearn_clone():
    est = Outer(inner=Inner(alpha=4.0), n=2)
    cl = clone(est)
    assert cl is not est
    assert cl.inner is not est.inner
    assert cl.get_params()["inner__alpha"] == 4.0
