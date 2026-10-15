from pathlib import Path

import pytest

from pydada2 import NotFittedError
from pydada2._validation import as_path_list, check_is_fitted


class Est:
    def __init__(self):
        self.param = 1


def test_not_fitted_error_hierarchy():
    assert issubclass(NotFittedError, ValueError)
    assert issubclass(NotFittedError, AttributeError)


def test_check_is_fitted_default():
    est = Est()
    with pytest.raises(NotFittedError, match="Est instance is not fitted"):
        check_is_fitted(est)
    est.result_ = 1
    check_is_fitted(est)


def test_check_is_fitted_attributes():
    est = Est()
    est.a_ = 1
    check_is_fitted(est, ["a_"])
    with pytest.raises(NotFittedError):
        check_is_fitted(est, ["a_", "b_"])


def test_as_path_list():
    assert as_path_list("x.fq") == [Path("x.fq")]
    assert as_path_list(Path("x.fq")) == [Path("x.fq")]
    assert as_path_list(["a", Path("b")]) == [Path("a"), Path("b")]
    with pytest.raises(TypeError, match="fls must be"):
        as_path_list(3, name="fls")  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="only str"):
        as_path_list(["a", 3])  # type: ignore[list-item]
