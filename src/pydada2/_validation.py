"""Input validation helpers."""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

__all__ = ["NotFittedError", "as_path_list", "check_is_fitted"]

PathLike = str | os.PathLike[str]


class NotFittedError(ValueError, AttributeError):
    """Raised when an estimator is used before ``fit`` has been called."""


def check_is_fitted(estimator: Any, attributes: Sequence[str] | None = None) -> None:
    """Raise :class:`NotFittedError` if ``estimator`` has not been fitted.

    Parameters
    ----------
    estimator : object
        The estimator to check.
    attributes : sequence of str, optional
        Attributes that must exist. By default, any public attribute ending in ``_``.

    Raises
    ------
    NotFittedError
        If the estimator is not fitted.
    """
    if attributes is None:
        fitted = [k for k in vars(estimator) if k.endswith("_") and not k.startswith("_")]
    else:
        fitted = [a for a in attributes if hasattr(estimator, a)]
        if len(fitted) != len(attributes):
            fitted = []
    if not fitted:
        raise NotFittedError(
            f"This {type(estimator).__name__} instance is not fitted yet. "
            "Call 'fit' with appropriate arguments before using this estimator."
        )


def as_path_list(paths: PathLike | Iterable[PathLike], *, name: str = "paths") -> list[Path]:
    """Normalise a path or an iterable of paths to a list of :class:`pathlib.Path`.

    Parameters
    ----------
    paths : str, os.PathLike or iterable of those
        One or more paths.
    name : str
        Argument name used in error messages.

    Returns
    -------
    list of pathlib.Path
        The paths, in input order.

    Raises
    ------
    TypeError
        If ``paths`` is not a path or an iterable of paths.
    """
    if isinstance(paths, str | os.PathLike):
        return [Path(paths)]
    try:
        items = list(paths)
    except TypeError:
        raise TypeError(f"{name} must be a path or an iterable of paths.") from None
    out = []
    for p in items:
        if not isinstance(p, str | os.PathLike):
            raise TypeError(f"{name} must contain only str or os.PathLike, got {type(p).__name__}.")
        out.append(Path(p))
    return out
