"""Minimal scikit-learn compatible estimator base class (no sklearn dependency)."""

from __future__ import annotations

import inspect
from typing import Any, Self

__all__ = ["BaseEstimator"]

_NO_DEFAULT = inspect.Parameter.empty


class BaseEstimator:
    """Base class providing ``get_params``/``set_params`` and an informative ``repr``.

    Subclasses must declare every hyper-parameter as an explicit keyword argument of
    ``__init__`` and store it unmodified under an attribute of the same name, exactly
    like scikit-learn estimators. No ``*args``/``**kwargs`` in ``__init__``.
    """

    @classmethod
    def _get_param_names(cls) -> list[str]:
        init = cls.__init__
        if init is object.__init__:
            return []
        sig = inspect.signature(init)
        names = []
        for p in sig.parameters.values():
            if p.name == "self":
                continue
            if p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
                raise TypeError(
                    f"{cls.__name__}.__init__ must not use *args/**kwargs (found {p.name!r})."
                )
            names.append(p.name)
        return sorted(names)

    def get_params(self, deep: bool = True) -> dict[str, Any]:
        """Get the hyper-parameters of this estimator.

        Parameters
        ----------
        deep : bool, default=True
            Also return the parameters of nested estimators (as ``name__param``).

        Returns
        -------
        dict
            Parameter names mapped to their values.
        """
        out: dict[str, Any] = {}
        for name in self._get_param_names():
            value = getattr(self, name)
            if deep and hasattr(value, "get_params") and not isinstance(value, type):
                for sub_name, sub_value in value.get_params(deep=True).items():
                    out[f"{name}__{sub_name}"] = sub_value
            out[name] = value
        return out

    def set_params(self, **params: Any) -> Self:
        """Set the hyper-parameters of this estimator.

        Parameters
        ----------
        **params : dict
            Parameter names and values; nested estimators use ``name__param``.

        Returns
        -------
        self
            The estimator itself.

        Raises
        ------
        ValueError
            If a parameter name is not valid for this estimator.
        """
        if not params:
            return self
        valid = self.get_params(deep=True)
        nested: dict[str, dict[str, Any]] = {}
        for key, value in params.items():
            name, sep, sub = key.partition("__")
            if name not in valid:
                raise ValueError(
                    f"Invalid parameter {name!r} for estimator {type(self).__name__}. "
                    f"Valid parameters are: {sorted(self._get_param_names())!r}."
                )
            if sep:
                nested.setdefault(name, {})[sub] = value
            else:
                setattr(self, name, value)
                valid[name] = value
        for name, sub_params in nested.items():
            valid[name].set_params(**sub_params)
        return self

    def __repr__(self) -> str:
        sig = inspect.signature(type(self).__init__)
        parts = []
        for name in self._get_param_names():
            value = getattr(self, name)
            default = sig.parameters[name].default
            if default is _NO_DEFAULT or not _equal(value, default):
                parts.append(f"{name}={value!r}")
        return f"{type(self).__name__}({', '.join(parts)})"


def _equal(a: Any, b: Any) -> bool:
    if a is b:
        return True
    try:
        result = a == b
    except Exception:  # pragma: no cover - exotic objects
        return False
    if isinstance(result, bool):
        return result
    return False
