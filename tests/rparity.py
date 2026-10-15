"""Loader for the R reference fixtures in tests/reference/data (see CLAUDE.md, "Fixtures").

R values are always JSON arrays (``auto_unbox=FALSE``); use :func:`scalar` to unwrap.
"""

from __future__ import annotations

import gzip
import json
from functools import cache
from importlib.resources import files
from pathlib import Path
from typing import Any

import numpy as np

DATA = Path(__file__).parent / "reference" / "data"
INPUTS = DATA / "inputs"
# Example data as installed with the package (so the suite also runs against a built wheel).
EXTDATA = Path(str(files("pydada2") / "data" / "extdata"))


@cache
def load(name: str) -> Any:
    """Load fixture ``name`` (without extension) from tests/reference/data."""
    plain = DATA / f"{name}.json"
    if plain.exists():
        return json.loads(plain.read_text())
    gz = DATA / f"{name}.json.gz"
    with gzip.open(gz, "rt") as fh:
        return json.load(fh)


def scalar(x: Any) -> Any:
    """Unwrap a length-1 JSON array written by R."""
    if isinstance(x, list):
        if len(x) != 1:
            raise ValueError(f"expected a length-1 array, got length {len(x)}")
        return x[0]
    return x


def vec(x: Any, dtype: Any = None) -> np.ndarray:
    """R vector -> 1-D numpy array (``null`` -> NaN for float dtypes)."""
    if x is None or x == {}:
        return np.array([], dtype=dtype)
    if not isinstance(x, list):
        x = [x]
    if dtype is not None and np.dtype(dtype).kind == "f":
        return np.array([np.nan if v is None else v for v in x], dtype=dtype)
    return np.array(x, dtype=dtype)


def matrix(obj: dict[str, Any], dtype: Any = float) -> np.ndarray:
    """R matrix fixture ``{rownames, colnames, values}`` -> 2-D array (null -> NaN)."""
    values = obj["values"]
    if len(values) == 0:
        ncol = len(obj.get("colnames") or [])
        return np.zeros((0, ncol), dtype=dtype)
    if np.dtype(dtype).kind == "f":
        return np.array([[np.nan if v is None else v for v in row] for row in values], dtype=dtype)
    return np.array(values, dtype=dtype)


def colnames(obj: dict[str, Any]) -> list[str]:
    return list(obj.get("colnames") or [])


def rownames(obj: dict[str, Any]) -> list[str]:
    return list(obj.get("rownames") or [])


# ---------------------------------------------------------------- domain helpers ---
def derep_fixture(name: str) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Return (sequences, abundances, mean quals [n x L], map 0-based/-1) for a derep fixture."""
    d = load("derep")[name]
    seqs = list(d["uniques"]["sequence"])
    abund = np.array(d["uniques"]["abundance"], dtype=np.int64)
    qsum = np.array(
        [[np.nan if v is None else v for v in row] for row in d["qual_sums"]], dtype=np.float64
    ).reshape(len(seqs), scalar(d["quals_ncol"]))
    quals = qsum / abund[:, None]
    m = np.array([-1 if v is None else v - 1 for v in d["map"]], dtype=np.int64)
    return seqs, abund, quals, m


def r_index(values: list[Any]) -> np.ndarray:
    """R 1-based integer vector with NA -> 0-based int64 array with -1."""
    return np.array([-1 if v is None else v - 1 for v in values], dtype=np.int64)


def r_int(values: list[Any]) -> np.ndarray:
    """R integer vector with NA -> int64 array with -1 for NA."""
    return np.array([-1 if v is None else v for v in values], dtype=np.int64)


def assert_close_r(actual: Any, r_values: Any, rtol: float = 1e-10, err_msg: str = "") -> None:
    """Compare doubles to an R vector where ``null`` stands for NA/NaN/Inf (jsonlite na="null").

    Positions that are null in R must be non-finite in ``actual``; the rest must match to rtol.
    """
    actual = np.asarray(actual, dtype=float)
    r_list = r_values if isinstance(r_values, list) else [r_values]
    mask = np.array([v is None for v in r_list], dtype=bool)
    expected = np.array([np.nan if v is None else v for v in r_list], dtype=float)
    assert actual.shape == expected.shape, (actual.shape, expected.shape, err_msg)
    assert not np.isfinite(actual[mask]).any(), f"{err_msg}: expected non-finite at R nulls"
    np.testing.assert_allclose(actual[~mask], expected[~mask], rtol=rtol, err_msg=err_msg)
