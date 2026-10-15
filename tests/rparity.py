"""Loader for the R reference fixtures in tests/reference/data (see CLAUDE.md, "Fixtures").

R values are always JSON arrays (``auto_unbox=FALSE``); use :func:`scalar` to unwrap.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import warnings
from collections.abc import Callable, Iterator
from contextlib import contextmanager
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
def derep_fixture(
    name: str, fixture: str = "derep"
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Return (sequences, abundances, mean quals [n x L], map 0-based/-1) for a derep fixture."""
    d = load(fixture)[name]
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


# ------------------------------------------------------------------- FASTQ outputs ---
def fastq_digest(path: Path) -> dict[str, Any]:
    """Python side of the generator's ``fq_digest``: ids, widths and md5 of "id\\nseq\\nqual\\n"."""
    import dnaio

    if not path.exists():
        return {"exists": False}
    ids, widths, txt = [], [], []
    with dnaio.open(path) as fh:
        for r in fh:
            ids.append(r.name)
            widths.append(len(r.sequence))
            txt.append(f"{r.name}\n{r.sequence}\n{r.qualities}\n")
    return {
        "exists": True,
        "n": len(ids),
        "ids": ids,
        "widths": widths,
        "md5": hashlib.md5("".join(txt).encode()).hexdigest(),
    }


def check_fastq_output(path: Path, exp: dict[str, Any]) -> None:
    """Assert that a FASTQ file has exactly the records R wrote (or is absent like R's)."""
    got = fastq_digest(path)
    assert got["exists"] == scalar(exp["exists"]), path
    if got["exists"]:
        assert got["n"] == scalar(exp["n"])
        assert got["ids"] == exp["ids"]
        assert got["widths"] == exp["widths"]
        assert got["md5"] == scalar(exp["md5"])


FILTER_ARGS = {
    "truncLen": "trunc_len", "maxN": "max_n", "maxEE": "max_ee", "truncQ": "trunc_q",
    "rm.phix": "rm_phix", "trimLeft": "trim_left", "trimRight": "trim_right", "minLen": "min_len",
    "maxLen": "max_len", "minQ": "min_q", "rm.lowcomplex": "rm_lowcomplex", "n": "n",
    "orient.fwd": "orient_fwd", "matchIDs": "match_ids",
}  # fmt: skip


def filter_kwargs(rparams: Any, paired: bool) -> dict[str, Any]:
    """filterAndTrim arguments recorded by the generator -> filter_and_trim keyword arguments."""
    out: dict[str, Any] = {}
    for k, v in (rparams if isinstance(rparams, dict) else {}).items():
        if isinstance(v, list):
            v = [math.inf if x == "Inf" else x for x in v]
            v = v[0] if (len(v) == 1 or not paired) else tuple(v)
        out[FILTER_ARGS[k]] = v  # KeyError for an argument the tests do not map yet
    if "n" in out:
        out["n"] = int(out["n"])
    return out


@contextmanager
def r_warnings(expected: Any) -> Iterator[None]:
    """Run a block and require exactly R's warnings (messages, in order) as UserWarnings.

    ``expected`` is the generator's recorded ``warnings`` (list of messages; [] for none).
    Warnings of other categories (e.g. ResourceWarning) are ignored.
    """
    exp = [] if expected is None else list(expected)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter("always")
        yield
    got = [str(w.message) for w in rec if issubclass(w.category, UserWarning)]
    assert got == exp, f"warnings differ from R: {got} != {exp}"


def call_with_r_warnings(expected: Any, fn: Callable[..., Any], *args: Any, **kw: Any) -> Any:
    """``fn(*args, **kw)`` inside :func:`r_warnings`."""
    with r_warnings(expected):
        return fn(*args, **kw)


def fixture_input(name: str) -> Path:
    """Input file of a fixture: bundled example data or tests/reference/data/inputs."""
    return EXTDATA / name if (EXTDATA / name).exists() else INPUTS / name


def run_filter_case(fx: dict[str, Any], tmp_path: Path) -> None:
    """Run one generator ``run_filter`` case; compare counts, output records and warnings."""
    import pydada2 as dada2

    paired = fx.get("inputs_rev") is not None
    fwd = [fixture_input(f) for f in fx["inputs"]]
    filt = [tmp_path / f"F_{i}_{p.name}" for i, p in enumerate(fwd)]
    kwargs = filter_kwargs(fx["params"], paired)
    rev = filt_rev = None
    if paired:
        rev = [fixture_input(f) for f in fx["inputs_rev"]]
        filt_rev = [tmp_path / f"R_{i}_{p.name}" for i, p in enumerate(rev)]
    with r_warnings(fx["warnings"]):
        df = dada2.filter_and_trim(fwd, filt, rev, filt_rev, compress=True, **kwargs)
    assert df["reads_in"].tolist() == fx["reads_in"]
    assert df["reads_out"].tolist() == fx["reads_out"]
    for p, exp in zip(filt, fx["outputs"], strict=True):
        check_fastq_output(p, exp)
    if paired:
        for p, exp in zip(filt_rev, fx["outputs_rev"], strict=True):
            check_fastq_output(p, exp)


# ---------------------------------------------------------------------- taxonomy ---
def r_taxonomy(obj: dict[str, list[Any]]) -> np.ndarray:
    """R taxonomy data.frame fixture (columns V1..Vk) -> object array (None for NA)."""
    # Built directly: pandas 3 would turn None into NaN in string columns.
    cols = sorted(obj, key=lambda c: int(c[1:]))
    return np.array([obj[c] for c in cols], dtype=object).T


def binomial_bound(p_a: Any, p_b: Any = None) -> np.ndarray:
    """3 sigma of the difference of two independent 100-draw bootstrap estimates (+2 for ties),
    using the pooled estimate of the underlying proportion (spec §8.3)."""
    pooled = p_a if p_b is None else (np.asarray(p_a) + np.asarray(p_b)) / 2
    p = np.clip(np.asarray(pooled, dtype=float) / 100.0, 0.0, 1.0)
    return np.asarray(3 * np.sqrt(2 * 100 * p * (1 - p)) + 2)


def assert_bootstrap_agreement(boot: Any, r_boot: Any) -> None:
    """Spec §8.3: median |boot_py - boot_R| <= 10 and every entry within 3 sigma."""
    boot, r_boot = np.asarray(boot, dtype=float), np.asarray(r_boot, dtype=float)
    assert boot.shape == r_boot.shape
    diff = np.abs(boot - r_boot)
    assert np.median(diff) <= 10
    assert (diff <= binomial_bound(boot, r_boot)).all(), np.argwhere(
        diff > binomial_bound(boot, r_boot)
    )


def taxonomy_disagreements(got: Any, exp: Any, r_boot: Any, min_boot: float = 50) -> list[Any]:
    """Spec §8.3 at ``min_boot``: (row, level, R bootstrap) of every first disagreement per row.

    Raises AssertionError for a disagreement whose R bootstrap is outside min_boot +- 3 sigma.
    """
    out = []
    for i in range(exp.shape[0]):
        for j in range(exp.shape[1]):
            if got[i, j] == exp[i, j] or (got[i, j] is None and exp[i, j] is None):
                continue
            assert abs(r_boot[i, j] - min_boot) <= binomial_bound(np.array(min_boot)), (i, j)
            out.append((i, j, int(r_boot[i, j])))
            break  # deeper levels follow from this one
    return out
