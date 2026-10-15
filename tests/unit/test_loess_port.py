"""The NumPy port of netlib dloess (1-D interpolate surface) against skmisc (same algorithm)."""

import numpy as np
import pytest

from pydada2._loess import loess_interpolate_1d
from pydada2.error_models import _skmisc_loess

from rparity import load, matrix


def _loess_inputs(trans, k=1, i=0):
    tot = trans[4 * i : 4 * i + 4].sum(axis=0)
    ok = tot > 0
    q = np.arange(trans.shape[1], dtype=float)
    return q[ok], np.log10((trans[k] + 1)[ok] / tot[ok]), tot[ok]


def test_loess_port_matches_skmisc_on_regular_data():
    """Where skmisc works, the NumPy port of the same algorithm agrees to rounding."""
    trans = matrix(load("error_models")["trans"])
    for k in (1, 2, 3):
        x, y, w = _loess_inputs(trans, k=k)
        xs = np.linspace(x.min(), x.max(), 57)
        np.testing.assert_allclose(
            loess_interpolate_1d(x, y, w, xs), _skmisc_loess(x, y, w, xs), rtol=1e-11
        )


def test_loess_port_random_regular_inputs():
    rng = np.random.default_rng(3)
    for _ in range(20):
        n = int(rng.integers(10, 45))
        x = np.sort(rng.choice(np.arange(60.0), n, replace=False))
        y = -0.05 * x + rng.normal(0, 0.3, n)
        w = 10 ** rng.uniform(0, 5, n)
        np.testing.assert_allclose(
            loess_interpolate_1d(x, y, w, x), _skmisc_loess(x, y, w, x), rtol=1e-10, atol=1e-12
        )


def test_loess_port_ties_and_errors():
    rng = np.random.default_rng(5)
    for rep in (2, 3):  # tied x values: the k-d tree splits at the nearest gap
        x = np.repeat(np.arange(0.0, 16.0), rep)
        y = np.sin(x / 3) + rng.normal(0, 0.1, x.size)
        w = 10 ** rng.uniform(0, 3, x.size)
        xs = np.unique(x)
        np.testing.assert_allclose(
            loess_interpolate_1d(x, y, w, xs), _skmisc_loess(x, y, w, xs), rtol=1e-11
        )
    with pytest.raises(ValueError, match="span"):
        loess_interpolate_1d([1.0], [1.0], [1.0], [1.0])
    with pytest.raises(ValueError, match="not finite"):
        loess_interpolate_1d([1.0, 2.0], [1.0, 2.0], [1.0, 1.0], [1.5])
