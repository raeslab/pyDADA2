"""One-dimensional loess with an interpolated (k-d tree) surface.

A clean-room re-implementation of the algorithm of the netlib ``dloess`` package
(W. S. Cleveland, E. Grosse and M.-J. Shyu, "Local regression models", 1992) for the single
case DADA2 needs: one predictor, ``family="gaussian"``, tricube weights, prior weights,
``surface="interpolate"``. It is the fallback of :class:`~pydada2.LoessErrorModel` for inputs on
which ``skmisc.loess`` (the same netlib code) aborts with what R's ``loess`` only reports as
warnings ("pseudoinverse used", "near singularities", "Chernobyl! trL>n"), i.e. transition
matrices with only a few observed quality scores. R returns a fit for these and DADA2's
``dada()`` suppresses the warnings, so pydada2 has to return the same fit.

The algorithm:

1. **k-d tree.** The bounding box of ``x`` is widened by 0.5 % of its range on both sides.
   A cell holding more than ``floor(n * span * cell)`` points is split at the median point
   (the ``floor((l + u) / 2)``-th order statistic of its points, moved to the nearest gap
   between distinct values), unless the split value is one of the cell's bounds.
2. **Vertex fits.** At every cell bound ``v`` a local polynomial of degree ``degree`` is fitted
   to the ``q = floor(n * span + 1e-5)`` nearest points with weights
   ``prior * (1 - (d / d_q) ** 3) ** 3``. The columns of the weighted design matrix are scaled
   to unit length and the least-squares problem is solved through the singular value
   decomposition, with singular values ``<= 100 * eps * s_max`` treated as zero (the
   "pseudoinverse"). The fitted value and slope at ``v`` are kept.
3. **Evaluation.** Inside a cell the surface is the cubic Hermite interpolant of the values and
   slopes at its two bounds.
"""

from __future__ import annotations

import math

import numpy as np
import numpy.typing as npt

__all__ = ["loess_interpolate_1d"]

FloatArray = npt.NDArray[np.float64]

_EPS = 2.0**-52  # netlib d1mach(4)
_BOX_MARGIN = 0.005


def _vertices(x: FloatArray, fc: int) -> FloatArray:
    """Cell bounds of the 1-D k-d tree (sorted)."""
    alpha, beta = float(x.min()), float(x.max())
    mu = _BOX_MARGIN * max(beta - alpha, 1e-10 * max(abs(alpha), abs(beta)) + 1e-30)
    lo, hi = alpha - mu, beta + mu
    xs = np.sort(x)
    verts = [lo, hi]
    cells = [(0, x.size - 1, lo, hi)]  # 0-based rank range [first, last] and cell bounds
    i = 0
    while i < len(cells):
        first, last, a, b = cells[i]
        i += 1
        if last - first + 1 <= fc or b - a <= 0:
            continue
        m = (first + last) // 2
        off = 0  # move the split to the nearest gap between distinct values: 0, +1, -1, +2, ...
        while first <= m + off < last:
            if xs[m + off] == xs[m + off + 1]:
                off = -off
                if off >= 0:
                    off += 1
            else:
                m += off
                break
        xi = float(xs[m])
        if xi in (a, b):
            continue
        cells.append((first, m, a, xi))
        cells.append((m + 1, last, xi, b))
        verts.append(xi)
    return np.array(sorted(verts), dtype=np.float64)


def _local_fit(
    x: FloatArray, y: FloatArray, w: FloatArray, at: float, nf: int, span: float, degree: int
) -> tuple[float, float]:
    """Value and slope at ``at`` of the local weighted polynomial fit."""
    dist = (x - at) ** 2
    near = np.argsort(dist, kind="stable")[:nf]
    rho = dist[near[-1]] * max(1.0, span)
    with np.errstate(divide="ignore", invalid="ignore"):
        u3 = np.sqrt(dist[near] / rho) ** 3
        wt = np.sqrt(w[near] * (1 - u3) ** 3)
    dx = x[near] - at
    design = np.column_stack([wt * dx**p for p in range(degree + 1)])
    eta = wt * y[near]
    if not (np.isfinite(design).all() and np.isfinite(eta).all()):
        return math.nan, math.nan  # empty neighbourhood (R's fit is NaN as well)
    k = design.shape[1]
    colnor = np.sqrt((design**2).sum(axis=0))
    colnor[~(colnor > 0)] = 1.0
    design = design / colnor
    qmat, r = np.linalg.qr(design, mode="reduced")
    rr = np.zeros((k, k))
    rr[: r.shape[0]] = r
    qty = np.zeros(k)
    qty[: qmat.shape[1]] = qmat.T @ eta
    uu, s, vt = np.linalg.svd(rr)
    keep = s > s[0] * (100 * _EPS)
    gamma = np.zeros(k)
    gamma[keep] = (uu[:, keep].T @ qty) / s[keep]
    coef = (vt.T @ gamma) / colnor
    return float(coef[0]), float(coef[1])


def loess_interpolate_1d(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    weights: npt.ArrayLike,
    newx: npt.ArrayLike,
    *,
    span: float = 0.75,
    degree: int = 2,
    cell: float = 0.2,
) -> FloatArray:
    """Fit a 1-D loess and evaluate its interpolated surface (see the module docstring).

    Parameters
    ----------
    x, y, weights : array-like of shape (n,)
        Predictor, response and prior weights (finite).
    newx : array-like
        Points to evaluate; they must lie within ``[min(x), max(x)]`` (R's ``predict`` returns
        NA outside, callers handle that).
    span, degree, cell : float, int, float
        As in R's ``loess`` / ``loess.control`` (``degree`` 1 or 2).

    Returns
    -------
    numpy.ndarray
        The fitted values at ``newx``.

    Raises
    ------
    ValueError
        If the span is too small or the fit is not finite (R's ``loess`` fails as well).
    """
    xv = np.asarray(x, dtype=np.float64)
    yv = np.asarray(y, dtype=np.float64)
    wv = np.asarray(weights, dtype=np.float64)
    n = xv.size
    nf = min(n, int(np.floor(n * span + 1e-5)))
    if nf <= 0:
        raise ValueError("span is too small")
    verts = _vertices(xv, int(np.floor(n * span * cell)))
    vals = np.array([_local_fit(xv, yv, wv, float(v), nf, span, degree) for v in verts])
    z = np.atleast_1d(np.asarray(newx, dtype=np.float64))
    c = np.clip(np.searchsorted(verts, z, side="left"), 1, verts.size - 1)  # z in (v[c-1], v[c]]
    a, b = verts[c - 1], verts[c]
    width = b - a
    h = (z - a) / width
    phi0 = (1 - h) ** 2 * (1 + 2 * h)
    phi1 = h**2 * (3 - 2 * h)
    psi0 = h * (1 - h) ** 2
    psi1 = -(h**2) * (1 - h)
    out = (
        phi0 * vals[c - 1, 0]
        + phi1 * vals[c, 0]
        + (psi0 * vals[c - 1, 1] + psi1 * vals[c, 1]) * width
    )
    if not np.all(np.isfinite(out)):
        raise ValueError("The loess fit is not finite.")
    return np.asarray(out, dtype=np.float64)
