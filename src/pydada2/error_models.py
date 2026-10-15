"""Error models: estimate per-quality substitution rates from transition counts.

R: ``loessErrfun``, ``noqualErrfun``, ``makeBinnedQualErrfun``, ``PacBioErrfun``, ``inflateErr``,
``getErrors`` and ``learnErrors``.

Error matrices and transition-count matrices have shape ``(16, n_qual)``: rows are the
transitions in :data:`TRANSITIONS` order (A2A, A2C, ..., T2T) and columns the quality scores
0, 1, ..., n_qual - 1.
"""

from __future__ import annotations

import itertools
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from ._base import BaseEstimator
from ._parallel import effective_n_jobs
from ._validation import check_is_fitted
from .io import QualityType

__all__ = [
    "TRANSITIONS",
    "BinnedQualErrorModel",
    "ErrorLearner",
    "ErrorModel",
    "LoessErrorModel",
    "NoQualErrorModel",
    "PacBioErrorModel",
    "accumulate_trans",
    "get_errors",
    "inflate_err",
    "learn_errors",
]

#: Row labels of error and transition matrices.
TRANSITIONS: tuple[str, ...] = tuple(f"{a}2{b}" for a in "ACGT" for b in "ACGT")

_MAX_ERROR_RATE = 0.25
_MIN_ERROR_RATE = 1e-7

FloatArray = npt.NDArray[np.float64]


def _as_trans(trans: Any) -> tuple[FloatArray, FloatArray]:
    """Return (16 x Q float matrix, quality scores) from an array or DataFrame."""
    if isinstance(trans, pd.DataFrame):
        q = np.asarray([float(c) for c in trans.columns], dtype=np.float64)
        mat = trans.to_numpy(dtype=np.float64)
    else:
        mat = np.asarray(trans, dtype=np.float64)
        q = np.arange(mat.shape[1] if mat.ndim == 2 else 0, dtype=np.float64)
    if mat.ndim != 2 or mat.shape[0] != 16:
        raise ValueError("The transition matrix must have 16 rows (A2A, A2C, ..., T2T).")
    return mat, q


def _expand(est: FloatArray) -> FloatArray:
    """12 substitution rows (A2C, A2G, A2T, C2A, ...) -> 16 rows with self-transitions."""
    e = est
    return np.vstack(
        [
            1 - e[0:3].sum(axis=0),
            e[0:3],
            e[3:4],
            1 - e[3:6].sum(axis=0),
            e[4:6],
            e[6:8],
            1 - e[6:9].sum(axis=0),
            e[8:9],
            e[9:12],
            1 - e[9:12].sum(axis=0),
        ]
    )


def _substitutions() -> list[tuple[int, int]]:
    return [(i, j) for i in range(4) for j in range(4) if i != j]


def _clamp(est: FloatArray) -> FloatArray:
    est = est.copy()
    est[est > _MAX_ERROR_RATE] = _MAX_ERROR_RATE
    est[est < _MIN_ERROR_RATE] = _MIN_ERROR_RATE
    return est


def _fill_edges(pred: FloatArray) -> FloatArray:
    ok = np.nonzero(~np.isnan(pred))[0]
    if ok.size == 0:
        raise ValueError("Error rates could not be estimated (no usable transitions).")
    pred = pred.copy()
    pred[: ok[0]] = pred[ok[0]]
    pred[ok[-1] + 1 :] = pred[ok[-1]]
    return pred


def _loess_fit(trans: FloatArray, q: FloatArray) -> FloatArray:
    from skmisc.loess import loess  # same netlib dloess code as R's stats::loess

    est = np.empty((12, trans.shape[1]), dtype=np.float64)
    for k, (i, j) in enumerate(_substitutions()):
        errs = trans[4 * i + j]
        tot = trans[4 * i : 4 * i + 4].sum(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            rlogp = np.log10((errs + 1) / tot)  # 1 pseudocount per error; tot=0 -> NA (dropped)
        ok = np.isfinite(rlogp)
        if ok.sum() < 3:
            raise ValueError("Error rates could not be estimated (too few observations).")
        lo = loess(
            q[ok],
            rlogp[ok],
            weights=tot[ok],
            span=0.75,
            degree=2,
            family="gaussian",
            surface="interpolate",
            statistics="approximate",
            trace_hat="exact",
            cell=0.2,
            iterations=4,
            normalize=True,
        )
        lo.fit()
        pred = np.full(q.size, np.nan)
        inside = (q >= q[ok].min()) & (q <= q[ok].max())  # R's predict() gives NA outside
        pred[inside] = lo.predict(q[inside]).values  # noqa: PD011 (skmisc object, not pandas)
        est[k] = 10 ** _fill_edges(pred)
    return _expand(_clamp(est))


class ErrorModel(BaseEstimator):
    """Base class of error models.

    An error model is fitted on a ``(16, n_qual)`` matrix of transition counts and estimates
    the error matrix ``err_``. Instances are also callable: ``model(trans)`` returns the
    estimated error matrix without modifying the model.
    """

    def _estimate(self, trans: FloatArray, q: FloatArray) -> FloatArray:
        raise NotImplementedError

    def fit(self, trans: Any) -> ErrorModel:
        """Estimate the error matrix from transition counts.

        Parameters
        ----------
        trans : array-like of shape (16, n_qual) or pandas.DataFrame
            Transition counts (rows A2A..T2T, columns quality scores).

        Returns
        -------
        self
            The fitted model, with the estimate in ``err_``.
        """
        mat, q = _as_trans(trans)
        self.err_ = self._estimate(mat, q)
        return self

    def __call__(self, trans: Any) -> FloatArray:
        """Estimate and return the error matrix (does not modify the model)."""
        mat, q = _as_trans(trans)
        return self._estimate(mat, q)

    def to_frame(self) -> pd.DataFrame:
        """Return the fitted error matrix as a DataFrame (rows: transitions, columns: quality)."""
        check_is_fitted(self, ["err_"])
        return pd.DataFrame(self.err_, index=list(TRANSITIONS))


class LoessErrorModel(ErrorModel):
    """Loess-smoothed error rates as a function of quality (R: ``loessErrfun``; the default).

    For each substitution type, ``log10((errors + 1) / total)`` is smoothed against the quality
    score with ``loess`` (span 0.75, degree 2, weights = total) and back-transformed. Rates are
    clamped to [1e-7, 0.25].

    Examples
    --------
    >>> import numpy as np
    >>> trans = np.zeros((16, 41)); trans[[0, 5, 10, 15]] = 1e4; trans[[1, 6, 11, 12]] = 10
    >>> LoessErrorModel().fit(trans).err_.shape
    (16, 41)
    """

    def __init__(self) -> None:
        pass

    def _estimate(self, trans: FloatArray, q: FloatArray) -> FloatArray:
        return _loess_fit(trans, q)


class NoQualErrorModel(ErrorModel):
    """Quality-independent error rates (R: ``noqualErrfun``).

    Parameters
    ----------
    pseudocount : float, default=1
        Added to every transition count before estimating the rates.
    """

    def __init__(self, pseudocount: float = 1) -> None:
        self.pseudocount = pseudocount

    def _estimate(self, trans: FloatArray, q: FloatArray) -> FloatArray:
        obs = trans.sum(axis=1) + self.pseudocount
        est = np.empty((12, trans.shape[1]), dtype=np.float64)
        for k, (i, j) in enumerate(_substitutions()):
            est[k] = obs[4 * i + j] / obs[4 * i : 4 * i + 4].sum()
        return _expand(est)


def _r_seq(lo: float, hi: float, n: int) -> FloatArray:
    """R's ``seq(from, to, length.out = n)``."""
    if n == 1:
        return np.array([lo])
    if n == 2:
        return np.array([lo, hi])
    if lo == hi:
        return np.full(n, lo)
    n1 = n - 1
    return np.concatenate([[lo], lo + np.arange(1, n1) * ((hi - lo) / n1), [hi]])


class BinnedQualErrorModel(ErrorModel):
    """Error rates for binned quality scores (R: ``makeBinnedQualErrfun``).

    Rates are observed at the binned scores and linearly interpolated in between.

    Parameters
    ----------
    binned_quals : sequence of int
        The quality scores present in the (binned) data, e.g. ``(2, 11, 25, 37)``.
    """

    def __init__(self, binned_quals: Sequence[int] = ()) -> None:
        self.binned_quals = binned_quals

    def _estimate(self, trans: FloatArray, q: FloatArray) -> FloatArray:
        bq = [int(b) for b in self.binned_quals]
        if not bq:
            raise ValueError("The quality scores used in your data must be provided.")
        used = trans.sum(axis=0) > 0
        qmax, qmin = q[used].max(), q[used].min()
        if qmax > max(bq):
            raise ValueError(
                "Input data contains a higher quality score than the provided binned values."
            )
        if qmin < min(bq):
            raise ValueError(
                "Input data contains a lower quality score than the provided binned values."
            )
        if qmax not in bq:
            warnings.warn(
                "Maximum observed quality score is not in the provided binned values.", stacklevel=3
            )
        if qmin not in bq:
            warnings.warn(
                "Minimum observed quality score is not in the provided binned values.", stacklevel=3
            )
        if not np.array_equal(q, np.arange(q.size)):
            raise ValueError("Unexpected Q score series.")
        est = np.empty((12, trans.shape[1]), dtype=np.float64)
        for k, (i, j) in enumerate(_substitutions()):
            errs = trans[4 * i + j]
            tot = trans[4 * i : 4 * i + 4].sum(axis=0)
            with np.errstate(divide="ignore", invalid="ignore"):
                p = errs / tot
            pred = np.full(q.size, np.nan)
            for lo_q, hi_q in itertools.pairwise(bq):
                lo_p = p[lo_q] if lo_q < p.size else np.nan
                hi_p = p[hi_q] if hi_q < p.size else np.nan
                if not np.isnan(lo_p) and not np.isnan(hi_p):
                    pred[lo_q : hi_q + 1] = _r_seq(lo_p, hi_p, hi_q - lo_q + 1)
            est[k] = _fill_edges(pred)
        return _expand(_clamp(est))


class PacBioErrorModel(ErrorModel):
    """Error rates for PacBio CCS data (R: ``PacBioErrfun``).

    Loess-smoothed rates for qualities below 93, and a separate maximum-likelihood estimate for
    the special maximum score 93. Without a Q93 column this equals :class:`LoessErrorModel`.
    """

    def __init__(self) -> None:
        pass

    def _estimate(self, trans: FloatArray, q: FloatArray) -> FloatArray:
        if 93.0 not in q:
            return _loess_fit(trans, q)
        i93 = int(np.nonzero(q == 93.0)[0][0])
        if i93 != q.size - 1:
            raise ValueError("Max qual score of 93 not the last column as expected.")
        err = _loess_fit(trans[:, :i93], q[:i93])
        tot93 = np.repeat([trans[4 * i : 4 * i + 4, i93].sum() for i in range(4)], 4)
        err93 = (trans[:, i93] + 1) / (tot93 + 4)
        return np.hstack([err, err93[:, None]])


def accumulate_trans(trans: Sequence[Any]) -> FloatArray:
    """Sum transition matrices with possibly different numbers of quality columns.

    Parameters
    ----------
    trans : sequence of array-like of shape (16, n_qual_i)
        Transition counts.

    Returns
    -------
    numpy.ndarray
        Shape (16, max n_qual_i).
    """
    mats = [np.asarray(t, dtype=np.float64) for t in trans]
    if not mats:
        raise ValueError("No transition matrices to accumulate.")
    out = np.zeros((16, max(m.shape[1] for m in mats)), dtype=np.float64)
    for m in mats:
        out[:, : m.shape[1]] += m
    return out


def get_errors(obj: Any, /, *, detailed: bool = False, enforce: bool = True) -> Any:
    """Extract the error matrix of an object (R: ``getErrors``).

    Parameters
    ----------
    obj : object
        An error matrix, a fitted :class:`ErrorModel` or learner (``err_``), a fitted ``DADA``
        (``err_out_``), a DADA result (``err_out``), or a list of DADA results with identical
        error rates.
    detailed : bool, default=False
        Return a dict ``{"err_out", "err_in", "trans"}`` instead of the matrix.
    enforce : bool, default=True
        Validate the matrix: 16 rows, values in [0, 1] (a zero issues a warning).

    Returns
    -------
    numpy.ndarray or dict
        The error matrix, or the detailed dict.

    Raises
    ------
    ValueError
        If the error matrix is missing or invalid.
    """
    rval: dict[str, Any] = {"err_out": None, "err_in": None, "trans": None}
    if isinstance(obj, np.ndarray | pd.DataFrame):
        rval["err_out"] = np.asarray(obj, dtype=np.float64)
    elif isinstance(obj, dict) and {"err_out", "err_in", "trans"} <= set(obj):
        rval = dict(obj)
    elif isinstance(obj, list | tuple) and obj and all(hasattr(o, "err_out") for o in obj):
        outs = [o.err_out for o in obj]
        if not all(_same(outs[0], o) for o in outs[1:]):
            raise ValueError(
                "If a list of DADA results is provided, all must have the same output error rates."
            )
        rval["err_out"] = outs[0]
        rval["err_in"] = obj[0].err_in
        rval["trans"] = accumulate_trans([o.trans for o in obj])
    elif hasattr(obj, "err_out_"):  # fitted DADA
        rval["err_out"] = obj.err_out_
        rval["err_in"] = getattr(obj, "err_in_", None)
        rval["trans"] = getattr(obj, "trans_", None)
    elif hasattr(obj, "err_out"):
        rval["err_out"] = obj.err_out
        rval["err_in"] = getattr(obj, "err_in", None)
        rval["trans"] = getattr(obj, "trans", None)
    elif hasattr(obj, "err_"):
        rval["err_out"] = obj.err_
        rval["err_in"] = getattr(obj, "err_in_", None)
        rval["trans"] = getattr(obj, "trans_", None)
    elif obj is not None:
        raise TypeError(f"Cannot extract an error matrix from {type(obj).__name__}.")
    if enforce:
        err = rval["err_out"]
        if err is None:
            raise ValueError("Error matrix is None.")
        err = np.asarray(err, dtype=np.float64)
        if err.ndim != 2 or err.shape[0] != 16:
            raise ValueError("Error matrix must have 16 rows (A2A, A2C, ...).")
        if not (err >= 0).all():
            raise ValueError("All error matrix entries must be >= 0.")
        if not (err <= 1).all():
            raise ValueError("All error matrix entries must be <= 1.")
        if (err == 0).any():
            warnings.warn("Zero in error matrix.", stacklevel=2)
        rval["err_out"] = err
    return rval if detailed else rval["err_out"]


def _same(a: Any, b: Any) -> bool:
    if a is None or b is None:
        return a is b
    return bool(np.array_equal(np.asarray(a), np.asarray(b)))


def inflate_err(
    err: Any, inflation: float, *, inflate_self_transitions: bool = False
) -> FloatArray:
    """Inflate error rates (R: ``inflateErr``).

    Each substitution rate ``e`` becomes ``e * f / (1 + (f - 1) * e)``, which keeps rates <= 1.

    Parameters
    ----------
    err : object
        Anything accepted by :func:`get_errors`.
    inflation : float
        Inflation factor ``f``.
    inflate_self_transitions : bool, default=False
        Also inflate A2A, C2C, G2G and T2T.

    Returns
    -------
    numpy.ndarray
        The inflated (16, n_qual) error matrix.
    """
    e = np.array(get_errors(err), dtype=np.float64, copy=True)
    rows = [k for k, t in enumerate(TRANSITIONS) if t[0] != t[2]]
    if inflate_self_transitions:
        rows = list(range(16))
    e[rows] = (e[rows] * inflation) / (1 + (inflation - 1) * e[rows])
    return e


class ErrorLearner(BaseEstimator):
    """Learn error rates from data (R: ``learnErrors``).

    Samples are dereplicated in order until more than ``n_bases`` bases have been read; then
    :class:`~pydada2.DADA` is run in self-consistent mode (starting from maximal error rates)
    to estimate the error matrix.

    Parameters
    ----------
    n_bases : float, default=1e8
        Minimum number of bases to use.
    error_model : ErrorModel, optional
        Error model (default :class:`LoessErrorModel`).
    randomize : bool, default=False
        Pick samples in random order (seeded by ``random_state``).
    random_state : int or numpy.random.Generator, optional
        Seed for ``randomize``.
    max_consist : int, default=10
        Maximum number of self-consistency rounds.
    omega_c : float, default=0
        ``omega_c`` used while learning (R default for ``learnErrors``).
    params : DadaParams or dict, optional
        Other DADA parameters (e.g. ``{"band_size": 32}`` for PacBio).
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        Quality encoding of FASTQ inputs.
    n_jobs : int, optional
        Threads for the C++ core.
    verbose : int, default=0
        Logging level.

    Attributes
    ----------
    err_ : numpy.ndarray
        Learned error matrix (16 x n_qual).
    err_in_ : list of numpy.ndarray
        Error matrices used in each self-consistency round.
    trans_ : numpy.ndarray
        Transition counts of the final round.
    n_iter_ : int
        Number of rounds.
    converged_ : bool
        Whether the error rates converged.
    n_bases_, n_reads_, n_samples_ : int
        Amount of data used.
    dada_ : DADA
        The fitted self-consistent DADA estimator.
    """

    def __init__(
        self,
        n_bases: float = 1e8,
        error_model: ErrorModel | None = None,
        randomize: bool = False,
        random_state: Any = None,
        max_consist: int = 10,
        omega_c: float = 0,
        params: Any = None,
        quality_type: QualityType = "auto",
        n_jobs: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.n_bases = n_bases
        self.error_model = error_model
        self.randomize = randomize
        self.random_state = random_state
        self.max_consist = max_consist
        self.omega_c = omega_c
        self.params = params
        self.quality_type = quality_type
        self.n_jobs = n_jobs
        self.verbose = verbose

    def fit(self, X: Any, y: None = None) -> ErrorLearner:
        """Learn the error rates.

        Parameters
        ----------
        X : path, Dereplication, or sequence of them
            Samples (FASTQ files are dereplicated as needed).
        y : None
            Ignored.

        Returns
        -------
        self
            The fitted learner.
        """
        import logging
        import os

        from .denoise import DADA, DadaParams
        from .derep import Dereplication, derep_fastq

        log = logging.getLogger("pydada2")
        items = [X] if isinstance(X, str | os.PathLike | Dereplication) else list(X)
        if not items:
            raise ValueError("No samples provided.")
        if self.randomize:
            rng = np.random.default_rng(self.random_state)
            items = [items[i] for i in rng.permutation(len(items))]
        dereps = []
        nbases = nreads = 0
        for item in items:
            drp = (
                item
                if isinstance(item, Dereplication)
                else derep_fastq(item, quality_type=self.quality_type)
            )
            dereps.append(drp)
            ab = drp.uniques.to_numpy(dtype=np.float64)
            nreads += int(ab.sum())
            nbases += int(
                (ab * np.fromiter((len(s) for s in drp.uniques.index), float, len(ab))).sum()
            )
            if nbases > self.n_bases:
                break
        if self.verbose:
            log.info(
                "%d total bases in %d reads from %d samples will be used for learning the error "
                "rates.",
                nbases,
                nreads,
                len(dereps),
            )
        if self.params is None:
            base = DadaParams()
        elif isinstance(self.params, DadaParams):
            base = self.params
        else:
            base = DadaParams().replace(**dict(self.params))
        params = base.replace(max_consist=int(self.max_consist), omega_c=float(self.omega_c))
        dd = DADA(
            err=None,
            error_model=self.error_model,
            self_consistent=True,
            params=params,
            n_jobs=effective_n_jobs(self.n_jobs),
            verbose=self.verbose,
        ).fit(dereps)
        self.dada_ = dd
        self.err_ = dd.err_out_
        self.err_in_ = dd.err_in_
        self.trans_ = dd.trans_
        self.n_iter_ = dd.n_iter_
        self.converged_ = dd.converged_
        self.n_bases_ = nbases
        self.n_reads_ = nreads
        self.n_samples_ = len(dereps)
        return self


def learn_errors(X: Any, /, **kwargs: Any) -> ErrorLearner:
    """Learn error rates from data (R: ``learnErrors``).

    Convenience wrapper for ``ErrorLearner(**kwargs).fit(X)``.

    Parameters
    ----------
    X : path, Dereplication, or sequence of them
        Samples.
    **kwargs
        Parameters of :class:`ErrorLearner`.

    Returns
    -------
    ErrorLearner
        Fitted learner; the error matrix is ``.err_``.
    """
    return ErrorLearner(**kwargs).fit(X)
