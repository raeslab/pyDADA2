"""Chimera (bimera) detection and removal.

R: ``isBimera``, ``isBimeraDenovo``, ``isBimeraDenovoTable``, ``removeBimeraDenovo``,
``isShiftDenovo``.
"""

from __future__ import annotations

import logging
import warnings
from collections.abc import Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd

from . import _core
from ._base import BaseEstimator
from ._parallel import effective_n_jobs
from ._validation import check_is_fitted
from .denoise import DadaResult
from .derep import Dereplication
from .sequences import get_sequences, get_uniques

__all__ = [
    "BimeraRemover",
    "is_bimera",
    "is_bimera_denovo",
    "is_bimera_denovo_table",
    "is_shift_denovo",
    "remove_bimera_denovo",
]

logger = logging.getLogger("pydada2")

# Alignment scores: DADA2's default options (getDadaOpt MATCH, MISMATCH, GAP_PENALTY).
_MATCH, _MISMATCH, _GAP = 5, -4, -8

Method = Literal["consensus", "pooled", "per-sample"]


def is_bimera(
    sq: str,
    parents: Sequence[str],
    /,
    *,
    allow_one_off: bool = False,
    min_one_off_parent_distance: int = 4,
    max_shift: int = 16,
) -> bool:
    """Whether a sequence is an exact bimera of two of the given parents (R: ``isBimera``).

    Parameters
    ----------
    sq : str
        Query sequence.
    parents : sequence of str
        Candidate parent sequences.
    allow_one_off : bool, default=False
        Also flag sequences one mismatch/indel away from an exact bimera.
    min_one_off_parent_distance : int, default=4
        Minimum distance of parents to the query to be considered for one-off bimeras.
    max_shift : int, default=16
        Maximum shift allowed when aligning to the parents.

    Returns
    -------
    bool
        True if ``sq`` is a bimera.

    Examples
    --------
    >>> a, b = "ACGTACGTAC" * 3, "TTGGCCAATT" * 3
    >>> is_bimera(a[:15] + b[15:], [a, b])
    True
    """
    return bool(
        _core.is_bimera(
            sq,
            list(parents),
            allow_one_off,
            int(min_one_off_parent_distance),
            _MATCH,
            _MISMATCH,
            _GAP,
            int(max_shift),
        )
    )


def is_bimera_denovo(
    unqs: Any,
    /,
    *,
    min_fold_parent_over_abundance: float = 2,
    min_parent_abundance: float = 8,
    allow_one_off: bool = False,
    min_one_off_parent_distance: int = 4,
    max_shift: int = 16,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.Series:
    """Identify bimeras among a set of sequences (R: ``isBimeraDenovo``).

    Each sequence is tested against the more abundant sequences (abundance >
    ``min_fold_parent_over_abundance`` times its own and > ``min_parent_abundance``).

    Parameters
    ----------
    unqs : object
        Anything accepted by :func:`~pydada2.get_uniques` (a sequence table is pooled).
    min_fold_parent_over_abundance : float, default=2
        Minimum fold over-abundance of parents.
    min_parent_abundance : float, default=8
        Minimum abundance of parents.
    allow_one_off, min_one_off_parent_distance, max_shift
        See :func:`is_bimera`.
    n_jobs : int, optional
        Threads.
    verbose : int, default=0
        Log a summary when > 0.

    Returns
    -------
    pandas.Series
        Boolean flag per input sequence (indexed by sequence; duplicates kept).
    """
    seqs_input = get_sequences(unqs)
    u = get_uniques(unqs)
    if len(seqs_input) != len(set(seqs_input)):
        logger.info("Duplicate sequences detected.")
    flags = np.asarray(
        _core.bimera_denovo(
            list(u.index),
            u.to_numpy(dtype=np.float64).tolist(),
            float(min_fold_parent_over_abundance),
            float(min_parent_abundance),
            bool(allow_one_off),
            int(min_one_off_parent_distance),
            _MATCH,
            _MISMATCH,
            _GAP,
            int(max_shift),
            effective_n_jobs(n_jobs),
        )
    )
    bims = {s for s, f in zip(u.index, flags, strict=True) if f}
    out = pd.Series([s in bims for s in seqs_input], index=pd.Index(seqs_input, dtype=object))
    if verbose:
        logger.info("Identified %d bimeras out of %d input sequences.", int(out.sum()), len(out))
    return out


def is_bimera_denovo_table(
    seqtab: pd.DataFrame,
    /,
    *,
    min_sample_fraction: float = 0.9,
    ignore_n_negatives: int = 1,
    min_fold_parent_over_abundance: float = 1.5,
    min_parent_abundance: float = 2,
    allow_one_off: bool = False,
    min_one_off_parent_distance: int = 4,
    max_shift: int = 16,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.Series:
    """Identify bimeras in a sequence table by sample consensus (R: ``isBimeraDenovoTable``).

    Sequences are tested in every sample in which they occur; a sequence is a bimera if it is
    flagged in all samples, or in a fraction of at least ``min_sample_fraction`` of them
    (ignoring up to ``ignore_n_negatives`` samples that do not flag it).

    Parameters
    ----------
    seqtab : pandas.DataFrame
        Samples x sequences table of non-negative integer counts.
    min_sample_fraction : float, default=0.9
        Required fraction of samples flagging a sequence.
    ignore_n_negatives : int, default=1
        Number of non-flagging samples to ignore.
    min_fold_parent_over_abundance : float, default=1.5
        Minimum fold over-abundance of parents (within a sample).
    min_parent_abundance : float, default=2
        Minimum abundance of parents (within a sample).
    allow_one_off, min_one_off_parent_distance, max_shift
        See :func:`is_bimera`.
    n_jobs : int, optional
        Threads.
    verbose : int, default=0
        Log a summary when > 0.

    Returns
    -------
    pandas.Series
        Boolean flag per sequence (column).

    Raises
    ------
    ValueError
        For duplicated sequences or non-integer tables.
    """
    if not isinstance(seqtab, pd.DataFrame):
        raise ValueError("Input must be a valid sequence table (a DataFrame).")
    sqs = [str(c) for c in seqtab.columns]
    if len(set(sqs)) != len(sqs):
        raise ValueError("Duplicate sequences detected in input.")
    vals = seqtab.to_numpy()
    if not np.issubdtype(vals.dtype, np.integer):
        raise ValueError("Input must be a valid sequence table (integer counts).")
    if vals.size and vals.max() > np.iinfo(np.int32).max:
        raise ValueError("Counts exceed the supported range (int32).")
    nflag, nsam = _core.table_bimera(
        np.ascontiguousarray(vals, dtype=np.int32),
        sqs,
        float(min_fold_parent_over_abundance),
        int(min_parent_abundance),
        bool(allow_one_off),
        int(min_one_off_parent_distance),
        _MATCH,
        _MISMATCH,
        _GAP,
        int(max_shift),
        effective_n_jobs(n_jobs),
    )
    nflag = np.asarray(nflag)
    nsam = np.asarray(nsam)
    bims = (nflag >= nsam) | (
        (nflag > 0) & (nflag >= (nsam - ignore_n_negatives) * min_sample_fraction)
    )
    out = pd.Series(bims, index=pd.Index(sqs, dtype=object))
    if verbose:
        logger.info("Identified %d bimeras out of %d input sequences.", int(out.sum()), len(out))
    return out


def is_shift_denovo(
    unqs: Any,
    /,
    *,
    min_overlap: int = 20,
    flag_subseqs: bool = False,
    n_jobs: int | None = None,
) -> pd.Series:
    """Identify sequences that are shifted copies of more abundant ones (R: ``isShiftDenovo``).

    Parameters
    ----------
    unqs : object
        Anything accepted by :func:`~pydada2.get_uniques`.
    min_overlap : int, default=20
        Minimum overlap to call a shift.
    flag_subseqs : bool, default=False
        Has no effect, as in R (``isShiftDenovo`` never passes ``flagSubseqs`` on); a warning is
        issued when True. Sequences contained in a more abundant one are not flagged.
    n_jobs : int, optional
        Threads.

    Returns
    -------
    pandas.Series
        Boolean flag per input sequence.
    """
    if flag_subseqs:
        warnings.warn(
            "flag_subseqs has no effect (as in R's isShiftDenovo); contained sequences are not "
            "flagged.",
            stacklevel=2,
        )
    seqs_input = get_sequences(unqs)
    u = get_uniques(unqs)
    flags = np.asarray(
        _core.is_shift_denovo(
            list(u.index),
            u.to_numpy(dtype=np.int64).tolist(),
            int(min_overlap),
            False,
            effective_n_jobs(n_jobs),
        )
    )
    shifts = {s for s, f in zip(u.index, flags, strict=True) if f}
    return pd.Series([s in shifts for s in seqs_input], index=pd.Index(seqs_input, dtype=object))


_DEFAULTS: dict[str, dict[str, float]] = {
    "consensus": {"min_fold_parent_over_abundance": 1.5, "min_parent_abundance": 2},
    "pooled": {"min_fold_parent_over_abundance": 2, "min_parent_abundance": 8},
    "per-sample": {"min_fold_parent_over_abundance": 2, "min_parent_abundance": 8},
}


class BimeraRemover(BaseEstimator):
    """Remove bimeras from a sequence table (R: ``removeBimeraDenovo``).

    Parameters
    ----------
    method : {"consensus", "pooled", "per-sample"}, default="consensus"
        ``"consensus"``: flag per sample and remove by consensus (:func:`is_bimera_denovo_table`);
        ``"pooled"``: pool all samples (:func:`is_bimera_denovo`); ``"per-sample"``: remove
        bimeras independently in each sample (counts set to 0).
    min_fold_parent_over_abundance : float, optional
        Default 1.5 for ``"consensus"``, 2 otherwise (as in R).
    min_parent_abundance : float, optional
        Default 2 for ``"consensus"``, 8 otherwise (as in R).
    min_sample_fraction : float, default=0.9
        ``"consensus"`` only.
    ignore_n_negatives : int, default=1
        ``"consensus"`` only.
    allow_one_off : bool, default=False
        See :func:`is_bimera`.
    min_one_off_parent_distance : int, default=4
        See :func:`is_bimera`.
    max_shift : int, default=16
        See :func:`is_bimera`.
    n_jobs : int, optional
        Threads.
    verbose : int, default=0
        Log summaries when > 0.

    Attributes
    ----------
    is_bimera_ : pandas.Series or pandas.DataFrame
        Bimera flags per sequence (per sample and sequence for ``"per-sample"``, where
        sequences absent from a sample are not tested and flagged False).

    Examples
    --------
    >>> import pandas as pd
    >>> a, b = "ACGTACGTAC" * 3, "TTGGCCAATT" * 3
    >>> tab = pd.DataFrame({a: [100, 80], b: [90, 70], a[:15] + b[15:]: [5, 4]})
    >>> BimeraRemover().fit_transform(tab).shape
    (2, 2)
    """

    def __init__(
        self,
        method: Method = "consensus",
        min_fold_parent_over_abundance: float | None = None,
        min_parent_abundance: float | None = None,
        min_sample_fraction: float = 0.9,
        ignore_n_negatives: int = 1,
        allow_one_off: bool = False,
        min_one_off_parent_distance: int = 4,
        max_shift: int = 16,
        n_jobs: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.method = method
        self.min_fold_parent_over_abundance = min_fold_parent_over_abundance
        self.min_parent_abundance = min_parent_abundance
        self.min_sample_fraction = min_sample_fraction
        self.ignore_n_negatives = ignore_n_negatives
        self.allow_one_off = allow_one_off
        self.min_one_off_parent_distance = min_one_off_parent_distance
        self.max_shift = max_shift
        self.n_jobs = n_jobs
        self.verbose = verbose

    def _common(self) -> dict[str, Any]:
        if self.method not in _DEFAULTS:
            raise ValueError("Valid values for method: 'pooled', 'consensus', or 'per-sample'.")
        d = _DEFAULTS[self.method]
        return {
            "min_fold_parent_over_abundance": d["min_fold_parent_over_abundance"]
            if self.min_fold_parent_over_abundance is None
            else self.min_fold_parent_over_abundance,
            "min_parent_abundance": d["min_parent_abundance"]
            if self.min_parent_abundance is None
            else self.min_parent_abundance,
            "allow_one_off": self.allow_one_off,
            "min_one_off_parent_distance": self.min_one_off_parent_distance,
            "max_shift": self.max_shift,
            "n_jobs": self.n_jobs,
        }

    def fit(self, seqtab: pd.DataFrame, y: None = None) -> BimeraRemover:
        """Identify the bimeras of a sequence table.

        Parameters
        ----------
        seqtab : pandas.DataFrame
            Samples x sequences table.
        y : None
            Ignored.

        Returns
        -------
        self
            The fitted estimator (flags in ``is_bimera_``).
        """
        kw = self._common()
        self.is_bimera_: pd.Series | pd.DataFrame
        if not isinstance(seqtab, pd.DataFrame):
            raise TypeError("seqtab must be a pandas DataFrame (samples x sequences).")
        if self.method == "consensus":
            self.is_bimera_ = is_bimera_denovo_table(
                seqtab,
                min_sample_fraction=self.min_sample_fraction,
                ignore_n_negatives=self.ignore_n_negatives,
                verbose=self.verbose,
                **kw,
            )
        elif self.method == "pooled":
            self.is_bimera_ = is_bimera_denovo(seqtab, verbose=self.verbose, **kw)
        else:
            # Only the sequences present in a sample are tested: absent ones cannot be parents
            # (abundance 0) and their counts stay 0, so the transformed table equals R's.
            cols = [str(c) for c in seqtab.columns]
            vals = seqtab.to_numpy(dtype=np.int64)
            flags = np.zeros(vals.shape, dtype=bool)
            for i in range(vals.shape[0]):
                present = np.flatnonzero(vals[i])
                if present.size:
                    u = pd.Series(vals[i, present], index=[cols[j] for j in present])
                    flags[i, present] = is_bimera_denovo(u, verbose=self.verbose, **kw).to_numpy()
            self.is_bimera_ = pd.DataFrame(flags, index=seqtab.index, columns=seqtab.columns)
        return self

    def transform(self, seqtab: pd.DataFrame) -> pd.DataFrame:
        """Remove the identified bimeras.

        Parameters
        ----------
        seqtab : pandas.DataFrame
            The table passed to :meth:`fit`.

        Returns
        -------
        pandas.DataFrame
            The table without bimeras (``"per-sample"``: bimeric counts set to 0 and all-zero
            columns dropped).
        """
        check_is_fitted(self, ["is_bimera_"])
        if isinstance(self.is_bimera_, pd.DataFrame):
            out: pd.DataFrame = seqtab.mask(self.is_bimera_.to_numpy(), 0)
            kept: pd.DataFrame = out.loc[:, out.sum(axis=0).to_numpy() != 0]
            return kept
        return seqtab.loc[:, ~self.is_bimera_.to_numpy()]

    def fit_transform(self, seqtab: pd.DataFrame, y: None = None) -> pd.DataFrame:
        """Fit and remove the bimeras in one step."""
        return self.fit(seqtab).transform(seqtab)


def remove_bimera_denovo(
    unqs: Any, /, *, method: Method = "consensus", verbose: int = 0, **kwargs: Any
) -> Any:
    """Remove bimeras (R: ``removeBimeraDenovo``).

    Parameters
    ----------
    unqs : object
        A sequence table (DataFrame of samples x sequences; uses ``method``), a uniques
        Series / mapping, a :class:`~pydada2.DadaResult`, a :class:`~pydada2.Dereplication`,
        a merged-pairs DataFrame, or a list of those.
    method : {"consensus", "pooled", "per-sample"}, default="consensus"
        For sequence tables, see :class:`BimeraRemover`.
    verbose : int, default=0
        Log summaries when > 0.
    **kwargs
        Passed to :class:`BimeraRemover` (tables) or :func:`is_bimera_denovo` (others).

    Returns
    -------
    object
        Same type as the input with bimeras removed (Series for uniques, DADA results and
        dereplications; a list for a list input).
    """
    if isinstance(unqs, list):
        return [remove_bimera_denovo(u, method=method, verbose=verbose, **kwargs) for u in unqs]
    if isinstance(unqs, pd.DataFrame) and not {"sequence", "abundance"} <= set(unqs.columns):
        return BimeraRemover(method=method, verbose=verbose, **kwargs).fit_transform(unqs)
    bim = is_bimera_denovo(unqs, verbose=verbose, **kwargs).to_numpy()
    if isinstance(unqs, pd.DataFrame):
        return unqs.loc[~bim].reset_index(drop=True)
    if isinstance(unqs, DadaResult):
        return unqs.denoised.loc[~bim]
    if isinstance(unqs, Dereplication):
        return unqs.uniques.loc[~bim]
    return get_uniques(unqs, collapse=False).loc[~bim]
