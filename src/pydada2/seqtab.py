"""Sequence tables (R: ``makeSequenceTable``, ``mergeSequenceTables``, ``collapseNoMismatch``).

A sequence table is a :class:`pandas.DataFrame` of int64 abundances with samples as rows and
sequences (ASVs) as columns.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping
from typing import Any, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import _core
from ._parallel import effective_n_jobs
from .denoise import DADA
from .sequences import get_uniques, reverse_complement

__all__ = ["collapse_no_mismatch", "make_sequence_table", "merge_sequence_tables"]

logger = logging.getLogger("pydada2")

OrderBy = Literal["abundance", "nsamples"] | None


def _order_columns(tab: pd.DataFrame, order_by: OrderBy) -> pd.DataFrame:
    if order_by is None:
        return tab
    if order_by == "abundance":
        key = tab.sum(axis=0).to_numpy()
    elif order_by == "nsamples":
        key = (tab > 0).sum(axis=0).to_numpy()
    else:
        raise ValueError("order_by must be 'abundance', 'nsamples' or None.")
    return tab.iloc[:, np.argsort(-key, kind="stable")]


def _uniques_arrays(obj: Any) -> tuple[list[str], npt.NDArray[np.int64]]:
    """Return ``get_uniques(obj)`` as (sequences, abundances).

    Avoids the per-sample pandas overhead for the common inputs (merged pairs, DADA results);
    the results are identical.
    """
    seqs: list[str] | None = None
    if isinstance(obj, pd.DataFrame) and {"sequence", "abundance"} <= set(obj.columns):
        seqs, ab = [str(s) for s in obj["sequence"].tolist()], obj["abundance"]
    elif isinstance(getattr(obj, "denoised", None), pd.Series):  # DadaResult
        seqs, ab = [str(s) for s in obj.denoised.index], obj.denoised
    if seqs is not None and len(set(seqs)) == len(seqs):  # duplicates: get_uniques merges them
        return seqs, ab.to_numpy(dtype=np.int64)
    u = get_uniques(obj)
    return [str(s) for s in u.index], u.to_numpy(dtype=np.int64)


def make_sequence_table(samples: Any, /, *, order_by: OrderBy = "abundance") -> pd.DataFrame:
    """Build a sequence table from samples (R: ``makeSequenceTable``).

    Parameters
    ----------
    samples : object
        A fitted :class:`~pydada2.DADA`, a list of samples, or a dict ``{name: sample}``. A
        sample is anything accepted by :func:`~pydada2.get_uniques` (DADA results, merged-pair
        frames, dereplications, ...).
    order_by : {"abundance", "nsamples", None}, default="abundance"
        Column order: decreasing total abundance, decreasing number of samples, or order of
        first appearance.

    Returns
    -------
    pandas.DataFrame
        Samples x sequences table of int64 abundances. The index holds the sample names (dict
        keys, or the names of a fitted DADA), otherwise 0, 1, ...
    """
    names: list[str] | None = None
    if isinstance(samples, DADA):
        items = list(samples.results_)
        names = samples.sample_names_
    elif isinstance(samples, Mapping):
        names = [str(k) for k in samples]
        items = list(samples.values())
    elif isinstance(samples, list | tuple):
        items = list(samples)
    else:
        items = [samples]
    unqs = [_uniques_arrays(s) for s in items]
    cols: dict[str, int] = {}
    for seqs, _ in unqs:
        for sq in seqs:
            cols.setdefault(sq, len(cols))
    mat = np.zeros((len(unqs), len(cols)), dtype=np.int64)
    for i, (seqs, ab) in enumerate(unqs):
        mat[i, [cols[sq] for sq in seqs]] = ab
    colnames = list(cols)
    if order_by is not None:  # as _order_columns, on the array (no pandas reindexing)
        if order_by == "abundance":
            key = mat.sum(axis=0)
        elif order_by == "nsamples":
            key = (mat > 0).sum(axis=0)
        else:
            raise ValueError("order_by must be 'abundance', 'nsamples' or None.")
        order = np.argsort(-key, kind="stable")
        mat = mat[:, order]
        colnames = [colnames[j] for j in order]
    return pd.DataFrame(mat, columns=colnames, index=pd.Index(names) if names is not None else None)


def _check_table(tab: Any, i: int) -> pd.DataFrame:
    if not isinstance(tab, pd.DataFrame):
        raise TypeError(f"Sequence table {i} is not a DataFrame.")
    if (tab.to_numpy() < 0).any():
        raise ValueError(f"Sequence table {i} has negative entries.")
    if any(len(str(c)) == 0 for c in tab.columns) or any(len(str(r)) == 0 for r in tab.index):
        raise ValueError(f"Sequence table {i} has blank sequence or sample names.")
    return tab


def merge_sequence_tables(
    *tables: pd.DataFrame,
    repeats: Literal["error", "sum"] = "error",
    order_by: OrderBy = "abundance",
    try_rc: bool = False,
) -> pd.DataFrame:
    """Merge sequence tables (R: ``mergeSequenceTables``).

    Parameters
    ----------
    *tables : pandas.DataFrame
        Two or more sequence tables.
    repeats : {"error", "sum"}, default="error"
        What to do with samples present in more than one table: raise, or sum their counts.
    order_by : {"abundance", "nsamples", None}, default="abundance"
        Column order of the result.
    try_rc : bool, default=False
        Re-orient sequences whose reverse complement occurs earlier in the tables.

    Returns
    -------
    pandas.DataFrame
        The merged table.

    Raises
    ------
    ValueError
        With fewer than two tables, duplicated sample names and ``repeats="error"``, or
        invalid tables.
    """
    if len(tables) == 1 and isinstance(tables[0], list | tuple):
        tables = tuple(tables[0])
    if len(tables) < 2:
        raise ValueError("At least two sequence tables are expected.")
    tabs = [_check_table(t, i) for i, t in enumerate(tables)]
    sample_names = [str(r) for t in tabs for r in t.index]
    dup = pd.Index(sample_names).duplicated()
    if dup.any():
        if repeats == "error":
            names = sorted({n for n, d in zip(sample_names, dup, strict=True) if d})
            raise ValueError(
                "Duplicated sample names detected in the sequence table row names: "
                + ", ".join(names)
            )
        logger.info("Duplicated sample names detected in the sequence table row names.")
        sample_names = list(dict.fromkeys(sample_names))
    seqs = list(dict.fromkeys(str(c) for t in tabs for c in t.columns))
    if try_rc and len(seqs) > 1:
        seen: set[str] = {seqs[0]}
        rc_seqs: set[str] = set()
        for s in seqs[1:]:
            if reverse_complement(s) in seen:
                rc_seqs.add(s)
            seen.add(s)
        if rc_seqs:
            logger.info("Reverse complemented sequences detected and re-oriented.")
            new = []
            for t in tabs:
                do_rc = [str(c) in rc_seqs for c in t.columns]
                if any(do_rc):
                    t = t.copy()
                    t.columns = [
                        reverse_complement(str(c)) if r else str(c)
                        for c, r in zip(t.columns, do_rc, strict=True)
                    ]
                    t = collapse_no_mismatch(t, identical_only=True)
                new.append(t)
            tabs = new
            seqs = [s for s in seqs if s not in rc_seqs]
    rval = pd.DataFrame(0, index=pd.Index(sample_names), columns=seqs, dtype=np.int64)
    for t in tabs:
        rows = [str(r) for r in t.index]
        cols = [str(c) for c in t.columns]
        rval.loc[rows, cols] = rval.loc[rows, cols].to_numpy() + t.to_numpy(dtype=np.int64)
    return _order_columns(rval.astype(np.int64), order_by)


def collapse_no_mismatch(
    seqtab: pd.DataFrame,
    /,
    *,
    min_overlap: int = 20,
    order_by: OrderBy = "abundance",
    identical_only: bool = False,
    vec: bool = True,
    band: int = -1,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.DataFrame:
    """Collapse sequences that are identical up to shifts or length (R: ``collapseNoMismatch``).

    Sequences are processed in order of decreasing total abundance; each is merged into the
    first more abundant sequence it matches with no mismatches or internal indels (end gaps
    allowed), otherwise it is kept.

    Parameters
    ----------
    seqtab : pandas.DataFrame
        Sequence table.
    min_overlap : int, default=20
        Length of the prefix used for the cheap pre-screen (and minimum overlap).
    order_by : {"abundance", "nsamples", None}
        Column order before the final ordering by abundance.
    identical_only : bool, default=False
        Only merge columns with identical sequences.
    vec : bool, default=True
        Use the vectorized aligner.
    band : int, default=-1
        Alignment band (-1: unbanded).
    n_jobs : int, optional
        Threads for the alignments.
    verbose : int, default=0
        Log a summary when > 0.

    Returns
    -------
    pandas.DataFrame
        The collapsed table, ordered by decreasing abundance.
    """
    tab = seqtab.copy()
    tab.columns = [str(c) for c in tab.columns]
    if tab.columns.duplicated().any():
        tab = tab.T.groupby(level=0, sort=False).sum().T
        logger.info("Duplicate sequences detected and merged.")
    if identical_only:
        return tab.astype(np.int64)
    totals = tab.sum(axis=0)
    order = np.argsort(-totals.to_numpy(), kind="stable")
    seqs = [tab.columns[i] for i in order]
    target = np.asarray(
        _core.collapse_no_mismatch(
            seqs, int(min_overlap), bool(vec), int(band), effective_n_jobs(n_jobs)
        )
    )
    # Column q of the abundance-ordered table is added to column target[q] (vectorised).
    vals = tab.to_numpy(dtype=np.int64)[:, order]
    summed = np.zeros_like(vals)
    np.add.at(summed.T, target.astype(np.intp), vals.T)
    kept = np.zeros(len(seqs), dtype=bool)
    kept[target] = True
    out_vals = np.zeros_like(vals)
    out_vals[:, order] = summed  # back to the input column order
    keep_cols = np.zeros(tab.shape[1], dtype=bool)
    keep_cols[order] = kept
    collapsed = pd.DataFrame(
        out_vals[:, keep_cols], index=tab.index, columns=tab.columns[keep_cols], dtype=np.int64
    )
    collapsed = _order_columns(collapsed, order_by)
    collapsed = _order_columns(collapsed, "abundance")
    if verbose:
        logger.info(
            "Output %d collapsed sequences out of %d input sequences.",
            collapsed.shape[1],
            seqtab.shape[1],
        )
    return collapsed
