"""Pairwise alignment helpers (R: ``nwalign``, ``nwhamming``)."""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple, overload

import numpy as np
import numpy.typing as npt

from . import _core

__all__ = ["AlignmentCounts", "eval_pair", "nwalign", "nwhamming", "pair_consensus"]


class AlignmentCounts(NamedTuple):
    """Matches, mismatches and internal indels of an alignment (end gaps excluded)."""

    match: int
    mismatch: int
    indel: int


def nwalign(
    s1: str,
    s2: str,
    /,
    *,
    match: int = 5,
    mismatch: int = -4,
    gap: int = -8,
    homo_gap: int | None = None,
    band: int = -1,
    endsfree: bool = True,
    vec: bool = False,
) -> tuple[str, str]:
    """Needleman-Wunsch alignment of two sequences (R: ``nwalign``).

    Parameters
    ----------
    s1, s2 : str
        Sequences (A/C/G/T only unless ``vec``).
    match, mismatch, gap : int
        Scores (DADA2 defaults 5, -4, -8).
    homo_gap : int, optional
        Gap penalty within homopolymers (default: ``gap``).
    band : int, default=-1
        Band width (-1: unbanded).
    endsfree : bool, default=True
        Free end gaps.
    vec : bool, default=False
        Use the vectorized aligner (no homopolymer gap penalty).

    Returns
    -------
    (str, str)
        The aligned sequences, with ``-`` for gaps.

    Raises
    ------
    ValueError
        For non-ACGT input (unless ``vec``) or ``homo_gap`` with ``vec``.

    Examples
    --------
    >>> nwalign("ACGTACGT", "ACGACGT")
    ('ACGTACGT', 'ACG-ACGT')
    """
    if not isinstance(s1, str) or not isinstance(s2, str):
        raise TypeError("Can only align character sequences.")
    homo = gap if homo_gap is None else homo_gap
    if vec:
        if homo != gap:
            raise ValueError(
                "Homopolymer gap penalties are not implemented in the vectorized aligner."
            )
        a, b = _core.nwvec([s1], [s2], match, mismatch, gap, band, endsfree)[0]
        return a, b
    if not all(_core.is_acgt([s1, s2])):
        raise ValueError("Sequences must contain only A/C/G/T characters.")
    a, b = _core.nwalign(s1, s2, match, mismatch, gap, homo, band, endsfree)
    return a, b


def eval_pair(al1: str, al2: str, /) -> AlignmentCounts:
    """Count matches, mismatches and internal indels of an alignment (R: ``C_eval_pair``).

    Parameters
    ----------
    al1, al2 : str
        Aligned sequences of equal length.

    Returns
    -------
    AlignmentCounts
        ``(match, mismatch, indel)``; gaps at either end are ignored.
    """
    m, mm, ind = _core.eval_pair(al1, al2)
    return AlignmentCounts(int(m), int(mm), int(ind))


def pair_consensus(al1: str, al2: str, /, *, prefer: int = 1, trim_overhang: bool = False) -> str:
    """Consensus of two aligned sequences, gaps removed (R: ``C_pair_consensus``).

    Parameters
    ----------
    al1, al2 : str
        Aligned sequences of equal length.
    prefer : {1, 2}, default=1
        Which sequence wins mismatches.
    trim_overhang : bool, default=False
        Remove the parts of ``al2`` overhanging the start of ``al1`` and vice versa.

    Returns
    -------
    str
        The consensus sequence.
    """
    return str(_core.pair_consensus(al1, al2, int(prefer), bool(trim_overhang)))


@overload
def nwhamming(s1: str, s2: str, /, **kwargs: object) -> int: ...  # type: ignore[overload-overlap]
@overload
def nwhamming(
    s1: Sequence[str], s2: Sequence[str], /, **kwargs: object
) -> npt.NDArray[np.int64]: ...


def nwhamming(s1: str | Sequence[str], s2: str | Sequence[str], /, **kwargs: object) -> object:
    """Mismatches plus internal indels after alignment (R: ``nwhamming``).

    Parameters
    ----------
    s1, s2 : str or sequence of str
        Sequences (vectorized pairwise over sequences of equal length).
    **kwargs
        Passed to :func:`nwalign`.

    Returns
    -------
    int or numpy.ndarray
        The distance(s).

    Examples
    --------
    >>> nwhamming("ACGTACGTTT", "ACGAACGTTT")
    1
    """
    if isinstance(s1, str) and isinstance(s2, str):
        c = eval_pair(*nwalign(s1, s2, **kwargs))  # type: ignore[arg-type]
        return c.mismatch + c.indel
    a = [s1] if isinstance(s1, str) else list(s1)
    b = [s2] if isinstance(s2, str) else list(s2)
    if len(a) != len(b):
        if len(a) == 1:
            a = a * len(b)
        elif len(b) == 1:
            b = b * len(a)
        else:
            raise ValueError("s1 and s2 must have the same length (or length 1).")
    return np.array([nwhamming(x, y, **kwargs) for x, y in zip(a, b, strict=True)], dtype=np.int64)
