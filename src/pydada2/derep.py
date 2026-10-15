"""Dereplication of FASTQ files (R: ``derepFastq`` / ``derep-class``)."""

from __future__ import annotations

import logging
import os
import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import overload

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import _core
from ._validation import PathLike
from .io import QualityType, read_fastq_chunks

__all__ = ["Dereplication", "derep_fastq"]

logger = logging.getLogger("pydada2")


@dataclass(frozen=True, repr=False, eq=False)
class Dereplication:
    """Dereplicated reads (R: ``derep-class``).

    Attributes
    ----------
    uniques : pandas.Series
        Abundance (int64) of each unique sequence, indexed by sequence, in decreasing
        abundance order.
    quals : numpy.ndarray
        Average quality per unique and position, shape (n_uniques, max_length); NaN beyond
        each sequence's length.
    map : numpy.ndarray
        For each input read, the 0-based index of its unique sequence (-1 for zero-length reads).
    """

    uniques: pd.Series
    quals: npt.NDArray[np.float64]
    map: npt.NDArray[np.int64]

    @property
    def sequences(self) -> list[str]:
        """The unique sequences, in order."""
        return [str(s) for s in self.uniques.index]

    @property
    def n_reads(self) -> int:
        """Total number of (non-empty) reads."""
        return int(self.uniques.sum())

    def __len__(self) -> int:
        return len(self.uniques)

    def __repr__(self) -> str:
        lens = [len(s) for s in self.uniques.index]
        rng = f"{min(lens)}-{max(lens)}" if lens else "0"
        return (
            f"Dereplication({self.n_reads} reads in {len(self.uniques)} unique sequences; "
            f"sequence lengths {rng})"
        )


@overload
def derep_fastq(  # type: ignore[overload-overlap]
    fls: PathLike, /, *, n: int = ..., quality_type: QualityType = ..., verbose: int = ...
) -> Dereplication: ...
@overload
def derep_fastq(
    fls: Sequence[PathLike], /, *, n: int = ..., quality_type: QualityType = ..., verbose: int = ...
) -> list[Dereplication]: ...


def derep_fastq(
    fls: PathLike | Sequence[PathLike],
    /,
    *,
    n: int = 1_000_000,
    quality_type: QualityType = "auto",
    verbose: int = 0,
) -> Dereplication | list[Dereplication]:
    """Dereplicate FASTQ file(s) (R: ``derepFastq``).

    Identical reads are collapsed into unique sequences with their abundance and average
    per-position quality. Files are streamed in chunks of ``n`` reads; as in R, the order of
    uniques with equal abundance depends on the chunking (lexical order within each chunk).

    Parameters
    ----------
    fls : path or sequence of paths
        FASTQ file(s), optionally compressed.
    n : int, default=1_000_000
        Number of reads per chunk.
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        Quality encoding (auto-detected per chunk by default).
    verbose : int, default=0
        Log progress at INFO level when > 0.

    Returns
    -------
    Dereplication or list of Dereplication
        One object per input file (a single object for a single path).

    Raises
    ------
    FileNotFoundError
        If a file does not exist.
    ValueError
        If a file has no reads, or a chunk contains only zero-length sequences (as in R).

    Examples
    --------
    >>> from pydada2.datasets import example_path
    >>> drp = derep_fastq(example_path("sam1F.fastq.gz"))
    >>> drp.n_reads
    1500
    """
    if isinstance(fls, str | os.PathLike):
        return _derep_one(Path(fls), n, quality_type, verbose)
    return [_derep_one(Path(f), n, quality_type, verbose) for f in fls]


def _derep_one(path: Path, n: int, quality_type: QualityType, verbose: int) -> Dereplication:
    if not path.is_file():
        raise FileNotFoundError(f"No such file: {path}")
    if verbose:
        logger.info("Dereplicating sequence entries in FASTQ file: %s", path)
    drp = _core.Dereplicator()
    n_chunks = 0
    for chunk in read_fastq_chunks(path, n=n, quality_type=quality_type):
        drp.add_chunk(chunk.seqs, chunk.quals, chunk.offset)
        n_chunks += 1
    if n_chunks == 0:  # R: the first (empty) chunk fails in qtables2
        raise ValueError("Only zero-length sequences detected during dereplication.")
    res = drp.finish()
    if res["n_zero_length"]:
        warnings.warn(
            "Zero-length sequences detected during dereplication. They will be ignored.",
            stacklevel=3,
        )
    uniques = pd.Series(
        np.asarray(res["abundance"], dtype=np.int64), index=pd.Index(res["uniques"], dtype=object)
    )
    uniques.name = "abundance"
    quals = np.asarray(res["quals"], dtype=np.float64)
    m = np.asarray(res["map"], dtype=np.int64)
    if verbose:
        logger.info(
            "Encountered %d unique sequences from %d total sequences read.",
            len(uniques),
            int(uniques.sum()),
        )
    return Dereplication(uniques=uniques, quals=quals, map=m)
