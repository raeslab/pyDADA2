"""Sequence utilities: reverse complement, complexity, phiX detection, uniques handling."""

from __future__ import annotations

import os
import warnings
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, overload

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import _core
from ._parallel import effective_n_jobs
from .datasets import phix_genome
from .io import read_fasta

if TYPE_CHECKING:
    from ._validation import PathLike

__all__ = [
    "get_sequences",
    "get_uniques",
    "is_phix",
    "reverse_complement",
    "seq_complexity",
    "uniques_to_fasta",
]

_COMPLEMENT = str.maketrans(
    "ACGTMRWSYKVHDBNacgtmrwsykvhdbn-+.",
    "TGCAKYWSRMBDHVNTGCAKYWSRMBDHVN-+.",
)


@overload
def reverse_complement(seqs: str) -> str: ...  # type: ignore[overload-overlap]
@overload
def reverse_complement(seqs: Iterable[str]) -> list[str]: ...


def reverse_complement(seqs: str | Iterable[str]) -> str | list[str]:
    """Reverse complement of DNA sequence(s), IUPAC aware (R: ``rc``).

    Parameters
    ----------
    seqs : str or iterable of str
        DNA sequence(s); lower case is converted to upper case, as in R.

    Returns
    -------
    str or list of str
        Same shape as the input.

    Examples
    --------
    >>> reverse_complement("ACGTN")
    'NACGT'
    >>> reverse_complement(["AAG", "rYk"])
    ['CTT', 'MRY']
    """
    if isinstance(seqs, str):
        return seqs.translate(_COMPLEMENT)[::-1]
    return [s.translate(_COMPLEMENT)[::-1] for s in seqs]


def seq_complexity(
    seqs: Any,
    /,
    *,
    kmer_size: int = 2,
    window: int | None = None,
    by: int = 5,
    n_jobs: int | None = None,
) -> npt.NDArray[np.float64]:
    """Sequence complexity: the effective number of k-mers (R: ``seqComplexity``).

    Complexity is ``exp`` of the Shannon entropy of the k-mer frequencies (only k-mers made of
    A/C/G/T are counted). Random sequences approach ``4**kmer_size``; homopolymers give 1.

    Parameters
    ----------
    seqs : sequences or object accepted by :func:`get_sequences`
        Input sequences.
    kmer_size : int, default=2
        k-mer size.
    window : int, optional
        If given, the minimum complexity over windows of this width (every ``by`` bases).
    by : int, default=5
        Step between windows.
    n_jobs : int, optional
        Number of threads.

    Returns
    -------
    numpy.ndarray
        Complexity per sequence (NaN when a sequence has no A/C/G/T k-mer).

    Examples
    --------
    >>> seq_complexity(["AAAAAAAA", "ACGTACGTTGCA"]).round(2).tolist()
    [1.0, 7.54]
    """
    s = get_sequences(seqs)
    if not all(_core.is_acgt(s)):
        warnings.warn(
            "Not all sequences were A/C/G/T only, which can interfere with the calculation "
            "of the Shannon richness.",
            stacklevel=2,
        )
    w = 0 if window is None else int(window)
    if window is not None and kmer_size >= w:
        raise ValueError(
            "The window over which kmer frequency is calculated must be larger than the kmer_size."
        )
    return np.asarray(_core.seq_complexity(s, int(kmer_size), w, int(by), effective_n_jobs(n_jobs)))


def is_phix(
    seqs: Any,
    /,
    *,
    word_size: int = 16,
    min_matches: int = 2,
    non_overlapping: bool = True,
) -> npt.NDArray[np.bool_]:
    """Flag sequences that match the phiX genome (R: ``isPhiX``).

    Parameters
    ----------
    seqs : sequences or object accepted by :func:`get_sequences`
        Sequences to test.
    word_size : int, default=16
        Size of the exact words compared with the genome.
    min_matches : int, default=2
        Minimum number of matching words (in either orientation) to flag a sequence.
    non_overlapping : bool, default=True
        Count only non-overlapping matching words.

    Returns
    -------
    numpy.ndarray of bool
        True where the sequence matches phiX.

    Examples
    --------
    >>> from pydada2.datasets import phix_genome
    >>> is_phix([phix_genome()[100:300], "ACGT" * 50]).tolist()
    [True, False]
    """
    s = get_sequences(seqs)
    return np.asarray(
        _core.is_phix(s, phix_genome(), int(word_size), int(min_matches), bool(non_overlapping))
    )


def _uniques_from_object(obj: Any) -> pd.Series:
    # Late imports avoid import cycles (these modules import sequences).
    from .derep import Dereplication

    if isinstance(obj, str | os.PathLike) and Path(obj).is_file():
        from .derep import derep_fastq

        derep = derep_fastq(obj)
        return derep.uniques.copy()
    if isinstance(obj, Dereplication):
        return obj.uniques.copy()
    denoised = getattr(obj, "denoised", None)  # DadaResult
    if isinstance(denoised, pd.Series):
        return denoised.copy()
    if isinstance(obj, pd.DataFrame) and {"sequence", "abundance"} <= set(obj.columns):
        return pd.Series(
            obj["abundance"].to_numpy(dtype=np.int64), index=obj["sequence"].astype(str).to_list()
        )
    if isinstance(obj, pd.DataFrame):  # sequence table: samples x sequences
        return obj.sum(axis=0).astype(np.int64)
    if isinstance(obj, pd.Series):
        return obj.astype(np.int64)
    if isinstance(obj, Mapping):
        return pd.Series({str(k): int(v) for k, v in obj.items()}, dtype=np.int64)
    raise TypeError(
        "Unrecognized format: requires a sequence->abundance mapping/Series, a FASTQ path, a "
        "Dereplication, a DADA result, a sequence table or a DataFrame with 'sequence' and "
        "'abundance' columns."
    )


def get_uniques(obj: Any, /, *, collapse: bool = True) -> pd.Series:
    """Uniques vector (sequence -> abundance) of an object (R: ``getUniques``).

    Parameters
    ----------
    obj : object
        A mapping or Series of sequence -> abundance, a FASTQ path (dereplicated), a
        :class:`~pydada2.Dereplication`, a DADA result, a sequence table, or a DataFrame with
        ``sequence`` and ``abundance`` columns (e.g. merged pairs).
    collapse : bool, default=True
        Sum the abundances of duplicated sequences.

    Returns
    -------
    pandas.Series
        int64 abundances indexed by sequence.

    Examples
    --------
    >>> get_uniques({"ACGT": 3, "AGGT": 1}).to_dict()
    {'ACGT': 3, 'AGGT': 1}
    """
    unqs = _uniques_from_object(obj)
    unqs.index = unqs.index.astype(str)
    if collapse and unqs.index.has_duplicates:
        unqs = unqs.groupby(level=0, sort=False).sum()
    unqs.name = "abundance"
    return unqs.astype(np.int64)


def get_sequences(obj: Any, /, *, collapse: bool = False) -> list[str]:
    """Sequences of an object (R: ``getSequences``).

    Parameters
    ----------
    obj : object
        A sequence or list of sequences, a FASTA/FASTQ path, or any object accepted by
        :func:`get_uniques`.
    collapse : bool, default=False
        Remove duplicated sequences (keeping the first occurrence).

    Returns
    -------
    list of str
        Upper-case sequences.

    Examples
    --------
    >>> get_sequences(["acgt", "ACGT", "GG"], collapse=True)
    ['ACGT', 'GG']
    """
    if isinstance(obj, str | os.PathLike) and Path(obj).is_file():
        seqs = [s for _, s in read_fasta(obj)]
    elif isinstance(obj, str):
        seqs = [obj]
    elif isinstance(obj, Sequence | np.ndarray | pd.Index) and not isinstance(obj, pd.Series):
        seqs = [str(s) for s in obj]
    else:
        seqs = list(get_uniques(obj, collapse=collapse).index)
    seqs = [s.upper() for s in seqs]
    if collapse:
        seqs = list(dict.fromkeys(seqs))
    return seqs


def uniques_to_fasta(
    unqs: Any,
    fout: PathLike,
    /,
    *,
    ids: Sequence[str] | None = None,
    mode: str = "w",
    width: int = 20000,
) -> None:
    """Write a uniques vector to FASTA (R: ``uniquesToFasta``; uchime-compatible ids).

    Parameters
    ----------
    unqs : object
        Anything accepted by :func:`get_uniques`.
    fout : str or os.PathLike
        Output path.
    ids : sequence of str, optional
        Record ids; default ``sq{i};size={abundance};``.
    mode : {"w", "a"}
        Write or append.
    width : int
        Maximum sequence line width.
    """
    from .io import write_fasta

    u = get_uniques(unqs, collapse=False)
    if ids is None:
        ids = [f"sq{i};size={a};" for i, a in enumerate(u.to_numpy(), start=1)]
    if len(ids) != len(u):
        raise ValueError("ids must have one entry per unique sequence.")
    write_fasta(zip(ids, u.index, strict=True), fout, mode=mode, width=width)
