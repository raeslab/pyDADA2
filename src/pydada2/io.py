"""FASTQ/FASTA input and output.

Reading uses :mod:`dnaio`. Sequences are upper-cased on input (as ShortRead/Biostrings do).
Quality strings are kept encoded; the quality offset (33 for Sanger/Illumina 1.8+, 64 for
Solexa/Illumina 1.3-1.7) is detected per chunk exactly like ShortRead's ``qualityType="Auto"``.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Literal

import dnaio
import numpy as np
import pandas as pd
from xopen import xopen

from . import _core
from ._validation import PathLike

__all__ = [
    "FastqChunk",
    "QualityType",
    "list_fastq",
    "quality_offset",
    "quality_profile",
    "read_fasta",
    "read_fastq_chunks",
    "write_fasta",
]

QualityType = Literal["auto", "FastqQuality", "SFastqQuality"]

_QUALITY_TYPES = {"auto": None, "fastqquality": 33, "sfastqquality": 64}


@dataclass(frozen=True)
class FastqChunk:
    """A chunk of FASTQ records (ShortRead ``FastqStreamer`` yield).

    Attributes
    ----------
    ids : list of str
        Full header lines without the leading ``@``.
    seqs : list of str
        Upper-case sequences.
    quals : list of str
        Encoded quality strings.
    offset : int
        Quality offset (33 or 64) of this chunk.
    """

    ids: list[str]
    seqs: list[str]
    quals: list[str]
    offset: int

    def __len__(self) -> int:
        return len(self.seqs)


def quality_offset(quals: Sequence[str], quality_type: QualityType = "auto") -> int:
    """Quality offset for a set of quality strings.

    Parameters
    ----------
    quals : sequence of str
        Encoded quality strings (only the first 10000 are inspected).
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        ``"auto"`` reproduces ShortRead's detection; the others force phred+33 / phred+64.

    Returns
    -------
    int
        33 or 64.

    Raises
    ------
    ValueError
        For an unknown ``quality_type``.
    """
    key = str(quality_type).lower()
    if key not in _QUALITY_TYPES:
        raise ValueError(
            f"quality_type must be 'auto', 'FastqQuality' or 'SFastqQuality', got {quality_type!r}."
        )
    forced = _QUALITY_TYPES[key]
    if forced is not None:
        return forced
    return int(_core.detect_quality_offset(list(quals)))


def read_fastq_chunks(
    path: PathLike, n: int = 1_000_000, quality_type: QualityType = "auto"
) -> Iterator[FastqChunk]:
    """Stream a FASTQ file in chunks of ``n`` records (like ShortRead's ``FastqStreamer``).

    Parameters
    ----------
    path : str or os.PathLike
        FASTQ file (optionally gzip/bzip2/xz compressed).
    n : int, default=1_000_000
        Records per chunk.
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}, default="auto"
        Quality encoding; ``"auto"`` is detected per chunk.

    Yields
    ------
    FastqChunk
        Consecutive chunks; the last one may be shorter.

    Raises
    ------
    ValueError
        If ``n`` is not positive.
    """
    n = int(n)
    if n < 1:
        raise ValueError("n must be a positive number of records.")
    with dnaio.open(os.fspath(path), mode="r") as reader:
        ids: list[str] = []
        seqs: list[str] = []
        quals: list[str] = []
        for rec in reader:
            if rec.qualities is None:  # pragma: no cover - dnaio yields FASTA without qualities
                raise ValueError(f"{path}: not a FASTQ file.")
            ids.append(rec.name)
            seqs.append(rec.sequence.upper())
            quals.append(rec.qualities)
            if len(seqs) == n:
                yield FastqChunk(ids, seqs, quals, quality_offset(quals, quality_type))
                ids, seqs, quals = [], [], []
        if seqs:
            yield FastqChunk(ids, seqs, quals, quality_offset(quals, quality_type))


def open_fastq_writer(path: PathLike, compress: bool) -> IO[bytes]:
    """Open a binary output stream (gzip-compressed when ``compress``)."""
    if compress:
        return xopen(os.fspath(path), "wb", compresslevel=6, threads=0, format="gz")
    return open(path, "wb")


def format_fastq(ids: Iterable[str], seqs: Iterable[str], quals: Iterable[str]) -> bytes:
    """Format records as FASTQ text (``@id``, sequence, ``+``, quality)."""
    return "".join(
        f"@{i}\n{s}\n+\n{q}\n" for i, s, q in zip(ids, seqs, quals, strict=True)
    ).encode()


def list_fastq(
    directory: PathLike,
    pattern: str | Sequence[str] = (r"\.fastq\.gz$", r"\.fastq\.bz2$", r"\.fastq$"),
) -> list[Path]:
    """List FASTQ files in a directory (R: ``parseFastqDirectory``).

    Parameters
    ----------
    directory : str or os.PathLike
        Directory to search (not recursive).
    pattern : str or sequence of str
        Regular expression(s) matched against the file names.

    Returns
    -------
    list of pathlib.Path
        Matching files, sorted by name.

    Raises
    ------
    NotADirectoryError
        If ``directory`` is not a directory.
    """
    d = Path(directory)
    if not d.is_dir():
        raise NotADirectoryError(f"{directory} is not a directory.")
    patterns = [pattern] if isinstance(pattern, str) else list(pattern)
    regs = [re.compile(p) for p in patterns]
    return sorted(p for p in d.iterdir() if p.is_file() and any(r.search(p.name) for r in regs))


def read_fasta(path: PathLike) -> list[tuple[str, str]]:
    """Read a FASTA (or FASTQ) file.

    Parameters
    ----------
    path : str or os.PathLike
        File path (optionally compressed).

    Returns
    -------
    list of (str, str)
        ``(id, sequence)`` pairs; sequences are upper-cased.
    """
    with dnaio.open(os.fspath(path), mode="r") as reader:
        return [(rec.name, rec.sequence.upper()) for rec in reader]


def write_fasta(
    records: Iterable[tuple[str, str]], path: PathLike, *, mode: str = "w", width: int = 20000
) -> None:
    """Write sequences to a FASTA file (R: ``writeFasta``).

    Parameters
    ----------
    records : iterable of (str, str)
        ``(id, sequence)`` pairs.
    path : str or os.PathLike
        Output path; written gzip-compressed if it ends in ``.gz``.
    mode : {"w", "a"}
        Write or append.
    width : int
        Maximum line width of the sequence lines.
    """
    if mode not in ("w", "a"):
        raise ValueError("mode must be 'w' or 'a'.")
    out = []
    for name, seq in records:
        out.append(f">{name}\n")
        for i in range(0, max(len(seq), 1), width):
            out.append(seq[i : i + width] + "\n")
    with xopen(os.fspath(path), "wt" if mode == "w" else "at") as fh:
        fh.write("".join(out))


def quality_profile(path: PathLike, n: int = 500_000) -> pd.DataFrame:
    """Per-cycle quality-score counts of a FASTQ file (ShortRead ``qa()``'s ``perCycle$quality``).

    Parameters
    ----------
    path : str or os.PathLike
        FASTQ file.
    n : int, default=500_000
        Number of reads to use. Unlike R (which samples ``n`` reads at random from larger
        files), the first ``n`` reads are used.

    Returns
    -------
    pandas.DataFrame
        Columns ``cycle`` (1-based), ``score`` and ``count``, sorted by cycle and score.
        ``df.attrs["read_count"]`` holds the number of reads used.
    """
    quals: list[bytes] = []
    offset = 33
    for chunk in read_fastq_chunks(path, n=n):
        quals = [q.encode() for q in chunk.quals]
        offset = chunk.offset
        break
    if not quals:
        df = pd.DataFrame({"cycle": [], "score": [], "count": []}, dtype="int64")
        df.attrs["read_count"] = 0
        return df
    lens = np.fromiter((len(q) for q in quals), dtype=np.int64, count=len(quals))
    flat = np.frombuffer(b"".join(quals), dtype=np.uint8).astype(np.int64) - offset
    starts = np.repeat(np.cumsum(lens) - lens, lens)
    cycle = np.arange(flat.size, dtype=np.int64) - starts + 1
    smin = int(flat.min())
    nscore = int(flat.max()) - smin + 1
    counts = np.bincount((cycle - 1) * nscore + (flat - smin), minlength=int(lens.max()) * nscore)
    cyc, sc = np.nonzero(counts.reshape(int(lens.max()), nscore))
    df = pd.DataFrame(
        {
            "cycle": cyc.astype(np.int64) + 1,
            "score": sc.astype(np.int64) + smin,
            "count": counts.reshape(int(lens.max()), nscore)[cyc, sc].astype(np.int64),
        }
    )
    df.attrs["read_count"] = len(quals)
    return df
