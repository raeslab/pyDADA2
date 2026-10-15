"""Read filtering and trimming.

R: ``filterAndTrim``, ``fastqFilter``, ``fastqPairedFilter``, ``removePrimers``.

The per-read work runs in the C++ kernels ``filter_mate`` (K1) and ``match_primer`` (K7); this
module streams the files in chunks of ``n`` reads exactly like the R functions (the chunking
matters for ``orient_fwd`` and ``match_ids``, whose effects are applied per chunk).
"""

from __future__ import annotations

import logging
import math
import numbers
import re
import warnings
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, NamedTuple

import pandas as pd

from . import _core
from ._parallel import effective_n_jobs
from ._validation import PathLike, as_path_list
from .datasets import phix_genome
from .io import (
    FastqChunk,
    QualityType,
    format_fastq,
    open_fastq_writer,
    read_fastq_chunks,
)
from .sequences import reverse_complement

__all__ = [
    "FilterStats",
    "fastq_filter",
    "fastq_paired_filter",
    "filter_and_trim",
    "remove_primers",
]

logger = logging.getLogger("pydada2")


class FilterStats(NamedTuple):
    """Number of reads read and written by a filtering function."""

    reads_in: int
    reads_out: int


@dataclass(frozen=True)
class _Mate:
    trunc_q: int
    trunc_len: int
    max_len: float
    min_len: int
    trim_left: int
    trim_right: int
    max_n: int
    min_q: int
    max_ee: float
    rm_phix: bool
    rm_lowcomplex: float


_MATE_PARAMS = (
    "trunc_q",
    "trunc_len",
    "max_len",
    "min_len",
    "trim_left",
    "trim_right",
    "max_n",
    "min_q",
    "max_ee",
    "rm_phix",
    "rm_lowcomplex",
)


def _as_pair(value: Any, name: str) -> tuple[Any, Any]:
    if isinstance(value, str) or not isinstance(value, Sequence):
        return (value, value)
    if len(value) == 1:
        return (value[0], value[0])
    if len(value) == 2:
        return (value[0], value[1])
    raise ValueError(f"{name} must be a scalar or a (forward, reverse) pair.")


_INT_MIN, _INT_MAX = -(2**31), 2**31 - 1
# Thresholds that are only compared with read properties: R accepts +-Inf for them.
_INF_OK = ("min_len", "max_n", "min_q")


def _make_mate(values: dict[str, Any]) -> _Mate:
    ints: dict[str, int] = {}
    for key in ("trunc_q", "trunc_len", "min_len", "trim_left", "trim_right", "max_n", "min_q"):
        v = values[key]
        inf_ok = key in _INF_OK and not isinstance(v, bool) and isinstance(v, numbers.Real)
        if inf_ok and math.isinf(float(v)):
            ints[key] = _INT_MAX if float(v) > 0 else _INT_MIN
            continue
        if isinstance(v, bool) or not math.isfinite(float(v)) or float(v) != int(v):
            raise ValueError(f"{key} must be an integer, got {v!r}.")
        ints[key] = int(v)
    return _Mate(
        trunc_q=ints["trunc_q"],
        trunc_len=ints["trunc_len"],
        max_len=float(values["max_len"]),
        min_len=ints["min_len"],
        trim_left=ints["trim_left"],
        trim_right=ints["trim_right"],
        max_n=ints["max_n"],
        min_q=ints["min_q"],
        max_ee=float(values["max_ee"]),
        rm_phix=bool(values["rm_phix"]),
        rm_lowcomplex=float(values["rm_lowcomplex"]),
    )


def _check_trunc_q(trunc_q: int, offset: int) -> None:
    lo, hi = (0, 93) if offset == 33 else (-5, 41)  # ShortRead FastqQuality / SFastqQuality
    if not lo <= trunc_q <= hi:
        raise ValueError("Encoding for this trunc_q value not found.")


def _run_mate(
    seqs: list[str], quals: list[str], offset: int, mate: _Mate, n_threads: int
) -> tuple[Any, Any, Any]:
    _check_trunc_q(mate.trunc_q, offset)
    return _core.filter_mate(
        seqs,
        quals,
        mate.trunc_q,
        mate.trunc_len,
        mate.max_len,
        mate.min_len,
        mate.trim_left,
        mate.trim_right,
        mate.max_n,
        mate.min_q,
        mate.max_ee,
        mate.rm_phix,
        mate.rm_lowcomplex,
        offset,
        phix_genome() if mate.rm_phix else "",
        n_threads,
    )


def _check_orient(orient_fwd: str | None) -> str | None:
    if orient_fwd is None:
        return None
    if not all(_core.is_acgt([orient_fwd])):
        raise ValueError("Non-ACGT characters detected in orient_fwd.")
    return orient_fwd


class _LazyWriter:
    """Opens the output file on the first write (R writes the first chunk even if empty)."""

    def __init__(self, path: Path, compress: bool) -> None:
        self.path = path
        self.compress = compress
        self.fh: IO[bytes] | None = None

    def write(self, data: bytes) -> None:
        if self.fh is None:
            self.fh = open_fastq_writer(self.path, self.compress)
        self.fh.write(data)

    def close(self) -> None:
        if self.fh is not None:
            self.fh.close()
            self.fh = None


def _remove_existing(path: Path) -> None:
    if path.exists():
        path.unlink()


def fastq_filter(
    fn: PathLike,
    fout: PathLike,
    /,
    *,
    trunc_q: int = 2,
    trunc_len: int = 0,
    max_len: float = math.inf,
    min_len: int = 20,
    trim_left: int = 0,
    trim_right: int = 0,
    max_n: int = 0,
    min_q: int = 0,
    max_ee: float = math.inf,
    rm_phix: bool = True,
    rm_lowcomplex: float = 0,
    orient_fwd: str | None = None,
    n: int = 1_000_000,
    quality_type: QualityType = "auto",
    compress: bool = True,
    verbose: int = 0,
    n_threads: int = 1,
) -> FilterStats:
    r"""Filter and trim one FASTQ file (R: ``fastqFilter``).

    See :func:`filter_and_trim` for the meaning of the filtering parameters. Steps are applied
    in R's order: ``orient_fwd``, ``max_len``, ``trim_left``, ``trim_right``, ``trunc_q``,
    ``trunc_len``, ``min_len``, ``max_n``, ``min_q``, ``max_ee``, ``rm_phix``, ``rm_lowcomplex``.

    Parameters
    ----------
    fn, fout : str or os.PathLike
        Input FASTQ and output path (gzip-compressed when ``compress``). The output is removed
        if no read passes.
    trunc_q, trunc_len, max_len, min_len, trim_left, trim_right, max_n, min_q, max_ee, \
rm_phix, rm_lowcomplex, orient_fwd
        Filtering parameters (see :func:`filter_and_trim`).
    n : int, default=1_000_000
        Reads processed per chunk.
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        Quality encoding.
    compress : bool, default=True
        Gzip the output.
    verbose : int, default=0
        Log a summary when > 0.
    n_threads : int, default=1
        Threads used by the C++ kernel.

    Returns
    -------
    FilterStats
        ``(reads_in, reads_out)``.

    Raises
    ------
    ValueError
        If input and output are the same file or a parameter is invalid.
    """
    fn_p, fout_p = Path(fn), Path(fout)
    if fn_p.resolve() == fout_p.resolve():
        raise ValueError("The output and input files must be different.")
    mate = _make_mate(
        dict(
            trunc_q=trunc_q,
            trunc_len=trunc_len,
            max_len=max_len,
            min_len=min_len,
            trim_left=trim_left,
            trim_right=trim_right,
            max_n=max_n,
            min_q=min_q,
            max_ee=max_ee,
            rm_phix=rm_phix,
            rm_lowcomplex=rm_lowcomplex,
        )
    )
    orient = _check_orient(orient_fwd)
    _remove_existing(fout_p)
    writer = _LazyWriter(fout_p, compress)
    inseqs = outseqs = 0
    try:
        for chunk in read_fastq_chunks(fn_p, n=n, quality_type=quality_type):
            inseqs += len(chunk)
            ids, seqs, quals = chunk.ids, chunk.seqs, chunk.quals
            if orient is not None:
                ids, seqs, quals = _orient_single(ids, seqs, quals, orient)
            passed, start, end = _run_mate(seqs, quals, chunk.offset, mate, n_threads)
            keep = passed.nonzero()[0]
            outseqs += len(keep)
            writer.write(
                format_fastq(
                    (ids[i] for i in keep),
                    (seqs[i][start[i] : end[i]] for i in keep),
                    (quals[i][start[i] : end[i]] for i in keep),
                )
            )
    finally:
        writer.close()
    if verbose:
        pct = round(outseqs * 100 / inseqs, 1) if inseqs else 0.0
        logger.info("Read in %d, output %d (%s%%) filtered sequences.", inseqs, outseqs, pct)
    if outseqs == 0:
        logger.info("The filter removed all reads: %s not written.", fout_p)
        _remove_existing(fout_p)
    return FilterStats(inseqs, outseqs)


def _orient_single(
    ids: list[str], seqs: list[str], quals: list[str], orient: str
) -> tuple[list[str], list[str], list[str]]:
    # R: keepF <- narrow(sread(fq),1,barlen) == orient.fwd
    #    keepR <- narrow(sread(fq.rc),1,barlen) == orient.fwd & !keepF ; fwd first, then rc
    keep_f = [s.startswith(orient) for s in seqs]
    out_ids = [i for i, k in zip(ids, keep_f, strict=True) if k]
    out_seqs = [s for s, k in zip(seqs, keep_f, strict=True) if k]
    out_quals = [q for q, k in zip(quals, keep_f, strict=True) if k]
    for i, k in enumerate(keep_f):
        if k:
            continue
        rc = reverse_complement(seqs[i])
        if rc.startswith(orient):
            out_ids.append(ids[i])
            out_seqs.append(rc)
            out_quals.append(quals[i][::-1])
    return out_ids, out_seqs, out_quals


def _id_field_auto(first_id: str, id_sep: str) -> tuple[int, str]:
    fields = re.split(id_sep, first_id)
    # R: sapply(gregexpr(":", fields), length) -> 1 when there is no ':' (gregexpr gives -1)
    ncolon = [max(1, f.count(":")) for f in fields]
    if max(ncolon) == 6 and ncolon.count(6) == 1:
        return ncolon.index(6), "Current"
    if max(ncolon) == 4 and ncolon.count(4) == 1:
        return ncolon.index(4), "Old"
    raise ValueError(
        "Couldn't automatically detect the sequence identifier field in the fastq id string."
    )


def _field(id_: str, id_sep: str, field: int) -> str | None:
    parts = re.split(id_sep, id_)
    return parts[field] if field < len(parts) else None


def fastq_paired_filter(
    fn: tuple[PathLike, PathLike],
    fout: tuple[PathLike, PathLike],
    /,
    *,
    trunc_q: Any = 2,
    trunc_len: Any = 0,
    max_len: Any = math.inf,
    min_len: Any = 20,
    trim_left: Any = 0,
    trim_right: Any = 0,
    max_n: Any = 0,
    min_q: Any = 0,
    max_ee: Any = math.inf,
    rm_phix: Any = True,
    rm_lowcomplex: Any = 0,
    match_ids: bool = False,
    orient_fwd: str | None = None,
    id_sep: str = r"\s",
    id_field: int | None = None,
    n: int = 1_000_000,
    quality_type: QualityType = "auto",
    compress: bool = True,
    verbose: int = 0,
    n_threads: int = 1,
) -> FilterStats:
    r"""Filter and trim a pair of FASTQ files (R: ``fastqPairedFilter``).

    A read pair is kept only if both mates pass. Each filtering parameter may be a scalar (used
    for both mates) or a ``(forward, reverse)`` pair. See :func:`filter_and_trim`.

    Parameters
    ----------
    fn, fout : (path, path)
        Input forward/reverse FASTQ files and output paths.
    trunc_q, trunc_len, max_len, min_len, trim_left, trim_right, max_n, min_q, max_ee, \
rm_phix, rm_lowcomplex
        Per-mate filtering parameters.
    match_ids : bool, default=False
        Keep only reads whose ids occur in both files (assumes the same relative order).
    orient_fwd : str, optional
        Swap forward/reverse mates when the reverse read starts with this sequence.
    id_sep : str, default=r"\s"
        Regular expression separating the fields of the read ids.
    id_field : int, optional
        0-based id field used by ``match_ids``; auto-detected (Illumina CASAVA formats) if None.
    n, quality_type, compress, verbose, n_threads
        As for :func:`fastq_filter`.

    Returns
    -------
    FilterStats
        ``(reads_in, reads_out)`` counted in read pairs (forward file).

    Raises
    ------
    ValueError
        On inconsistent inputs, e.g. different numbers of reads without ``match_ids``.
    """
    if len(fn) != 2 or len(fout) != 2:
        raise ValueError("Two paired input and two paired output file names are required.")
    fin = [Path(p) for p in fn]
    fo = [Path(p) for p in fout]
    allp = [p.resolve() for p in fin + fo]
    if len(set(allp)) != 4:
        raise ValueError("The output and input file names must be different.")
    pairs = {
        k: _as_pair(v, k)
        for k, v in dict(
            trunc_q=trunc_q, trunc_len=trunc_len, max_len=max_len, min_len=min_len,
            trim_left=trim_left, trim_right=trim_right, max_n=max_n, min_q=min_q,
            max_ee=max_ee, rm_phix=rm_phix, rm_lowcomplex=rm_lowcomplex,
        ).items()
    }  # fmt: skip
    mate_f = _make_mate({k: v[0] for k, v in pairs.items()})
    mate_r = _make_mate({k: v[1] for k, v in pairs.items()})
    orient = _check_orient(orient_fwd)
    for p in fo:
        _remove_existing(p)
    w_f, w_r = _LazyWriter(fo[0], compress), _LazyWriter(fo[1], compress)

    it_f = read_fastq_chunks(fin[0], n=n, quality_type=quality_type)
    it_r = read_fastq_chunks(fin[1], n=n, quality_type=quality_type)
    empty = FastqChunk([], [], [], 33)
    rem_f: FastqChunk = empty
    rem_r: FastqChunk = empty
    first = True
    casava = "Undetermined"
    field = id_field
    inseqs = outseqs = 0
    try:
        while True:
            cf = next(it_f, empty)
            cr = next(it_r, empty)
            if len(cf) == 0 and len(cr) == 0:
                break
            inseqs += len(cf)
            if match_ids:
                if first:
                    if field is None:
                        field, casava = _id_field_auto(cf.ids[0], id_sep)
                else:
                    cf = _prepend(rem_f, cf)
                    cr = _prepend(rem_r, cr)
                assert field is not None
                ids_f = [_field(i, id_sep, field) for i in cf.ids]
                ids_r = [_field(i, id_sep, field) for i in cr.ids]
                if casava == "Old":  # drop the "#index/pair" suffix (forward reads only, as in R)
                    ids_f = [None if i is None else i.split("#")[0] for i in ids_f]
                set_f, set_r = set(ids_f), set(ids_r)
                in_r = [i in set_r for i in ids_f]
                in_f = [i in set_f for i in ids_r]
                last_f = max((k + 1 for k, v in enumerate(in_r) if v), default=0)
                last_r = max((k + 1 for k, v in enumerate(in_f) if v), default=0)
                rem_f = _slice(cf, range(last_f, len(cf)))
                rem_r = _slice(cr, range(last_r, len(cr)))
                cf = _slice(cf, [k for k, v in enumerate(in_r) if v])
                cr = _slice(cr, [k for k, v in enumerate(in_f) if v])
            elif len(cf) != len(cr):
                raise ValueError(
                    f"Mismatched forward and reverse sequence files: {len(cf)}, {len(cr)}."
                )
            if len(cf) != len(cr):
                raise ValueError("Forward and reverse reads could not be matched by id.")
            if orient is not None:
                cf, cr = _orient_paired(cf, cr, orient)
            pass_f, s_f, e_f = _run_mate(cf.seqs, cf.quals, cf.offset, mate_f, n_threads)
            pass_r, s_r, e_r = _run_mate(cr.seqs, cr.quals, cr.offset, mate_r, n_threads)
            keep = (pass_f & pass_r).nonzero()[0]
            outseqs += len(keep)
            w_f.write(
                format_fastq(
                    (cf.ids[i] for i in keep),
                    (cf.seqs[i][s_f[i] : e_f[i]] for i in keep),
                    (cf.quals[i][s_f[i] : e_f[i]] for i in keep),
                )
            )
            w_r.write(
                format_fastq(
                    (cr.ids[i] for i in keep),
                    (cr.seqs[i][s_r[i] : e_r[i]] for i in keep),
                    (cr.quals[i][s_r[i] : e_r[i]] for i in keep),
                )
            )
            first = False
    finally:
        w_f.close()
        w_r.close()
    if verbose:
        pct = round(outseqs * 100 / inseqs, 1) if inseqs else 0.0
        logger.info(
            "Read in %d paired-sequences, output %d (%s%%) filtered paired-sequences.",
            inseqs,
            outseqs,
            pct,
        )
    if outseqs == 0:
        logger.info("The filter removed all reads: %s and %s not written.", fo[0], fo[1])
        for p in fo:
            _remove_existing(p)
    return FilterStats(inseqs, outseqs)


def _prepend(rem: FastqChunk, chunk: FastqChunk) -> FastqChunk:
    if len(rem) == 0:
        return chunk
    if len(chunk) == 0:
        return rem
    return FastqChunk(
        rem.ids + chunk.ids, rem.seqs + chunk.seqs, rem.quals + chunk.quals, chunk.offset
    )


def _slice(chunk: FastqChunk, idx: Sequence[int] | range) -> FastqChunk:
    return FastqChunk(
        [chunk.ids[i] for i in idx],
        [chunk.seqs[i] for i in idx],
        [chunk.quals[i] for i in idx],
        chunk.offset,
    )


def _orient_paired(cf: FastqChunk, cr: FastqChunk, orient: str) -> tuple[FastqChunk, FastqChunk]:
    keep_f = [s.startswith(orient) for s in cf.seqs]
    keep_r = [s.startswith(orient) and not kf for s, kf in zip(cr.seqs, keep_f, strict=True)]
    i_f = [k for k, v in enumerate(keep_f) if v]
    i_r = [k for k, v in enumerate(keep_r) if v]

    def cat(a: FastqChunk, ia: list[int], b: FastqChunk, ib: list[int]) -> FastqChunk:
        return FastqChunk(
            [a.ids[k] for k in ia] + [b.ids[k] for k in ib],
            [a.seqs[k] for k in ia] + [b.seqs[k] for k in ib],
            [a.quals[k] for k in ia] + [b.quals[k] for k in ib],
            a.offset,
        )

    return cat(cf, i_f, cr, i_r), cat(cr, i_f, cf, i_r)


def _prepare_outputs(inputs: list[Path], outputs: list[Path], what: str) -> None:
    if len(inputs) != len(outputs):
        raise ValueError(f"Every input file ({what}) must have a corresponding output file.")
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        raise FileNotFoundError(f"Some input files ({what}) do not exist: {', '.join(missing)}")
    for d in {p.parent for p in outputs}:
        if not d.exists():
            logger.info("Creating output directory: %s", d)
            d.mkdir(parents=True, exist_ok=True)


def filter_and_trim(
    fwd: PathLike | Sequence[PathLike],
    filt: PathLike | Sequence[PathLike],
    /,
    rev: PathLike | Sequence[PathLike] | None = None,
    filt_rev: PathLike | Sequence[PathLike] | None = None,
    *,
    compress: bool = True,
    trunc_q: Any = 2,
    trunc_len: Any = 0,
    trim_left: Any = 0,
    trim_right: Any = 0,
    max_len: Any = math.inf,
    min_len: Any = 20,
    max_n: Any = 0,
    min_q: Any = 0,
    max_ee: Any = math.inf,
    rm_phix: Any = True,
    rm_lowcomplex: Any = 0,
    orient_fwd: str | None = None,
    match_ids: bool = False,
    id_sep: str = r"\s",
    id_field: int | None = None,
    n_jobs: int | None = None,
    n: int = 100_000,
    quality_type: QualityType = "auto",
    verbose: int = 0,
) -> pd.DataFrame:
    r"""Filter and trim FASTQ files (R: ``filterAndTrim``).

    Parameters
    ----------
    fwd : path or sequence of paths
        Input FASTQ files (forward reads when ``rev`` is given).
    filt : path or sequence of paths
        Output files, one per input (gzip-compressed when ``compress``). Outputs of samples
        in which no read passes are not written.
    rev, filt_rev : path or sequence of paths, optional
        Reverse-read inputs and outputs for paired-end data.
    compress : bool, default=True
        Gzip the output files.
    trunc_q : int, default=2
        Truncate reads at the first base with quality <= ``trunc_q``.
    trunc_len : int, default=0
        Truncate reads after ``trunc_len`` bases; shorter reads are discarded (0: off).
    trim_left, trim_right : int, default=0
        Bases removed from the start / end of each read.
    max_len : float, default=inf
        Remove reads longer than this (before trimming).
    min_len : int, default=20
        Remove reads shorter than this (after trimming). ``-inf``/``inf`` are accepted (as in R).
    max_n : int, default=0
        Maximum number of non-A/C/G/T bases (``inf``: no limit).
    min_q : int, default=0
        Remove reads with a quality <= ``min_q`` after truncation (only used when
        ``min_q > trunc_q``, as in R). ``-inf``/``inf`` are accepted.
    max_ee : float, default=inf
        Maximum expected errors, ``sum(10**(-Q/10))``.
    rm_phix : bool, default=True
        Remove reads matching the phiX genome.
    rm_lowcomplex : float, default=0
        Remove reads with :func:`~pydada2.seq_complexity` below this value (0: off).
    orient_fwd : str, optional
        Keep only reads starting with this sequence, reverse-complementing (single-end) or
        swapping mates (paired-end) when needed.
    match_ids, id_sep, id_field
        Paired-end id matching (see :func:`fastq_paired_filter`).
    n_jobs : int, optional
        Number of files processed in parallel (-1: all CPUs).
    n : int, default=100_000
        Reads per chunk.
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        Quality encoding.
    verbose : int, default=0
        Log per-file summaries when > 0.

    Paired-end parameters (``trunc_q`` … ``rm_lowcomplex``) may be ``(forward, reverse)``
    pairs.

    Returns
    -------
    pandas.DataFrame
        ``reads_in`` and ``reads_out`` per input file (index: file name).

    Raises
    ------
    ValueError
        If the inputs and outputs are inconsistent.
    FileNotFoundError
        If an input file does not exist.

    Examples
    --------
    >>> import tempfile, pathlib
    >>> from pydada2.datasets import example_path
    >>> out = pathlib.Path(tempfile.mkdtemp()) / "sam1F.fastq.gz"
    >>> filter_and_trim(example_path("sam1F.fastq.gz"), out, trunc_len=240, max_ee=2)
                    reads_in  reads_out
    sam1F.fastq.gz      1500       1087
    """
    fwd_l = as_path_list(fwd, name="fwd")
    filt_l = as_path_list(filt, name="filt")
    _prepare_outputs(fwd_l, filt_l, "fwd")
    paired = rev is not None
    rev_l: list[Path] = []
    filt_rev_l: list[Path] = []
    if paired:
        if filt_rev is None:
            raise ValueError("Output files for the reverse reads are required.")
        assert rev is not None
        rev_l = as_path_list(rev, name="rev")
        filt_rev_l = as_path_list(filt_rev, name="filt_rev")
        if len(rev_l) != len(fwd_l):
            raise ValueError("Paired forward and reverse input files must correspond.")
        _prepare_outputs(rev_l, filt_rev_l, "rev")
    outs = [p.resolve() for p in filt_l + filt_rev_l]
    if len(set(outs)) != len(outs):
        raise ValueError("All output files must be distinct.")
    if set(outs) & {p.resolve() for p in fwd_l + rev_l}:
        raise ValueError("Output files must be distinct from the input files.")

    params = dict(
        trunc_q=trunc_q,
        trunc_len=trunc_len,
        max_len=max_len,
        min_len=min_len,
        trim_left=trim_left,
        trim_right=trim_right,
        max_n=max_n,
        min_q=min_q,
        max_ee=max_ee,
        rm_phix=rm_phix,
        rm_lowcomplex=rm_lowcomplex,
    )
    if not paired:
        for k, v in params.items():
            if not isinstance(v, str) and isinstance(v, Sequence) and len(v) != 1:
                raise ValueError(f"{k}: filtering arguments must be scalars for single-end data.")
            params[k] = _as_pair(v, k)[0]
    workers = effective_n_jobs(n_jobs)
    nfiles = len(fwd_l)
    kernel_threads = max(1, workers // max(1, nfiles))

    def run(i: int) -> FilterStats:
        if paired:
            return fastq_paired_filter(
                (fwd_l[i], rev_l[i]),
                (filt_l[i], filt_rev_l[i]),
                **params,
                match_ids=match_ids,
                orient_fwd=orient_fwd,
                id_sep=id_sep,
                id_field=id_field,
                n=n,
                quality_type=quality_type,
                compress=compress,
                verbose=verbose,
                n_threads=kernel_threads,
            )
        return fastq_filter(
            fwd_l[i],
            filt_l[i],
            **params,
            orient_fwd=orient_fwd,
            n=n,
            quality_type=quality_type,
            compress=compress,
            verbose=verbose,
            n_threads=kernel_threads,
        )

    if workers > 1 and nfiles > 1:
        with ThreadPoolExecutor(max_workers=min(workers, nfiles)) as ex:
            stats = list(ex.map(run, range(nfiles)))
    else:
        stats = [run(i) for i in range(nfiles)]
    df = pd.DataFrame(
        {"reads_in": [s.reads_in for s in stats], "reads_out": [s.reads_out for s in stats]},
        index=pd.Index([p.name for p in fwd_l]),
    )
    if (df["reads_out"] == 0).all():
        warnings.warn(
            "No reads passed the filter. Please revisit your filtering parameters.", stacklevel=2
        )
    elif (df["reads_out"] == 0).any():
        logger.info("Some input samples had no reads pass the filter.")
    return df


def remove_primers(
    fn: PathLike | Sequence[PathLike],
    fout: PathLike | Sequence[PathLike],
    /,
    primer_fwd: str,
    primer_rev: str | None = None,
    *,
    max_mismatch: int = 2,
    allow_indels: bool = False,
    trim_fwd: bool = True,
    trim_rev: bool = True,
    orient: bool = True,
    compress: bool = True,
    quality_type: QualityType = "auto",
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.DataFrame:
    """Remove primers and orient reads (R: ``removePrimers``; typically for PacBio CCS reads).

    Reads must contain the forward primer and, if given, the reverse primer (in the read's
    orientation, i.e. usually the reverse complement of the reverse PCR primer). Primers may
    contain IUPAC codes (then ambiguity codes in reads match too, as in Biostrings).

    Parameters
    ----------
    fn, fout : path or sequence of paths
        Input FASTQ file(s) and outputs (gzip-compressed when ``compress``).
    primer_fwd : str
        Forward primer sequence.
    primer_rev : str, optional
        Reverse primer, in read orientation.
    max_mismatch : int, default=2
        Mismatches allowed (edit distance when ``allow_indels``).
    allow_indels : bool, default=False
        Allow insertions/deletions in primer matches.
    trim_fwd, trim_rev : bool, default=True
        Trim the primers (and anything outside them) from the reads.
    orient : bool, default=True
        Reverse-complement reads in which only the reverse complement matches the forward primer.
    compress : bool, default=True
        Gzip the outputs.
    quality_type : {"auto", "FastqQuality", "SFastqQuality"}
        Quality encoding (only used for auto-detection; qualities are written unchanged).
    n_jobs : int, optional
        Threads for primer matching.
    verbose : int, default=0
        Log per-file summaries when > 0.

    Returns
    -------
    pandas.DataFrame
        ``reads_in`` and ``reads_out`` per input file.
    """
    fn_l = as_path_list(fn, name="fn")
    fout_l = as_path_list(fout, name="fout")
    _prepare_outputs(fn_l, fout_l, "fn")
    outs = [p.resolve() for p in fout_l]
    if len(set(outs)) != len(outs):
        raise ValueError("All output files must be distinct.")
    if set(outs) & {p.resolve() for p in fn_l}:
        raise ValueError("Output files must be distinct from the input files.")
    if not isinstance(primer_fwd, str) or (
        primer_rev is not None and not isinstance(primer_rev, str)
    ):
        raise TypeError("Primer sequences must be provided as strings.")
    if allow_indels:
        logger.info("Primer matching with indels allowed is slower.")
    nt = effective_n_jobs(n_jobs)
    fixed_fwd = bool(_core.is_acgt([primer_fwd])[0])
    pf = primer_fwd.upper()
    has_rev = primer_rev is not None
    fixed_rev = bool(_core.is_acgt([primer_rev])[0]) if primer_rev is not None else True
    pr = primer_rev.upper() if primer_rev is not None else ""

    def match(seqs: list[str], primer: str, fixed: bool) -> dict[str, Any]:
        return _core.match_primer(seqs, primer, int(max_mismatch), bool(allow_indels), fixed, nt)

    rows = []
    for f_in, f_out in zip(fn_l, fout_l, strict=True):
        chunks = list(read_fastq_chunks(f_in, n=2**62, quality_type=quality_type))
        ids = [i for c in chunks for i in c.ids]
        seqs = [s for c in chunks for s in c.seqs]
        quals = [q for c in chunks for q in c.quals]
        inseqs = len(seqs)
        hf = match(seqs, pf, fixed_fwd)
        hr = match(seqs, pr, fixed_rev) if has_rev else None
        # R: files without any forward (or reverse) hit are skipped entirely (no output)
        if int(hf["nhits"].sum()) == 0 or (hr is not None and int(hr["nhits"].sum()) == 0):
            rows.append((inseqs, 0))
            if verbose:
                logger.info("%s: no primer matches; nothing written.", f_in.name)
            continue
        if orient:
            rcs = reverse_complement(seqs)
            hf_rc = match(rcs, pf, fixed_fwd)
            hr_rc = match(rcs, pr, fixed_rev) if has_rev else None
            flip = (hf["nhits"] == 0) & (hf_rc["nhits"] > 0)
            if verbose and flip.any():
                logger.info(
                    "%d sequences out of %d are being reverse-complemented.",
                    int(flip.sum()),
                    inseqs,
                )
            for k in flip.nonzero()[0]:
                seqs[k] = rcs[k]
                quals[k] = quals[k][::-1]
                for key in hf:
                    hf[key][k] = hf_rc[key][k]
                    if hr is not None and hr_rc is not None:
                        hr[key][k] = hr_rc[key][k]
        out_ids, out_seqs, out_quals = [], [], []
        for k in range(inseqs):
            if hf["nhits"][k] == 0 or (hr is not None and hr["nhits"][k] == 0):
                continue
            first = int(hf["first_end"][k]) if trim_fwd else 0  # R: end(match)+1 (1-based)
            width = len(seqs[k])
            last = int(hr["last_start"][k]) if (hr is not None and trim_rev) else width
            # R keeps reads with last > first (1-based), i.e. at least two bases
            if not last - 1 > first:
                continue
            out_ids.append(ids[k])
            out_seqs.append(seqs[k][first:last])
            out_quals.append(quals[k][first:last])
        _remove_existing(f_out)
        with open_fastq_writer(f_out, compress) as fh:
            fh.write(format_fastq(out_ids, out_seqs, out_quals))
        rows.append((inseqs, len(out_seqs)))
        if verbose:
            logger.info("%s: read in %d, output %d sequences.", f_in.name, inseqs, len(out_seqs))
    df = pd.DataFrame(rows, columns=["reads_in", "reads_out"], index=[p.name for p in fn_l])
    if (df["reads_out"] == 0).all():
        warnings.warn("No reads passed the primer detection.", stacklevel=2)
    elif (df["reads_out"] == 0).any():
        logger.info("Some input samples had no reads pass the primer detection.")
    return df
