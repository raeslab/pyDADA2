"""File formats of the command-line interface (:mod:`pydada2.cli`).

* ``pydada2.dada/1``: one denoised mate of one sample (a :class:`~pydada2.DadaResult` and its
  :class:`~pydada2.Dereplication`) as (gzipped) JSON, with everything
  :func:`~pydada2.merge_pairs` needs (:func:`write_dada_json`, :func:`read_dada_json`).
* Sequence tables as TSV, one row per ASV: ``asv_id``, ``sequence``, then one column of counts
  per sample (:func:`write_seqtab`, :func:`read_seqtab`). This is the transpose of the in-memory
  (R) orientation, which keeps the files readable with ``cut``/``awk``.

Files whose name ends in ``.gz`` are gzip-compressed (with a zero timestamp, so identical
content gives identical files).
"""

from __future__ import annotations

import dataclasses
import gzip
import json
import math
import os
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import IO, Any, NamedTuple

import numpy as np
import pandas as pd

from ._validation import PathLike
from .denoise import DadaParams, DadaResult, _clusters_frame, _uniques_of
from .derep import Dereplication

__all__ = [
    "DADA_FORMAT",
    "DadaSample",
    "read_dada_json",
    "read_seqtab",
    "write_dada_json",
    "write_seqtab",
]

#: Format identifier stored in (and required by) the denoising JSON files.
DADA_FORMAT = "pydada2.dada/1"

_CLUSTER_COLUMNS = (
    "sequence",
    "abundance",
    "n0",
    "n1",
    "nunq",
    "pval",
    "birth_from",
    "birth_pval",
    "birth_fold",
    "birth_ham",
    "birth_qave",
)
_BIRTH_SUB_COLUMNS = ("pos", "ref", "sub", "qual", "clust")


class DadaSample(NamedTuple):
    """Contents of a ``pydada2.dada/1`` file (see :func:`read_dada_json`).

    Attributes
    ----------
    sample : str
        Sample name.
    mate : str
        ``"fwd"`` or ``"rev"``.
    result : DadaResult
        The denoising result. ``quality`` has zero columns and ``err_in``/``err_out`` are None
        (they are not stored).
    derep : Dereplication
        The dereplicated reads (``quals`` has zero columns: not stored).
    source : str or None
        Name of the FASTQ file that was denoised.
    pydada2_version : str
        Version of pydada2 that wrote the file.
    """

    sample: str
    mate: str
    result: DadaResult
    derep: Dereplication
    source: str | None
    pydada2_version: str


@contextmanager
def _open_write(path: PathLike) -> Iterator[IO[bytes]]:
    """Open ``path`` for binary writing, gzip-compressed if it ends in ``.gz``.

    The gzip header holds no file name and a zero timestamp: identical content gives identical
    files (useful for workflow caching).
    """
    with open(path, "wb") as raw:
        if not os.fspath(path).endswith(".gz"):
            yield raw
            return
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as gz:
            yield gz  # type: ignore[misc]


def _jsonable(values: Any) -> list[Any]:
    """Convert a column to a JSON list: NumPy scalars to Python, NaN and +-inf to None."""
    out = np.asarray(values).tolist()
    return [None if isinstance(v, float) and not math.isfinite(v) else v for v in out]


def _floats(values: Sequence[Any]) -> np.ndarray[Any, np.dtype[np.float64]]:
    """Convert a JSON list back to float64 (None -> NaN)."""
    return np.array([np.nan if v is None else v for v in values], dtype=np.float64)


def write_dada_json(
    path: PathLike,
    result: DadaResult,
    derep: Dereplication,
    /,
    *,
    sample: str,
    mate: str,
    source: str | None = None,
) -> None:
    """Write one denoised mate of one sample in the ``pydada2.dada/1`` format.

    The file holds the clusters (R ``$clustering``), the birth substitutions, the per-unique
    cluster map and p-values, the transition counts, the dereplicated reads (sequences,
    abundances and the per-read map) and the DADA parameters: everything
    :func:`~pydada2.merge_pairs` needs. Per-position qualities (``result.quality``,
    ``derep.quals``) are not stored. Columns are stored as lists; NaN and +-inf as ``null``.

    Parameters
    ----------
    path : str or os.PathLike
        Output file; gzip-compressed if it ends in ``.gz`` (e.g. ``sam1_F.dada.json.gz``).
    result : DadaResult
        Denoising result of the sample (e.g. ``DADA(...).fit(derep).results_[0]``).
    derep : Dereplication
        The dereplicated reads that were denoised.
    sample : str
        Sample name.
    mate : {"fwd", "rev"}
        Which mate of a read pair this is.
    source : str, optional
        Name of the FASTQ file the reads came from (informational).

    Raises
    ------
    ValueError
        If ``mate`` is invalid or ``result`` does not correspond to ``derep``.

    See Also
    --------
    read_dada_json : Read the file back.
    """
    from . import __version__

    if mate not in ("fwd", "rev"):
        raise ValueError("mate must be 'fwd' or 'rev'.")
    if len(result.map) != len(derep):
        raise ValueError("The DADA result does not correspond to the dereplicated reads.")
    doc = {
        "format": DADA_FORMAT,
        "pydada2_version": __version__,
        "sample": str(sample),
        "mate": mate,
        "source": None if source is None else str(source),
        "clusters": {c: _jsonable(result.clusters[c].to_numpy()) for c in _CLUSTER_COLUMNS},
        "birth_subs": {c: _jsonable(result.birth_subs[c].to_numpy()) for c in _BIRTH_SUB_COLUMNS},
        "dada_map": _jsonable(result.map),
        "pval": _jsonable(result.pval),
        "trans": np.asarray(result.trans).tolist(),
        "derep": {
            "sequence": derep.sequences,
            "abundance": _jsonable(derep.uniques.to_numpy(dtype=np.int64)),
            "map": _jsonable(derep.map),
        },
        "params": {
            k: (None if isinstance(v, float) and not math.isfinite(v) else v)
            for k, v in dataclasses.asdict(result.params).items()
        },
    }
    with _open_write(path) as fh:
        fh.write(json.dumps(doc, separators=(",", ":")).encode())


def read_dada_json(path: PathLike) -> DadaSample:
    """Read a ``pydada2.dada/1`` file written by :func:`write_dada_json`.

    Parameters
    ----------
    path : str or os.PathLike
        The file (gzip-compressed or not).

    Returns
    -------
    DadaSample
        ``(sample, mate, result, derep, source, pydada2_version)``. ``result`` and ``derep``
        can be passed to :func:`~pydada2.merge_pairs` directly.

    Raises
    ------
    ValueError
        If the file is not in the ``pydada2.dada/1`` format.
    """
    with open(path, "rb") as raw:
        gz = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if gz else open
    with opener(path, "rb") as fh:
        try:
            doc = json.loads(fh.read())
        except ValueError as err:
            raise ValueError(f"{path}: not a {DADA_FORMAT} file ({err}).") from None
    if not isinstance(doc, dict) or doc.get("format") != DADA_FORMAT:
        found = doc.get("format") if isinstance(doc, dict) else None
        raise ValueError(f"{path}: not a {DADA_FORMAT} file (format: {found!r}).")

    clusters = _clusters_frame(
        {c: _floats(v) if "pval" in c or c in ("birth_fold", "birth_qave") else v
         for c, v in doc["clusters"].items()}
    )  # fmt: skip
    bs = doc["birth_subs"]
    birth_subs = pd.DataFrame(
        {
            "pos": np.asarray(bs["pos"], dtype=np.int64),
            "ref": list(bs["ref"]),
            "sub": list(bs["sub"]),
            "qual": _floats(bs["qual"]),
            "clust": np.asarray(bs["clust"], dtype=np.int64),
        }
    )
    defaults = DadaParams()
    params = DadaParams(
        **{
            f.name: (getattr(defaults, f.name) if doc["params"].get(f.name) is None
                     else doc["params"][f.name])
            for f in dataclasses.fields(DadaParams)
            if f.name in doc["params"]
        }
    )  # fmt: skip
    trans = np.asarray(doc["trans"], dtype=np.int32)
    result = DadaResult(
        denoised=_uniques_of(clusters),
        clusters=clusters,
        quality=np.empty((len(clusters), 0), dtype=np.float64),
        birth_subs=birth_subs,
        trans=trans.reshape(16, -1) if trans.size else np.zeros((16, 0), dtype=np.int32),
        map=np.asarray(doc["dada_map"], dtype=np.int64),
        pval=_floats(doc["pval"]),
        err_in=None,
        err_out=None,
        params=params,
    )
    d = doc["derep"]
    uniques = pd.Series(
        np.asarray(d["abundance"], dtype=np.int64),
        index=pd.Index(list(d["sequence"]), dtype=object),
        name="abundance",
    )
    derep = Dereplication(
        uniques=uniques,
        quals=np.empty((len(uniques), 0), dtype=np.float64),
        map=np.asarray(d["map"], dtype=np.int64),
    )
    return DadaSample(
        sample=str(doc["sample"]),
        mate=str(doc["mate"]),
        result=result,
        derep=derep,
        source=doc.get("source"),
        pydada2_version=str(doc.get("pydada2_version", "")),
    )


def _empty_dada(params: DadaParams) -> tuple[DadaResult, Dereplication]:
    """Return the DADA result and dereplication of a sample without reads (``infer-asv``)."""
    clusters = _clusters_frame({c: [] for c in _CLUSTER_COLUMNS})
    result = DadaResult(
        denoised=_uniques_of(clusters),
        clusters=clusters,
        quality=np.empty((0, 0), dtype=np.float64),
        birth_subs=pd.DataFrame(
            {
                "pos": np.empty(0, dtype=np.int64),
                "ref": pd.Series(dtype=object),
                "sub": pd.Series(dtype=object),
                "qual": np.empty(0, dtype=np.float64),
                "clust": np.empty(0, dtype=np.int64),
            }
        ),
        trans=np.zeros((16, 0), dtype=np.int32),
        map=np.empty(0, dtype=np.int64),
        pval=np.empty(0, dtype=np.float64),
        err_in=None,
        err_out=None,
        params=params,
    )
    uniques = pd.Series(np.empty(0, dtype=np.int64), index=pd.Index([], dtype=object))
    uniques.name = "abundance"
    derep = Dereplication(
        uniques=uniques,
        quals=np.empty((0, 0), dtype=np.float64),
        map=np.empty(0, dtype=np.int64),
    )
    return result, derep


def write_seqtab(
    seqtab: pd.DataFrame, path: PathLike, /, *, asv_ids: Sequence[str] | None = None
) -> None:
    """Write a sequence table as TSV, one row per ASV (the transpose of ``seqtab``).

    Columns: ``asv_id``, ``sequence``, then one column per sample (in ``seqtab``'s row order)::

        asv_id  sequence     sam1  sam2
        ASV1    TACGGAGG...  579   345

    Parameters
    ----------
    seqtab : pandas.DataFrame
        Samples x sequences table of counts (e.g. from :func:`~pydada2.make_sequence_table`).
    path : str or os.PathLike
        Output file; gzip-compressed if it ends in ``.gz``.
    asv_ids : sequence of str, optional
        One id per column of ``seqtab`` (default ``ASV1``, ``ASV2``, ... in column order).

    Raises
    ------
    ValueError
        If ``asv_ids`` has the wrong length.
    """
    seqs = [str(c) for c in seqtab.columns]
    ids = [f"ASV{i + 1}" for i in range(len(seqs))] if asv_ids is None else list(asv_ids)
    if len(ids) != len(seqs):
        raise ValueError("asv_ids must have one id per sequence (column) of seqtab.")
    out = pd.DataFrame(seqtab.to_numpy(dtype=np.int64).T, columns=[str(s) for s in seqtab.index])
    out.insert(0, "sequence", seqs)
    out.insert(0, "asv_id", ids)
    with _open_write(path) as fh:
        out.to_csv(fh, sep="\t", index=False)


def read_seqtab(path: PathLike, /, *, return_ids: bool = False) -> Any:
    """Read a sequence table written by :func:`write_seqtab` (or ``pydada2 merge``).

    Parameters
    ----------
    path : str or os.PathLike
        TSV file (gzip-compressed or not) with ``asv_id`` and ``sequence`` columns followed by
        one column of counts per sample.
    return_ids : bool, default=False
        Also return the ASV ids.

    Returns
    -------
    pandas.DataFrame or (pandas.DataFrame, pandas.Series)
        The table in the in-memory (R) orientation: samples x sequences, int64. With
        ``return_ids``, also the ``asv_id`` of each sequence (a Series indexed by sequence).

    Raises
    ------
    ValueError
        If the ``asv_id`` and ``sequence`` columns are missing.
    """
    df = pd.read_csv(
        Path(path), sep="\t", dtype={"asv_id": object, "sequence": object}, keep_default_na=False
    )
    if list(df.columns[:2]) != ["asv_id", "sequence"]:
        raise ValueError(f"{path}: expected the columns 'asv_id' and 'sequence' first.")
    seqs = [str(s) for s in df["sequence"]]
    counts = df.iloc[:, 2:]
    tab = pd.DataFrame(
        counts.to_numpy(dtype=np.int64).T,
        index=pd.Index([str(c) for c in counts.columns]),
        columns=pd.Index(seqs, dtype=object),
    )
    if not return_ids:
        return tab
    ids = pd.Series([str(i) for i in df["asv_id"]], index=pd.Index(seqs, dtype=object))
    ids.name = "asv_id"
    return tab, ids
