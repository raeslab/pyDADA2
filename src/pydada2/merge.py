"""Merging of denoised forward and reverse reads (R: ``mergePairs``)."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from . import _core
from ._parallel import effective_n_jobs
from .denoise import DADA, DadaResult
from .derep import Dereplication
from .sequences import reverse_complement

__all__ = ["merge_pairs"]

logger = logging.getLogger("pydada2")

_COLUMNS = [
    "sequence",
    "abundance",
    "forward",
    "reverse",
    "nmatch",
    "nmismatch",
    "nindel",
    "prefer",
    "accept",
]


def _as_results(x: Any, what: str) -> tuple[list[DadaResult], bool]:
    if isinstance(x, DADA):
        return list(x.results_), len(x.results_) == 1
    if isinstance(x, DadaResult):
        return [x], True
    items = list(x)
    if not all(isinstance(i, DadaResult) for i in items):
        raise TypeError(f"{what} must be a fitted DADA, a DadaResult or a list of DadaResult.")
    return items, False


def _as_dereps(x: Any, what: str) -> list[Dereplication]:
    if isinstance(x, Dereplication):
        return [x]
    items = list(x)
    if not all(isinstance(i, Dereplication) for i in items):
        raise TypeError(f"{what} must be a Dereplication or a list of Dereplication.")
    return items


def merge_pairs(
    dada_f: Any,
    derep_f: Any,
    dada_r: Any,
    derep_r: Any,
    /,
    *,
    min_overlap: int = 12,
    max_mismatch: int = 0,
    return_rejects: bool = False,
    propagate_col: Sequence[str] = (),
    just_concatenate: bool = False,
    trim_overhang: bool = False,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.DataFrame | list[pd.DataFrame]:
    """Merge denoised forward and reverse reads (R: ``mergePairs``).

    Each read pair is assigned the forward and reverse ASVs of its reads; every distinct
    (forward, reverse) ASV pair is aligned (forward against the reverse complement of the
    reverse ASV) and merged if the overlap is long enough and has few enough mismatches.

    Parameters
    ----------
    dada_f, dada_r : DADA, DadaResult or list of DadaResult
        Denoised forward / reverse reads.
    derep_f, derep_r : Dereplication or list of Dereplication
        The corresponding dereplicated reads (same samples, same order).
    min_overlap : int, default=12
        Minimum overlap (matching positions) to merge.
    max_mismatch : int, default=0
        Maximum mismatches plus indels in the overlap.
    return_rejects : bool, default=False
        Also return the pairs that could not be merged (with ``accept`` False and an empty
        sequence).
    propagate_col : sequence of str, default=()
        Columns of the DADA ``clusters`` frames to copy (as ``F.<col>`` and ``R.<col>``).
    just_concatenate : bool, default=False
        Concatenate the reads (forward + ``NNNNNNNNNN`` + reverse complement of reverse)
        instead of merging.
    trim_overhang : bool, default=False
        Trim overhanging ends of the merged sequences.
    n_jobs : int, optional
        Threads for the alignments.
    verbose : int, default=0
        Log a summary per sample when > 0.

    Returns
    -------
    pandas.DataFrame or list of pandas.DataFrame
        One frame per sample with columns ``sequence``, ``abundance``, ``forward`` and
        ``reverse`` (0-based ASV indices), ``nmatch``, ``nmismatch``, ``nindel``, ``prefer``
        and ``accept``, sorted by decreasing abundance. A single frame for a single sample.

    Raises
    ------
    ValueError
        If the DADA and dereplication inputs do not correspond.
    """
    rf, single = _as_results(dada_f, "dada_f")
    rr, _ = _as_results(dada_r, "dada_r")
    df_ = _as_dereps(derep_f, "derep_f")
    dr_ = _as_dereps(derep_r, "derep_r")
    if not (len(rf) == len(rr) == len(df_) == len(dr_)):
        raise ValueError("dada_f, derep_f, dada_r and derep_r must have the same length.")
    nt = effective_n_jobs(n_jobs)
    out = [
        _merge_one(
            rf[i],
            df_[i],
            rr[i],
            dr_[i],
            min_overlap,
            max_mismatch,
            return_rejects,
            propagate_col,
            just_concatenate,
            trim_overhang,
            nt,
            verbose,
        )
        for i in range(len(rf))
    ]
    return out[0] if (single and len(out) == 1) else out


def _merge_one(
    dd_f: DadaResult,
    drp_f: Dereplication,
    dd_r: DadaResult,
    drp_r: Dereplication,
    min_overlap: int,
    max_mismatch: int,
    return_rejects: bool,
    propagate_col: Sequence[str],
    just_concatenate: bool,
    trim_overhang: bool,
    n_threads: int,
    verbose: int,
) -> pd.DataFrame:
    map_f, map_r = drp_f.map, drp_r.map
    if not (
        len(map_f) == len(map_r)
        and int(map_f.max(initial=-1)) + 1 == len(dd_f.map)
        and int(map_r.max(initial=-1)) + 1 == len(dd_r.map)
    ):
        raise ValueError("Non-corresponding Dereplication and DADA objects.")
    r_f = np.where(map_f >= 0, dd_f.map[np.clip(map_f, 0, None)], -1)
    r_r = np.where(map_r >= 0, dd_r.map[np.clip(map_r, 0, None)], -1)
    valid = (r_f >= 0) & (r_r >= 0)
    # unique (forward, reverse) pairs in order of first appearance (R: unique(data.frame))
    pairs = pd.DataFrame({"forward": r_f[valid], "reverse": r_r[valid]})
    if pairs.empty:
        if verbose:
            logger.info(
                "No paired-reads (in ZERO unique pairings) successfully merged out of %d "
                "pairings) input.",
                len(r_f),
            )
        return pd.DataFrame({c: pd.Series(dtype=object) for c in _COLUMNS})
    counts = pairs.value_counts(sort=False)
    ups = pairs.drop_duplicates().reset_index(drop=True)
    fseq = dd_f.clusters["sequence"].to_numpy()[ups["forward"].to_numpy()].tolist()
    rseq = reverse_complement(
        dd_r.clusters["sequence"].to_numpy()[ups["reverse"].to_numpy()].tolist()
    )
    n = len(ups)
    if just_concatenate:
        ups["sequence"] = [f + "NNNNNNNNNN" + r for f, r in zip(fseq, rseq, strict=True)]
        ups["nmatch"] = 0
        ups["nmismatch"] = 0
        ups["nindel"] = 0
        ups["prefer"] = pd.array([pd.NA] * n, dtype="Int64")
        ups["accept"] = True
    else:
        mm, gap = (-64, -64) if max_mismatch == 0 else (-8, -8)
        n0_f = dd_f.clusters["n0"].to_numpy()[ups["forward"].to_numpy()]
        n0_r = dd_r.clusters["n0"].to_numpy()[ups["reverse"].to_numpy()]
        prefer = 1 + (n0_r > n0_f).astype(np.int64)
        res = _core.merge_pairs_core(
            fseq, rseq, prefer.tolist(), 1, mm, gap, bool(trim_overhang), n_threads
        )
        ups["nmatch"] = np.asarray(res["nmatch"], dtype=np.int64)
        ups["nmismatch"] = np.asarray(res["nmismatch"], dtype=np.int64)
        ups["nindel"] = np.asarray(res["nindel"], dtype=np.int64)
        ups["prefer"] = prefer
        ups["accept"] = (ups["nmatch"] >= min_overlap) & (
            (ups["nmismatch"] + ups["nindel"]) <= max_mismatch
        )
        ups["sequence"] = res["sequence"]
    ups["abundance"] = [
        int(counts[(f, r)]) for f, r in zip(ups["forward"], ups["reverse"], strict=True)
    ]
    ups.loc[~ups["accept"], "sequence"] = ""
    cols = [c for c in propagate_col if c in dd_f.clusters.columns]
    for col in cols:
        ups[f"F.{col}"] = dd_f.clusters[col].to_numpy()[ups["forward"].to_numpy()]
        ups[f"R.{col}"] = dd_r.clusters[col].to_numpy()[ups["reverse"].to_numpy()]
    ups = ups[_COLUMNS + [c for c in ups.columns if c not in _COLUMNS]]
    order = np.argsort(-ups["abundance"].to_numpy(), kind="stable")
    ups = ups.iloc[order].reset_index(drop=True)
    if verbose:
        acc = ups["accept"]
        logger.info(
            "%d paired-reads (in %d unique pairings) successfully merged out of %d (in %d "
            "pairings) input.",
            int(ups.loc[acc, "abundance"].sum()),
            int(acc.sum()),
            int(ups["abundance"].sum()),
            len(ups),
        )
    if not return_rejects:
        ups = ups.loc[ups["accept"]].reset_index(drop=True)
    if ups["sequence"].duplicated().any():
        logger.info("Duplicate sequences in merged output.")
    return ups
