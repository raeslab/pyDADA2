"""Diagnostic plots (R: ``plotQualityProfile``, ``plotErrors``, ``plotComplexity``).

Requires matplotlib (``pip install pydada2[plot]``). The data-preparation functions
(:func:`quality_profile_stats`, :func:`error_plot_data`) do not need matplotlib.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from ._validation import PathLike, as_path_list
from .error_models import TRANSITIONS, get_errors
from .io import quality_profile, read_fastq_chunks
from .sequences import seq_complexity

if TYPE_CHECKING:  # pragma: no cover
    from matplotlib.figure import Figure

__all__ = [
    "error_plot_data",
    "plot_complexity",
    "plot_errors",
    "plot_quality_profile",
    "quality_profile_stats",
]


def _plt() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as err:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Plotting requires matplotlib; install it with `pip install pydada2[plot]`."
        ) from err
    return plt


def quality_profile_stats(profile: pd.DataFrame, n_reads: int | None = None) -> pd.DataFrame:
    """Per-cycle summary of a quality profile (as drawn by R's ``plotQualityProfile``).

    Parameters
    ----------
    profile : pandas.DataFrame
        Output of :func:`pydada2.io.quality_profile` (columns cycle, score, count).
    n_reads : int, optional
        Number of reads for the cumulative-length line (default: ``profile.attrs["read_count"]``).

    Returns
    -------
    pandas.DataFrame
        Per cycle: ``mean`` quality, quartiles ``q25``, ``q50``, ``q75`` (smallest score whose
        cumulative fraction reaches the quantile) and ``cum`` (10 x fraction of reads reaching
        the cycle).
    """
    n = n_reads if n_reads is not None else int(profile.attrs.get("read_count", 0))
    rows = []
    for cyc, df in profile.sort_values(["cycle", "score"]).groupby("cycle", sort=True):
        sc = df["score"].to_numpy(dtype=float)
        ct = df["count"].to_numpy(dtype=float)
        tot = ct.sum()
        frac = np.cumsum(ct) / tot

        def q(p: float, sc: Any = sc, frac: Any = frac) -> float:
            return float(sc[np.nonzero(frac >= p)[0][0]])

        rows.append(
            {
                "cycle": int(str(cyc)),
                "mean": float((sc * ct).sum() / tot),
                "q25": q(0.25),
                "q50": q(0.5),
                "q75": q(0.75),
                "cum": 10 * tot / n if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_quality_profile(
    fl: PathLike | Sequence[PathLike], /, *, n: int = 500_000, aggregate: bool = False
) -> Figure:
    """Plot per-cycle quality-score distributions (R: ``plotQualityProfile``).

    Grey-scale heat map of score frequencies, mean (green) and quartiles (orange) per cycle and,
    for variable-length reads, the fraction of reads reaching each cycle (red).

    Parameters
    ----------
    fl : path or sequence of paths
        FASTQ file(s).
    n : int, default=500_000
        Reads used per file (the first ``n``).
    aggregate : bool, default=False
        Combine all files in one panel.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt = _plt()
    paths = as_path_list(fl, name="fl")
    profiles = [(p.name, quality_profile(p, n=n)) for p in paths]
    if aggregate:
        allp = pd.DataFrame(
            pd.concat([d for _, d in profiles])
            .groupby(["cycle", "score"], as_index=False)["count"]
            .sum()
        )
        allp.attrs["read_count"] = sum(int(d.attrs["read_count"]) for _, d in profiles)
        profiles = [(f"{len(paths)} files (aggregated)", allp)]
    fig, axes = plt.subplots(1, len(profiles), figsize=(5 * len(profiles), 4), squeeze=False)
    for ax, (label, prof) in zip(axes[0], profiles, strict=True):
        if len(prof):
            ax.scatter(prof["cycle"], prof["score"], c=prof["count"], cmap="Greys", marker="s", s=6)
            st = quality_profile_stats(prof)
            ax.plot(st["cycle"], st["mean"], color="#66C2A5")
            ax.plot(st["cycle"], st["q50"], color="#FC8D62", lw=0.5)
            ax.plot(st["cycle"], st["q25"], color="#FC8D62", lw=0.5, ls="--")
            ax.plot(st["cycle"], st["q75"], color="#FC8D62", lw=0.5, ls="--")
            if (st["cum"] != st["cum"].iloc[0]).any():
                ax.plot(st["cycle"], st["cum"], color="red", lw=0.5)
        rc = int(prof.attrs.get("read_count", 0))
        ax.set_title(label)
        ax.text(0, 0, f"Reads: {rc}", color="red", ha="left", va="bottom")
        ax.set_xlabel("Cycle")
        ax.set_ylabel("Quality Score")
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    out: Figure = fig
    return out


def error_plot_data(dq: Any, /) -> pd.DataFrame:
    """Long table of observed/estimated error rates (the data of R's ``plotErrors``).

    Parameters
    ----------
    dq : object
        Anything accepted by :func:`pydada2.get_errors` with transitions and/or estimates
        (e.g. a fitted :class:`~pydada2.ErrorLearner` or ``DADA``).

    Returns
    -------
    pandas.DataFrame
        Columns ``transition``, ``qual``, ``count``, ``tot``, ``observed``, ``estimated``,
        ``input`` and ``nominal`` (self-transitions blanked as in R).
    """
    d = get_errors(dq, detailed=True, enforce=False)
    trans = None if d["trans"] is None else np.asarray(d["trans"], dtype=float)
    err_out = None if d["err_out"] is None else np.asarray(d["err_out"], dtype=float)
    err_in = d["err_in"]
    if isinstance(err_in, list):
        err_in = err_in[0] if err_in else None
    err_in = None if err_in is None else np.asarray(err_in, dtype=float)
    base = trans if trans is not None else err_out
    if base is None:
        raise ValueError("Observed and/or estimated error rates must be available.")
    if base.shape[1] <= 1:
        raise ValueError("Error plots are only supported when quality scores are used.")
    nq = base.shape[1]
    rows = []
    for t, name in enumerate(TRANSITIONS):
        frm, to = name[0], name[2]
        tot_rows = [4 * "ACGT".index(frm) + k for k in range(4)]
        for qv in range(nq):
            count = trans[t, qv] if trans is not None else np.nan
            tot = trans[tot_rows, qv].sum() if trans is not None else np.nan
            nominal = 1 - 10 ** (-qv / 10) if frm == to else (1 / 3) * 10 ** (-qv / 10)
            rows.append(
                {
                    "transition": name,
                    "qual": qv,
                    "count": count,
                    "tot": tot,
                    "observed": count / tot if (trans is not None and tot > 0) else np.nan,
                    "estimated": err_out[t, qv]
                    if err_out is not None and qv < err_out.shape[1]
                    else np.nan,
                    "input": err_in[t, qv]
                    if err_in is not None and qv < err_in.shape[1]
                    else np.nan,
                    "nominal": nominal,
                }
            )
    df = pd.DataFrame(rows)
    self_t = df["transition"].str[0] == df["transition"].str[2]
    df.loc[self_t, ["count", "tot", "observed", "estimated", "input"]] = np.nan
    return df


def plot_errors(
    dq: Any,
    /,
    *,
    nti: Sequence[str] = ("A", "C", "G", "T"),
    ntj: Sequence[str] = ("A", "C", "G", "T"),
    obs: bool = True,
    err_out: bool = True,
    err_in: bool = False,
    nominal_q: bool = False,
) -> Figure:
    """Plot observed and estimated error rates per transition (R: ``plotErrors``).

    Parameters
    ----------
    dq : object
        Fitted :class:`~pydada2.ErrorLearner`/``DADA`` or anything accepted by
        :func:`pydada2.get_errors`.
    nti, ntj : sequence of {"A", "C", "G", "T"}
        Source and target nucleotides to show.
    obs : bool, default=True
        Show observed frequencies (points, size ~ count).
    err_out : bool, default=True
        Show the estimated error rates.
    err_in : bool, default=False
        Show the input error rates (dashed).
    nominal_q : bool, default=False
        Show the error rates implied by the quality scores (red).

    Returns
    -------
    matplotlib.figure.Figure
        A grid of panels, one per transition.

    Raises
    ------
    ValueError
        For invalid nucleotides.
    """
    acgt = ("A", "C", "G", "T")
    if (
        not (all(x in acgt for x in nti) and all(x in acgt for x in ntj))
        or len(set(nti)) != len(nti)
        or len(set(ntj)) != len(ntj)
    ):
        raise ValueError("nti and ntj must be nucleotide(s): A/C/G/T.")
    plt = _plt()
    df = error_plot_data(dq)
    fig, axes = plt.subplots(
        len(nti), len(ntj), figsize=(3 * len(ntj), 2.5 * len(nti)), squeeze=False
    )
    for i, a in enumerate(nti):
        for j, b in enumerate(ntj):
            ax = axes[i][j]
            sub = df[df["transition"] == f"{a}2{b}"]
            if obs and sub["observed"].notna().any():
                cnt = sub["count"].fillna(0).to_numpy()
                size = 16 * np.sqrt(cnt / cnt.max()) if cnt.max() > 0 else 1
                ax.scatter(sub["qual"], sub["observed"], s=size, color="0.4")
            if err_in and sub["input"].notna().any():
                ax.plot(sub["qual"], sub["input"], ls="--", color="black")
            if err_out and sub["estimated"].notna().any():
                ax.plot(sub["qual"], sub["estimated"], color="black")
            if nominal_q:
                ax.plot(sub["qual"], sub["nominal"], color="red")
            ax.set_yscale("log")
            ax.set_title(f"{a}2{b}")
            if i == len(nti) - 1:
                ax.set_xlabel("Consensus quality score")
            if j == 0:
                ax.set_ylabel("Error frequency (log10)")
    fig.tight_layout()
    out: Figure = fig
    return out


def plot_complexity(
    fl: PathLike | Sequence[PathLike],
    /,
    *,
    kmer_size: int = 2,
    window: int | None = None,
    by: int = 5,
    n: int = 100_000,
    bins: int = 100,
    aggregate: bool = False,
) -> Figure:
    """Histogram of sequence complexity (R: ``plotComplexity``).

    Parameters
    ----------
    fl : path or sequence of paths
        FASTQ file(s).
    kmer_size, window, by
        See :func:`pydada2.seq_complexity`.
    n : int, default=100_000
        Reads used per file (the first ``n``; R samples them at random).
    bins : int, default=100
        Number of histogram bins.
    aggregate : bool, default=False
        Combine all files in one histogram.

    Returns
    -------
    matplotlib.figure.Figure
        The figure.
    """
    plt = _plt()
    paths = as_path_list(fl, name="fl")
    data: list[tuple[str, np.ndarray]] = []
    for p in paths:
        seqs: list[str] = []
        for chunk in read_fastq_chunks(p, n=n):
            seqs = chunk.seqs
            break
        data.append((Path(p).name, seq_complexity(seqs, kmer_size=kmer_size, window=window, by=by)))
    if aggregate:
        data = [(f"{len(paths)} files (aggregated)", np.concatenate([d for _, d in data]))]
    fig, axes = plt.subplots(1, len(data), figsize=(5 * len(data), 3.5), squeeze=False)
    for ax, (label, cx) in zip(axes[0], data, strict=True):
        ax.hist(cx[np.isfinite(cx)], bins=bins, range=(0, 4**kmer_size), color="0.3")
        ax.set_title(label)
        ax.set_xlabel("Effective Oligonucleotide Number")
        ax.set_ylabel("Count")
    fig.tight_layout()
    out: Figure = fig
    return out
