"""DADA2 sample inference (R: ``dada``, ``setDadaOpt``/``getDadaOpt``)."""

from __future__ import annotations

import dataclasses
import logging
import math
import os
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import numpy.typing as npt
import pandas as pd

from . import _core
from ._base import BaseEstimator
from ._parallel import effective_n_jobs
from .derep import Dereplication, derep_fastq
from .error_models import (
    ErrorModel,
    LoessErrorModel,
    NoQualErrorModel,
    accumulate_trans,
    get_errors,
)
from .sequences import get_sequences

__all__ = ["DADA", "DadaParams", "DadaResult", "combine_dereps"]

logger = logging.getLogger("pydada2")

# Memory budget of the alignment cache used across the rounds of one fit (all samples).
_ALIGNMENT_CACHE_BYTES = 1024**3

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]


@dataclass(frozen=True)
class DadaParams:
    """Algorithm parameters of :class:`DADA` (R: ``setDadaOpt``/``getDadaOpt``).

    Defaults equal those of DADA2 1.36. Field names are the lower-case R option names.

    Attributes
    ----------
    omega_a : float
        Threshold for the abundance p-value at which a new cluster is formed.
    omega_p : float
        Threshold for the abundance p-value of sequences with prior evidence.
    omega_c : float
        Threshold below which unique sequences are not corrected (left unassigned).
    detect_singletons : bool
        Allow singletons to form new clusters.
    use_kmers : bool
        Use a k-mer distance screen before aligning.
    kdist_cutoff : float
        k-mer distance above which sequences are not aligned.
    max_consist : int
        Maximum number of self-consistency iterations.
    match, mismatch, gap_penalty : int
        Alignment scores.
    band_size : int
        Alignment band (-1: unbanded).
    vectorized_alignment : bool
        Use the vectorized aligner.
    max_clust : int
        Maximum number of clusters (0: no limit).
    min_fold, min_hamming, min_abundance
        Minimum fold over-abundance, Hamming distance and abundance for new clusters.
    use_quals : bool
        Use quality scores (False: quality-independent error model).
    homopolymer_gap_penalty : int or None
        Gap penalty in homopolymers (None: same as ``gap_penalty``).
    sse : int
        SSE mode of the k-mer routines (0, 1, 2).
    gapless : bool
        Use the gapless k-mer screen.
    greedy : bool
        Greedy cluster recruitment.
    pseudo_prevalence, pseudo_abundance
        Thresholds defining pseudo-priors in ``pool="pseudo"`` mode.
    """

    omega_a: float = 1e-40
    omega_p: float = 1e-4
    omega_c: float = 1e-40
    detect_singletons: bool = False
    use_kmers: bool = True
    kdist_cutoff: float = 0.42
    max_consist: int = 10
    match: int = 5
    mismatch: int = -4
    gap_penalty: int = -8
    band_size: int = 16
    vectorized_alignment: bool = True
    max_clust: int = 0
    min_fold: float = 1
    min_hamming: int = 1
    min_abundance: int = 1
    use_quals: bool = True
    homopolymer_gap_penalty: int | None = None
    sse: int = 2
    gapless: bool = True
    greedy: bool = True
    pseudo_prevalence: int = 2
    pseudo_abundance: float = math.inf

    def replace(self, **changes: Any) -> DadaParams:
        """Return a copy with some fields replaced.

        Raises
        ------
        TypeError
            For unknown field names.
        """
        names = {f.name for f in dataclasses.fields(self)}
        unknown = set(changes) - names
        if unknown:
            raise TypeError(f"Unknown DADA option(s): {', '.join(sorted(unknown))}.")
        return dataclasses.replace(self, **changes)

    def resolved(self) -> DadaParams:
        """Apply R's normalisations (penalty signs, homopolymer gapping, band 0)."""
        p = self
        if p.gap_penalty > 0:
            p = p.replace(gap_penalty=-p.gap_penalty)
        homo = p.gap_penalty if p.homopolymer_gap_penalty is None else p.homopolymer_gap_penalty
        if homo > 0:
            homo = -homo
        p = p.replace(homopolymer_gap_penalty=homo)
        if homo != p.gap_penalty:
            p = p.replace(vectorized_alignment=False)
        if p.vectorized_alignment and p.band_size == 0:
            p = p.replace(vectorized_alignment=False)
        return p


@dataclass(frozen=True, repr=False, eq=False)
class DadaResult:
    """Denoising result for one sample (R: ``dada-class``).

    Attributes
    ----------
    denoised : pandas.Series
        Abundance of each inferred sequence (ASV), indexed by sequence, in cluster order.
    clusters : pandas.DataFrame
        Per-cluster information (R ``$clustering``): sequence, abundance, n0, n1, nunq, pval,
        birth_from (0-based, -1 for the first cluster), birth_pval, birth_fold, birth_ham
        (-1 for the first cluster), birth_qave.
    quality : numpy.ndarray
        Average quality per cluster and position (n_clusters x max_length; NaN padded).
    birth_subs : pandas.DataFrame
        Substitutions that drove the birth of each cluster: pos (0-based), ref, sub, qual,
        clust (0-based).
    trans : numpy.ndarray
        Transition counts (16 x n_qual) of this sample.
    map : numpy.ndarray
        For each unique sequence of the input, the 0-based cluster index (-1 if uncorrected).
    pval : numpy.ndarray
        Abundance p-value of each unique sequence.
    err_in : numpy.ndarray or list of numpy.ndarray
        Error matrix used (list of matrices in self-consistent mode).
    err_out : numpy.ndarray or None
        Error matrix estimated from this run's transitions.
    params : DadaParams
        Parameters used.
    """

    denoised: pd.Series
    clusters: pd.DataFrame
    quality: FloatArray
    birth_subs: pd.DataFrame
    trans: npt.NDArray[Any]
    map: IntArray
    pval: FloatArray
    err_in: Any
    err_out: FloatArray | None
    params: DadaParams = field(default_factory=DadaParams)

    @property
    def sequence(self) -> list[str]:
        """The inferred sequences, in cluster order."""
        return [str(s) for s in self.denoised.index]

    def __repr__(self) -> str:
        return (
            f"DadaResult({len(self.denoised)} sequence variants were inferred from "
            f"{len(self.map)} input unique sequences; "
            f"omega_a={self.params.omega_a}, omega_c={self.params.omega_c}, "
            f"band_size={self.params.band_size})"
        )


def combine_dereps(dereps: Sequence[Dereplication]) -> Dereplication:
    """Pool dereplicated samples into one (R: ``combineDereps2``).

    Parameters
    ----------
    dereps : sequence of Dereplication
        Samples to pool.

    Returns
    -------
    Dereplication
        Pooled uniques (decreasing abundance), abundance-weighted average qualities and the
        concatenated read maps.
    """
    if not dereps:
        raise ValueError("No dereplicated samples to combine.")
    maxlen = max(d.quals.shape[1] for d in dereps)
    seqs: dict[str, int] = {}
    for d in dereps:
        for s in d.uniques.index:
            seqs.setdefault(str(s), len(seqs))
    counts = np.zeros(len(seqs), dtype=np.int64)
    quals = np.zeros((len(seqs), maxlen), dtype=np.float64)
    maps = []
    for d in dereps:
        idx = np.fromiter((seqs[str(s)] for s in d.uniques.index), dtype=np.int64, count=len(d))
        ab = d.uniques.to_numpy(dtype=np.int64)
        counts[idx] += ab
        q = np.full((len(d), maxlen), np.nan)
        q[:, : d.quals.shape[1]] = d.quals
        quals[idx] = quals[idx] + q * ab[:, None].astype(np.float64)
        m = d.map
        maps.append(np.where(m < 0, -1, idx[np.clip(m, 0, None)]))
    quals = quals / counts[:, None].astype(np.float64)
    order = np.argsort(-counts, kind="stable")
    rank = np.empty_like(order)
    rank[order] = np.arange(order.size)
    allmap = np.concatenate(maps) if maps else np.zeros(0, dtype=np.int64)
    names = list(seqs)
    uniques = pd.Series(counts[order], index=pd.Index([names[i] for i in order], dtype=object))
    uniques.name = "abundance"
    return Dereplication(
        uniques=uniques,
        quals=quals[order],
        map=np.where(allmap < 0, -1, rank[np.clip(allmap, 0, None)]).astype(np.int64),
    )


def _as_dereps(X: Any) -> tuple[list[Dereplication], list[str] | None]:
    if isinstance(X, Dereplication):
        return [X], None
    if isinstance(X, str | os.PathLike):
        return [derep_fastq(X)], [Path(X).name]
    if isinstance(X, Mapping):
        names = [str(k) for k in X]
        items = list(X.values())
    else:
        items = list(X)
        names = None
        if items and all(isinstance(x, str | os.PathLike) for x in items):
            names = [Path(x).name for x in items]
    out = []
    for x in items:
        if isinstance(x, Dereplication):
            out.append(x)
        elif isinstance(x, str | os.PathLike):
            out.append(derep_fastq(x))
        else:
            raise TypeError(
                "DADA.fit requires Dereplication objects, FASTQ paths, or a list/dict of them."
            )
    if not out:
        raise ValueError("No samples provided.")
    return out, names


def _validate_derep(d: Dereplication, i: int, use_quals: bool) -> int:
    seqs = d.sequences
    if not all(_core.is_acgt(seqs)):
        raise ValueError(
            "Invalid Dereplication uniques. Sequences must be made up only of A/C/G/T."
        )
    if not use_quals:
        return 0
    if d.quals is None or d.quals.shape[0] != len(seqs):
        raise ValueError("Dereplication quals must have one row for each unique sequence.")
    lens = np.fromiter((len(s) for s in seqs), dtype=np.int64, count=len(seqs))
    if (lens > d.quals.shape[1]).any():
        raise ValueError(
            "Dereplication quals must have as many columns as the length of the sequences."
        )
    inside = np.arange(d.quals.shape[1])[None, :] < lens[:, None]
    if np.isnan(d.quals[inside]).any():
        raise ValueError(
            "NAs in Dereplication quals. Check that all input sequences had valid qualities."
        )
    if np.nanmin(d.quals) < 0:
        raise ValueError("Invalid Dereplication quals. Quality values must be positive.")
    qmax = math.ceil(np.nanmax(d.quals))
    if qmax > 250:
        raise ValueError(f"Sample {i} has an invalid maximum Phred quality score of {qmax}.")
    return qmax


def _clusters_frame(out: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sequence": out["sequence"],
            "abundance": np.asarray(out["abundance"], dtype=np.int64),
            "n0": np.asarray(out["n0"], dtype=np.int64),
            "n1": np.asarray(out["n1"], dtype=np.int64),
            "nunq": np.asarray(out["nunq"], dtype=np.int64),
            "pval": np.asarray(out["pval"], dtype=np.float64),
            "birth_from": np.asarray(out["birth_from"], dtype=np.int64),
            "birth_pval": np.asarray(out["birth_pval"], dtype=np.float64),
            "birth_fold": np.asarray(out["birth_fold"], dtype=np.float64),
            "birth_ham": np.asarray(out["birth_ham"], dtype=np.int64),
            "birth_qave": np.asarray(out["birth_qave"], dtype=np.float64),
        }
    )


def _birth_subs_frame(out: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "pos": np.asarray(out["bs_pos"], dtype=np.int64),
            "ref": out["bs_ref"],
            "sub": out["bs_sub"],
            "qual": np.asarray(out["bs_qual"], dtype=np.float64),
            "clust": np.asarray(out["bs_clust"], dtype=np.int64),
        }
    )


def _uniques_of(clusters: pd.DataFrame) -> pd.Series:
    s = pd.Series(
        clusters["abundance"].to_numpy(dtype=np.int64),
        index=pd.Index(clusters["sequence"].to_list(), dtype=object),
    )
    s.name = "abundance"
    return s


class DADA(BaseEstimator):
    """DADA2 sample inference: infer exact sequence variants from dereplicated reads.

    R: ``dada()``. With ``self_consistent=True`` the error rates are re-estimated from the
    inferred transitions until convergence (``err=None`` starts from the maximum possible error
    rates), which is what :func:`~pydada2.learn_errors` does.

    Parameters
    ----------
    err : array-like of shape (16, n_qual), optional
        Error matrix (rows A2A..T2T). Required unless ``self_consistent=True``.
    error_model : ErrorModel, optional
        Estimates error rates from transitions (default :class:`LoessErrorModel`;
        :class:`NoQualErrorModel` is used when ``params.use_quals`` is False).
    self_consistent : bool, default=False
        Alternate sample inference and error-rate estimation until convergence.
    pool : {False, True, "pseudo"}, default=False
        Process samples independently, pooled, or pseudo-pooled.
    priors : sequence of str, default=()
        Sequences with prior evidence of existence (tested against ``omega_p``).
    params : DadaParams or dict, optional
        Algorithm parameters (a dict updates the defaults).
    n_jobs : int, optional
        Threads used by the C++ core (-1: all CPUs).
    verbose : int, default=0
        Logging level (1: progress, 2: detailed core output).

    Attributes
    ----------
    results_ : list of DadaResult
        One result per input sample.
    sample_names_ : list of str or None
        Sample names (dict keys or file names), if known.
    err_in_ : numpy.ndarray or list of numpy.ndarray
        Error matrix used (all iterations in self-consistent mode).
    err_out_ : numpy.ndarray or None
        Final estimated error matrix.
    trans_ : numpy.ndarray
        Accumulated transition counts of the final iteration.
    n_iter_ : int
        Number of sample-inference rounds.
    converged_ : bool
        Whether the self-consistency loop converged (always True otherwise).
    params_ : DadaParams
        Resolved parameters.

    Examples
    --------
    >>> from pydada2 import derep_fastq
    >>> from pydada2.datasets import example_path, err_balanced_f
    >>> drp = derep_fastq(example_path("sam1F.fastq.gz"))
    >>> dd = DADA(err=err_balanced_f()).fit(drp)
    >>> len(dd.results_[0].denoised) > 0
    True
    """

    def __init__(
        self,
        err: Any = None,
        error_model: ErrorModel | None = None,
        self_consistent: bool = False,
        pool: bool | str = False,
        priors: Sequence[str] = (),
        params: DadaParams | Mapping[str, Any] | None = None,
        n_jobs: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.err = err
        self.error_model = error_model
        self.self_consistent = self_consistent
        self.pool = pool
        self.priors = priors
        self.params = params
        self.n_jobs = n_jobs
        self.verbose = verbose

    # ------------------------------------------------------------------ helpers ---
    def _resolve_params(self) -> DadaParams:
        if self.params is None:
            p = DadaParams()
        elif isinstance(self.params, DadaParams):
            p = self.params
        elif isinstance(self.params, Mapping):
            p = DadaParams().replace(**dict(self.params))
        else:
            raise TypeError("params must be a DadaParams, a dict or None.")
        return p.resolved()

    def fit(self, X: Any, y: None = None) -> DADA:
        """Run sample inference.

        Parameters
        ----------
        X : Dereplication, path, or list/dict of them
            Dereplicated samples (FASTQ paths are dereplicated first).
        y : None
            Ignored.

        Returns
        -------
        self
            The fitted estimator (results in ``results_``).

        Raises
        ------
        ValueError
            For invalid inputs or parameters.
        """
        p = self._resolve_params()
        dereps, names = _as_dereps(X)
        n_threads = effective_n_jobs(self.n_jobs)
        verbose = int(self.verbose)
        priors = set(get_sequences(list(self.priors))) if len(self.priors) else set()

        pool, pseudo = self.pool, False
        if len(dereps) <= 1:
            pool = False
        if isinstance(pool, str):
            if pool != "pseudo":
                raise ValueError("Invalid pool argument (use False, True or 'pseudo').")
            pool, pseudo = False, True
        pool = bool(pool)
        dereps_in = dereps
        if pool:
            dereps = [combine_dereps(dereps_in)]

        initialize_err = False
        err: FloatArray | None
        if self.self_consistent and self.err is None:
            err = None
            initialize_err = True
        else:
            err = np.array(get_errors(self.err), dtype=np.float64)

        if not 0 <= p.omega_a < 1:
            raise ValueError("omega_a must be between zero and one.")
        if not 0 <= p.omega_p < 1:
            raise ValueError("omega_p must be between zero and one.")
        if p.omega_p < p.omega_a and priors:
            warnings.warn("omega_p should generally be larger than omega_a.", stacklevel=2)
        if p.omega_c > 1e-10 and self.self_consistent:
            warnings.warn(
                "Strict error correction (omega_c < 1e-10) is not recommended when learning "
                "error rates.",
                stacklevel=2,
            )
        if p.omega_c >= 1 and self.self_consistent:
            raise ValueError("Some error correction required when learning error rates.")
        model: ErrorModel = self.error_model if self.error_model is not None else LoessErrorModel()
        if not p.use_quals:
            model = NoQualErrorModel()

        homo_gap = int(
            p.gap_penalty if p.homopolymer_gap_penalty is None else p.homopolymer_gap_penalty
        )
        core_log = (lambda m: logger.debug("%s", m)) if verbose >= 2 else None
        nconsist = 0 if initialize_err else 1
        errs: list[FloatArray | None] = []
        pseudo_priors: set[str] = set()
        converged = True
        qmaxes = [_validate_derep(d, i, p.use_quals) for i, d in enumerate(dereps)]
        seqs_all = [d.sequences for d in dereps]
        abunds_all = [d.uniques.to_numpy(dtype=np.int64) for d in dereps]
        for ab in abunds_all:
            if ab.size and ab.max() > np.iinfo(np.int32).max:
                raise ValueError("Abundances exceed the supported range (int32).")
        abund_lists = [ab.tolist() for ab in abunds_all]
        quals_all = [np.ascontiguousarray(d.quals, dtype=np.float64) for d in dereps]
        # Several rounds on the same samples (self-consistency, pseudo-pooling): reuse the
        # alignments (they do not depend on the error rates), within a memory budget.
        cache = (
            _core.AlignmentCache(len(dereps), _ALIGNMENT_CACHE_BYTES)
            if p.use_kmers and (self.self_consistent or pseudo)
            else None
        )
        # Samples run concurrently, each may use all n_jobs threads (one shared pool). Core log
        # messages only reach the calling thread, so verbose >= 2 runs samples one at a time.
        sample_threads = 1 if verbose >= 2 else min(n_threads, len(dereps))
        while True:
            if nconsist > 0:
                errs.append(err)
            err_list = []
            prior_list = []
            for i, d in enumerate(dereps):
                qmax = qmaxes[i]
                if initialize_err:
                    erri = np.ones((16, max(41, qmax + 1)), dtype=np.float64)
                else:
                    assert err is not None
                    erri = err
                if erri.shape[1] < qmax + 1:
                    if verbose:
                        logger.info(
                            "The supplied error matrix does not extend to the maximum observed "
                            "quality score in sample %d (q=%d); repeating its last column.",
                            i,
                            qmax,
                        )
                    erri = np.hstack([erri] + [erri[:, -1:]] * (qmax + 1 - erri.shape[1]))
                if verbose and nconsist == 1:
                    if self.self_consistent:
                        logger.info("selfConsist step 1: sample %d", i)
                    elif pool:
                        logger.info(
                            "%d samples were pooled: %d reads in %d unique sequences.",
                            len(dereps_in),
                            d.n_reads,
                            len(d),
                        )
                    else:
                        logger.info(
                            "Sample %d - %d reads in %d unique sequences.", i, d.n_reads, len(d)
                        )
                err_list.append(np.ascontiguousarray(erri, dtype=np.float64))
                prior_list.append([(s in priors) or (s in pseudo_priors) for s in seqs_all[i]])
            outs: list[dict[str, Any]] = _core.dada_uniques_batch(
                seqs_all,
                abund_lists,
                prior_list,
                err_list,
                quals_all,
                match=p.match,
                mismatch=p.mismatch,
                gap=p.gap_penalty,
                use_kmers=p.use_kmers,
                kdist_cutoff=p.kdist_cutoff,
                band_size=int(p.band_size),
                omega_a=p.omega_a,
                omega_p=p.omega_p,
                omega_c=p.omega_c,
                detect_singletons=p.detect_singletons,
                max_clust=1 if initialize_err else int(p.max_clust),
                min_fold=p.min_fold,
                min_hamming=int(p.min_hamming),
                min_abund=int(p.min_abundance),
                use_quals=True,
                final_consensus=False,
                vectorized_alignment=p.vectorized_alignment,
                homo_gap=homo_gap,
                sse=int(p.sse),
                gapless=p.gapless,
                greedy=p.greedy,
                n_threads=n_threads,
                sample_threads=sample_threads,
                verbose=verbose >= 2,
                log=core_log,
                cache=cache,
            )
            cur = accumulate_trans([o["trans"] for o in outs])
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    err = model(cur)
            except Exception:
                if self.self_consistent or verbose >= 2:
                    logger.warning(
                        "Error rates could not be estimated (this is usually because of very "
                        "few reads)."
                    )
                err = None
            if self.self_consistent:
                get_errors(err)  # validates; raises if the estimate failed
            if initialize_err:
                initialize_err = False
                assert err is not None
                err = err.copy()
                err[[0, 5, 10, 15], :] = 1.0
            done = (
                not self.self_consistent
                or any(e is not None and err is not None and np.array_equal(e, err) for e in errs)
                or nconsist >= p.max_consist
            )
            if done and (not pseudo or nconsist >= 2):
                if self.self_consistent and nconsist >= p.max_consist:
                    converged = any(
                        e is not None and err is not None and np.array_equal(e, err) for e in errs
                    )
                break
            if pseudo and nconsist >= 1:
                st = make_table([_uniques_of(_clusters_frame(o)) for o in outs])
                prev = (st > 0).sum(axis=0)
                tot = st.sum(axis=0)
                keep = (prev >= p.pseudo_prevalence) | (tot >= p.pseudo_abundance)
                pseudo_priors = set(st.columns[keep.to_numpy()])
            nconsist += 1

        if self.self_consistent and verbose:
            if converged:
                logger.info("Convergence after %d rounds.", nconsist)
            else:
                logger.info("Self-consistency loop terminated before convergence.")

        err_in: Any = list(errs) if self.self_consistent else errs[0]
        results: list[DadaResult] = []
        for i in range(len(dereps)):
            o = outs[i]
            clusters = _clusters_frame(o)
            results.append(
                DadaResult(
                    denoised=_uniques_of(clusters),
                    clusters=clusters,
                    quality=np.asarray(o["cluster_quals"], dtype=np.float64),
                    birth_subs=_birth_subs_frame(o),
                    trans=np.asarray(o["trans"]),
                    map=np.asarray(o["map"], dtype=np.int64),
                    pval=np.asarray(o["pval_raw"], dtype=np.float64),
                    err_in=err_in,
                    err_out=err,
                    params=p,
                )
            )
        if pool:
            results = _expand_pooled(results[0], dereps[0], dereps_in)
        self.results_ = results
        self.sample_names_ = names
        self.err_in_ = err_in
        self.err_out_ = err
        self.trans_ = cur
        self.n_iter_ = nconsist
        self.converged_ = converged
        self.params_ = p
        return self

    def __getitem__(self, i: int) -> DadaResult:
        """Return the result of sample ``i``."""
        return self.results_[i]


def make_table(uniques: Sequence[pd.Series]) -> pd.DataFrame:
    """Build a samples x sequences abundance table (R: ``makeSequenceTable``, unordered)."""
    cols: dict[str, int] = {}
    for u in uniques:
        for s in u.index:
            cols.setdefault(str(s), len(cols))
    mat = np.zeros((len(uniques), len(cols)), dtype=np.int64)
    for i, u in enumerate(uniques):
        mat[i, [cols[str(s)] for s in u.index]] = u.to_numpy(dtype=np.int64)
    return pd.DataFrame(mat, columns=list(cols))


def _expand_pooled(
    pooled: DadaResult, pooled_derep: Dereplication, dereps_in: Sequence[Dereplication]
) -> list[DadaResult]:
    """Split a pooled result into per-sample results (R's ``if(pool)`` block in ``dada``)."""
    pmap = pd.Series(pooled.map, index=pooled_derep.sequences)
    out = []
    nclust = len(pooled.denoised)
    for d in dereps_in:
        seqs = d.sequences
        member = pmap.loc[seqs].to_numpy()
        keep_idx = np.unique(member[member >= 0])
        keep = np.zeros(nclust, dtype=bool)
        keep[keep_idx] = True
        new_bi = np.cumsum(keep) - 1  # pooled cluster index -> sample cluster index (0-based)
        clusters = pooled.clusters.loc[keep].reset_index(drop=True)
        bs = pooled.birth_subs
        bs = bs.loc[keep[bs["clust"].to_numpy()]].copy()
        bs["clust"] = new_bi[bs["clust"].to_numpy()]
        bs = bs.reset_index(drop=True)
        smap = np.where(member >= 0, new_bi[np.clip(member, 0, None)], -1).astype(np.int64)
        ab = np.bincount(smap[smap >= 0], weights=d.uniques.to_numpy()[smap >= 0],
                         minlength=int(keep.sum())).astype(np.int64)  # fmt: skip
        clusters["abundance"] = ab
        denoised = pd.Series(ab, index=pd.Index(clusters["sequence"].to_list(), dtype=object))
        denoised.name = "abundance"
        out.append(
            dataclasses.replace(
                pooled,
                denoised=denoised,
                clusters=clusters,
                quality=pooled.quality[keep],
                birth_subs=bs,
                map=smap,
            )
        )
    return out
