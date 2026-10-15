"""Taxonomic assignment (R: ``assignTaxonomy``, ``assignSpecies``, ``addSpecies``).

The genus-level classifier is the naive Bayesian classifier of Wang et al. (2007), as in DADA2:
8-mer profiles of the reference sequences per taxonomy ("genus"), best hit by summed log
probabilities, and 100 bootstrap replicates (each using 1/8 of the query's k-mers) for
confidence. Species are assigned by exact matching.
"""

from __future__ import annotations

import logging
import re
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd

from . import _core
from ._base import BaseEstimator
from ._parallel import effective_n_jobs
from ._validation import PathLike, check_is_fitted
from .io import read_fasta
from .sequences import get_sequences, reverse_complement

__all__ = [
    "DEFAULT_TAX_LEVELS",
    "SpeciesAssigner",
    "TaxonomyClassifier",
    "add_species",
    "assign_species",
    "assign_taxonomy",
]

logger = logging.getLogger("pydada2")

DEFAULT_TAX_LEVELS = ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
_UNSPECIFIED = "_DADA2_UNSPECIFIED"
_MIN_REF_LEN = 20
_MIN_TAX_LEN = 50


def _r_strsplit(s: str, sep: str) -> list[str]:
    """R's ``strsplit(s, sep, fixed)``: a single trailing empty field is dropped."""
    parts = s.split(sep)
    if s.endswith(sep):
        parts = parts[:-1]
    return parts


def _seed(random_state: Any) -> int:
    if isinstance(random_state, int | np.integer) and not isinstance(random_state, bool):
        return int(random_state) & 0xFFFFFFFFFFFFFFFF
    rng = (
        random_state
        if isinstance(random_state, np.random.Generator)
        else np.random.default_rng(random_state)
    )
    return int(rng.integers(0, 2**63 - 1))


def _read_reference(ref: PathLike | Sequence[tuple[str, str]]) -> list[tuple[str, str]]:
    if isinstance(ref, list | tuple):
        return [(str(i), str(s).upper()) for i, s in ref]
    return read_fasta(ref)  # type: ignore[arg-type]


class TaxonomyClassifier(BaseEstimator):
    """Naive Bayesian taxonomic classifier (R: ``assignTaxonomy``).

    Parameters
    ----------
    min_boot : int, default=50
        Minimum bootstrap confidence (0-100) for assigning a taxonomic level.
    try_rc : bool, default=False
        Also try the reverse complement of each query and keep the better orientation.
    tax_levels : sequence of str
        Names of the taxonomic levels (columns of the output).
    random_state : int, numpy.random.Generator or None
        Seed of the bootstrap. Results are identical for the same ``random_state`` regardless
        of ``n_jobs`` (R uses its own RNG, so bootstrap values differ slightly from R).
    n_jobs : int, optional
        Threads.
    verbose : int, default=0
        Logging level.

    Attributes
    ----------
    n_levels_ : int
        Depth of the reference taxonomy.
    genera_ : list of str
        The distinct reference taxonomy strings (";"-separated, padded to ``n_levels_``).
    n_references_ : int
        Number of reference sequences used.
    """

    def __init__(
        self,
        min_boot: int = 50,
        try_rc: bool = False,
        tax_levels: Sequence[str] = DEFAULT_TAX_LEVELS,
        random_state: Any = None,
        n_jobs: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.min_boot = min_boot
        self.try_rc = try_rc
        self.tax_levels = tax_levels
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.verbose = verbose

    def fit(self, ref: PathLike | Sequence[tuple[str, str]], y: None = None) -> TaxonomyClassifier:
        """Build the reference k-mer model from a training FASTA.

        Parameters
        ----------
        ref : path or sequence of (id, sequence)
            Training set: FASTA whose ids are ";"-separated taxonomies (e.g. the DADA2-formatted
            SILVA/RDP files), or UNITE general FASTA releases.
        y : None
            Ignored.

        Returns
        -------
        self
            The fitted classifier.

        Raises
        ------
        ValueError
            If the reference is not in a supported format.
        """
        records = _read_reference(ref)
        short = [len(s) < _MIN_REF_LEN for _, s in records]
        if any(short):
            records = [r for r, sh in zip(records, short, strict=True) if not sh]
            warnings.warn(
                f"Some reference sequences were too short (<{_MIN_REF_LEN}nts) and were excluded.",
                stacklevel=2,
            )
        if not records:
            raise ValueError("No usable reference sequences.")
        refs = [s for _, s in records]
        tax = [i.strip() for i, _ in records]
        if all(re.search(r"FU\|re[pf]s", t) for t in tax[:10]):  # UNITE general FASTA
            logger.info("UNITE fungal taxonomic reference detected.")
            tax = [t.split("|")[4] if len(t.split("|")) > 4 else "" for t in tax]
            tax = [re.sub(r"[pcofg]__unidentified;", _UNSPECIFIED + ";", t) for t in tax]
            tax = [re.sub(r";s__(\w+)_", ";s__", t) for t in tax]
            tax = [re.sub(r";s__sp$", ";" + _UNSPECIFIED, t) for t in tax]
        if ";" not in tax[0]:
            if len(re.split(r"\s", tax[0])) == 3:
                raise ValueError(
                    "Incorrect reference file format for assign_taxonomy (this looks like a "
                    "file formatted for assign_species)."
                )
            raise ValueError("Incorrect reference file format for assign_taxonomy.")
        depth = [len(_r_strsplit(t, ";")) for t in tax]
        td = max(depth)
        tax = [t + (_UNSPECIFIED + ";") * (td - d) for t, d in zip(tax, depth, strict=True)]
        genera = list(dict.fromkeys(tax))
        index = {g: k for k, g in enumerate(genera)}
        levels = [_r_strsplit(g, ";") for g in genera]
        if any(len(lv) != td for lv in levels):
            raise ValueError("Inconsistent taxonomy depths in the reference.")
        codes = np.empty((len(genera), td), dtype=np.int32)
        for j in range(td):
            col: dict[str, int] = {}
            for k, lv in enumerate(levels):
                codes[k, j] = col.setdefault(lv[j], len(col))
        self.model_ = _core.TaxonomyModel(refs, [index[t] for t in tax], codes)
        self.genera_ = genera
        self.levels_ = levels
        self.n_levels_ = td
        self.n_references_ = len(refs)
        if self.verbose:
            logger.info("Finished processing reference fasta (%d sequences).", len(refs))
        return self

    def _columns(self) -> list[str]:
        names = list(self.tax_levels)
        return [names[j] if j < len(names) else f"Level{j + 1}" for j in range(self.n_levels_)]

    def _classify(self, seqs: Any) -> tuple[list[str], dict[str, Any]]:
        check_is_fitted(self, ["model_"])
        s = get_sequences(seqs)
        if not s:
            raise ValueError("No sequences to classify.")
        if min(len(x) for x in s) < _MIN_TAX_LEN:
            warnings.warn(
                f"Some sequences were shorter than {_MIN_TAX_LEN} nts and will not receive a "
                "taxonomic classification.",
                stacklevel=3,
            )
        res = _core.assign_taxonomy(
            self.model_,
            s,
            reverse_complement(s),
            bool(self.try_rc),
            _seed(self.random_state),
            effective_n_jobs(self.n_jobs),
        )
        return s, res

    def _tax_frame(self, s: list[str], res: dict[str, Any]) -> pd.DataFrame:
        td = self.n_levels_
        out: list[list[str | None]] = []
        boots = np.asarray(res["boot"])
        for i, g in enumerate(np.asarray(res["tax"])):
            row: list[str | None] = [None] * td
            if g >= 0:
                lv = self.levels_[int(g)]
                keep = [lv[j] for j in range(td) if boots[i, j] >= self.min_boot]
                for j, name in enumerate(keep):
                    row[j] = None if name == _UNSPECIFIED else name
            out.append(row)
        return pd.DataFrame(
            out, index=pd.Index(s, dtype=object), columns=self._columns(), dtype=object
        )

    def predict(self, seqs: Any, /, *, return_bootstrap: bool = False) -> Any:
        """Assign taxonomy.

        Parameters
        ----------
        seqs : object
            Query sequences (or any object accepted by :func:`~pydada2.get_sequences`, e.g. the
            columns of a sequence table).
        return_bootstrap : bool, default=False
            Also return the bootstrap confidences.

        Returns
        -------
        pandas.DataFrame or (pandas.DataFrame, pandas.DataFrame)
            Taxonomy (sequences x levels; None where unassigned) and, optionally, the bootstrap
            confidences (int, 0-100) of the best-hit taxonomy at every level.
        """
        s, res = self._classify(seqs)
        tax = self._tax_frame(s, res)
        if not return_bootstrap:
            return tax
        return tax, self._boot_frame(s, res)

    def _boot_frame(self, s: list[str], res: dict[str, Any]) -> pd.DataFrame:
        return pd.DataFrame(
            np.asarray(res["boot"], dtype=np.int64),
            index=pd.Index(s, dtype=object),
            columns=self._columns(),
        )

    def predict_bootstrap(self, seqs: Any, /) -> pd.DataFrame:
        """Bootstrap confidences (0-100) of the best-hit taxonomy at every level.

        Parameters
        ----------
        seqs : object
            Query sequences.

        Returns
        -------
        pandas.DataFrame
            Sequences x levels, int64.
        """
        s, res = self._classify(seqs)
        return self._boot_frame(s, res)


def assign_taxonomy(
    seqs: Any,
    ref_fasta: PathLike | Sequence[tuple[str, str]],
    /,
    *,
    min_boot: int = 50,
    try_rc: bool = False,
    output_bootstraps: bool = False,
    tax_levels: Sequence[str] = DEFAULT_TAX_LEVELS,
    random_state: Any = None,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> Any:
    """Assign taxonomy with the naive Bayesian classifier (R: ``assignTaxonomy``).

    Convenience wrapper around :class:`TaxonomyClassifier`.

    Parameters
    ----------
    seqs : object
        Query sequences.
    ref_fasta : path or sequence of (id, sequence)
        Training set.
    min_boot, try_rc, tax_levels, random_state, n_jobs, verbose
        See :class:`TaxonomyClassifier`.
    output_bootstraps : bool, default=False
        Also return the bootstrap confidences.

    Returns
    -------
    pandas.DataFrame or (pandas.DataFrame, pandas.DataFrame)
        Taxonomy table, and bootstraps when ``output_bootstraps``.

    Examples
    --------
    >>> from pydada2.datasets import example_path
    >>> from pydada2.io import read_fasta
    >>> q = [s for _, s in read_fasta(example_path("ten_16s.100.fa.gz"))[:2]]
    >>> tax = assign_taxonomy(q, example_path("example_train_set.fa.gz"), random_state=1)
    >>> tax.shape
    (2, 6)
    """
    clf = TaxonomyClassifier(
        min_boot=min_boot,
        try_rc=try_rc,
        tax_levels=tax_levels,
        random_state=random_state,
        n_jobs=n_jobs,
        verbose=verbose,
    ).fit(ref_fasta)
    return clf.predict(seqs, return_bootstrap=output_bootstraps)


def _map_hits(hits: list[int], names: list[str | None], keep: float, sep: str = "/") -> str | None:
    """R's ``mapHits``: NA counts as a distinct value but is dropped when pasting."""
    vals = [names[h] for h in hits]
    vals = [
        "Escherichia/Shigella" if (v is not None and ("Escherichia" in v or "Shigella" in v)) else v
        for v in vals
    ]
    distinct = {v for v in vals}
    if len(distinct) > keep:
        return None
    joined = sep.join(sorted(v for v in distinct if v is not None))
    return joined or None


def _match_genera(gen_tax: Any, gen_binom: Any, split_glyph: str = "/") -> bool:
    if gen_tax is None or gen_binom is None or pd.isna(gen_tax) or pd.isna(gen_binom):
        return False
    gen_tax, gen_binom = str(gen_tax), str(gen_binom)
    if not gen_tax or not gen_binom:
        return False
    if gen_tax == gen_binom:
        return True
    try:  # R builds unescaped regular expressions from the names
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            return bool(
                re.search(f"^{gen_binom}[ _{split_glyph}]", gen_tax)
                or re.search(f"{split_glyph}{gen_binom}$", gen_tax)
            )
    except re.error:
        prefixes = tuple(gen_binom + c for c in " _" + split_glyph)
        return bool(gen_tax.startswith(prefixes) or gen_tax.endswith(split_glyph + gen_binom))


class SpeciesAssigner(BaseEstimator):
    """Species assignment by exact matching (R: ``assignSpecies``).

    Parameters
    ----------
    allow_multiple : bool or int, default=False
        Report up to this many equally matching species, separated by "/" (True: all).
    try_rc : bool, default=False
        Also match the reverse complement of the references.
    n_jobs : int, optional
        Threads.
    verbose : int, default=0
        Logging level.
    """

    def __init__(
        self,
        allow_multiple: bool | int = False,
        try_rc: bool = False,
        n_jobs: int | None = None,
        verbose: int = 0,
    ) -> None:
        self.allow_multiple = allow_multiple
        self.try_rc = try_rc
        self.n_jobs = n_jobs
        self.verbose = verbose

    def fit(self, ref: PathLike | Sequence[tuple[str, str]], y: None = None) -> SpeciesAssigner:
        """Load the species reference (ids ``"<id> <Genus> <species> ..."``).

        Parameters
        ----------
        ref : path or sequence of (id, sequence)
            Species reference FASTA.
        y : None
            Ignored.

        Returns
        -------
        self
            The fitted assigner.

        Raises
        ------
        ValueError
            If the file is not in the species-reference format.
        """
        records = _read_reference(ref)
        if not records:
            raise ValueError("Empty species reference.")
        ids = [i for i, _ in records]
        if len(re.split(r"\s", ids[0])) < 3:
            if ids[0].count(";") >= 3:
                raise ValueError(
                    "Incorrect reference file format for assign_species (this looks like a file "
                    "formatted for assign_taxonomy)."
                )
            raise ValueError("Incorrect reference file format for assign_species.")
        fields = [re.split(r"\s", i) for i in ids]
        self.genus_ = [f[1] if len(f) > 1 else None for f in fields]
        self.species_ = [f[2] if len(f) > 2 else None for f in fields]
        self.refs_ = [s for _, s in records]
        return self

    def predict(self, seqs: Any, /) -> pd.DataFrame:
        """Assign genus and species.

        Parameters
        ----------
        seqs : object
            Query sequences (A/C/G/T only).

        Returns
        -------
        pandas.DataFrame
            Columns ``Genus`` and ``Species`` (None when unassigned), indexed by sequence.

        Raises
        ------
        ValueError
            If a query contains non-ACGT characters.
        """
        check_is_fitted(self, ["refs_"])
        s = get_sequences(seqs)
        if not all(_core.is_acgt(s)):
            raise ValueError("Non-ACGT characters present in the query sequences.")
        if isinstance(self.allow_multiple, bool):
            keep = float("inf") if self.allow_multiple else 1.0
        else:
            keep = float(int(self.allow_multiple))
        hits = _core.species_hits(s, self.refs_, bool(self.try_rc), effective_n_jobs(self.n_jobs))
        genus = [_map_hits(h, self.genus_, 1) for h in hits]
        species = [_map_hits(h, self.species_, keep) for h in hits]
        out = pd.DataFrame(
            {"Genus": genus, "Species": species}, index=pd.Index(s, dtype=object), dtype=object
        )
        if self.verbose:
            logger.info(
                "%d out of %d were assigned to the species level.",
                int(out["Species"].notna().sum()),
                len(s),
            )
        return out


def assign_species(
    seqs: Any,
    ref_fasta: PathLike | Sequence[tuple[str, str]],
    /,
    *,
    allow_multiple: bool | int = False,
    try_rc: bool = False,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.DataFrame:
    """Assign genus and species by exact matching (R: ``assignSpecies``).

    Parameters
    ----------
    seqs : object
        Query sequences (A/C/G/T only).
    ref_fasta : path or sequence of (id, sequence)
        Species reference.
    allow_multiple, try_rc, n_jobs, verbose
        See :class:`SpeciesAssigner`.

    Returns
    -------
    pandas.DataFrame
        ``Genus`` and ``Species`` per sequence.
    """
    return (
        SpeciesAssigner(
            allow_multiple=allow_multiple, try_rc=try_rc, n_jobs=n_jobs, verbose=verbose
        )
        .fit(ref_fasta)
        .predict(seqs)
    )


def add_species(
    taxtab: pd.DataFrame,
    ref_fasta: PathLike | Sequence[tuple[str, str]],
    /,
    *,
    allow_multiple: bool | int = False,
    try_rc: bool = False,
    n_jobs: int | None = None,
    verbose: int = 0,
) -> pd.DataFrame:
    """Add a ``Species`` column to a taxonomy table (R: ``addSpecies``).

    Species are kept only when the genus of the exact match agrees with the genus of the table.

    Parameters
    ----------
    taxtab : pandas.DataFrame
        Taxonomy table indexed by sequence (e.g. from :func:`assign_taxonomy`).
    ref_fasta : path or sequence of (id, sequence)
        Species reference.
    allow_multiple, try_rc, n_jobs, verbose
        See :class:`SpeciesAssigner`.

    Returns
    -------
    pandas.DataFrame
        ``taxtab`` with an additional ``Species`` column.
    """
    seqs = [str(i) for i in taxtab.index]
    binom = assign_species(
        seqs,
        ref_fasta,
        allow_multiple=allow_multiple,
        try_rc=try_rc,
        n_jobs=n_jobs,
        verbose=verbose,
    )
    gcol = "Genus" if "Genus" in taxtab.columns else taxtab.columns[-1]
    gen_match = [
        _match_genera(g, b)
        for g, b in zip(taxtab[gcol].tolist(), binom["Genus"].tolist(), strict=True)
    ]
    out = taxtab.copy()
    out["Species"] = [
        sp if ok else None for sp, ok in zip(binom["Species"].tolist(), gen_match, strict=True)
    ]
    if verbose:
        logger.info(
            "Of which %d had genera consistent with the input table.",
            int(out["Species"].notna().sum()),
        )
    return out
