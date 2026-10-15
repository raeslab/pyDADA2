# mypy: disable-error-code="type-arg"
"""Private C++ core of pydada2 (ported from DADA2). Not a public API."""

from collections.abc import Callable, Sequence
from typing import Annotated

import numpy
from numpy.typing import NDArray

def core_version() -> str:
    """Version of the pydada2 C++ core."""

def upstream_version() -> str:
    """Version of upstream DADA2 the core was ported from."""

def ppois_upper(n: int, lam: float) -> float:
    """
    Upper Poisson tail P(X > n), X ~ Poisson(lam) (R: ppois(n, lam, lower.tail=FALSE)).
    """

def dada_uniques(
    seqs: Sequence[str],
    abundances: Sequence[int],
    priors: Sequence[bool],
    err: Annotated[
        NDArray[numpy.float64], dict(shape=(None, None), order="C", device="cpu", writable=False)
    ],
    quals: Annotated[
        NDArray[numpy.float64], dict(shape=(None, None), order="C", device="cpu", writable=False)
    ]
    | None,
    match: int = 5,
    mismatch: int = -4,
    gap: int = -8,
    use_kmers: bool = True,
    kdist_cutoff: float = 0.42,
    band_size: int = 16,
    omega_a: float = 1e-40,
    omega_p: float = 0.0001,
    omega_c: float = 1e-40,
    detect_singletons: bool = False,
    max_clust: int = 0,
    min_fold: float = 1.0,
    min_hamming: int = 1,
    min_abund: int = 1,
    use_quals: bool = True,
    final_consensus: bool = False,
    vectorized_alignment: bool = True,
    homo_gap: int = -8,
    sse: int = 2,
    gapless: bool = True,
    greedy: bool = True,
    n_threads: int = 1,
    verbose: bool = False,
    log: Callable | None = None,
) -> dict:
    """
    Run the DADA2 sample-inference algorithm on dereplicated sequences.

    err is 16 x Q (rows A2A..T2T); quals is n_uniques x max_len (or None). Returns a dict of
    arrays (clustering columns, birth_subs columns, trans, cluster_quals, map, pval_raw);
    indices are 0-based, -1/NaN mean NA.
    """

class AlignmentCache:
    """
    Cache of centre-vs-raw comparisons for repeated dada runs on the same samples (one slot
    per sample); comparisons beyond budget_bytes are computed without being cached.
    """

    def __init__(self, n_slots: int, budget_bytes: int) -> None: ...
    @property
    def n_slots(self) -> int: ...
    @property
    def budget_bytes(self) -> int: ...
    @property
    def bytes_used(self) -> int: ...

def dada_uniques_batch(
    seqs: Sequence[Sequence[str]],
    abundances: Sequence[Sequence[int]],
    priors: Sequence[Sequence[bool]],
    errs: Sequence[
        Annotated[
            NDArray[numpy.float64],
            dict(shape=(None, None), order="C", device="cpu", writable=False),
        ]
    ],
    quals: Sequence[
        Annotated[
            NDArray[numpy.float64],
            dict(shape=(None, None), order="C", device="cpu", writable=False),
        ]
        | None
    ],
    match: int = 5,
    mismatch: int = -4,
    gap: int = -8,
    use_kmers: bool = True,
    kdist_cutoff: float = 0.42,
    band_size: int = 16,
    omega_a: float = 1e-40,
    omega_p: float = 0.0001,
    omega_c: float = 1e-40,
    detect_singletons: bool = False,
    max_clust: int = 0,
    min_fold: float = 1.0,
    min_hamming: int = 1,
    min_abund: int = 1,
    use_quals: bool = True,
    final_consensus: bool = False,
    vectorized_alignment: bool = True,
    homo_gap: int = -8,
    sse: int = 2,
    gapless: bool = True,
    greedy: bool = True,
    n_threads: int = 1,
    sample_threads: int = 1,
    verbose: bool = False,
    log: Callable | None = None,
    cache: AlignmentCache | None = None,
) -> list:
    """
    dada_uniques on several samples (lists with one entry per sample; sample i uses cache
    slot i). Up to sample_threads samples run concurrently, each with up to n_threads
    threads, on one shared pool. Returns one dict per sample, as dada_uniques; the results
    do not depend on the thread counts. Log messages only come from the calling thread.
    """

def nwalign(
    s1: str,
    s2: str,
    match: int,
    mismatch: int,
    gap_p: int,
    homo_gap_p: int,
    band: int,
    endsfree: bool,
) -> tuple[str, str]:
    """
    Needleman-Wunsch alignment (ACGT only). Returns the two aligned strings.
    """

def nwvec(
    s1: Sequence[str],
    s2: Sequence[str],
    match: int,
    mismatch: int,
    gap_p: int,
    band: int,
    endsfree: bool,
) -> list[tuple[str, str]]:
    """Vectorized Needleman-Wunsch alignment of s1[i] with s2[i]."""

def eval_pair(s1: str, s2: str) -> list[int]:
    """Count (match, mismatch, indel) in the internal part of an alignment."""

def pair_consensus(s1: str, s2: str, prefer: int, trim_overhang: bool) -> str:
    """
    Consensus of two aligned sequences (prefer=1|2 wins mismatches), gaps removed.
    """

def is_acgt(seqs: Sequence[str]) -> NDArray[numpy.bool_]:
    """Whether each sequence consists only of A/C/G/T."""

def kmer_dist(s1: Sequence[str], s2: Sequence[str], kmer_size: int) -> list[float]:
    """k-mer distance between s1[i] and s2[i]."""

def kord_dist(s1: Sequence[str], s2: Sequence[str], kmer_size: int, sse: int) -> list[float]:
    """Ordered k-mer distance between s1[i] and s2[i]."""

def kmer_matches(s1: Sequence[str], s2: Sequence[str], kmer_size: int) -> list[int]:
    """Positional k-mer matches between s1[i] and s2[i]."""

def kdist_matches(s1: Sequence[str], s2: Sequence[str], kmer_size: int) -> list[int]:
    """Shared k-mer counts between s1[i] and s2[i]."""

def match_ref(seqs: Sequence[str], ref: str, word_size: int, non_overlapping: bool) -> list[int]:
    """
    Number of word_size-mers of each sequence found in the (circular) reference.
    """

NA_QUAL: int = -2147483648

def matrix_ee(
    quals: Annotated[
        NDArray[numpy.int32], dict(shape=(None, None), order="C", device="cpu", writable=False)
    ],
) -> NDArray[numpy.float64]:
    """
    Expected errors per row of an int32 quality matrix (NA_QUAL ends a row).
    """

def is_bimera(
    sq: str,
    pars: Sequence[str],
    allow_one_off: bool,
    min_one_off_par_dist: int,
    match: int,
    mismatch: int,
    gap_p: int,
    max_shift: int,
) -> bool:
    """Whether sq is a bimera of two of the given parents."""

def table_bimera(
    mat: Annotated[
        NDArray[numpy.int32], dict(shape=(None, None), order="C", device="cpu", writable=False)
    ],
    seqs: Sequence[str],
    min_fold: float,
    min_abund: int,
    allow_one_off: bool,
    min_one_off_par_dist: int,
    match: int,
    mismatch: int,
    gap_p: int,
    max_shift: int,
    n_threads: int = 1,
) -> tuple[NDArray[numpy.int32], NDArray[numpy.int32]]:
    """
    Per-sequence (nflag, nsam) bimera votes for a samples x sequences int32 table.
    """

class TaxonomyModel:
    """Reference k-mer model for the naive Bayesian classifier."""

    def __init__(
        self,
        refs: Sequence[str],
        ref_to_genus: Sequence[int],
        genusmat: Annotated[
            NDArray[numpy.int32], dict(shape=(None, None), order="C", device="cpu", writable=False)
        ],
        n_threads: int = 1,
    ) -> None: ...
    @property
    def ngenus(self) -> int: ...
    @property
    def nlevel(self) -> int: ...
    def __getstate__(self) -> tuple: ...
    def __setstate__(self, arg: tuple, /) -> None: ...

NBOOT: int = 100

def assign_taxonomy(
    model: TaxonomyModel,
    seqs: Sequence[str],
    rcs: Sequence[str],
    try_rc: bool,
    seed: int,
    n_threads: int = 1,
) -> dict:
    """
    Classify sequences; returns dict(tax, boot, boot_tax) with 0-based genus indices.
    """

def detect_quality_offset(quals: Sequence[str]) -> int:
    """
    ShortRead quality-type auto detection on a chunk: 64 (Solexa) or 33 (Sanger).
    """

def filter_mate(
    seqs: Sequence[str],
    quals: Sequence[str],
    trunc_q: int,
    trunc_len: int,
    max_len: float,
    min_len: int,
    trim_left: int,
    trim_right: int,
    max_n: int,
    min_q: int,
    max_ee: float,
    rm_phix: bool,
    rm_lowcomplex: float,
    qual_offset: int,
    phix: str,
    n_threads: int = 1,
) -> tuple:
    """
    Per-mate filter/trim kernel (K1). Returns (pass, start, end) with 0-based [start, end).
    """

def is_phix(
    seqs: Sequence[str],
    phix: str,
    word_size: int = 16,
    min_matches: int = 2,
    non_overlapping: bool = True,
) -> NDArray[numpy.bool_]:
    """
    Whether each sequence matches the phiX genome (or its reverse complement).
    """

def seq_complexity(
    seqs: Sequence[str], kmer_size: int = 2, window: int = 0, by: int = 5, n_threads: int = 1
) -> NDArray[numpy.float64]:
    """
    Effective number of k-mers (exp Shannon entropy); window <= 0 uses the whole sequence (K6).
    """

def dna_less(a: str, b: str) -> bool:
    """Biostrings lexical order of two DNA strings (IUPAC bit codes)."""

class Dereplicator:
    """Streaming dereplication (K2)."""

    def __init__(self) -> None: ...
    def add_chunk(self, seqs: Sequence[str], quals: Sequence[str], qual_offset: int) -> None: ...
    def finish(self) -> dict: ...

def match_primer(
    seqs: Sequence[str],
    primer: str,
    max_mismatch: int,
    with_indels: bool,
    fixed: bool,
    n_threads: int = 1,
) -> dict:
    """
    Primer matches per sequence (K7): counts and first/last match (0-based start, end excl.).
    """

def match_pattern(
    subject: str, pattern: str, max_mismatch: int, with_indels: bool, fixed: bool
) -> list[tuple[int, int]]:
    """
    All matches (0-based start, exclusive end) of pattern in subject (Biostrings semantics).
    """

def merge_pairs_core(
    fwd: Sequence[str],
    rev_rc: Sequence[str],
    prefer: Sequence[int],
    match: int,
    mismatch: int,
    gap: int,
    trim_overhang: bool,
    n_threads: int = 1,
) -> dict:
    """Align, evaluate and merge read pairs (K3)."""

def collapse_no_mismatch(
    seqs: Sequence[str], min_overlap: int, vec: bool, band: int, n_threads: int = 1
) -> NDArray[numpy.int32]:
    """Greedy no-mismatch collapsing (K5): target index per sequence."""

def is_shift_denovo(
    seqs: Sequence[str],
    abunds: Sequence[int],
    min_overlap: int,
    flag_subseqs: bool,
    n_threads: int = 1,
) -> NDArray[numpy.bool_]:
    """Shifted copies of more abundant sequences (K5)."""

def bimera_denovo(
    seqs: Sequence[str],
    abunds: Sequence[float],
    min_fold: float,
    min_parent_abund: float,
    allow_one_off: bool,
    min_one_off_par_dist: int,
    match: int,
    mismatch: int,
    gap_p: int,
    max_shift: int,
    n_threads: int = 1,
) -> NDArray[numpy.bool_]:
    """De novo bimera flags for a set of uniques (K4)."""

def species_hits(
    queries: Sequence[str], refs: Sequence[str], try_rc: bool, n_threads: int = 1
) -> list[list[int]]:
    """Indices of the references containing each query exactly (K8)."""
