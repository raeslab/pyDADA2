"""pydada2: exact amplicon sequence variant inference (DADA2) for Python.

An unofficial Python port of the DADA2 R package
(https://github.com/benjjneb/dada2, Callahan et al. 2016, Nature Methods).
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from ._validation import NotFittedError
from .chimeras import (
    BimeraRemover,
    is_bimera,
    is_bimera_denovo,
    is_bimera_denovo_table,
    is_shift_denovo,
    remove_bimera_denovo,
)
from .denoise import DADA, DadaParams, DadaResult
from .derep import Dereplication, derep_fastq
from .error_models import (
    BinnedQualErrorModel,
    ErrorLearner,
    LoessErrorModel,
    NoQualErrorModel,
    PacBioErrorModel,
    get_errors,
    inflate_err,
    learn_errors,
)
from .filtering import (
    FilterStats,
    fastq_filter,
    fastq_paired_filter,
    filter_and_trim,
    remove_primers,
)
from .merge import merge_pairs
from .seqtab import collapse_no_mismatch, make_sequence_table, merge_sequence_tables
from .sequences import (
    get_sequences,
    get_uniques,
    is_phix,
    reverse_complement,
    seq_complexity,
    uniques_to_fasta,
)
from .taxonomy import (
    SpeciesAssigner,
    TaxonomyClassifier,
    add_species,
    assign_species,
    assign_taxonomy,
)

try:
    __version__ = version("pydada2")
except PackageNotFoundError:  # pragma: no cover - running from a source tree
    __version__ = "0.0.0+unknown"

__all__ = [
    "DADA",
    "BimeraRemover",
    "BinnedQualErrorModel",
    "DadaParams",
    "DadaResult",
    "Dereplication",
    "ErrorLearner",
    "FilterStats",
    "LoessErrorModel",
    "NoQualErrorModel",
    "NotFittedError",
    "PacBioErrorModel",
    "SpeciesAssigner",
    "TaxonomyClassifier",
    "__version__",
    "add_species",
    "assign_species",
    "assign_taxonomy",
    "collapse_no_mismatch",
    "derep_fastq",
    "fastq_filter",
    "fastq_paired_filter",
    "filter_and_trim",
    "get_errors",
    "get_sequences",
    "get_uniques",
    "inflate_err",
    "is_bimera",
    "is_bimera_denovo",
    "is_bimera_denovo_table",
    "is_phix",
    "is_shift_denovo",
    "learn_errors",
    "make_sequence_table",
    "merge_pairs",
    "merge_sequence_tables",
    "remove_bimera_denovo",
    "remove_primers",
    "reverse_complement",
    "seq_complexity",
    "uniques_to_fasta",
]
