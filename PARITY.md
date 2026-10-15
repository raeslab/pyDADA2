# PARITY.md — pydada2 vs R dada2 1.36

The phase-10 sign-off was updated after the REVIEW.md fixes (2026-10-15).

**Scope.** Every R NAMESPACE export maps to a Python object, as checked by
`tests/unit/test_api_coverage.py`. The parity tests (`uv run pytest -m parity`, 874 tests, of
which 600 are Biostrings pattern-matching cases) compare pydada2 against
fixtures produced by `tests/reference/generate_reference.R` with R 4.6.1 / dada2 1.36.0. The
whole suite has 1067 tests.

**Where it has been verified.**
- **Linux x86_64, Python 3.13** (editable install): verified locally, for the full suite and for
  the dataset benchmarks.
- **sdist and wheels:** the phase-10 checks on Python 3.11 (sdist) and 3.12 (abi3 wheel) predate
  these changes.
- **macOS (arm64, x86_64) and Windows:** CI (`.github/workflows/ci.yml`) runs the same suite,
  incl. tests of the built wheels, once the repository is pushed to raeslab/pyDADA2. That run
  has not happened yet.

Status legend:
- **exact:** identical sequences, counts, maps, flags and orders, with doubles within rtol 1e-10
  (most are bit-identical);
- **tolerance:** within the documented spec §8.3 tolerance;
- **deviation #n:** see `DEVIATIONS.md`.

| R function | Python | Status | Tests |
|---|---|---|---|
| `filterAndTrim` (single, paired, chunked, orient.fwd, matchIDs, all filters) | `filter_and_trim` | exact (output records via canonical md5); maxN, maxLen, minLen, rm.phix and rm.lowcomplex are each shown to remove reads on synthetic inputs (the generator asserts it); R's warnings compared | `parity/test_filtering.py`, `parity/test_filter_derep_coverage.py` |
| `fastqFilter`, `fastqPairedFilter` | `fastq_filter`, `fastq_paired_filter` | exact (via the above) | idem |
| `removePrimers` (mismatch and indel modes, orient, trim options) | `remove_primers` | exact, incl. `max.mismatch` ≥ primer length | `parity/test_filtering.py`, `parity/test_match_pattern.py` (300 Biostrings cases), `parity/test_core_lowlevel.py` (300 degenerate cases) |
| `derepFastq` (chunk-dependent tie order, IUPAC/N/lowercase sort order, NaN padding) | `derep_fastq` | exact (quals rtol 1e-12) | `parity/test_derep_sequences.py`, `parity/test_filter_derep_coverage.py` |
| `seqComplexity`, `isPhiX`, `rc` | `seq_complexity`, `is_phix`, `reverse_complement` | exact (complexity bit-identical, incl. uniform-spectrum threshold cases) | idem, `parity/test_core_lowlevel.py` |
| `qa` / `plotQualityProfile` summaries | `io.quality_profile`, `plotting.quality_profile_stats` | exact (first n reads, deviation #3) | `parity/test_derep_sequences.py`, `unit/test_plotting_align.py` |
| `nwalign`, `nwhamming`, `C_eval_pair`, `C_pair_consensus`, k-mer distances | `align.*`, `_core.*` | exact (through the public `align` API) | `parity/test_core_lowlevel.py`, `parity/test_align_api.py` |
| `C_matrixEE`, `C_matchRef`, `ppois` | `_core.matrix_ee`, `match_ref`, `ppois_upper` | exact; ppois rel ≤ 1e-12 (≤ 2.5e-12 below 1e-292, deviation #1), zero/non-zero pattern identical | `parity/test_core_lowlevel.py` |
| `loessErrfun`, `noqualErrfun`, `makeBinnedQualErrfun`, `PacBioErrfun`, `inflateErr`, `getErrors` | error model classes, `inflate_err`, `get_errors` | rtol 1e-12 (loess: deviation #4, incl. the few-quality-scores fallback) | `parity/test_error_models.py`, `parity/test_loess_fallback.py` |
| `learnErrors` (Illumina F/R, noqual, nbases, PacBio, binned-quality singletons) | `learn_errors` / `ErrorLearner` | exact transitions, error rates rtol 1e-10 | `parity/test_dada.py`, `parity/test_loess_fallback.py` |
| `dada` (default, pooled, pseudo-pooled, priors, 15 option variants, PacBio, selfConsist) | `DADA` | exact, incl. `selfConsist=TRUE` (cluster order, every column, `birth_subs`) | `parity/test_dada.py`, `parity/test_core_dada.py`, `parity/test_review_python.py` (priors as string/path) |
| `mergePairs` (rejects, mismatches, minOverlap, concatenate, trimOverhang, propagateCol) | `merge_pairs` | exact | `parity/test_merge_seqtab.py` |
| `makeSequenceTable`, `mergeSequenceTables` (repeats, tryRC), `collapseNoMismatch` | same names, snake_case | exact, incl. tie order with duplicate sequences | `parity/test_merge_seqtab.py`, `parity/test_review_python.py` |
| `isBimera`, `isBimeraDenovo`, `isBimeraDenovoTable`, `removeBimeraDenovo` (3 methods), `isShiftDenovo` | `is_bimera*`, `BimeraRemover`, `remove_bimera_denovo`, `is_shift_denovo` | exact, on tables where methods and options disagree | `parity/test_chimeras.py`, `parity/test_chimera_coverage.py` |
| `assignTaxonomy` | `TaxonomyClassifier` / `assign_taxonomy` | top hits exact; bootstraps within binomial noise (deviation #2) | `parity/test_taxonomy.py` |
| `assignSpecies`, `addSpecies` | `assign_species` / `SpeciesAssigner`, `add_species` | exact (`allow_multiple` join order: deviation #7) | `parity/test_taxonomy.py`, `parity/test_review_python.py` |
| End-to-end vignette pipeline with learned errors (track table, seqtab.nochim, taxonomy) | — | exact; taxonomy per §8.3 at every level | `e2e/test_pipeline.py` |

**Thread-count independence** (1 vs 2/4 threads) is tested for:

- `dada_uniques` and DADA (incl. pseudo- and full pooling);
- `table_bimera`, `is_bimera_denovo`, `is_shift_denovo` and every removal method;
- `collapse_no_mismatch`;
- `filter_and_trim` and the K1 kernel;
- `seq_complexity`, `remove_primers` (incl. indels), `merge_pairs`, `assign_species` and
  taxonomy.

**Memory and thread safety.**
- AddressSanitizer: `scripts/asan_check.sh`.
- ThreadSanitizer: concurrency, batch, cache, nested-pool and taxonomy tests.
- Both report no errors in the C++ code.

**Full data sets** (`benchmarks/dataset/REPORT.md`, `REPORT_stability.md`):
- **Outputs vs R:**
  - MiSeq SOP (20 samples) and the full Stability data set (362 samples, 3.65M read pairs): every
    exact output identical to R, incl. per-sample denoising and merged pairs;
  - error rates agree to rel ≤ 1.4e-14;
  - taxonomy meets the §8.3 criterion.
- **Determinism:** the same runs reproduce bit-identically across repeats and thread counts.

**Performance** (README, `benchmarks/`; numbers at `5e6453a`):

| data set | threads | pydada2 vs R | pydada2 vs previous version (`dff04da`) |
|---|---|---|---|
| MiSeq SOP | 1 | 6.0x faster | 3.6x faster |
| MiSeq SOP | 8 | 12.1x faster | 6.6x faster |
| Stability | 8 | 4.5x faster | 2.7x faster |

**Micro-benchmarks:**
- **vs R:** every stage is faster except `make_sequence_table`, which is level (~0.1 ms per
  call).
- **vs the previous version:** no stage is more than marginally slower.
