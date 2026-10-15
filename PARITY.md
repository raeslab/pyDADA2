# PARITY.md — pydada2 vs R dada2 1.36 on the bundled example data

Sign-off for spec phase 10. Every R NAMESPACE export maps to a Python object (checked by
`tests/unit/test_api_coverage.py`). The parity tests (`uv run pytest -m parity`, 141 tests) compare
pydada2 against fixtures produced by `tests/reference/generate_reference.R` with R 4.6.1 /
dada2 1.36.0.

Verified locally on **Linux x86_64** (Python 3.11 from the sdist, 3.12 and 3.13 from the abi3
wheel). macOS (arm64, x86_64) and Windows are covered by the CI matrix
(`.github/workflows/ci.yml`), which runs the same parity suite once the repository is pushed to
raeslab/pyDADA2. That run has not happened yet.

Status legend: **exact** = identical sequences, counts, maps, flags and orders, with doubles
within rtol 1e-10 (most are bit-identical); **tolerance** = within the documented spec §8.3
tolerance; **deviation #n** = see `DEVIATIONS.md`.

| R function | Python | Status | Tests |
|---|---|---|---|
| `filterAndTrim` (single, paired, chunked, orient.fwd, matchIDs, all filters) | `filter_and_trim` | exact (output records via canonical md5) | `parity/test_filtering.py` (14 configurations) |
| `fastqFilter`, `fastqPairedFilter` | `fastq_filter`, `fastq_paired_filter` | exact (via the above) | idem |
| `removePrimers` (mismatch and indel modes, orient, trim options) | `remove_primers` | exact | `parity/test_filtering.py`; Biostrings matching on 300 random cases in `match_pattern.json` |
| `derepFastq` (incl. chunk-dependent tie order, IUPAC sort order, NaN padding) | `derep_fastq` | exact (quals rtol 1e-12) | `parity/test_derep_sequences.py` |
| `seqComplexity`, `isPhiX`, `rc` | `seq_complexity`, `is_phix`, `reverse_complement` | exact / rtol 1e-12 | idem |
| `qa` / `plotQualityProfile` summaries | `io.quality_profile`, `plotting.quality_profile_stats` | exact (first n reads, deviation #3) | `parity/test_derep_sequences.py`, `unit/test_plotting_align.py` |
| `nwalign`, `nwhamming`, `C_eval_pair`, `C_pair_consensus`, k-mer distances | `align.*`, `_core.*` | exact | `parity/test_core_lowlevel.py` |
| `C_matrixEE`, `C_matchRef`, `ppois` | `_core.matrix_ee`, `match_ref`, `ppois_upper` | exact / rel 1e-12 (deviation #1 for subnormals) | idem |
| `loessErrfun`, `noqualErrfun`, `makeBinnedQualErrfun`, `PacBioErrfun`, `inflateErr`, `getErrors` | error model classes, `inflate_err`, `get_errors` | rtol 1e-12 (loess: deviation #4) | `parity/test_error_models.py` |
| `learnErrors` (Illumina F/R, noqual, nbases, PacBio) | `learn_errors` / `ErrorLearner` | exact transitions, error rates rtol 1e-10 | `parity/test_dada.py` |
| `dada` (default, pooled, pseudo-pooled, priors, 15 option variants, PacBio) | `DADA` | exact | `parity/test_dada.py`, `parity/test_core_dada.py` |
| `dada(selfConsist=TRUE)` with a supplied error matrix | `DADA(self_consistent=True)` | deviation #4 (cluster creation order; same ASVs, abundances, membership) | `parity/test_dada.py` |
| `mergePairs` (rejects, mismatches, minOverlap, concatenate, trimOverhang, propagateCol) | `merge_pairs` | exact | `parity/test_merge_seqtab.py` |
| `makeSequenceTable`, `mergeSequenceTables` (repeats, tryRC), `collapseNoMismatch` | same names, snake_case | exact | idem |
| `isBimera`, `isBimeraDenovo`, `isBimeraDenovoTable`, `removeBimeraDenovo` (3 methods), `isShiftDenovo` | `is_bimera*`, `BimeraRemover`, `remove_bimera_denovo`, `is_shift_denovo` | exact | `parity/test_chimeras.py`, `parity/test_core_lowlevel.py` |
| `assignTaxonomy` | `TaxonomyClassifier` / `assign_taxonomy` | top hits exact; bootstraps within binomial noise (deviation #2) | `parity/test_taxonomy.py` |
| `assignSpecies`, `addSpecies` | `assign_species` / `SpeciesAssigner`, `add_species` | exact | idem |
| End-to-end vignette pipeline (track table, seqtab.nochim) | — | exact | `e2e/test_pipeline.py` |

Thread-count independence (1 vs 2/4 threads) is tested for `dada_uniques`, DADA (pseudo-pooling),
`table_bimera`, `collapse_no_mismatch`, `filter_and_trim` and taxonomy. Memory safety is checked
with AddressSanitizer/LeakSanitizer (`scripts/asan_check.sh`): no errors, no leaks in the C++ code.

Phase 11 repeats the comparison on the full MiSeq SOP dataset with SILVA 138.2
(`benchmarks/dataset/REPORT.md`). All exact outputs are identical, error rates agree to rel 1.4e-14,
and taxonomy meets the §8.3 criterion (229/232 ASVs identical; the other 3 are borderline calls).

Performance on the same data: see `benchmarks/micro/REPORT.md`. Every stage is as fast as R or
faster, so no stage exceeds the 1.5x budget.
