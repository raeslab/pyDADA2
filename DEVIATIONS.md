# DEVIATIONS.md — known differences between pydada2 and R dada2 1.36

Every known difference from the R reference is listed here (spec §8.3). Each entry gives:

- the root cause;
- the outputs it affects;
- its measured magnitude;
- the test that pins that magnitude.

| # | Area | Difference | Cause | Magnitude | Pinned by |
|---|---|---|---|---|---|
| 1 | ppois (abundance p-values) | **Below the mode (λ < n+1):** clean-room P(X ≥ m) = dpois(m; λ)·Σⱼ Πᵢ≤ⱼ λ/(m+i), compensated sum. dpois uses Loader's saddle-point form with a series `bd0`, is computed directly as λᵐe^−λ/m! for λ < 1, and is exponentiated once in log space below `DBL_MIN`. **At or above the mode:** Boost.Math `gamma_p`. | no GPL code (spec §5.4); plain `gamma_p` underflowed to 0 where R returns a subnormal and lost 1e-10 at large n | **Vs a 45-digit reference:** ≤ 5.4e-13 for n ∈ [0, 2e6] and λ down to 1e-320; zero/non-zero pattern identical to R. **Vs R:** ≤ 1e-12 for p ≥ 1e-292. Below that, R itself (exp(log p)) is up to 1.9e-12 off the reference, so the tolerance is 2.5e-12 there, and at least one subnormal ulp. | `tests/parity/test_core_lowlevel.py::test_ppois_upper`, `::test_ppois_upper_extended` (4,358-point R grid incl. large n and the underflow boundary) |
| 2 | taxonomy RNG | Bootstrap uniforms come from `std::mt19937_64(random_state)` instead of R's RNG; ties are broken deterministically. | spec §5.4 | **Top hits:** identical. **Bootstrap \|Δ\|:** median 0, max 15. **At `min_boot=50`:** 8/109 sequences differ at one level (7–9 across R regenerations, because R's own tie-breaking is not reproducible), all with an R bootstrap of 42–55. | `tests/parity/test_taxonomy.py::test_deviation_2_magnitude` (median 0, max ≤ 20, ≤ 12 sequences, \|R boot − 50\| ≤ 10) |
| 3 | `quality_profile` | Uses the first `n` reads. | R's `qa()` samples `n` reads at random from larger files | identical on files with ≤ n reads (all fixtures) | `test_derep_sequences.py::test_quality_profile_parity` |
| 4 | loess error model | `skmisc.loess` agrees with R's `loess` to a few ULP. When skmisc aborts on R's loess *warnings* (≤ ~8 observed quality scores per nucleotide, e.g. binned qualities with mostly singletons), pydada2 uses its own NumPy port of the same netlib dloess 1-D algorithm (`pydada2._loess`). | different builds of the netlib dloess code; skmisc turns warnings into errors | **skmisc path:** error rates rel ≤ 1e-12 (1.4e-14 on MiSeq SOP). **Fallback:** rel ≤ 1.2e-12 vs R on 56 degenerate transition matrices. **Where R itself fails** (< 3 scores), pydada2 raises too. | `test_error_models.py::test_loess`, `tests/parity/test_loess_fallback.py` |
| 5 | DADA cluster order with `use_quals=False` | R's `ppois` can return exactly equal p for λ values one ulp apart; pydada2 separates them. In `b_bud`, R then takes the raw it iterates first, while pydada2 takes the smaller p. | last-ulp differences of the p-value algorithms (#1) | **Ties:** 498 of 9,000 ulp-adjacent λ pairs tie in R but not in pydada2; never the reverse. **Effect when a tie is hit:** cluster *creation order* can change, i.e. the row order of `clusters`/`map`/`birth_subs` and the post-hoc `pval` of the swapped clusters (~1 %). ASVs, abundances and membership stay identical. **Seen:** 2–4 rows reordered on synthetic sets of 63–283 clusters through the private `use_quals=False` path only. **Not seen:** with quality-aware λ (24 synthetic runs, all fixtures, MiSeq SOP, Stability). | differential runs in REVIEW.md §1.2; quality-aware paths pinned exactly by `test_dada.py` |
| 6 | `learn_errors(randomize=True)` | Sample order comes from `numpy.random.default_rng(random_state).permutation`, not R's `sample()`. | different RNG | Different samples are used than in R. Identical with `randomize=False` (the default) or when all samples are needed. | `test_dada.py::test_learn_errors_*` (randomize=False) |
| 7 | `assign_species`/`add_species` with `allow_multiple` | Multiple names are joined in code-point order (`sorted`). | R's `sort()` uses the locale's collation | differs only for names that differ in case or punctuation, e.g. `Anthracis/anthracis` vs R (en_US) `anthracis/Anthracis` | `tests/parity/test_review_python.py` |
| 8 | `BimeraRemover(method="per-sample").is_bimera_` | Sequences absent from a sample are not tested and are reported False. | speed (89 s → 1 s on 30 × 3000) | R computes flags for them internally but never returns them; the output table is identical | `tests/unit/test_perf_equivalence.py::test_per_sample_bimeras_skip_absent_sequences` |

The performance work of REVIEW.md §5 changes no output. Each item below was checked bit-identical
on MiSeq SOP (every output, taxonomy included at the same seed), at 1 and 8 threads, and on
the Stability data set:

- taxonomy model and classifier;
- alignment cache;
- sample concurrency;
- persistent thread pool;
- phiX matching;
- merge batching.

## Upstream non-determinism (not a pydada2 deviation)

- **Taxonomy tie-breaking.** R dada2 breaks ties between equally likely genera with
  `std::random_device`, so `set.seed()` does not make `assignTaxonomy` reproducible.
  - Measured R-vs-R on the fixture data: bootstrap |Δ| ≤ 2 (median 0), and `minBoot=50`
    assignments 99/100 identical.
  - Reruns can also change the `minBoot=50` disagreement count (deviation #2).
  - pydada2 is reproducible for a given `random_state`.

## Documented R quirks preserved on purpose

- **`minQ` in `fastqFilter`/`filterAndTrim`:** it only applies when `minQ > truncQ`, and reads
  are kept when `min(qual) > minQ` (the docs say reads *below* minQ are discarded).
- **`maxN`** counts every non-ACGT letter (ShortRead `nFilter`). Pinned by
  `filter_coverage.json`.
- **`derep_fastq(n=...)`** raises "Only zero-length sequences detected during dereplication."
  when any chunk of `n` reads contains only zero-length reads, as R does. Pinned by
  `test_properties.py::test_derep_zero_length_chunk_is_an_error`.
- **`is_shift_denovo(flag_subseqs=True)`** has no effect, because R's `isShiftDenovo` never
  passes `flagSubseqs` on to `isShift`. pydada2 warns when the flag is set. Pinned by
  `test_chimera_coverage.py::test_is_shift_denovo[subseqs]`.

## R quirks not reproduced

- **Quality encoding.** pydada2 detects the encoding once per input chunk with ShortRead's rule,
  which gives identical results for phred+33 data. R's output can instead depend on `n`:
  - `filterAndTrim(orient.fwd=...)` re-detects the encoding on the re-oriented subset only;
  - `fastqPairedFilter(matchIDs=TRUE)` forces phred+33 on chunks after the first.

## Inputs R rejects but pydada2 accepts

- **`seq_complexity` with `window` ≥ the sequence length.** R fails with "wrong sign in 'by'
  argument"; pydada2 returns the whole-sequence value.
- **Reads containing U, X or `*`.** R fails with "_DNAencode(): invalid DNAString input
  character" in `derepFastq`/`filterAndTrim`. pydada2 counts them as non-ACGT letters, e.g. for
  `max_n`.
- **Reads shorter than the orientation sequence in `orient_fwd`.** In paired mode they are
  treated as not matching, as in single-end mode (see below). R errors in both modes.

## Behaviour changes at the C++ level (no effect on valid inputs)

All entries are pinned by `tests/unit/test_behaviour_changes.py`.

**Raise instead of misbehaving**

- **Unequal-length aligned strings** in `eval_pair`/`pair_consensus` raise `ValueError`
  (R: warning + `NULL`). Empty strings give zero counts (R: out-of-bounds read).
- **`remove_primers` with an empty primer** raises `ValueError`.
- **DADA on sequences ≥ 9999 nt** raises `ValueError`, as upstream.
- **DADA on a sample without sequences** raises "Sample i has zero input sequences."
  (R fails with "subscript out of bounds").
- **`dada_uniques` input checks:** it rejects error matrices that are not finite or not in
  [0, 1] (R enforces the same in `getErrors(enforce=TRUE)`). It also rejects
  `use_quals=True` without qualities (upstream read a NULL array).
- **`_core.nwalign` on non-ACGT input** raises `ValueError` (upstream read out of bounds).
  `nwvec` (R's `nwalign(vec=TRUE)`) still accepts any letters, as R does.
- **`table_bimera`/`matrix_ee`** reject values outside int32 instead of wrapping.
- **k-mer functions** reject k-mer sizes outside their range before allocating.
- **Transition range check in `error.cpp`:** it now throws where upstream wrote out of bounds.

**Upstream bugs fixed**

- **Prior-born clusters:** clusters born from priors (`birth_type` "P", used with `priors` or
  `pool="pseudo"`) now report `birth_from` as the cluster they were born from. R/upstream leave
  it uninitialised (random values), so there is no parity to keep.
- **`C_is_bimera`:** the upstream memory leak (alignment not freed on `continue`) is fixed.
- **Taxonomy:**
  - sequences longer than 9999 nt no longer overflow fixed-size stack buffers;
  - queries shorter than the k-mer size no longer blow up the uniform-buffer size.
- **`rm_lowcomplex`:** reads without any A/C/G/T k-mer (complexity NA in R) are removed
  (R would index with NA).
- **Single-end `orient_fwd`:** reads shorter than the orientation sequence are treated as not
  matching (R: a `narrow()` error aborts the whole file).

**Platform consistency**

- **SSE modes off x86** are emulated exactly. R on ARM forces `SSE=0`, which can change the
  k-mer/gapless screens. pydada2 matches R-on-x86 everywhere, including `kord_dist(SSE=1)`.
- **`seq_complexity` sums** reproduce R's long-double `sum()`. MSVC has no extended long
  double, so it uses a double-double accumulator instead. Compiled on Linux, that branch is
  bit-identical to R on all 1,179 fixture values. A last-bit difference from R-on-Windows
  (80-bit mingw) is possible in principle but was not observed.
