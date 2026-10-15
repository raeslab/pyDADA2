# REVIEW.md: review of the pyDADA2 port (R dada2 1.36 → C++/Python)

Reviewed: commit `57811e2` (phases 0–11), 2026-10-15. Scope: the de-R'd C++ core (`cpp/core`),
the R→C++ kernels (`cpp/core/ported`), the nanobind bindings, the Python package, tests and
fixtures, packaging/CI and both benchmark suites. The goal it is measured against: **faster than R,
functionally equivalent, (near-)exact output, with differences only from random numbers**.

## How this was reviewed

- Every `cpp/core` file was diffed against `reference/r-dada2/src`. Every kernel and Python driver
  was read against the R function it replaces.
- **Differential testing against R.** R dada2 1.36 is installed in the dev container, so the R and
  Python versions were run on the same inputs:
  - several hundred synthetic edge-case inputs, covering filtering, derep, primers, merging,
    sequence tables, chimeras, species, error models and `dada_uniques` (17 option variants × 5
    data sets, plus 24 binned-quality runs);
  - the bundled example data;
  - the MiSeq SOP data.
- **Sanitizers and stress.** The C++ core was built separately under ASan/UBSan (-O0) and TSan, and
  run through threaded stress tests.
- **Profiling.** The MiSeq SOP pipeline was profiled stage by stage at 1 and 8 threads.
- **Prototypes.** Two prototypes were built in scratch copies, outside the repo and the venv:
  - the ppois fix (§1.1), run against the whole test suite;
  - a combined performance prototype (§5.2), whose MiSeq SOP outputs were compared with HEAD's.
- **Gates at HEAD.** The full test suite, ruff, ruff-format and mypy were run.

The review itself changed nothing in the repository except adding this file. The benchmark
scripts were extended afterwards for the larger-dataset run (README, "Larger dataset"). Reproducer
scripts and prototype diffs were kept in a temporary session scratchpad, not in the repo. The
key ppois patch is reproduced below so that this document stands alone.

## Status after the fixes (2026-10-15)

The findings below were addressed on branches `fix/numerics`, `fix/python`, `fix/tests` and
`fix/perf`, merged into `pydada2-main`, followed by a few direct commits. The rest of this
document is the original review, unchanged.

**Verification of the merged code:**
- 1067 tests pass (874 parity), coverage 97 %, ruff/mypy clean.
- The fixture generator reproduces every committed fixture byte-for-byte, apart from the two
  documented nondeterministic files.
- **MiSeq SOP:** every exact output is identical to R, taxonomy within §8.3, at 1 and 8 threads.
  All outputs, taxonomy included, are bit-identical to the previous pydada2.
- **Stability data set (362 samples):** every exact output is identical to R, **including the
  per-sample denoising and merged pairs that differed before §1.1 was fixed**.

| Finding | Status |
|---|---|
| §1.1 ppois underflow → cluster order | **Fixed.** Clean-room log-space tail. The `self_consist` parity test is now strict and ordered, and the Stability mismatches are gone. |
| §1.2a ppois accuracy at large n | **Fixed:** max error vs a 45-digit reference 2.0e-10 → 2.7e-13 (Loader saddle-point + series). |
| §1.2b ulp ties with `use_quals=False` | **Documented** (DEVIATIONS #5). Not seen with quality-aware λ. |
| §1.3 loess failures | **Fixed:** NumPy port of netlib's 1-D loess used when skmisc aborts; ≤ 1.2e-12 vs R. |
| §1.4 seqtab tie order; §1.5 priors | **Fixed**, pinned by R fixtures. |
| §1.6 a, b, c, f, g, h, i, j, k; pickling; numpy `n_jobs` | **Fixed**, pinned by tests. b matches R: `flag_subseqs` has no effect and warns. h: `ErrorLearner` raises on conflicting `params`. |
| §1.6 d, e, l, m | **Documented** in DEVIATIONS. |
| §1.7 `birth_from` | **Fixed** (upstream bug). Other carry-overs unchanged. |
| §2.1 thread-unsafe handlers | **Fixed:** thread-local with save/restore; TSan clean. |
| §2.2 safety items | **Fixed:** err/quals validation, NaN-safe checks, `B` leak, thread-spawn failure, `kord_dist` off x86, k-mer size checks, `nwalign` ACGT check, int32 range checks, Dereplicator lock. |
| §3 T1–T8 | **Fixed:** discriminating filter/derep/chimera fixtures, `match_pattern` test, Hypothesis and determinism tests, behaviour-change tests, R warnings compared, coverage gate on the full run only, DEVIATIONS/PARITY/CLAUDE.md rewritten. |
| §4.1 CI/packaging | **Fixed in the workflow, but not yet run on GitHub.** Covers the slow step, `macos-15-intel`, a working sdist job, wheels tested with the suite incl. aarch64, SHA-pinned actions, `scikit-misc>=0.5.1`, `xopen`, C++ tests, offline Boost. **Open, for the maintainers:** licence version (spec Q9) and the PyPI name (Q1). |
| §4.2 benchmark methodology | **Fixed:** looped sub-ms stages, equal work and threads, gc caveats. Reports regenerated (README). |
| §5 P1–P7 and most smaller items | **Implemented**, bit-identical outputs; numbers in README. P8 (AVX2) was measured with no gain and dropped. Reusing the per-round raw/k-mer setup was replaced by parallelising it. |

The final old-vs-new benchmarks found regressions introduced by the fixes themselves. All are
fixed, with identical outputs:

| Regression | Cause | Fix |
|---|---|---|
| `add_species` slower | the verbatim-header change read the 530 MB SILVA species file twice; the binding copied ~0.5 GB of reference strings on every call | one dnaio-equivalent FASTA pass that keeps verbatim headers; references passed to C++ as string views |
| `assign_species` (micro) 0.8 → 12 ms | xopen started a decompression thread on every open | threaded decompression only for files ≥ 4 MB |
| `make_sequence_table` slower than R | pandas overhead | columns ordered in NumPy before building the DataFrame; now level with R |
| `derep_fastq` ~6 % slower on Stability | an `s if s.isupper() else s.upper()` guard meant to skip copies is 3–4x slower than `str.upper()`, which has an ASCII fast path that `isupper()` lacks | reverted to plain `upper()` |

A further apparent regression in the micro suite was an artefact, not a code change. The main
checkout sits on a slow WSL bind mount, where each `stat`/`open` costs ~1 ms. Both pydada2
versions are therefore now benchmarked from `/tmp` worktrees.

Still open:
- the MANIFEST hard-codes `upstream_commit` even for Bioconductor installs;
- CI has to run once on GitHub;
- macOS/Windows parity is therefore still unverified.

## Verdict

**The conversion works, and the core is faithful.** The C++ diffs are glue apart from the items
below:

- **dada_uniques:** clustering, birth_subs, subqual, clusterquals, map and p-values were
  bit-identical to R across 17 option variants on 5 synthetic data sets. The two exceptions are
  §1.1 and §1.2.
- **Kernels:** filtering (45 configurations, incl. chunking, matchIDs, orient, phred64), derep,
  primers, merging, `table_bimera`, alignment and k-mer functions, and taxonomy top hits all
  matched R exactly.
- **Gates at HEAD:** all 257 tests pass with 96.6 % coverage; ruff, ruff-format and mypy
  (strict) are clean.

**But the claim "every exact-category output is identical to R on the example data" does not
hold yet.**

- One real parity bug in the core changes cluster order and `birth_*` values on the bundled data
  (§1.1). It is documented in DEVIATIONS.md under the wrong cause. A ~10-line clean-room fix,
  prototyped and tested here, removes it.
- A handful of medium edge-case divergences remain (§1.3–§1.6).
- The verification story has gaps: several filters are never exercised by any fixture (§3), and
  CI cannot go green as written (§4).

**On speed**, pydada2 is already 1.7–1.8x faster than R end to end on MiSeq SOP. It could
be about **11x faster than R at 8 threads, and 6x at 1 thread, with bit-identical output**:

- **Prototype:** combining the changes in §5.2 brings MiSeq SOP from 78.4 s to 12.7 s at 8
  threads and from 151.5 s to 40.6 s at 1 thread. Every output value in its `results.json` was
  checked to be identical to HEAD's.
- **Taxonomy fit:** 62 % of today's 8-thread runtime is a single-threaded model build that does
  ~59 billion redundant byte operations. Fixing it alone saves 45 s.
- **learn_errors:** 82 % of its alignments are repeats across self-consistency rounds.
- **Blocker:** running samples concurrently, one of the largest wins, first needs the
  thread-safety bug in the runtime handlers fixed (§2.1).

Severity scale: **High** = wrong results or crashes in plausible use · **Medium** = divergence from
R or a crash in less common but legitimate use · **Low** = edge cases, API traps,
documentation.

---

## 1. Correctness and parity with R

### 1.1 HIGH: `dada()` creates clusters in a different order than R (Poisson tail underflow)

**Mechanism**
- `ppois_upper` (`cpp/core/stats.cpp:25-32`) uses Boost `gamma_p`, which returns **exactly 0**
  where R's `ppois(..., lower.tail=FALSE)` returns a subnormal (e.g. 4.9e-324).
- `b_bud` picks the most significant raw with `p < min || (p == min && reads > min_reads)`
  (`cluster.cpp:185-186`).
- So a strict "R raw has the smaller p" becomes a tie in the port, and the tie goes to the raw
  with more reads.
- The final per-unique p-values agree to 2e-15, so the difference only shows up inside the
  clustering loop.

**Verified:**
- **Underflow boundary.** On 133,437 (n, λ) points whose upper tail lies in 1e-330…1e-295, there
  are 154 cases where R gives > 0 and Boost gives 0, and none the other way.
- **Example data.** R's own final self-consistent error matrix was fed to plain, non-self-consistent
  `dada()` in both languages, using the identical derep of `sam1F`.
  - R produced cluster abundances `386 291 48 37 115 45 7 57 4`, the same order as its
    self-consistent run.
  - pydada2 produced `386 291 48 115 37 …`.
  - The swap moves `map`, `birth_from` and `birth_subs.clust`, and changes `birth_fold` by up to
    24 % relative (1.548e6 vs 1.664e6).

**Consequences:**
- DEVIATIONS.md #4 blames skmisc's loess for this. That is wrong: the error matrices are
  identical, and the difference comes from the core.
- The deviation is not limited to `self_consistent=True`. Any `dada()` call can hit it,
  including the standard vignette path with a learned error matrix. Phase 11 happened not to.
- It is pinned only by `_compare_unordered` (`tests/parity/test_dada.py:96-117`). That function
  skips `pval`, `birth_pval`, `birth_fold`, `birth_qave`, `quality` and `birth_subs`.
- §8.3 condition 4 ("does not change any exact-category output on the example data") is
  violated for cluster order and `map`.

**Fix: prototyped and tested.**
- When `gamma_p` returns a value below `DBL_MIN` (and λ < n+1), recompute the tail in log space
  from the textbook series P(X ≥ m) = dpois(m; λ) · Σⱼ λʲ/((m+1)…(m+j)), then exponentiate.
- This is a clean-room formula, not R code.
- Results:
  - The 154 zero-vs-subnormal disagreements drop to 0.
  - Subnormal mismatches drop from 2498 to 343; the remaining ones are last-ulp differences
    among subnormals.
  - pydada2 then reproduces R's cluster order exactly.
  - With `UNORDERED = set()`, i.e. the self-consistent case compared **strictly and in order**,
    the whole parity + unit + e2e suite passes: 239 passed.
- The prototype, in `ppois_upper`, after the `gamma_p` call:

  ```cpp
  const double p = boost::math::gamma_p(static_cast<double>(n) + 1.0, lambda, Policy());
  if (p >= std::numeric_limits<double>::min() || lambda >= n + 1.0) return p;
  // Tiny tail in log space: P(X >= m) = dpois(m; lambda) * sum_j lambda^j m!/(m+j)!
  const double m = static_cast<double>(n) + 1.0;
  const double ld = m * std::log(lambda) - lambda - boost::math::lgamma(m + 1.0, Policy());
  double s = 1.0, t = 1.0;
  for (int j = 1; j < 100000; ++j) { t *= lambda / (m + j); s += t; if (t < 1e-17 * s) break; }
  return std::exp(ld + std::log(s));
  ```
- Two details for the real implementation:
  - use `boost::math::lgamma`, not `std::lgamma` (not thread-safe on some libcs);
  - keep the normal range on the existing `gamma_p` path, so validated results don't move.
- After the fix:
  - make the self-consistent parity test ordered;
  - add a dense near-underflow fixture that asserts zero/non-zero agreement;
  - rewrite DEVIATIONS #1 and #4.

### 1.2 MEDIUM: `ppois_upper` accuracy for large counts; ulp-level ties in `b_bud`

- **(a) Accuracy for large counts.** Boost `gamma_p` misses spec §5.4's rel ≤ 1e-12 for large
  n, and Boost is the inaccurate side (checked with 80-digit mpmath):

  | n | Boost vs R, relative | Boost error vs mpmath | R error vs mpmath |
  |---|---|---|---|
  | 1e4 | 1.2e-12 | | |
  | 1e5 | 1.0e-11 | 9.9e-12 | 1.1e-13 |
  | 1e6 | 1.4e-10 | 1.4e-10 | 7e-14 |

  The fixture grid (`generate_reference.R:400-401`) never reaches this region (λ/n ∈ [0.6, 1.4],
  n ≥ 1e4). So DEVIATIONS #1 ("rel ≤ 1e-12") is inaccurate for uniques with ≥ 1e4 reads.
- **(b) Ulp-level ties.** With `use_quals=False`, two λ values one ulp apart get *equal* p in R
  but different p in Boost, so a different raw buds first. Observed effects:
  - cluster rows 108/109 and 23/24 swapped on two synthetic data sets;
  - same ASVs and abundances, and post-hoc `pval` ~1 % different.
  - On ulp-adjacent λ pairs, 496 of 9000 comparisons are ties in R but not in Boost.
  - Not seen with quality-aware λ (24 runs).
- **Fix:**
  - Extend the fixture grid to n ∈ [1e3, 2e6] with λ/n ∈ [0.5, 1.5].
  - Consider a clean-room large-n path (Loader's saddle-point dpois plus a continued fraction)
    so that p-values are as accurate as R's.
  - Document the residual tie sensitivity in DEVIATIONS.

### 1.3 MEDIUM: loess error model fails where R returns an estimate

- **Where:** `error_models.py:103-134` (skmisc loess).
- **Cause:** when only ~3–7 quality columns contain counts, skmisc raises `ValueError` on
  netlib's *warnings* ("near singularities", "Chernobyl! trL>n"). R only warns, and `dada()`
  suppresses it.
- **Symptoms:**
  - `learn_errors` fails with a misleading "Error matrix is None".
  - `DADA(...).err_out_` is silently `None`.
- **Repro (re-run here):** 5 transition matrices with binned qualities (2/12/23/37 etc.); R
  returns values for all 5, pydada2 raises for all 5. In a random sweep, R succeeded on 298/300
  and pydada2 failed on 83 of those.
- **Realistic trigger:** binned-quality NovaSeq-style data with mostly singleton uniques.
- **Fix:** vendor the netlib dloess that skmisc wraps (BSD-style licence, *not* R's GPL copy) with
  warnings made non-fatal, or patch skmisc's error handling. Pin the result with the 5 configs.

### 1.4 MEDIUM: `make_sequence_table` column order differs from R when a sample has duplicate sequences

- **Where:** `sequences.py:215` collapses duplicates with `groupby(level=0, sort=False)`
  (first-appearance order).
- **R behaviour:** `getUniques` uses `tapply(..., names, sum)`, which sorts by sequence. The
  later stable sort by total then breaks ties differently.
- **Synthetic repro (re-run here):** R `GGGG ACAC CCCC TTTT AAAA` vs pydada2
  `GGGG ACAC TTTT CCCC AAAA`.
- **Real data:** MiSeq SOP with `merge_pairs(max_mismatch=2)` (64 duplicate merged sequences):
  42/398 column positions differ, with identical counts.
- **Fix:** `sort=True`, which makes the SOP case identical to R. In `get_sequences` on a
  string list, de-duplicate before upper-casing as R does. Add a fixture with duplicate mergers.

### 1.5 MEDIUM: `DADA(priors=...)` ignores a single sequence string or a FASTA path

- **Where:** `denoise.py:451` calls `get_sequences(list(self.priors))`. A `str` becomes a list
  of characters, so nothing matches.
- **R behaviour:** `getSequences(priors)` reads a file path or treats the string as one sequence.
- **Repro:** sam1F with `OMEGA_A=1e-60`.

  | priors | R | pydada2 |
  |---|---|---|
  | one string | 8 ASVs | 7 ASVs |
  | FASTA path | 8 ASVs | 7 ASVs |

- **Fix:** `get_sequences(self.priors)`. This also makes Dereplication/DadaResult priors work.

### 1.6 LOW–MEDIUM: smaller divergences from R (all reproduced against R)

| # | Where | R | pydada2 | Suggested action |
|---|---|---|---|---|
| a | `complexity.cpp:50-54` | `sum()` accumulates in long double | double accumulation; only 2516/3116 values bit-identical. `NGAGTACCNN` → 5.999999999999998 vs R 6.0, so `rm_lowcomplex=6` drops a read R keeps | accumulate in `long double` (reproduces R on 3116/3116); note the MSVC gap |
| b | `chimeras.py:262-274` | `isShiftDenovo` never passes `flagSubseqs` to `isShift` (`chimeras.R:391`), so the flag has no effect | flag honoured: `F F T T T` vs R `F F F F T` | ignore it to match R, or list it in DEVIATIONS |
| c | `taxonomy.py:334, 415` | empty species field → `""`; `strsplit` drops trailing empty token | `None`; extra empty token gives `"/coli"` vs `"coli"` with `allow_multiple` | mirror R's tokenisation |
| d | `taxonomy.py:333` | `sort()` uses locale collation | code-point order (`Anthracis/anthracis` vs `anthracis/Anthracis`) | document (R is locale-dependent) |
| e | `filtering.py` orient/matchIDs | re-detects quality encoding on the oriented subset; forces phred33 on later chunks with `matchIDs` | detects once per chunk (more sensible) | document as R bugs (R's output depends on `n`) |
| f | `filtering.py:98` | `maxN=Inf` accepted | `OverflowError` from `int(inf)` | map ±inf to INT limits |
| g | `derep.py:133-136`, `denoise.py:292` | empty file → "Only zero-length sequences detected" | empty Dereplication, then numpy "zero-size array … fmin" in DADA | raise R's error |
| h | `error_models.py:526, 556-562` | `learnErrors` takes a named list; `OMEGA_C`/`MAX_CONSIST` are formals | a dict gives `FileNotFoundError: s1`; `params={"omega_c":…, "max_consist":…}` silently overridden | accept Mappings; let explicit params win or raise |
| i | `sequences.py:241-248` | `getSequences(taxa)` returns row names | tries to sum string columns → `ValueError` | detect non-numeric tables |
| j | `merge.py:209` | — | `propagate_col="n0"` iterates characters and is silently ignored | accept `str` |
| k | `primers.cpp:66` | `max.mismatch ≥ primer length`: Biostrings returns every overlapping window | different window set (308 degenerate fuzz cases; 0 mismatches in 2692 normal ones) | special-case `k ≥ m` |
| l | `error_models.py:529-531` | `randomize=TRUE` uses R's `sample()` | NumPy permutation | unavoidable; add to DEVIATIONS |
| m | various | errors on: `removePrimers` with no hits; paired `orient.fwd` with short reads; `seqComplexity(window ≥ len)`; U/X/`*` read letters | lenient | DEVIATIONS lists only the single-end orient case |

Also: a fitted `TaxonomyClassifier` cannot be pickled or deep-copied
(`cannot pickle 'pydada2._core.TaxonomyModel'`). That blocks joblib caching and process pools;
add `__getstate__`/`__setstate__`. `n_jobs=np.int64(2)` is rejected (`_parallel.py:31`).

### 1.7 Upstream bugs carried over (fix and document)

- **`birth_from` uninitialised for prior-born clusters** (`cluster.cpp:224-233`, upstream too).
  The "P" branch never sets it and `new Bi` default-initialises PODs. The same input in one
  process returned `birth_from` ∈ {0, 2, 2949164, 9175179, 408680624}. This affects every run
  with priors, including `pool="pseudo"`. R returns garbage too, so there is no parity to keep.
  Set it to `mini_prior` and use `new Bi()`.
- **Taxonomy uniforms:** windows of consecutive sequences overlap (offset `j*max_arraylen`), so
  neighbouring bootstraps are correlated. This is inherent to "same count/usage as upstream".
- Narrower carry-overs:
  - `get_lr` reads `al[0][pos]` before its `pos >= 0` check;
  - the int16 vectorised aligner can overflow above ~6.5 kb;
  - there is no interrupt check inside the first all-vs-centre pass.

---

## 2. Robustness and memory safety (C++/bindings)

### 2.1 HIGH: process-global log/interrupt handlers race and segfault under concurrent calls

- **Where:** `runtime.cpp:11-12, 33-41` and `module.cpp:65-84`. `g_log` and `g_interrupt` are
  plain globals.
- **Mechanism:** every `dada_uniques`/`assign_taxonomy` call installs them, then nulls them in
  `~ScopedHandlers`. Other Python threads, with the GIL released, are calling them at the same
  time.
- **Reproduced:**
  - SIGSEGV in 3/3 runs with verbose `dada_uniques` from 3 threads (`dada2_log ← b_bud`);
  - 2/3 runs with 8 threads calling `assign_taxonomy` and no logging (null `std::function` in
    `dada2_check_interrupt`);
  - a concurrent call silently removes another call's interrupt handler.
- **Who is affected:** the library itself calls the core serially. A user who parallelises
  samples with a `ThreadPoolExecutor` crashes, and that is natural since the GIL is released. It
  also blocks performance item P1 (§5).
- **Fix:** make both handlers `thread_local`, and have `ScopedHandlers` save and restore the
  previous values. They are only invoked on the calling thread, never inside `parallel_for`
  bodies.

### 2.2 MEDIUM–LOW: smaller safety items

All of these are only reachable through the private `_core` API, or after an exception.

- **Missing qualities:** `dada_uniques` with `use_quals=True` and no quals segfaults
  (`pval.cpp:122`, `error.cpp:173` loads before its `has_quals` check). Validate in the core.
- **NaN in the error matrix:** a NaN λ skips the `comp` push, and then `b_shuffle2` reads past
  `bi[0]->comp` (ASan heap-buffer-overflow, `cluster.cpp:117`). Validate that err is finite and
  in [0, 1], and write the range checks as `!(λ >= 0 && λ <= 1)`.
- **Leak on exception:** `run_dada` (`dada_main.cpp:248`) leaks the whole `B` when an exception
  escapes (interrupt, raising log callback, `dada2_stop`): +14.6 MiB per 60 aborted calls. Store
  it in `DadaState::bb` immediately.
- **Thread-spawn failure:** if `parallel_for` (`parallel.h:59`) fails to spawn a thread
  part-way, the process aborts via `std::terminate`. `n_threads` is uncapped. Catch the failure,
  run with the threads that did start, and cap at the number of chunks.
- **`kord_dist` off x86:** it uses scalar code when `SSE==1` (`evaluate.cpp:225` `&& X64`),
  contradicting "matches x86 everywhere". The exact emulation exists; drop `&& X64`.
- **Footguns:**
  - `kmer_dist(k ≥ 16)` is UB, and k=15 allocates 2×2 GiB before validating;
  - `_core.nwalign` on non-ACGT input reads out of bounds;
  - nanobind narrows int64 to int32 silently;
  - `Dereplicator.add_chunk` releases the GIL on a mutable object.

**Checked and fine:**
- row-major/column-major conversions of err, quals, trans and clusterquals;
- 0-based indexing;
- Fortran-ordered and non-contiguous arrays (converted, not misread);
- GIL released around all heavy calls;
- exception mapping;
- SSE2 emulation (bit-identical on 20,000 random cases, including wrap and saturation);
- MSVC x64 gets the SSE2 path;
- no GPL nmath code remains;
- ASan/UBSan clean on all dada variants apart from the two items above.

---

## 3. Verification (tests, fixtures, DEVIATIONS/PARITY)

What is solid:
- No skip or xfail anywhere in the parity suite.
- Exact outputs are compared as ordered lists.
- filterAndTrim compares ids, widths and the md5 of every record.
- Defaults match R: 171 default arguments across 24 exported functions.
- Fixtures regenerate byte-identically, except taxonomy and the MANIFEST timestamp.

Gaps:

| # | Sev. | Finding |
|---|---|---|
| T1 | High | **No fixture exercises `maxN`, `maxLen`, `rm.phix` or `rm.lowcomplex`.** Not one read in any input FASTQ contains a non-ACGT letter, and disabling each of these filters in the fixture configurations gives identical counts. PARITY.md's "all filters" and "IUPAC sort order" rest on stand-alone `isPhiX`/`seqComplexity` checks only. (Differential tests in this review covered them by hand and found them correct apart from §1.6a.) Add synthetic inputs with N/IUPAC bases, phiX reads, low-complexity reads, and PacBio configs where `minLen`/`maxLen` bite. |
| T2 | Medium | `match_pattern.json` (300 Biostrings cases, cited by PARITY.md and CLAUDE.md) is loaded by no test. All 600 checks pass today; add the test. |
| T3 | Medium | DEVIATIONS.md is inaccurate. #4 has the wrong cause and scope (§1.1). #1's magnitude is wrong: there are 0-vs-subnormal cases (§1.1) and 1e-10 errors at large n (§1.2), and `abs=1e-305` in `test_ppois_upper` accepts all of them. #2 says 7/109 but the committed fixture gives 8/109, and it is not pinned. Rows are ordered 1, 3, 4, 2. |
| T4 | Medium | Chimera fixtures can't tell methods or options apart: consensus, pooled and per-sample keep the same 9 columns; `allow_one_off` and `min_sample_fraction` give the defaults' calls; `is_bimera_denovo` on the DADA object has 0 positives. Build tables where they disagree. |
| T5 | Medium | No Hypothesis property tests, although spec §8.4 requires them and `hypothesis` is a dependency. No n_jobs-determinism tests for `seq_complexity`, `remove_primers`, `assign_species`, `is_bimera_denovo`, `is_shift_denovo` or pooled/per-sample removal (all identical when checked here). None of the "behaviour changes" in DEVIATIONS are tested. |
| T6 | Medium | `uv run pytest -m parity` (the documented command) passes 141 tests but **exits 1**. The 90 % coverage gate applies to every marker subset (71.8 %), and the same breaks `-m slow`. |
| T7 | Low | `test_get_errors_matrix` feeds R's output back through `get_errors` (tautological). `test_nwhamming` re-implements nwhamming instead of calling the public API. The e2e test compares Kingdom/Phylum only. Unused fixture keys: `get_errors_detailed`, `get_uniques_seqtab_row`, … Warnings are swallowed everywhere, so R's warnings are never compared. No test for "fit returns self" or "inputs not mutated". |
| T8 | Low | Stale docs: CLAUDE.md says fixture regeneration takes "~20 s" (it is ~3.5 min) and that `error_plot_data` is tested against R (it is not). `generate_reference.R:10` cites a non-existent `tests/reference/fixtures.py`. PARITY.md says "every stage as fast as R or faster" (false, §4.2). |

---

## 4. Packaging, CI and benchmark methodology

### 4.1 CI cannot produce a green run (not yet run anywhere)

| Sev. | Finding |
|---|---|
| High | `ci.yml:41` `uv run pytest -m slow` exits 1 on the coverage gate (43.9 %), so the test job fails and `wheels` never runs. Add `--no-cov`, or set `--cov-fail-under` only on the full run. |
| High | `macos-13` (`ci.yml:30, 48`) is a retired runner label, so the spec-required macOS x86_64 leg is lost. |
| High | The `sdist` job (`ci.yml:69`) runs `uv build --sdist` without `uv sync`. Under `no-build-isolation-package` it fails with "No module named scikit_build_core" (reproduced). |
| Medium | cibuildwheel only tests `import pydada2`. The shipped wheels are never parity-tested, and linux-aarch64 has no test leg at all. |
| Medium | `pip install pydada2` (README) installs a *different* project. |
| Medium | Every source build downloads Boost from GitHub, even with `DADA2_FETCH_BOOST=OFF` when headers are missing, so offline and distro builds fail. |
| Medium | Licence: metadata `LGPL-2.0-or-later` vs LGPL-3 text in LICENSE; no GPL-3 text; no nanobind/robin_map notices in the wheel (spec Q9 still open). |
| Low | `scikit-misc>=0.5` admits the yanked 0.5.0. `xopen` is imported but not declared. `DADA2_BUILD_CPP_TESTS=ON` fails (no `cpp/tests`). `.benchmarks/` and `.hypothesis/` are not git-ignored. Actions are not SHA-pinned. `dist/` predates HEAD. |

**Solid:** clean-clone `uv sync` builds; the sdist excludes the right files; real abi3 wheel;
fp-contract off and `/fp:precise`; LF enforced; README quickstart runs; no build artefacts
tracked.

### 4.2 Benchmark reports overstate some ratios

Both suites run the same functions with the same parameters, and the phase-11 dataset benchmark
is sound (interleaved runs, medians, outputs compared). The micro report misleads in places:

- **`assign_species` 53.7 s vs 17.9 ms:** R's time is ~129 `gc()` calls (98.98 % in gc per
  Rprof), on 325 queries that hit nothing in a 14-sequence reference. The report also predates the
  phase-11 species rewrite, which was 3x *slower* than R before it.
- **`remove_bimera_pooled` 410 ms vs 0.5 ms:** one `gc()` in R on a 2×9 table with no chimeras.
- **Timer resolution:** R's 1 ms timer hides that `make_sequence_table` is ~6x *slower* in Python
  (0.66 vs 0.11 ms in tight loops; the dataset benchmark shows 2x).
- **`seq_complexity`:** R's timing includes reading the FASTQ.
- **Unequal thread counts:** micro `filter_and_trim@4` is R ≤ 2 processes vs Python 4 threads.
  Dataset `add_species@8` is R 1 thread vs Python 8.
- **Peak RSS:** `wait4` gives the largest process in the tree, not the sum (run.py docstring).
- **Dataset report, taxonomy rows:** ~98 % of the stage, in both tools, is the single-threaded
  reference-model build. The 8-thread column therefore says nothing about classifier scaling.
- **Dataset report, merge_pairs:** "no thread scaling" only because `pipeline.py` never passes
  `n_jobs`.
- **Dataset report, `make_sequence_table` "0.5x":** 2 vs 4 ms.
- **Narrow 8-thread margins on dada/learn_errors:** ~20 % of pydada2's time there is thread
  creation per parallel region (§5.1). That cost is not intrinsic to the algorithm.

Recommendation: report these rows with a caveat or drop them, re-run the micro suite at HEAD, and
fix the PARITY.md sentence.

---

## 5. Performance: making pydada2 clearly faster than R

### 5.1 Where the time goes today (MiSeq SOP, pydada2 at HEAD)

| stage | 1 thread | 8 threads | what dominates |
|---|---|---|---|
| assign_taxonomy | 58.4 s | 48.6 s | **`TaxonomyModel` constructor: 47.9 s, single-threaded** (98.5 % of the stage at 8 threads). Classification is only 12.4 s / 2.5 s. R has the same serial dense build (~70 of its 76 s). |
| learn_errors F+R | 68.3 s | 21.0 s | 97.5 % in `dada_uniques`, of which `nwalign_vectorized2` is 73 % (the scalar diagonal fill at `nwalign_vectorized.cpp:160-161` alone is 22.7 %). Python glue < 2.5 %. |
| dada F+R | 11.8 s | 3.5 s | same kernel as learn_errors |
| filter_and_trim | 7.0 s | 1.6 s | phiX `match_ref` 44.5 % (serial, `substr` + `unordered_set<string>` per window), zlib-ng level 6 20 %, `pow` 5.8 % |
| add_species | 4.9 s | 2.2 s | reading 1.3 s, Python prep 0.7 s, matching 0.6 s |
| derep, merge, seqtab, chimeras | 1.7 s | 1.6 s | merge_pairs runs on 1 thread because `pipeline.py` never passes `n_jobs` |

**Why 8 threads give only ~3.5x on learn_errors:**
- `parallel_for` spawns and joins `std::thread`s on every call. learn_errors_F has **12,776
  parallel regions** at ~200 µs each, ≈ 2.5 s. R uses TBB's persistent pool.
- ~1.9 s per run is serial: raw/k-mer setup, shuffle, bud and output.
- `final_subs` is parallelised over clusters with grain 10, and cluster 0 holds most raws.
- Samples run one after another, in R and in the port alike.

**Build flags:** `-O3` at the x86-64 SSE2 baseline. The SSE2 aligner and `kmer_dist_SSEi_8`
really are the default path at runtime. LTO covers only `module.cpp`; enabling it for
`dada2_core` made no measurable difference.

### 5.2 Ranked, output-preserving optimisations

All of these were prototyped together in a scratch build: `prototype.diff`, about 900 lines over
14 files, in the review scratchpad. **I checked independently that its MiSeq SOP `results.json`
is bit-identical to the HEAD pydada2 run for every output, at 1 and 8 threads, and therefore
still exact against R.** The outputs compared were track, error rates, derep, dada F/R, mergers,
sequence table, chimera calls, taxonomy and bootstraps.

| | HEAD pydada2 | prototype | R dada2 | prototype vs R |
|---|---|---|---|---|
| 8 threads | 78.4 s | **12.7 s** | 142 s | **11x faster** |
| 1 thread | 151.5 s | **40.6 s** | 250 s | **6.2x faster** |

The machine was shared during these measurements, so treat them as ±5–10 %. AVX2 was off, and
merge_pairs still ran on 1 thread.

| # | Change | Where | Measured effect | Why outputs stay identical |
|---|---|---|---|---|
| P1 | **Taxonomy:** build 8-mers with a rolling 2-bit code and a per-reference stamp array (no 65,536-entry zero+scan per reference); count in parallel over genus blocks (integer counts); store the log-probability table k-mer-major `[kmer][genus padded to 16]` and sum the query's rows for all genera at once, then run upstream's sequential arg-max and tie draws | `taxonomy.cpp:51-68, 95-130, 146-204, 242-301` | model build 47.9 → 4.45 s (1 thr) / 0.98 s (8 thr); classification 12.4 → 0.42 s; **stage 58.4 → 7.4 s (1 thr), 48.6 → 3.6 s (8 thr)** | counts < 2^24 are exact in float; same float expressions for prior and `logf`; per-genus sums keep their order; all log-probabilities ≤ 0, so the early exit never changes the arg-max or tie set. `tax`, `boot` and all 100 `boot_tax` per query are bit-identical (232 and 664 queries, try_rc on/off, 1/8 threads). |
| P2 | **Reuse alignments across self-consistency rounds:** cache, per sample and centre, the k-mer distance and the compact substitutions; re-apply the k-mer cutoff and recompute λ exactly as `compute_lambda` does | `cluster.cpp:23-71`, `dada_main.cpp` | learn_errors_F performs 5.70M alignments, **of which only 1.04M are distinct (82 % repeats)**. learn_F 47.5 → 16.6 s (1 thr). | alignments and k-mer distances do not depend on the error matrix; λ is recomputed from the same transitions in the same order |
| P3 | **Run samples concurrently** in `DADA.fit`/`learn_errors` (thread pool over samples, or a C++ batch call) | `denoise.py:476-560` | with P2: learn_F at 8 threads 6.8 → 2.5–2.9 s; without P2: 13.0 → 8.0 s, dada_F 2.2 → 1.1 s | samples are independent. **Requires the §2.1 fix first** (thread-local handlers). |
| P4 | **Persistent thread pool** in `parallel_for`, with a fallback to spawning when the pool is busy (`filter_and_trim` calls kernels from several Python threads) | `parallel.h` | −2.5 s per learn_errors_F at 8 threads | same static chunking contract |
| P5 | **Aligner kernel:** use a reversed copy of `s1`, so the diagonal fill reads contiguous memory and auto-vectorises; skip `kord_dist` when the k-mer screen has already rejected the pair (74 % of comparisons) | `nwalign_vectorized.cpp:160-162`, `nwalign_endsfree.cpp:40-46` | uncached learn_F −22 % (43.4 → 33.9 s), dada_F 7.6 → 5.9 s | integer arithmetic; the skipped value is unused. (The prototype's SEQLEN stack buffer must become dynamic: `nwvec` accepts longer strings.) |
| P6 | **filter_and_trim:** phiX words as 2-bit codes in a 2 MB bitmap with binary-search verification (fallback to the old path for non-ACGT references), run in parallel; a 256-entry `pow` lookup table built from the same expression | `filter.cpp:13-41`, `ported/fastq_filter.cpp` | 7.2 → 3.7 s (1 thr), 1.76 → 1.25 s (8 thr) | output records identical; `match_ref` identical on 960k randomised cases |
| P7 | **merge_pairs:** pass `n_jobs` in the pipeline (0.69 → 0.17 s); batch all samples into one kernel call and de-duplicate (F, R) pairs (3,663 per-sample pairs are only 1,446 distinct) | `merge.py`, `pipeline.py` | ~4x and ~2.5x | same alignments |
| P8 | **AVX2 by runtime dispatch** (`target_clones` or per-file `-mavx2` behind a CPU check; never a global flag, since wheels must run on CPUs without AVX2) | build | learn_F −12 %, dada_F 8.0 → 6.85 s | integer code; keep `-ffp-contract=off` on any float code |

### 5.3 Smaller items, same outputs

**Pipeline stages**
- **`BimeraRemover(method="per-sample")`** checks every sequence in every sample, including
  zero-count cells (`chimeras.py:397-405`). On 30 × 3000: 83.5 s on full rows vs 1.17 s on
  non-zero cells.
  - The transformed table is unchanged. Only `is_bimera_` entries for absent sequences would
    differ, so either keep computing those lazily or document the change.
  - Also hoist `[str(c) for c in columns]` out of the row loop.
- **`collapse_no_mismatch`** assigns one pandas column per ASV in a Python loop
  (`seqtab.py:228-233`): 3.4 s → 0.05 s on 200 × 20k with `np.add.at`. In C++
  (`collapse.cpp:44-61`), precompute the prefixes instead of allocating a string per pair, and
  stop spawning threads per block.
- **`is_shift_denovo`, per-sample chimera removal, `assign_species`** default to 1 thread:
  `is_shift_denovo` takes 14.8 s at 1 thread vs 2.0 s at 8 (R: 16 s).
- **`derep_fastq`** threaded over files: 0.51 → 0.21 s.

**C++ core**
- `final_subs` parallelised over raws instead of clusters.
- Thread-local DP and alignment buffers instead of ~15 mallocs per comparison (~3 %).
- Reuse the per-round raw/k-mer setup (0.75 s per learn_F).
- Parallelise `raw_new` and the k-mer passes.
- `bimera_denovo` copies parent strings for every query; pass indices.

**Memory**
- **`C_unifs`** (`taxonomy.cpp:231`) draws nseq·100·(max_arraylen/8) doubles, but ~92 % are
  never read. Truncate to the used prefix, which keeps the same stream: 100k ASVs × 450 bp drops
  from 4.4 GB to 0.36 GB.
- **Taxonomy fit:** the nanobind copy of 452k reference strings costs ~0.66 GB of peak RSS.
  Pass one buffer. Compute taxonomy depths on the 5,129 unique strings, not on all 452k.

**Output (FASTQ writing)**
- **gzip output** runs in the calling thread, with the two mates written one after another.
  Use a threaded writer and run the mates concurrently.

### 5.4 Caveats for implementing §5.2

- **P2 scope.** In the prototype the cache is process-wide and keyed by sample content. That is
  why `dada_F` after `learn_errors` drops to 0.18 s. Scope it per fit, or attach it to the
  `Dereplication`, and bound its memory: +300 MiB on MiSeq SOP, and it grows with the data.
  The cache key must include every alignment parameter.
- **Interrupts with P3.** Ctrl-C is then only detected on the main thread. A C++ batch call
  that schedules samples itself is the cleaner end state.
- **Verification.** Gate every change on the strict parity suite, a byte-compare of the MiSeq
  SOP `results.json` against HEAD, and the larger-dataset run (README).

---

## 6. Recommended order of work

1. **Parity:** land the ppois log-tail fix (§1.1). Make the self-consistent parity test ordered
   and strict, and add the near-underflow fixture. Rewrite DEVIATIONS #1/#2/#4.
2. **Safety:** make the runtime handlers thread-local (§2.1) and fix `birth_from` (§1.7).
   Validate err/quals in `dada_uniques` and fix the `B` leak (§2.2).
3. **Medium divergences:** fix §1.3–§1.5 and §1.6a/b/f/g/h, each with an R-generated fixture.
   Fill the fixture gaps T1/T2/T4.
4. **CI:** fix it (§4.1) so that macOS/Windows parity is actually established. Parity is required
   on all three OSes and has only been verified on Linux.
5. **Performance:** apply the §5 items in this order:
   - P1 taxonomy (largest and self-contained);
   - P4 thread pool and P5 aligner (local);
   - P2 alignment cache and P3 sample concurrency (P3 after the §2.1 fix);
   - P6–P8.

   Verify each against the full parity suite, the MiSeq SOP outputs (bit-identical) and the
   larger Stability dataset (see README). Then refresh both benchmark reports, fixing the
   methodology issues in §4.2.
