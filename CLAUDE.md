# CLAUDE.md — working notes for pyDADA2

`spec.md` is the work order; this file holds conventions, commands and decisions/gotchas
discovered while executing it. Keep both up to date at the end of every phase.

## Project

- Python port of DADA2 (R, `benjjneb/dada2@b61fe2e`, v1.36.0), published as **pydada2**
  (import `pydada2`) at https://github.com/raeslab/pyDADA2 (fresh git history).
- C++ core reused from upstream (de-R'd) + hot R logic ported to C++ (spec §6), bound with
  nanobind; Pythonic sklearn-style API on top.
- License: LGPL (continuing DADA2 code). Exact SPDX pending (spec Q9): upstream `DESCRIPTION`
  says `LGPL-2`, upstream `LICENSE` contains LGPL-3 text. Placeholder `LGPL-2.0-or-later`.

## Hard rules

- **Git:** never commit, push, open PRs/issues, or touch the `origin` remote (benjjneb/dada2).
  End each phase with a summary + proposed commit message; the human commits.
- **Gated phases:** never start phase 11 (large dataset) or phase 12 (remove R reference)
  without an explicit instruction from the human.
- `reference/r-dada2/` is read-only (the R reference implementation).
- Never `#include <Rcpp.h>` anywhere; never Python/nanobind headers in `cpp/core/`.
- **No code copied from R itself (GPL-2+).** Note: upstream `src/pval.cpp` carries a
  commented-out copy of R's `nmath` ppois/pgamma — delete that block when porting.
  Use Boost.Math for `ppois_upper` (spec §5.4).
- Never loosen a tolerance or skip a parity test to get green; document real deviations in
  `DEVIATIONS.md` (spec §8.3).
- R behaviour (code, not docs) is the spec.
- Parity must hold on Linux, macOS and Windows.

## Commands

```bash
# R reference (installed by scripts/devcontainer-post-create.sh into $R_LIBS_USER)
R CMD INSTALL --no-docs --library="$R_LIBS_USER" reference/r-dada2
Rscript tests/reference/generate_reference.R      # regenerate fixtures (~3.5 min; then `git checkout --`
                                                  # taxonomy.json.gz + MANIFEST.json, nondeterministic)
```

```bash
uv sync --all-groups                    # env + editable build (C++ rebuilds on import)
uv run pytest -n auto --cov-fail-under=90   # full run with the coverage gate (CI); marker subsets
                                        # (-m parity, -m "not slow", -m slow --no-cov) don't gate
uv run ruff format . && uv run ruff check . && uv run mypy src     # gates (spec §10.1)
./scripts/gen_stubs.sh                  # regenerate src/pydada2/_core.pyi after binding changes
uv sync --all-groups --reinstall-package pydada2   # after deleting build/ (editable hook needs it)
Rscript benchmarks/micro/run_r.R 7 /tmp/r.json && uv run python benchmarks/micro/run_py.py 7 /tmp/py.json \
  && uv run python benchmarks/micro/report.py /tmp/r.json /tmp/py.json   # micro-benchmarks
benchmarks/dataset/download.sh && uv run python benchmarks/dataset/run.py 3 1 8 \
  && uv run python benchmarks/dataset/report.py           # phase 11 MiSeq SOP benchmark (~1 h)
benchmarks/dataset/download.sh --stability && uv run python benchmarks/dataset/run.py \
  --fastq-dir Stability --results benchmarks/dataset/results/stability 2 8   # 362 samples (~30 min)
# old vs new pydada2: run.py --tool py --python <other venv>/bin/python --results <dir_before>, then
# report.py <results> <REPORT.md> <dataset> <dir_before>; micro: report.py r.json py.json out.md py_before.json
cmake -S . -B build/cpp-tests -DDADA2_BUILD_CPP_TESTS=ON -Dnanobind_DIR=... && ctest --test-dir build/cpp-tests
uv build --no-build-isolation --python "$UV_PROJECT_ENVIRONMENT/bin/python"   # sdist + wheel
```

Build notes:
- `[tool.uv] no-build-isolation-package = ["pydada2"]` + scikit-build-core/nanobind in the dev
  group: required for `editable.rebuild`. Consequence: plain `uv build` fails (isolated build
  env lacks scikit-build-core); use the command above. cibuildwheel/pip use normal isolation.
- `[tool.uv] exclude-newer` pins resolution to releases before 2026-10-01 (supply-chain hygiene).
- Wheels: `wheel.py-api = "cp312"` → one abi3 wheel for CPython ≥ 3.12, per-version for 3.11.
- clangd reads `build/cp312-abi3-linux_x86_64/compile_commands.json` (see `.clangd`).
- Boost.Math: pinned standalone 1.86 fetched at configure time by default (`DADA2_FETCH_BOOST=ON`,
  identical on every platform); `DADA2_FETCH_BOOST=OFF` never downloads: system headers or
  `-DDADA2_BOOST_MATH_INCLUDE=...` (FATAL_ERROR if missing).
- Memory check: `./scripts/asan_check.sh` (ASan build of `_core`, throwaway venv, preloads
  libasan + libstdc++, fails on ASan errors or leaks with frames in cpp/core|cpp/bindings).
- Fast C++ iteration without the editable hook: `cmake --build build/dev` (configure once as in
  `scripts/asan_check.sh` without the sanitizer flags).
- Test helpers for fixtures live in `tests/rparity.py` (`pythonpath = ["tests"]`).

## Fixtures (`tests/reference/data/`)

- Written by `tests/reference/generate_reference.R` with `jsonlite`, `digits=I(17)` (exact
  double round-trip), `auto_unbox=FALSE` (**every R value is a JSON array, scalars too**),
  `NA → null` (**also NaN/±Inf for doubles**: compare with `rparity.assert_close_r`), matrices row-major as `{rownames, colnames, values}`, data.frames as
  `{column: [...]}`. Files > 100 kB are `.json.gz`.
- R indices are stored 1-based as-is (`map`, `pairs`, `forward`/`reverse` in mergers);
  Python tests convert.
- Dereplication quals stored as integer **sums** (`qual_sums`); expected mean = sum / abundance.
- FASTQ outputs stored as `ids`, `widths` and `md5` of the canonical text
  `"id\nseq\nqual\n"` per record (file format/compression independent).
- `err/*.npy` are the bundled `errBalancedF/R`, `tperr1` (16 × 41, float64, C order).
- `inputs/` holds small synthetic FASTQ inputs (orientation mixes, matchIDs membership).
- Fixtures are deterministic across runs **except** `taxonomy.json.gz` (and the MANIFEST
  timestamp): upstream breaks taxonomy ties with `std::random_device`. Measured R-vs-R on the
  fixture data: bootstrap |Δ| ≤ 2 (median 0); `minBoot=50` assignments 99/100 identical on
  full-length queries. Python-vs-R taxonomy tolerances must account for this.
- Review-fix fixtures (sections at the end of the generator, each self-seeded):
  - numerics: `ppois_extended`, `complexity_thresholds`, `match_pattern_degenerate`;
  - python layer: `loess_degenerate`, `review_python`;
  - coverage: `filter_coverage` (every filter bites; the generator asserts it), `derep_coverage`,
    `chimera_coverage` (methods/options disagree);
  - inputs: `cov_F/R`, `derep_iupac`, `chimeric_sam1F`, `binned_singletons`, `priors.fa`,
    `species_tokens.fa`, `with_n`.
- `run_filter` records R's warnings, and the parity tests require exactly those warnings
  (`rparity.r_warnings`).
- `e2e.json` uses learned errors and holds no taxonomy, because R's taxonomy is not reproducible
  even with a fixed seed.
- `rparity` takes `EXTDATA` from the installed package, so the suite also runs against wheels.

## C++ core decisions (phase 2)

- `cpp/core` = upstream `src/*.cpp` de-R'd; public C++ API in `api.h` (row-major buffers,
  0-based, -1/NaN for NA). Each file's header lists its changes vs upstream.
- `b_compare`/`b_compare_parallel` merged into one `parallel_for` implementation; results are
  identical for any `n_threads` (tested 1/2/4).
- Non-x86: the `*_SSEi` k-mer routines are exact scalar emulations of the SSE2 code (upstream
  compiled dummies and forced SSE=0, so R-on-ARM can differ from R-on-x86; we match x86).
- `ppois_upper` = Boost `gamma_p(n+1, λ)` with `promote_double<false>` (no long double ⇒ same
  result on all platforms). Matches R to rel 1e-12 except subnormal outputs (~1e-318).
- Taxonomy: `TaxonomyModel` (fit once) + `assign_taxonomy`; uniforms from `mt19937_64(seed)` in
  upstream's count/order; tie-breaks from a per-sequence generator ⇒ independent of n_threads.
- Errors: `std::invalid_argument` → ValueError, `dada2::Error` → RuntimeError,
  `dada2::Interrupted` → KeyboardInterrupt (interrupt checks only on the main thread).

## Phase 3 decisions

- IO: dnaio streaming in chunks of `n` records (`io.read_fastq_chunks`), sequences upper-cased
  (ShortRead does), qualities kept encoded; quality offset auto-detected per chunk exactly like
  ShortRead `.qualityTypeAuto` (first 10000 strings: min > ':' and max > 'J' → 64).
- K1 `filter_mate` applies R's per-read steps to one mate; paired keep = pass_F & pass_R (each R
  step only looks at its own read, so this equals R's sequential whole-chunk filtering).
  Orientation and `match_ids` (incl. remainder carry-over and CASAVA auto-detection) stay in
  Python because they reorder chunks.
- K2 `Dereplicator`: Biostrings order of DNA letters (IUPAC bit codes, `N` after `T`) decides tie
  order; quality sums are exact so the averages match R bit-for-bit.
- K7 primer matching: mismatch mode = Biostrings naive-inexact incl. out-of-limits matches;
  indel mode re-implements Biostrings' actual `match_pattern_indels` algorithm (its documented
  "best local match" definition is *not* what it computes for large max.mismatch); validated on
  300 random R cases (`match_pattern.json`).
- K9 quality profile is plain NumPy (`io.quality_profile`), not C++: bincount is as fast.
- Output FASTQ: `@id\nseq\n+\nqual\n`, gzip via xopen; parity compares records (canonical md5).

## Phase 4-7 decisions

- `DADA(err, error_model, self_consistent, pool, priors, params, n_jobs, verbose)`: options go in
  `params` (`DadaParams` or dict) instead of `**kwargs` (sklearn's `get_params`/`clone` need
  explicit `__init__` arguments). Deviation from spec §7.2 wording, same capability.
- `learn_errors`/`ErrorLearner` live in `error_models.py` and wrap `DADA(self_consistent=True)`.
- `get_errors` accepts matrices, fitted error models/learners (`err_`), fitted DADA (`err_out_`),
  DadaResults and lists of them.
- Pooled DADA results keep the pooled `pval` and `trans` per sample (as R).
- `BimeraRemover` per-method defaults follow R (consensus 1.5/2, pooled & per-sample 2/8).

## Phase 8-9 decisions

- Taxonomy parity criterion refined (spec §8.3 updated): exact top hits; bootstrap differences
  within 3σ of binomial noise (pooled p); at min_boot=50 disagreements only where R's bootstrap
  is within 50 ± 3σ. The original "≥95 % identical" is not attainable on the small fixture sets.
- `TaxonomyClassifier.predict(seqs, return_bootstrap=True)` returns both in one pass;
  `predict_bootstrap` exists for the spec'd API (re-runs with the same seed ⇒ same values).
- Plot functions return matplotlib Figures; `quality_profile_stats` is tested exactly against
  R's plotQualityProfile summary (`error_plot_data` is not compared with R).

## Phase 11 decisions

- Dataset agreed with the human: MiSeq SOP (20 samples) + SILVA 138.2 (toGenus trainset +
  assignSpecies, Zenodo 14169026); 1 and 8 threads; 3 repeats; spec §9.2 acceptance criteria.
- `pipeline.R`/`pipeline.py` time each stage in-process and dump all compared outputs to
  `results/<tool>_t<thr>_r<rep>/results.json`; `run.py` adds wall time + peak RSS (`wait4`
  rusage, includes R's forked `mclapply` workers; no `/usr/bin/time` in the container).
  Chimeras are timed via `isBimeraDenovoTable` (= `removeBimeraDenovo(method="consensus")`).
- Outcome: all exact outputs identical; error rates rel 1.4e-14; taxonomy 229/232 identical at
  every level, the 3 differences have R bootstrap 45–48. R-vs-R taxonomy was identical here.
- `species_hits` (K8) was 3x slower than R on SILVA (`substr` + full-string hash per window):
  now a rolling hash + bitmap pre-filter + memcmp verification, 18x faster, identical hits.
- `scripts/asan_check.sh` venv now installs matplotlib (needed since the phase-9 plotting tests;
  delete `/tmp/pydada2-asan/venv` to recreate an old one).

## REVIEW.md fixes (2026-10-15) decisions

- Development of large changes in parallel used git worktrees with their **own venvs**
  (`UV_PROJECT_ENVIRONMENT=<wt>/.venv uv sync --all-groups`, then `<wt>/.venv/bin/python -m ...`);
  never `uv run` there (the global UV_PROJECT_ENVIRONMENT points at the main tree's venv).
- `ppois_upper`: Loader saddle-point dpois + series below the mode (clean-room, stats.cpp header),
  direct λ^m e^-λ/m! for λ < 1, log space below DBL_MIN; Boost `gamma_p` at/above the mode. Don't
  revert to plain `gamma_p`: it underflows to 0 where R returns a subnormal (changed DADA cluster
  order vs R, REVIEW.md §1.1) and loses 1e-10 at large n.
- `seq_complexity` sums in long double like R's `sum()` (double-double on MSVC).
- loess: skmisc on the normal path; when skmisc raises on R *warnings* (few quality scores),
  `pydada2._loess` (NumPy port of the netlib 1-D interpolate algorithm) gives R's fit (≤ 1.2e-12).
- `get_uniques` sorts when merging duplicates (R `tapply`); species references are read with
  `read_fasta(raw_ids=True)` (dnaio strips trailing whitespace); `ErrorLearner` raises when
  `params` conflicts with its own `max_consist`/`omega_c`.
- Threading: `parallel_for` runs on a process-wide persistent pool (`parallel.cpp`) serving
  concurrent and nested regions (the caller works on its own region; capped at chunks/hardware;
  thread-creation failure non-fatal; reset after fork). Log/interrupt handlers are thread-local
  (`ScopedHandlers` restores); `ScopedCancel` stops concurrent work.
- `DADA.fit` uses `_core.dada_uniques_batch`: samples concurrently (largest first), n_threads per
  sample on the shared pool; with verbose >= 2 one at a time. Multi-round fits (self-consistency,
  pseudo-pooling) own an `AlignmentCache` (1 GiB budget, `denoise._ALIGNMENT_CACHE_BYTES`): k-mer
  distances + packed substitutions per (centre, raw); λ recomputed with the same multiplication
  order (prefix products). No global state: `dada()` after `learn_errors` reuses nothing.
- Taxonomy table is k-mer-major `[4^8][ngenus padded to 16]`, built sparsely in parallel; the
  classifier sums query rows for all genera at once (exact: per-genus order kept, log-probs < 0).
  Fitted models pickle (the table, ~1.3 GB for SILVA). Only the used prefix of the uniform stream
  is drawn.
- No `-march`/`target_clones`: AVX2 was measured with no gain (rows ~9 cells; taxonomy is
  memory-bound).
- FASTQ writing runs in a background thread (≤ 1 chunk pending); paired mates filtered
  concurrently; phiX words matched as 2-bit codes with a bitmap, in parallel.
- Every performance change was verified bit-identical on MiSeq SOP (all outputs incl. taxonomy at
  the same seed, 1 and 8 threads) and on Stability.

## Gotchas found so far

- **Benchmark from comparable storage.** `/workspaces` (the dev-container bind mount) costs
  ~1 ms per `stat`/`open`. Example data read from a checkout there makes tiny micro stages look
  10x slower than the same code under `/tmp`. Compare pydada2 versions from equivalent locations
  (the final runs used `/tmp` worktrees with their own venvs).
- **`s if s.isupper() else s.upper()` is 3–4x slower than `s.upper()`:** `upper()` has an ASCII
  fast path and `isupper()` does not. Measure "copy-avoiding" micro-optimisations before keeping
  them.
- **xopen starts a decompression thread per open by default** (~5 ms of start-up/shutdown).
  `read_fasta` uses `threads=0` for files under 4 MB.

- `generate_reference.R` sections run at top level: a loop variable named like a helper
  (`mm <- matchPattern(...)`) silently replaced the `mm()` mergePairs helper used by a later
  section. Use distinct names (or `local()`) in new sections.
- The sdist exclude pattern `"tests/"` is gitignore-style and also matched `cpp/tests/`; use
  `"/tests/"`.

- R partial argument matching: a helper with a `name` formal swallowed `n = 300`
  (`filterAndTrim`'s chunk size). Use dot-prefixed formals in R helpers.
- `filterAndTrim(..., n=)`: chunking changes output **order** when `orient.fwd` is used
  (within each chunk: forward-matching reads first, then the reverse-complemented ones) and
  `derepFastq(n=)` tie ordering (new uniques per chunk appended in lexical order, then a
  stable sort by abundance). Both are pinned by fixtures (`F7_orient_chunked`,
  `filtF_sam1_n200`, `raw_samPB_n100`).
- Single-end filtering drops zero-length reads inside `trimTails`; paired filtering uses
  `trimTails(ranges=TRUE)` and drops them afterwards — same effect, different code path.
- `fastqFilter` applies `minQ` only when `minQ > truncQ` and keeps reads with
  `min(qual) > minQ` (strict), contrary to the docs.
- `assignSpecies` rejects non-ACGT queries; fixtures use the ACGT-only subset (`queries_species`).
- Inside `namespace dada2`, unqualified `nwalign(...)`/`kmer_dist(...)` resolve to the API
  functions; call the upstream C functions as `::nwalign`, `::kmer_dist`.
- abi3 wheels: `find_package(Python ... ${SKBUILD_SABI_COMPONENT})` is required, otherwise
  nanobind silently builds a version-specific `.so` inside a `cp312-abi3` wheel (found by
  installing the wheel on 3.12).
- `uv pip install <sdist>` from inside the repo applies the repo's `[tool.uv]` settings
  (no build isolation) — test installs from another directory.
- pandas 3 infers `str` dtype: `None` in string columns becomes NaN; build object arrays directly
  when None must be preserved.
- ASan + nanobind: preload libstdc++ as well as libasan, otherwise `__cxa_throw` interception
  fails on the first C++ exception.
- The example data are low-diversity (8–9 merged ASVs per sample); synthetic columns are added
  for collapse/chimera/shift fixtures.
- `.devcontainer/devcontainer.json` originally had no `postCreateCommand` and the script
  `cd`'d to `/workspace`; fixed to run from the repo root (`/workspaces/dada2` here).

## Phase log (proposed commit messages; the human commits)

Fresh-history setup (human, once): `git switch --orphan pydada2-main`, then
`git remote add pydada2 https://github.com/raeslab/pyDADA2.git`, review, commit, and
`git push pydada2 pydada2-main:main`.

- phase 0: move R package to reference/r-dada2, R fixture generator + fixtures, CLAUDE.md, DEVIATIONS.md, devcontainer post-create fix
- phase 1: pyproject/CMake/nanobind skeleton, estimator base + validation, tests, pre-commit, CI
- phase 2: de-R'd C++ core (cpp/core), nanobind bindings, low-level + dada_uniques parity, ASan check
- phase 3: IO, derep (K2), filter_and_trim (K1), remove_primers (K7), complexity (K6), phiX, quality profile; parity exact
- phase 4: error models (loess via skmisc, noqual, binned, PacBio), inflate_err, get_errors
- phase 5: DADA estimator (pool/pseudo/priors/self-consistent), DadaParams, DadaResult, ErrorLearner/learn_errors
- phase 6: merge_pairs (K3), make/merge_sequence_tables, collapse_no_mismatch (K5); parity exact
- phase 7: chimeras (K4 de novo, table consensus, per-sample), BimeraRemover, is_shift_denovo (K5); parity exact
- phase 8: TaxonomyClassifier/assign_taxonomy (fit-once model), SpeciesAssigner/assign_species (K8), add_species
- phase 9: plotting (+ exact quality-profile stats), align module, README, CITATION.cff, e2e test, abi3 wheel fix
- phase 10: PARITY.md sign-off, micro-benchmarks (benchmarks/micro: run_r.R, run_py.py, report.py, REPORT.md)
- phase 11: MiSeq SOP + SILVA 138.2 benchmark vs R (benchmarks/dataset), faster species matching, README benchmark section
- review: REVIEW.md + Stability (362-sample) benchmark (dff04da)
- review fixes: branches fix/numerics, fix/python, fix/tests, fix/perf merged (ppois underflow/accuracy, loess fallback, R-parity fixes, discriminating fixtures, CI, ~6-12x faster with identical output)

