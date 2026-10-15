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
Rscript tests/reference/generate_reference.R      # regenerate fixtures (~20 s)
```

```bash
uv sync --all-groups                    # env + editable build (C++ rebuilds on import)
uv run pytest                           # tests + coverage;  -m "not slow" / -m parity / -n auto
uv run ruff format . && uv run ruff check . && uv run mypy src     # gates (spec §10.1)
./scripts/gen_stubs.sh                  # regenerate src/pydada2/_core.pyi after binding changes
uv sync --all-groups --reinstall-package pydada2   # after deleting build/ (editable hook needs it)
Rscript benchmarks/micro/run_r.R 7 /tmp/r.json && uv run python benchmarks/micro/run_py.py 7 /tmp/py.json \
  && uv run python benchmarks/micro/report.py /tmp/r.json /tmp/py.json   # micro-benchmarks
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
  identical on every platform); `DADA2_FETCH_BOOST=OFF` uses system headers (offline builds).
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
- Plot functions return matplotlib Figures; data helpers (`quality_profile_stats`,
  `error_plot_data`) are tested exactly against R's plotQualityProfile summary.

## Gotchas found so far

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

