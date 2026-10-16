# CLAUDE.md — working notes for pyDADA2

Conventions, commands and decisions/gotchas for working on this repo. Background on the port, the
parity status and every known difference from R: `docs/porting-notes.md`.

## Project

- Python port of DADA2 (R, `benjjneb/dada2@b61fe2e`, v1.36.0), published as **pydada2**
  (import `pydada2`) at https://github.com/raeslab/pyDADA2 (fresh git history).
- C++ core reused from upstream (de-R'd) + hot R logic ported to C++ (kernels K1–K9, see
  `docs/porting-notes.md`), bound with nanobind; Pythonic sklearn-style API on top.
- License: `LGPL-3.0-or-later` (continues DADA2's LGPL code; `LICENSE` = LGPL-3, `COPYING` =
  GPL-3, `THIRD_PARTY_NOTICES.md` = nanobind/robin-map/Boost). Upstream `DESCRIPTION` says
  `LGPL-2`, upstream `LICENSE` holds the LGPL-3 text.
- **Credit:** users are asked to cite the DADA2 papers, never the port (README "Citation";
  CITATION.cff `preferred-citation` = DADA2 paper, upstream package + PacBio/Wang et al. papers
  under `references`). `authors` in pyproject/CITATION.cff is the author of the port only. Keep
  it that way in any new docs or metadata.

## Hard rules

- **Git:** never commit, push, open PRs/issues unless asked. The only remote is `pydada2`
  (raeslab/pyDADA2). **Never push to, open PRs/issues on, or add a remote for benjjneb/dada2:**
  this project has diverged and is independent (references to it are for credit only). End work
  with a summary + proposed commit message.
- Never `#include <Rcpp.h>` anywhere; never Python/nanobind headers in `cpp/core/`.
- **No code copied from R itself (GPL-2+).** Use Boost.Math / clean-room code for numerics.
- Never loosen a tolerance or skip a parity test to get green; document real deviations in
  `docs/porting-notes.md` §3 (cause, affected outputs, magnitude, pinning test).
- R behaviour (code, not docs) is the spec. The upstream source is at the commit above.
- Parity must hold on Linux, macOS and Windows.

## Commands

```bash
uv sync --all-groups                    # env + editable build (C++ rebuilds on import)
uv run pytest -n auto --cov-fail-under=90   # full run with the coverage gate (CI); marker subsets
                                        # (-m parity, -m "not slow", -m slow --no-cov) don't gate
uv run ruff format . && uv run ruff check . && uv run mypy src     # gates
./scripts/gen_stubs.sh                  # regenerate src/pydada2/_core.pyi after binding changes
uv sync --all-groups --reinstall-package pydada2   # after deleting build/ (editable hook needs it)
cmake -S . -B build/cpp-tests -DDADA2_BUILD_CPP_TESTS=ON -Dnanobind_DIR=... && ctest --test-dir build/cpp-tests
uv build --no-build-isolation --python "$UV_PROJECT_ENVIRONMENT/bin/python"   # sdist + wheel
```

R (only for fixtures and R benchmarks; dev container built with `INSTALL_R=true`):

```bash
bash scripts/install_r_dada2.sh                   # dada2 @ b61fe2e into $R_LIBS_USER
Rscript tests/reference/generate_reference.R      # regenerate fixtures (~3.5 min; then `git checkout --`
                                                  # taxonomy.json.gz + MANIFEST.json, nondeterministic)
Rscript benchmarks/micro/run_r.R 7 /tmp/r.json && uv run python benchmarks/micro/run_py.py 7 /tmp/py.json \
  && uv run python benchmarks/micro/report.py /tmp/r.json /tmp/py.json   # micro-benchmarks
benchmarks/dataset/download.sh && uv run python benchmarks/dataset/run.py 3 1 8 \
  && uv run python benchmarks/dataset/report.py           # MiSeq SOP benchmark (~1 h)
benchmarks/dataset/download.sh --stability && uv run python benchmarks/dataset/run.py \
  --fastq-dir Stability --results benchmarks/dataset/results/stability 2 8   # 362 samples (~30 min)
# old vs new pydada2: run.py --tool py --python <other venv>/bin/python --results <dir_before>, then
# report.py <results> <REPORT.md> <dataset> <dir_before>; micro: report.py r.json py.json out.md py_before.json
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

- Format conventions: `tests/reference/README.md` (every R value is a JSON array, `NA`/NaN/Inf
  → `null`, 1-based indices stored as-is, FASTQ outputs as canonical md5, ...).
- Fixtures are deterministic across runs **except** `taxonomy.json.gz` (and the MANIFEST
  timestamp): upstream breaks taxonomy ties with `std::random_device`. Measured R-vs-R on the
  fixture data: bootstrap |Δ| ≤ 2 (median 0); `minBoot=50` assignments 99/100 identical on
  full-length queries. Python-vs-R taxonomy tolerances must account for this.
- Discriminating fixtures (each section of the generator self-seeded):
  - numerics: `ppois_extended`, `complexity_thresholds`, `match_pattern_degenerate`;
  - python layer: `loess_degenerate`, `review_python`;
  - coverage: `filter_coverage` (every filter bites; the generator asserts it), `derep_coverage`,
    `chimera_coverage` (methods/options disagree);
  - inputs: `cov_F/R`, `derep_iupac`, `chimeric_sam1F`, `binned_singletons`, `priors.fa`,
    `species_tokens.fa`, `with_n`.
- `e2e.json` uses learned errors and holds no taxonomy, because R's taxonomy is not reproducible
  even with a fixed seed.
- `rparity` takes `EXTDATA` from the installed package, so the suite also runs against wheels.

## Design decisions

C++ core:
- `cpp/core` = upstream `src/*.cpp` de-R'd; public C++ API in `api.h` (row-major buffers,
  0-based, -1/NaN for NA). Each file's header lists its changes vs upstream.
- `b_compare`/`b_compare_parallel` merged into one `parallel_for` implementation; results are
  identical for any `n_threads` (tested 1/2/4).
- Non-x86: the `*_SSEi` k-mer routines are exact scalar emulations of the SSE2 code (upstream
  compiled dummies and forced SSE=0, so R-on-ARM can differ from R-on-x86; we match x86).
- `ppois_upper`: Loader saddle-point dpois + series below the mode (clean-room, stats.cpp header),
  direct λ^m e^-λ/m! for λ < 1, log space below DBL_MIN; Boost `gamma_p` (`promote_double<false>`)
  at/above the mode. Don't revert to plain `gamma_p`: it underflows to 0 where R returns a
  subnormal (changes DADA cluster order vs R) and loses 1e-10 at large n.
- `seq_complexity` sums in long double like R's `sum()`; double-double where `LDBL_MANT_DIG < 64`
  (MSVC, Apple arm64: `long double` is plain double there).
- Taxonomy: `TaxonomyModel` (fit once) + `assign_taxonomy`; uniforms from `mt19937_64(seed)` in
  upstream's count/order (only the used prefix is drawn); tie-breaks from a per-sequence
  generator ⇒ independent of n_threads. Table is k-mer-major `[4^8][ngenus padded to 16]`, built
  sparsely in parallel; the classifier sums query rows for all genera at once (exact: per-genus
  order kept). Fitted models pickle (~1.3 GB for SILVA).
- Errors: `std::invalid_argument` → ValueError, `dada2::Error` → RuntimeError,
  `dada2::Interrupted` → KeyboardInterrupt (interrupt checks only on the main thread).
- Threading: `parallel_for` runs on a process-wide persistent pool (`parallel.cpp`) serving
  concurrent and nested regions (the caller works on its own region; capped at chunks/hardware;
  thread-creation failure non-fatal; reset after fork). Log/interrupt handlers are thread-local
  (`ScopedHandlers` restores); `ScopedCancel` stops concurrent work.
- No `-march`/`target_clones`: AVX2 was measured with no gain (rows ~9 cells; taxonomy is
  memory-bound).

IO, filtering, dereplication, primers:
- IO: dnaio streaming in chunks of `n` records (`io.read_fastq_chunks`), sequences upper-cased
  (ShortRead does), qualities kept encoded; quality offset auto-detected per chunk exactly like
  ShortRead `.qualityTypeAuto` (first 10000 strings: min > ':' and max > 'J' → 64).
- K1 `filter_mate` applies R's per-read steps to one mate; paired keep = pass_F & pass_R (each R
  step only looks at its own read, so this equals R's sequential whole-chunk filtering).
  Orientation and `match_ids` (incl. remainder carry-over and CASAVA auto-detection) stay in
  Python because they reorder chunks. Paired mates are filtered concurrently; phiX words are
  matched as 2-bit codes with a bitmap, in parallel; FASTQ writing runs in a background thread.
- K2 `Dereplicator`: Biostrings order of DNA letters (IUPAC bit codes, `N` after `T`) decides tie
  order; quality sums are exact so the averages match R bit-for-bit.
- K7 primer matching: mismatch mode = Biostrings naive-inexact incl. out-of-limits matches;
  indel mode re-implements Biostrings' actual `match_pattern_indels` algorithm (its documented
  "best local match" definition is *not* what it computes for large max.mismatch).
- K8 `species_hits`: rolling hash + bitmap pre-filter + memcmp verification.
- K9 quality profile is plain NumPy (`io.quality_profile`): bincount is as fast as C++.
- Output FASTQ: `@id\nseq\n+\nqual\n`, gzip via xopen; parity compares records (canonical md5).

Python API:
- `DADA(err, error_model, self_consistent, pool, priors, params, n_jobs, verbose)`: options go in
  `params` (`DadaParams` or dict) instead of `**kwargs` (sklearn's `get_params`/`clone` need
  explicit `__init__` arguments).
- `DADA.fit` uses `_core.dada_uniques_batch`: samples concurrently (largest first), n_threads per
  sample on the shared pool; with verbose >= 2 one at a time. Multi-round fits (self-consistency,
  pseudo-pooling) own an `AlignmentCache` (1 GiB budget, `denoise._ALIGNMENT_CACHE_BYTES`): k-mer
  distances + packed substitutions per (centre, raw); λ recomputed with the same multiplication
  order (prefix products). No global state: `dada()` after `learn_errors` reuses nothing.
- `learn_errors`/`ErrorLearner` live in `error_models.py` and wrap `DADA(self_consistent=True)`;
  `ErrorLearner` raises when `params` conflicts with its own `max_consist`/`omega_c`.
- `get_errors` accepts matrices, fitted error models/learners (`err_`), fitted DADA (`err_out_`),
  DadaResults and lists of them.
- Pooled DADA results keep the pooled `pval` and `trans` per sample (as R).
- loess: skmisc on the normal path; when skmisc raises on R *warnings* (few quality scores),
  `pydada2._loess` (NumPy port of the netlib 1-D interpolate algorithm) gives R's fit (≤ 1.2e-12).
- `BimeraRemover` per-method defaults follow R (consensus 1.5/2, pooled & per-sample 2/8).
- `get_uniques` sorts when merging duplicates (R `tapply`); species references are read with
  `read_fasta(raw_ids=True)` (dnaio strips trailing whitespace).
- Taxonomy parity criterion: exact top hits; bootstrap differences within 3σ of binomial noise
  (pooled p); at min_boot=50 disagreements only where R's bootstrap is within 50 ± 3σ.
- `TaxonomyClassifier.predict(seqs, return_bootstrap=True)` returns both in one pass;
  `predict_bootstrap` re-runs with the same seed ⇒ same values.
- Plot functions return matplotlib Figures; `quality_profile_stats` is tested exactly against
  R's plotQualityProfile summary (`error_plot_data` is not compared with R).

CLI (`cli.py`, spec in `CLI_SPEC.md`; formats in `formats.py`):
- `formats.py` holds `read/write_dada_json` + `read/write_seqtab` (not `io.py`: `io` is imported
  by `derep`/`denoise`, circular). `pydada2.dada/1` omits per-position quals (merge doesn't need
  them); restored `DadaResult`s have `quality` with 0 columns and `err_in/err_out=None`.
- `GreedyOption` wraps click's private parser `_Option.process` (click 8.5); re-check on click
  upgrades. Per-mate options accept a 2nd value only if it parses as a number (`-inf`, `-4`).
- Gzipped outputs: `_open_write` writes no file name and mtime 0 (byte-identical reruns).
- Float TSV columns are written with `repr` (round-trip exact); read them back with
  `float_precision="round_trip"` (pandas' default parser is not exact).

Benchmarks:
- `pipeline.R`/`pipeline.py` time each stage in-process and dump all compared outputs to
  `results/<tool>_t<thr>_r<rep>/results.json`; `run.py` adds wall time + peak RSS (`wait4`
  rusage, includes R's forked `mclapply` workers). Chimeras are timed via
  `isBimeraDenovoTable` (= `removeBimeraDenovo(method="consensus")`).
- Every performance change must be verified bit-identical on MiSeq SOP (all outputs incl.
  taxonomy at the same seed, 1 and 8 threads), ideally also on Stability.
- Large parallel changes: use git worktrees with their **own venvs**
  (`UV_PROJECT_ENVIRONMENT=<wt>/.venv uv sync --all-groups`, then `<wt>/.venv/bin/python -m ...`);
  never `uv run` there (the global UV_PROJECT_ENVIRONMENT points at the main tree's venv).

## Gotchas

- **Benchmark from comparable storage.** `/workspaces` (the dev-container bind mount) costs
  ~1 ms per `stat`/`open`. Example data read from a checkout there makes tiny micro stages look
  10x slower than the same code under `/tmp`. Compare pydada2 versions from equivalent locations.
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
  nanobind silently builds a version-specific `.so` inside a `cp312-abi3` wheel.
- `uv pip install <sdist>` from inside the repo applies the repo's `[tool.uv]` settings
  (no build isolation) — test installs from another directory.
- pandas 3 infers `str` dtype: `None` in string columns becomes NaN; build object arrays directly
  when None must be preserved.
- ASan + nanobind: preload libstdc++ as well as libasan, otherwise `__cxa_throw` interception
  fails on the first C++ exception.
- The example data are low-diversity (8–9 merged ASVs per sample); synthetic columns are added
  for collapse/chimera/shift fixtures.
- scikit-misc ≥ 0.5.3 ships no macOS x86_64 wheels (CI's cibuildwheel tests use
  `PIP_PREFER_BINARY=1`); the lock file pins 0.5.2.
- Repo scripts are stored without the executable bit unless `git add --chmod=+x` is used
  (`core.fileMode=false` in the container); call them with `bash` from automation.
