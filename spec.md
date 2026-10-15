# spec.md — Porting DADA2 from R to a modern Python package (pyDADA2)

This spec is the work order for Claude Code running inside this repo's dev container.
Read it fully before starting, then work through the phases in §10 in order.

---

## 0. How to use this spec (for the human)

1. Open the repo in VS Code → **Dev Containers: Reopen in Container**. The first build takes
   a while (R + Bioconductor binaries, compilers, Python 3.13, Claude Code).
2. In the container terminal run `claude`, log in once (the login is kept in a named volume).
3. Prompt: *"Read spec.md and execute phase 0. Stop at the end of the phase and summarise."*
   Then continue phase by phase. Each phase has acceptance criteria; review the diff and
   commit/push yourself (see §1.1) before moving on.
4. Phases marked **GATED** (§10, phases 11 and 12) are never started by Claude on its own. They
   run only when the human explicitly says e.g. *"start phase 12 (cleanup)"*.

---

## 1. Goal

Produce a Python package that provides DADA2's amplicon denoising pipeline (filtering,
dereplication, error learning, sample inference, read merging, chimera removal, taxonomy
assignment) by wrapping **DADA2's existing C++ core** (`src/*.cpp`) via **nanobind**. It also
**moves the hot R-level logic into C++** (§6), and exposes a **Pythonic, scikit-learn-style API**.

- Home of the Python package: **https://github.com/raeslab/pyDADA2**
- Distribution name `pydada2`, import name `pydada2` (decided).
- This working copy is a fork of upstream `benjjneb/dada2` (`origin`). It is only the
  workspace and the source of the R reference.
- **pyDADA2 starts with a fresh git history** (decided). Upstream is credited by commit SHA
  (`benjjneb/dada2@b61fe2e`) in file headers, README and CITATION.cff, not by carrying its history.

Hard requirements:

- Modern packaging: `pyproject.toml`, **uv** for env/lock/build, **scikit-build-core** + CMake
  as build backend, **nanobind** bindings, `src/` layout, `py.typed`.
- Tests with **pytest** + **pytest-cov** (branch coverage). Results match the R package on the
  bundled example data (§8) before any larger dataset is attempted.
- The upstream C++ algorithms are reused, not re-implemented. R/Rcpp glue is replaced (§5).
  Performance-relevant R logic is ported to C++ (§6), staying as close to R semantics as possible.
- Python ≥ 3.11 supported; development on 3.13.

Non-goals:

- 1:1 API compatibility with R. Names and shapes follow Python conventions (§7).
- Re-implementing Bioconductor (Biostrings/ShortRead) in general. Only the narrow pieces
  DADA2 uses are re-implemented (§4.3, §6).
- GPU. (Linux, macOS **and Windows** wheels are all required; see Q6.)

### 1.1 Git & publishing rules (strict)

- **Never push, never open PRs or issues, and never change anything on `origin`
  (`benjjneb/dada2`).** Do not run `git push` to any remote.
- **Do not commit.** At the end of each phase, leave changes in the working tree and print a
  proposed commit message (`phase N: <summary>`). The human reviews, commits and pushes to
  `raeslab/pyDADA2`.
- Fresh-history mechanics: in phase 0, Claude prints (does not run) the commands for the human
  to run:
  - `git switch --orphan pydada2-main` (fresh history in the same working tree; untracked and
    working-tree files are kept);
  - `git remote add pydada2 https://github.com/raeslab/pyDADA2.git`;
  - after review, the first commit, then `git push pydada2 pydada2-main:main`.
  From then on the working tree is pyDADA2. The old `master` branch / `origin` remote stay
  untouched.
- Never publish to PyPI/TestPyPI. Never upload benchmark data or results anywhere.

---

## 2. Environment (already provided by the dev container)

| Tool | Notes |
|---|---|
| Ubuntu 24.04, gcc/g++ 13, clang/clangd, cmake, ninja, ccache | `CMAKE_GENERATOR=Ninja`, ccache enabled via env |
| gdb, valgrind, lcov, gcovr | C++ debugging and coverage |
| libboost-dev | header-only Boost.Math is available for local builds (wheels vendor it; §5.4) |
| uv 0.12.x, Python 3.13 (uv-managed) | venv at `/home/vscode/.venv` (outside the bind mount) |
| R (CRAN current) + Rcpp, RcppParallel, Biostrings, ShortRead, … | to build the **reference** R package and generate fixtures |
| Claude Code | `claude` |

`scripts/devcontainer-post-create.sh` runs on container creation. It installs the R package
from `reference/r-dada2/` (or the repo root before phase 0) into `$R_LIBS_USER`. It also runs
`uv sync --all-groups` once `pyproject.toml` exists.

---

## 3. Target repository layout

The Python package becomes the root project. Until the gated cleanup (phase 12), the R package
is kept unchanged under `reference/` so it can be built for parity testing.

```
.
├── pyproject.toml
├── CMakeLists.txt                 # top-level, used by scikit-build-core
├── uv.lock                        # committed
├── README.md                      # rewritten for the Python package (credit upstream)
├── LICENSE                        # LGPL-2 (unchanged; see §11)
├── CITATION.cff
├── CLAUDE.md                      # working notes/conventions for Claude Code (create in phase 0)
├── DEVIATIONS.md                  # every known difference vs R, with cause + magnitude (§8.3)
├── spec.md
├── Dockerfile, .dockerignore, .devcontainer/, .vscode/, scripts/
├── cpp/
│   ├── core/                      # pure C++17, NO Python or R headers
│   │   ├── dada.h                 # adapted from src/dada.h
│   │   ├── containers.cpp cluster.cpp pval.cpp kmers.cpp misc.cpp
│   │   ├── nwalign_endsfree.cpp nwalign_vectorized.cpp
│   │   ├── chimera.cpp taxonomy.cpp filter.cpp evaluate.cpp error.cpp
│   │   ├── dada_main.cpp          # was Rmain.cpp: dada_uniques() returning plain structs
│   │   ├── matrix.h               # tiny row-/col-major view types replacing Rcpp matrices
│   │   ├── parallel.h             # parallel_for replacing RcppParallel
│   │   ├── stats.h                # ppois upper tail via Boost.Math (§5.4)
│   │   ├── rng.h                  # seeded std::mt19937_64 streams (§5.4)
│   │   └── ported/                # NEW C++ code ported from R-level logic (§6)
│   │       ├── fastq_filter.cpp   # per-chunk filter/trim kernel
│   │       ├── derep.cpp          # streaming dereplication accumulator
│   │       ├── merge.cpp          # batch pair merging
│   │       ├── bimera_denovo.cpp  # per-sample de novo bimera loop
│   │       ├── collapse.cpp       # collapse_no_mismatch, is_shift_denovo
│   │       ├── complexity.cpp     # seq_complexity (k-mer Shannon richness)
│   │       ├── primers.cpp        # IUPAC/mismatch(/indel) primer matching
│   │       ├── species.cpp        # exact-substring species matching (PDict replacement)
│   │       └── qprofile.cpp       # per-cycle quality histograms for plotting
│   ├── bindings/
│   │   └── module.cpp             # nanobind module -> pydada2._core
│   └── tests/                     # optional C++ unit tests (Catch2 via FetchContent), off by default
├── src/pydada2/
│   ├── __init__.py                # public API re-exports, __version__
│   ├── _core.pyi                  # type stubs for the extension (generate with nanobind stubgen)
│   ├── py.typed
│   ├── _base.py                   # BaseEstimator-like mixin (get_params/set_params/repr)
│   ├── _validation.py             # input checking, check_is_fitted, NotFittedError
│   ├── _parallel.py               # n_jobs handling
│   ├── io.py                      # FASTQ/FASTA reading/writing, quality encoding detection
│   ├── derep.py                   # Dereplication / derep_fastq
│   ├── filtering.py               # filter_and_trim, fastq_filter, remove_primers
│   ├── error_models.py            # error-model estimators + learn_errors
│   ├── denoise.py                 # DADA estimator, DadaParams, DadaResult
│   ├── merge.py                   # merge_pairs
│   ├── seqtab.py                  # make_sequence_table, merge_sequence_tables, collapse_no_mismatch
│   ├── chimeras.py                # BimeraRemover, is_bimera*, remove_bimera_denovo
│   ├── taxonomy.py                # TaxonomyClassifier, SpeciesAssigner, add_species
│   ├── align.py                   # nwalign, nwhamming, eval_pair, kmer distances
│   ├── sequences.py               # reverse_complement, seq_complexity, is_phix, get_sequences
│   ├── plotting.py                # optional matplotlib plots
│   ├── datasets.py                # example file paths + bundled error matrices
│   └── data/                      # extdata copies, errBalancedF/R & tperr1 as .npy
├── tests/
│   ├── conftest.py
│   ├── reference/
│   │   ├── generate_reference.R   # produces fixtures from the R package
│   │   └── data/                  # committed fixtures (small JSON / CSV / .npy) + MANIFEST.json
│   ├── unit/                      # fast, per-module tests
│   ├── parity/                    # compare against R fixtures
│   └── e2e/                       # full pipeline on example data
├── benchmarks/
│   ├── micro/                     # per-kernel timings, Python vs R, example data (phase 10)
│   └── dataset/                   # GATED phase 11: download scripts, R & Python pipelines, report
├── docs/                          # optional: mkdocs-material + mkdocstrings
└── reference/
    └── r-dada2/                   # the original R package @ b61fe2e, copied in phase 0 (removed in phase 12)
        ├── DESCRIPTION NAMESPACE R/ src/ man/ inst/ data/ vignettes/ ...
```

Rules:

- `reference/r-dada2/` is **read-only** after phase 0 (except to fix an R build issue, documented).
- `cpp/core/` starts as a copy of `reference/r-dada2/src/*.cpp|h`, then gets de-R'd. Each file
  keeps a header comment naming the upstream file (`benjjneb/dada2@b61fe2e:src/<file>`) and what
  was changed. Each file in `cpp/core/ported/` names the R function(s) it was ported from
  (`benjjneb/dada2@b61fe2e:R/<file>.R::<function>`).
- Never `#include <Python.h>`/nanobind in `cpp/core/`, and never `#include <Rcpp.h>` anywhere.
- Upstream provenance references use the **commit SHA**, never the local `reference/` path.
  This keeps them valid after cleanup.

---

## 4. Build & packaging

### 4.1 pyproject.toml (shape, not verbatim)

```toml
[build-system]
requires = ["scikit-build-core>=0.10", "nanobind>=2.4"]
build-backend = "scikit_build_core.build"

[project]
name = "pydada2"                    # verify PyPI availability (Q1)
version = "0.1.0"
description = "Exact amplicon sequence variant inference (DADA2) for Python"
readme = "README.md"
license = "LGPL-2.0-or-later"       # placeholder, see Q9
requires-python = ">=3.11"
dependencies = [
  "numpy>=1.26",
  "pandas>=2.2",
  "dnaio>=1.2",                     # fast FASTQ/FASTA IO (gz/bz2/xz)
  "scikit-misc>=0.5",               # skmisc.loess == R's loess (netlib dloess) -> parity for LoessErrorModel
]

[project.urls]
Homepage = "https://github.com/raeslab/pyDADA2"
Source = "https://github.com/raeslab/pyDADA2"
Upstream = "https://github.com/benjjneb/dada2"

[project.optional-dependencies]
plot = ["matplotlib>=3.8"]

[dependency-groups]
dev = [
  "pytest>=8", "pytest-cov>=5", "pytest-xdist>=3", "hypothesis>=6", "pytest-benchmark>=4",
  "ruff>=0.6", "mypy>=1.11", "scikit-learn>=1.5",   # sklearn only to test clone()/get_params compat
  "nanobind>=2.4", "scikit-build-core>=0.10", "matplotlib>=3.8", "pre-commit>=3",
]

[tool.scikit-build]
minimum-version = "build-system.requires"
build-dir = "build/{wheel_tag}"
cmake.build-type = "Release"
wheel.packages = ["src/pydada2"]
editable.rebuild = true            # uv sync gives an editable install that rebuilds C++ on import

[tool.scikit-build.cmake.define]
DADA2_ENABLE_COVERAGE = {env = "DADA2_ENABLE_COVERAGE", default = "OFF"}

[tool.pytest.ini_options]
addopts = "-ra --strict-markers --cov=pydada2 --cov-branch --cov-report=term-missing"
testpaths = ["tests"]
markers = [
  "parity: compares against fixtures generated from the R package",
  "slow: long-running (full pipeline)",
  "bench: micro-benchmarks (not run by default)",
]

[tool.coverage.report]
fail_under = 0                     # ramped: 60 after phase 3, 80 after phase 7, 90 from phase 9
exclude_also = ["if TYPE_CHECKING:", "raise NotImplementedError"]

[tool.ruff]
line-length = 100
target-version = "py311"
[tool.ruff.lint]
select = ["E", "F", "W", "I", "B", "UP", "N", "SIM", "RUF", "NPY", "PD", "D"]
[tool.ruff.lint.pydocstyle]
convention = "numpy"

[tool.mypy]
strict = true
```

Notes:

- With `editable.rebuild = true`, scikit-build-core needs build deps in the env. Install with
  `uv sync --all-groups` and set `[tool.uv] no-build-isolation-package = ["pydada2"]` if
  rebuilds fail. Document the final choice in CLAUDE.md.
- Commit `uv.lock`.

### 4.2 CMake

- `CMakeLists.txt` uses:
  - C++17;
  - `find_package(Python 3.11 COMPONENTS Interpreter Development.Module REQUIRED)`;
  - `find_package(nanobind CONFIG REQUIRED)`, via `python -m nanobind --cmake_dir`.
- It builds a static lib `dada2_core` from `cpp/core/**/*.cpp` and the module `_core` via
  `nanobind_add_module(... NB_STATIC STABLE_ABI LTO)`, linking `dada2_core`. Install with
  `install(TARGETS _core LIBRARY DESTINATION pydada2)`.
- **Floating-point reproducibility:** compile `dada2_core` with `-ffp-contract=off` and no
  `-ffast-math`, so FMA contraction can't change results vs R (R builds with `-O2`, no contraction
  on x86_64). On MSVC use `/fp:precise` and `/fp:contract-` is the default (do not enable
  `/fp:fast`). Parity is required on Linux, macOS and Windows (Q6). Any last-bit platform
  difference must be investigated, kept below §8.3 tolerances, and listed in DEVIATIONS.md.
- Export `compile_commands.json` (`CMAKE_EXPORT_COMPILE_COMMANDS ON`) and symlink/copy it to
  `build/` so clangd works.
- SSE2: upstream uses `emmintrin.h` when `__x86_64`. Keep that, and also enable it for MSVC x64
  (`_M_X64`), which does not define `__x86_64`. On aarch64 the scalar path is used.
- Option `DADA2_ENABLE_COVERAGE` adds `--coverage -O0 -g` for gcov/gcovr.
- Option `DADA2_BUILD_CPP_TESTS` (default OFF) builds Catch2 tests from `cpp/tests/`.

### 4.3 Replacing Bioconductor / R dependencies

| R dependency | Used for | Python/C++ replacement |
|---|---|---|
| ShortRead (`FastqStreamer`, `readFastq`, `writeFastq`) | streaming FASTQ IO | `dnaio` (streaming, gz-aware); records are passed to C++ kernels per chunk |
| ShortRead (`trimTails`, `narrow`, `srsort`, `srrank`, `qa`) | trimming, dereplication, quality profiles | C++ kernels `fastq_filter`, `derep`, `qprofile` (§6) |
| Biostrings `reverseComplement`, `alphabetFrequency`, `oligonucleotideFrequency` | seq ops, N filter, complexity | C++ (`misc`, `fastq_filter`, `complexity`) |
| Biostrings `vmatchPattern`/`matchPattern` | primer matching | C++ `primers` (§6) |
| Biostrings `PDict`/`vcountPDict` | species assignment | C++ `species` (§6) |
| `stats::loess` | `loessErrfun` | `skmisc.loess.loess` (same Fortran) — required for parity |
| `parallel::mclapply` | per-file parallelism | `concurrent.futures` (threads around GIL-free C++; processes only if needed) |
| RcppParallel | intra-sample parallel loops in C++ | `cpp/core/parallel.h` (§5.3) |
| ggplot2 / reshape2 | plots | matplotlib (optional extra) |
| `.rda` data (`errBalancedF/R`, `tperr1`) | example error matrices | convert once with R to `.npy` in `src/pydada2/data/` |

---

## 5. C++ core: de-R'ing the code

Inventory (from `src/`): Rcpp is used in every file. RcppParallel `Worker`/`parallelFor` appears
in `cluster.cpp` (CompareParallel), `chimera.cpp` (BimeraTableParallel), `taxonomy.cpp`
(AssignParallel) and `Rmain.cpp` (FinalSubsParallel).

Randomness: `Rcpp::runif` in `taxonomy.cpp` draws the bootstrap uniforms. In addition,
**`std::random_device` + `std::mt19937` is used in `taxonomy.cpp` for tie-breaking between
equally likely genera**, which makes even R non-reproducible on ties. `Rcpp::ppois` is used in
`pval.cpp`.

Outputs are built as `Rcpp::DataFrame/List/Matrix` in `error.cpp`, `Rmain.cpp`, `chimera.cpp`,
`evaluate.cpp` and `taxonomy.cpp`. Exported entry points (see `RcppExports.cpp`):
`dada_uniques`, `C_is_bimera`, `C_table_bimera2`, `C_nwalign`, `C_eval_pair`,
`C_pair_consensus`, `C_isACGT`, `kmer_dist`, `kord_dist`, `kmer_matches`, `kdist_matches`,
`C_matchRef`, `C_matrixEE`, `C_nwvec`, `C_assign_taxonomy2`.

### 5.1 Principles

- **Minimal algorithmic diff.** Change only types/IO/glue. Any change to numerics must be
  justified in a comment, covered by a parity test, and listed in DEVIATIONS.md.
- Inputs: `std::vector`, `std::string`, and a `MatrixView<T>` (pointer + nrow + ncol + layout).
  R matrices are **column-major**, so be explicit. The 16×Q error matrix and quality matrices
  must be indexed exactly as upstream expects. In Python, `err` is stored as shape
  `(16, n_qual)` float64 (rows = transitions A2A, A2C, … T2T, same order as R) and passed as a
  contiguous array. The binding layer converts the layout if needed.
- Outputs: plain structs (e.g. `DadaOutput { std::vector<Cluster> clusters; std::vector<int> map; ... }`)
  with `std::vector`s. nanobind converts them to NumPy arrays / dicts; Python assembles pandas
  DataFrames.
- Errors: `Rcpp::stop(msg)` → `throw std::invalid_argument` (bad input) or `std::runtime_error`
  (internal). nanobind maps them to `ValueError`/`RuntimeError`.
- Logging: `Rprintf` → an optional `std::function<void(const std::string&)>` log callback
  (called only from the main thread), or removed if purely debug.
- Interrupts: replace `Rcpp::checkUserInterrupt()` with a callback `bool should_stop()`, checked
  between parallel chunks. The binding implements it by acquiring the GIL and calling
  `PyErr_CheckSignals()`, which raises `KeyboardInterrupt`.
- Release the GIL around all heavy calls.
- Keep the fixed constants (`SEQLEN`, `KMER_SIZE`, …) identical.

### 5.2 Memory

Upstream uses `malloc/free` and raw pointers. Keep them (minimal diff). Run the parity suite
under valgrind or an ASan build once per phase that touches C++. Fix leaks on error paths
introduced by replacing `Rcpp::stop` with exceptions, using RAII guards where a throw can skip
a free. New code in `cpp/core/ported/` uses RAII/`std::vector` only, with no raw `malloc`.

### 5.3 Parallelism (`parallel.h`)

Provide `parallel_for(size_t begin, size_t end, F&& body, size_t grain, int n_threads)`
implemented with `std::thread`, using simple static chunking or a small pool. `n_threads <= 1`
runs serially in the calling thread. Convert each `RcppParallel::Worker` into a lambda/functor.
Results must be **deterministic and independent of `n_threads`**. Add a test that runs with 1 and
4 threads and asserts identical outputs. No TBB dependency, to keep wheels simple.

### 5.4 Stats & RNG

- **No code copied from R itself** (R's `nmath`/RNG are GPL-2+; the package must stay LGPL; §11).
- `Rcpp::ppois(n, lambda, lower=false)` → `ppois_upper(n, lambda)` = P(X > n) =
  `boost::math::gamma_p(n + 1, lambda)`. Vendor header-only Boost.Math "standalone" via CMake
  `FetchContent`, pinned to a release tag, so wheels don't need system Boost. Required accuracy:
  rel. err ≤ 1e-12 vs R across the tested range, including tiny tails (~1e-300) and the
  `TAIL_APPROX_CUTOFF` branch. If that fails, write a **clean-room** implementation from the
  literature (e.g. Temme / continued fraction), never a copy of R's source, and record it in
  DEVIATIONS.md.
- `Rcpp::runif` (bootstrap uniforms) → `std::mt19937_64` seeded from `random_state` (an int, or
  a `numpy.random.Generator` that draws a seed). Pre-generate the uniforms in C++ exactly as
  upstream does (same count and usage), so only the generator differs. **Bit-identical
  bootstrap values vs R are not a goal** (decided); parity is statistical (§8.3).
- Tie-breaking (`std::random_device` in upstream `taxonomy.cpp`) → a second generator seeded
  deterministically from `random_state`, so Python is reproducible. R is non-deterministic here,
  so parity on sequences with tied top hits is checked only as "assigned genus ∈ set of tied
  genera".
- Determinism of the generator stream must not depend on `n_jobs`: the uniforms are drawn up front
  and the tie-breaking stream is per-sequence (seeded from `(random_state, sequence index)`).

### 5.5 Binding surface (`pydada2._core`, private)

There is one nanobind function per upstream export, snake_case, documented in `_core.pyi`:
`dada_uniques`, `is_bimera`, `table_bimera`, `nwalign`, `eval_pair`, `pair_consensus`,
`is_acgt`, `kmer_dist`, `kord_dist`, `kmer_matches`, `kdist_matches`, `match_ref`,
`matrix_ee`, `nwvec`, `assign_taxonomy`. If the refactor is clean, add `TaxonomyIndex`, a
fit-once reference k-mer model, so that `TaxonomyClassifier.fit/predict` doesn't rebuild the
reference model per call. The §6 kernels come on top. `_core` is private; the public API is Python.

---

## 6. Porting R-level logic to C++ (performance)

### 6.1 Why

A large share of DADA2's pipeline logic is written in R, not C++. It falls into two groups:

- **(A) R code that is fast only because it calls Bioconductor C code** (ShortRead/Biostrings).
  A naive Python port with per-read loops would be *much slower* than R, so these **must** become
  C++ kernels to avoid a regression.
- **(B) R-level loops over sequences/pairs** (`lapply`/`mapply`/`for` around `nwalign`,
  `C_is_bimera`, `C_eval_pair`, `grepl`). Moving the whole loop into C++ (with `parallel_for`)
  removes per-call overhead and enables threading: a genuine **speed-up over R**.

What stays in Python/NumPy/pandas: orchestration, file handling, the self-consistency loop in
`learnErrors`/`dada`, pseudo-pooling bookkeeping, loess error models (`skmisc`), sequence-table
construction/merging, taxonomy-file parsing and formatting helpers, and plotting.

### 6.2 Kernels

| # | Kernel (`cpp/core/ported/`) | Ported from (R) | Group | What moves to C++ | Parity target |
|---|---|---|---|---|---|
| K1 | `fastq_filter` | `fastqFilter`, `fastqPairedFilter` (filter.R), `.nFilter`, `trimTails` | A | One pass per chunk of records: orient (`orient.fwd`), `maxLen`, `trimLeft`, `trimRight`, `truncQ` (first base with Q ≤ truncQ), `truncLen`, `minLen`, `maxN`, `minQ`, `maxEE`, `rm.phix` (existing `match_ref`), `rm.lowcomplex` (K6). Returns keep mask + trimmed (start, end) per read + reverse-complement flag. The paired variant ANDs both masks; ID matching (`matchIDs`) stays in Python. | exact (counts, read IDs, sequences, qualities) |
| K2 | `derep` | `derepFastq`, `qtables2` (sequenceIO.R) | A | Streaming accumulator over chunks: hash map seq → (count, Σquals as double, first-seen rank), plus read→unique map. It reproduces R's ordering exactly: **within each chunk new uniques are appended in lexical order of that chunk; final order is a stable sort by decreasing abundance**, so ties depend on chunk size `n` as in R. Zero-length reads → map = −1 + warning. Quality columns padded with NaN. | exact (uniques, order, map); quals rtol 1e-12 |
| K3 | `merge` | `mergePairs` (paired.R) | B | Given unique (F,R) denoised pairs for a sample, run per pair in parallel: rc(R), unbanded `nwalign` with merge scoring (`MATCH=1`, `MISMATCH=-64/-8`, `GAP=-64/-8` per `maxMismatch`), `eval_pair`, `pair_consensus(prefer, trimOverhang)`, `justConcatenate`. Pair tabulation, `propagateCol`, sorting and `returnRejects` stay in pandas. | exact |
| K4 | `bimera_denovo` | `isBimeraDenovo` (chimeras.R) | B | Loop over uniques: select parents (`abund > minFold·a & abund > minParentAbundance`), call `is_bimera` and `parallel_for` over sequences. Serves `remove_bimera_denovo(method="pooled"|"per-sample")`. `"consensus"` already uses `table_bimera`. | exact |
| K5 | `collapse` | `collapseNoMismatch` (multiSample.R), `isShiftDenovo`/`isShift`/`isShiftedPair` (chimeras.R) | B | Greedy collapse in decreasing-abundance order with prefix prescreen + `nwhamming == 0`. The outer loop stays sequential (order-dependent). The inner search over accepted refs may be parallel with a **min-index reduction**, so the first matching ref wins exactly as in R. Shift detection: all-vs-more-abundant `nwalign` + `eval_pair` in parallel. | exact |
| K6 | `complexity` | `seqComplexity`, `sindex` (filter.R) | A | k-mer counts (`oligonucleotideFrequency`, non-ACGT k-mers skipped as Biostrings does) → exp(Shannon). Includes the windowed variant (`window`, `by`, running min). | rtol 1e-12 |
| K7 | `primers` | `removePrimers` (filter.R) | A | Primer search with `max.mismatch`, IUPAC semantics matching Biostrings `fixed=FALSE` when the primer is non-ACGT, first fwd hit / last rev hit, orientation (`orient=TRUE`), trimming. `allow.indels=TRUE`: banded edit-distance search mirroring `matchPattern(with.indels=TRUE)` start/end reporting. | exact without indels; with indels "very similar" allowed (≥ 99.9 % reads identical, deviations listed) |
| K8 | `species` | `assignSpecies` (taxonomy.R) | A | Exact-substring matching of each query against all references (and their rc if `tryRC`): hash all length-L windows of the references for each distinct query length L. Returns the hit lists; `mapHits`/`matchGenera`/`allowMultiple` stay in Python. | exact |
| K9 | `qprofile` | `plotQualityProfile` (`ShortRead::qa`) | A | Per-cycle quality histograms + read-length distribution for a chunk; Python aggregates for plotting. | exact counts |

Notes:

- `matrix_ee`, `match_ref` (phiX), `table_bimera` and `assign_taxonomy` are already C++ upstream;
  K1/K4/K8 only call them in a loop C++-side.
- Each kernel gets (i) a parity test vs an R fixture and (ii) a micro-benchmark in
  `benchmarks/micro/` vs the R function on the example data (§9). Expected gains are modest on the
  tiny example files. Real speed-ups are confirmed on the larger dataset (phase 11).
- **R behaviour, not R docs, is the spec.** Example: `fastqFilter` keeps reads with
  `min(qual) > minQ`, although the docs say "less than minQ is discarded". Port the code
  behaviour and note it in DEVIATIONS.md under "documented R quirks preserved".
- Known ordering trap: R sorts DNA with Biostrings/ShortRead ordering. This equals ASCII order for
  A, C, G, N, T but differs for other IUPAC letters. Verify the lexical order used by K2 against R
  on a fixture containing IUPAC codes.

---

## 7. Public Python API

### 7.1 Conventions

- snake_case everywhere; classes in CapWords. Positional-only for the main data argument,
  keyword-only for options (`def f(x, /, *, opt=...)`).
- Estimator pattern (scikit-learn-like), **without** a hard sklearn dependency:
  - All hyper-parameters are set in `__init__` and stored unmodified under the same names. No
    logic in `__init__`; validation happens in `fit`.
  - `fit(...)` returns `self`. Learned state ends with a trailing underscore (`err_`, `results_`).
  - Implement `get_params(deep=True)`, `set_params(**params)` and an informative `__repr__`
    showing non-default params in `_base.py`. Test that `sklearn.base.clone` works on every
    estimator.
  - Calling a method that needs fitted state before `fit` raises `pydada2.NotFittedError`
    (subclass of `ValueError` and `AttributeError`, like sklearn).
- Common parameter names: `n_jobs` (None → 1, -1 → all cores), `random_state`, `verbose` (int, 0 = silent).
  Use `logging.getLogger("pydada2")`; `verbose` only adjusts what's emitted/progress output.
- No global mutable options. R's `setDadaOpt/getDadaOpt` become a frozen dataclass
  `DadaParams` that `DADA` accepts. Its fields can also be passed directly as kwargs.
- Paths: accept `str | os.PathLike`, and lists thereof. Directories are **not** silently
  expanded; provide `pydada2.io.list_fastq(directory, pattern=...)` instead.
- Return types: pandas `DataFrame`/`Series` for tabular data with meaningful index/columns;
  NumPy arrays for matrices; small frozen dataclasses for structured results. Full type hints
  and numpy-style docstrings, each with an `Examples` section that runs as a doctest on bundled
  data.
- Sequence tables: `pandas.DataFrame`, index = sample names, columns = ASV sequences,
  dtype `int64`. Columns are ordered by decreasing total abundance, like R.

### 7.2 Proposed API (adapt where the implementation suggests something cleaner; record deviations in CLAUDE.md)

```python
import pydada2 as dada2
from pydada2 import (
    filter_and_trim,
    Dereplication,
    derep_fastq,
    LoessErrorModel,
    NoQualErrorModel,
    BinnedQualErrorModel,
    PacBioErrorModel,
    learn_errors,
    DADA,
    DadaParams,
    merge_pairs,
    make_sequence_table,
    BimeraRemover,
    TaxonomyClassifier,
)

stats = filter_and_trim(
    fwd,
    filt_fwd,
    rev=rev,
    filt_rev=filt_rev,
    trunc_len=(240, 160),
    max_ee=(2, 2),
    trunc_q=2,
    rm_phix=True,
    n_jobs=-1,
)
# -> DataFrame[index=input file, columns=["reads_in", "reads_out"]]

err_f = learn_errors(
    filt_fwd, n_bases=1e8, error_model=LoessErrorModel(), n_jobs=-1, random_state=0
)
# -> fitted ErrorLearner-like result: .err_ (16 x Q), .err_in_, .trans_, .n_iter_, .converged_

drp_f = [derep_fastq(p) for p in filt_fwd]  # Dereplication objects
dada_f = DADA(err=err_f.err_, pool=False, n_jobs=-1).fit(drp_f)
dada_f.results_  # list[DadaResult], one per sample (pool="pseudo" | True | False)
dada_f.results_[0].denoised  # pd.Series: ASV sequence -> abundance
dada_f.results_[0].clusters  # pd.DataFrame (R's $clustering)
dada_f.results_[0].map  # np.ndarray: unique index -> ASV index (R's $map, 0-based, -1 for NA)

mergers = merge_pairs(dada_f, drp_f, dada_r, drp_r, min_overlap=12, max_mismatch=0)
seqtab = make_sequence_table(mergers)  # DataFrame samples x ASVs
seqtab_nochim = BimeraRemover(method="consensus", n_jobs=-1).fit_transform(seqtab)

clf = TaxonomyClassifier(min_boot=50, try_rc=False, random_state=0, n_jobs=-1).fit(
    "silva_train.fa.gz"
)
taxa = clf.predict(seqtab_nochim.columns)  # DataFrame seq x levels (None where unassigned)
boot = clf.predict_bootstrap(seqtab_nochim.columns)  # DataFrame seq x levels, int 0..100
```

R → Python name map. **All 38 NAMESPACE exports must be covered.** A test parses
`tests/reference/data/NAMESPACE_exports.txt` (a copy of the export list, so it survives cleanup)
and asserts that each export maps to a Python name or appears in an explicit `NOT_PORTED` dict
with a reason.

| R | Python |
|---|---|
| `filterAndTrim`, `fastqFilter`, `fastqPairedFilter` | `filter_and_trim` (paired via `rev=`), `fastq_filter` (low level) |
| `removePrimers` | `remove_primers` |
| `derepFastq` (+ `derep` class) | `derep_fastq` → `Dereplication` (frozen dataclass: `uniques: pd.Series`, `quals: np.ndarray`, `map: np.ndarray`) |
| `learnErrors` | `learn_errors` (function) returning a fitted `ErrorLearner`; `ErrorLearner` also usable directly |
| `loessErrfun`, `noqualErrfun`, `makeBinnedQualErrfun`, `PacBioErrfun` | `LoessErrorModel`, `NoQualErrorModel`, `BinnedQualErrorModel(bins)`, `PacBioErrorModel`. These are objects with `fit(trans) -> self` and `err_`, and are also callable as `model(trans) -> err`. |
| `inflateErr`, `getErrors` | `inflate_err`, `get_errors` |
| `dada` (+ `selfConsist`, `pool`, `priors`) | `DADA(err=None, error_model=..., self_consistent=False, pool=False, priors=(), params=DadaParams(), n_jobs=None, verbose=0, **param_overrides).fit(derep_or_list)` |
| `setDadaOpt`, `getDadaOpt` | `DadaParams` (frozen dataclass; fields lower-case: `omega_a=1e-40`, `omega_p=1e-4`, `omega_c=1e-40`, `detect_singletons=False`, `use_kmers=True`, `kdist_cutoff=0.42`, `max_consist=10`, `match=5`, `mismatch=-4`, `gap_penalty=-8`, `band_size=16`, `vectorized_alignment=True`, `max_clust=0`, `min_fold=1`, `min_hamming=1`, `min_abundance=1`, `use_quals=True`, `homopolymer_gap_penalty=None`, `sse=2`, `gapless=True`, `greedy=True`, `pseudo_prevalence=2`, `pseudo_abundance=inf`). Defaults must equal `reference/r-dada2/R/dada.R`. |
| `mergePairs` | `merge_pairs` → DataFrame (one) or list[DataFrame] |
| `makeSequenceTable`, `mergeSequenceTables`, `collapseNoMismatch` | `make_sequence_table`, `merge_sequence_tables`, `collapse_no_mismatch` |
| `isBimera`, `isBimeraDenovo`, `isBimeraDenovoTable`, `removeBimeraDenovo` | `is_bimera`, `is_bimera_denovo`, `is_bimera_denovo_table`, and `BimeraRemover(method="consensus"|"pooled"|"per-sample").fit(seqtab)` with `is_bimera_` mask and `.transform()`; convenience `remove_bimera_denovo` |
| `isShiftDenovo` | `is_shift_denovo` |
| `assignTaxonomy` | `TaxonomyClassifier(...).fit(ref).predict(seqs)`; convenience `assign_taxonomy(seqs, ref, ...)` |
| `assignSpecies`, `addSpecies` | `SpeciesAssigner`, `assign_species`, `add_species` |
| `nwalign`, `nwhamming` | `pydada2.align.nwalign`, `pydada2.align.nwhamming` |
| `rc` | `reverse_complement` |
| `seqComplexity`, `isPhiX` | `seq_complexity`, `is_phix` |
| `getSequences`, `getUniques`, `uniquesToFasta`, `getErrors` | `get_sequences`, `get_uniques`, `uniques_to_fasta`, `get_errors` |
| `plotQualityProfile`, `plotErrors`, `plotComplexity` | `pydada2.plotting.plot_quality_profile`, `plot_errors`, `plot_complexity`. They return `matplotlib.figure.Figure`, import matplotlib lazily, and raise a clear error if the extra is missing. |
| internal `makeTaxonomyFasta_*`, `makeSpeciesFasta_*` (not exported) | optional `pydada2.taxonomy.reference_builders` (low priority) |
| `show` methods | `__repr__` on result classes |
| `c()` methods on dada/derep | not ported (use lists) |

Indexing: all public indices are **0-based**; NA → `-1` in integer arrays (document it).

---

## 8. Testing strategy

### 8.1 Layout & commands

```bash
uv run pytest                          # all fast tests with coverage
uv run pytest -m "not slow"            # quick loop
uv run pytest -m parity                # R-parity tests
uv run pytest -n auto                  # parallel (xdist)
uv run pytest -m bench -p no:cacheprovider --benchmark-only   # micro-benchmarks (phase 10)
DADA2_ENABLE_COVERAGE=ON uv sync --reinstall-package pydada2 && uv run pytest && \
  gcovr -r . --filter cpp/core --print-summary   # C++ line coverage
```

Targets: Python branch coverage ≥ 90 % from phase 9 (ramped before that, see §4.1). C++ core line
coverage is reported (aim ≥ 75 %).

### 8.2 Reference fixtures from R

`tests/reference/generate_reference.R` runs with `Rscript` against the installed reference
package. It writes **small, committed** fixtures to `tests/reference/data/` using `jsonlite`
(and `.csv` for tables).

- Seed R with `set.seed(100)`.
- Use the bundled `inst/extdata` files: `sam1F/R`, `sam2F/R`, `samPB*`, `example_train_set.fa.gz`,
  `example_species_assignment.fa.gz`, `ten_16s.100.fa.gz`, `phix_genome.fa`, `example_seqs.fa`.
- Record `packageVersion("dada2")`, the R version, `sessionInfo()` package versions and the
  upstream commit SHA in `tests/reference/data/MANIFEST.json`.
- The script must also run against a **Bioconductor-installed dada2 1.36**, not only
  `reference/r-dada2`, so fixtures can be regenerated after cleanup.

Cover at least:

- `filterAndTrim` read counts + output read IDs/seqs/quals for several parameter combos (single &
  paired; `truncQ`, `truncLen`, `trimLeft/Right`, `maxN`, `minQ`, `maxEE`, `rm.phix`,
  `rm.lowcomplex`, `orient.fwd`, `matchIDs`).
- `removePrimers` on `samPBprimers` (with and without `allow.indels`, `orient=TRUE`).
- `derepFastq` uniques, order, average quals, map — with default `n` **and** a small `n`
  (e.g. 500), to pin down the chunk-dependent tie ordering.
- `nwalign`, `nwhamming`, `C_eval_pair`, `C_pair_consensus`, `kmer_dist` on a seq panel.
- `seqComplexity` (default, `kmerSize=3`, windowed) and `isPhiX`, `rc`.
- `loessErrfun`, `noqualErrfun`, `makeBinnedQualErrfun`, `PacBioErrfun` output for fixed `trans`.
- `learnErrors` on the example forward/reverse files (err_out, trans, convergence).
- `dada` with fixed `err` (e.g. `errBalancedF`/R) for: default, `pool=TRUE`, `pool="pseudo"`,
  `USE_QUALS=FALSE`, `GAPLESS=FALSE`, `BAND_SIZE=32`, `priors`, `selfConsist=TRUE`. Store
  `denoised`, `map`, `clustering` (birth_type, birth_pval, birth_fold, …) and `pval`.
- `mergePairs` incl. `maxMismatch>0`, `returnRejects`, `justConcatenate`, `trimOverhang`,
  `propagateCol`.
- `makeSequenceTable`, `mergeSequenceTables`, `collapseNoMismatch`, `isShiftDenovo`,
  `removeBimeraDenovo` (all methods), `isBimera`, `isBimeraDenovo`.
- `assignTaxonomy` (`set.seed(100)`, `minBoot=0`, `outputBootstraps=TRUE`, `tryRC` both ways),
  `assignSpecies` (`allowMultiple`, `tryRC`), `addSpecies`. For taxonomy, also store per sequence
  whether the top hit was tied (needed for §5.4).
- Quality-profile summary stats used by `plotQualityProfile`.
- Export `errBalancedF`, `errBalancedR`, `tperr1` to `.npy` for `src/pydada2/data/`.
- Copy the NAMESPACE export list to `NAMESPACE_exports.txt`.

Parity tests load fixtures. They must **not** require R at test time (R is only needed to
regenerate).

### 8.3 Tolerances & deviation budget

| Output | Comparison |
|---|---|
| sequences, counts, maps, ordering, flags, filter counts, merge results, bimera calls, ASV sets, species hits | exact |
| C++ doubles (lambda, p-values, birth stats, alignment scores, averaged quals) | `rtol=1e-10` |
| loess error rates | `rtol=1e-6` |
| taxonomy | top-hit assignment (`min_boot=0`) exact for untied sequences, tied cases genus ∈ tied set; bootstrap values: median `|boot_py − boot_R|` ≤ 10 and every entry within 3σ of binomial noise; at `min_boot=50`, assignments identical at every level whose R bootstrap is outside 50 ± 3σ (≈ ±23). *(Refined in phase 8: the original "≥ 95 % identical" is not attainable on small fixture sets with many confidences near 50, because the RNG differs by design.)* |
| `remove_primers(allow_indels=True)` | ≥ 99.9 % reads identical |

"Negligible differences" are acceptable **only** if they meet all of these conditions:

1. The root cause is understood.
2. They are listed in `DEVIATIONS.md` with cause, affected outputs and magnitude.
3. They have a test pinning the current magnitude, so it can't silently grow.
4. They do not change any exact-category output on the example data (ASV set, counts, chimera
   calls, merged sequences).

Never loosen a tolerance or skip a parity test to get green.

### 8.4 Other tests

- Unit tests for every public function and every §6 kernel cover:
  - the happy path;
  - edge cases: empty input, N's, lowercase, sequence length ≥ SEQLEN, mismatched lengths,
    zero-length reads, all reads filtered;
  - error messages, via `pytest.raises(..., match=)`.
- Estimator compliance: `get_params`/`set_params`/`clone` round-trip, `NotFittedError`,
  `fit` returns `self`, no mutation of inputs.
- Determinism: same output for `n_jobs=1` vs `n_jobs=4` (every kernel that uses
  `parallel_for`), and the same `random_state` gives the same output.
- Hypothesis property tests: `reverse_complement` is an involution; `nwalign`
  symmetry/score bounds; dereplication preserves total read count; the derep result is
  independent of chunk size *except* for tie order; filtering never lengthens reads.
- e2e (`slow`): the MiSeq-style vignette pipeline on `sam1/sam2` matches the R pipeline's
  final `seqtab.nochim` and taxonomy within §8.3.

---

## 9. Benchmarking

### 9.1 Micro-benchmarks on example data (phase 10, in scope now)

- `benchmarks/micro/run_r.R` and `benchmarks/micro/test_bench_*.py` (pytest-benchmark) time each
  stage and each §6 kernel on the bundled example data: filter, derep, learn_errors, dada
  (`pool` False/True/pseudo), merge, seqtab, chimeras (all methods), collapse, taxonomy, species.
- Run single-threaded and with 4 threads, ≥ 5 repeats, and report the median.
- Output: `benchmarks/micro/REPORT.md` with a Python/R ratio per stage. Raw JSON stays out of git.
- Expectation (not a gate): Python ≤ 1.5× R for every stage single-threaded, and faster than R
  for K3–K5 multithreaded. Investigate and explain any stage > 1.5×.

### 9.2 Larger dataset (phase 11, GATED, details to be agreed with the human)

Default proposal, to be confirmed in the phase 11 interview:

- **Dataset:** MiSeq SOP (Kozich et al., 20 samples, 16S V4, paired 2×250) + SILVA training set
  (the version the DADA2 tutorial currently uses), and optionally a second dataset of a different
  type (ITS with variable lengths + primer removal, or PacBio full-length 16S).
- Data is fetched by `benchmarks/dataset/download.sh` (SHA-256 checked) into a git-ignored
  directory. It is never committed and never part of CI.
- Identical pipelines: `benchmarks/dataset/pipeline.R` and `pipeline.py` with the same parameters.
- Results compared:
  - the per-step read-tracking table: exact;
  - the ASV set and sequence table: exact;
  - chimera calls: exact;
  - taxonomy per §8.3.
- Performance: wall time per stage (hyperfine or repeated runs), peak RSS
  (`/usr/bin/time -v`), at 1 and N threads.
- Output: `benchmarks/dataset/REPORT.md` (+ optional HTML).

---

## 10. Phased plan

Work strictly phase by phase. At the end of each phase:

- all gates are green (§10.1);
- `CLAUDE.md` and `DEVIATIONS.md` are updated;
- print a summary and a proposed commit message, then **stop** (no commit, no push; §1.1).

**Phase 0 — Restructure & reference**
- Move (plain `mv`; history is fresh, §1.1) the R package into `reference/r-dada2/`: DESCRIPTION, NAMESPACE, NEWS, R/, src/,
  man/, inst/, data/, vignettes/, tests/ (if present), .Rbuildignore, .travis.yml, README.md as
  `reference/r-dada2/README.md`. Keep LICENSE at the root and copy it into reference too.
- Re-run `scripts/devcontainer-post-create.sh` → the R package installs from `reference/r-dada2`.
- Write `tests/reference/generate_reference.R`, generate fixtures (§8.2), and leave them in the
  working tree for the human to commit.
- Create `CLAUDE.md` (project conventions, commands, git rules §1.1, key rules of this spec)
  and an empty `DEVIATIONS.md`.
- ✅ `Rscript -e 'library(dada2); packageVersion("dada2")'` prints 1.36.0; fixtures generated;
  MANIFEST.json complete.

**Phase 1 — Skeleton & build**
- Add `pyproject.toml`, `CMakeLists.txt`, `src/pydada2/__init__.py`, `_core` with a trivial
  function, `_base.py`, `_validation.py`, ruff/mypy/pytest config, pre-commit, and CI workflow
  (`.github/workflows/ci.yml`, written but not pushed).
- ✅ `uv sync` builds the extension; `uv run python -c "import pydada2, pydada2._core"` works;
  `uv build` produces sdist + wheel; `uv run pytest` passes; editing a `.cpp` triggers a rebuild.

**Phase 2a — C++ core compiles without R**
- Copy `reference/r-dada2/src` → `cpp/core` and remove Rcpp/RcppParallel per §5.
- Add `matrix.h`, `parallel.h`, `stats.h` (Boost.Math `ppois_upper`) and `rng.h`, each with
  its own unit tests. `ppois_upper` is checked against R values over a grid (fixture).
- ✅ `dada2_core` builds with zero R headers; `ppois_upper` within rel. 1e-12 of R.

**Phase 2b — Bindings & low-level parity**
- Bind all upstream functions (§5.5) and write `_core.pyi`.
- ✅ Low-level parity tests pass for nwalign/nwhamming/eval_pair/pair_consensus/kmer_dist/
  matrix_ee/is_bimera/table_bimera and `dada_uniques` with a fixed err matrix. The
  thread-count determinism test passes. No valgrind/ASan errors on the parity tests.

**Phase 3 — IO, dereplication, filtering (K1, K2, K6, K7, K9)**
- `io.py` (dnaio, phred33; `quality_type="auto"|"FastqQuality"|"SFastqQuality"` as upstream),
  `derep.py` (K2), `filtering.py` (`filter_and_trim`, `fastq_filter`, paired, `match_ids`,
  `rm_phix`, `rm_lowcomplex`, `orient_fwd`, chunked streaming of `n` reads, `n_jobs` over files,
  using K1/K6), `remove_primers` (K7), `sequences.py`, `datasets.py`, and the quality-profile
  data function (K9).
- ✅ Filter/derep/complexity/primer parity per §8.3. Output FASTQ gz files are equivalent in
  content (records), not necessarily in compression.

**Phase 4 — Error models**
- Error-model classes, `ErrorLearner`/`learn_errors`, `inflate_err`, `get_errors`.
- ✅ Error-model and learn_errors parity within §8.3.

**Phase 5 — DADA estimator**
- `DadaParams`, `DADA` (single/list input, `pool=False|True|"pseudo"`, `priors`,
  `self_consistent`), `DadaResult`.
- ✅ All `dada` parity fixtures pass.

**Phase 6 — Merging & sequence tables (K3, K5-collapse)**
- `merge_pairs` (K3), `make_sequence_table`, `merge_sequence_tables`, `collapse_no_mismatch` (K5).
- ✅ Parity exact; n_jobs determinism.

**Phase 7 — Chimeras (K4, K5-shift)**
- `is_bimera`, `is_bimera_denovo` (K4), `is_bimera_denovo_table`, `is_shift_denovo` (K5),
  `BimeraRemover`, `remove_bimera_denovo`.
- ✅ Parity exact for all methods.

**Phase 8 — Taxonomy (K8)**
- Split `C_assign_taxonomy2` into fit (reference k-mer model) / predict if feasible. Implement
  `TaxonomyClassifier`, `assign_taxonomy`, `SpeciesAssigner`, `assign_species` (K8),
  `add_species`.
- ✅ Parity per §8.3; `random_state` determinism test (incl. `n_jobs=1` vs `4`).

**Phase 9 — Plotting, docs, e2e**
- `plotting.py` (optional extra).
- README rewrite: install, quickstart pipeline, R→Python table, link to raeslab/pyDADA2, and
  credit/citation of Callahan et al. 2016 + "unofficial port".
- `CITATION.cff`, docstring examples as doctests, optional mkdocs site, and the e2e slow test.
- ✅ Coverage ≥ 90 %, e2e passes, `uv build` + `twine check` clean, wheel installs in a fresh
  venv and runs the quickstart.

**Phase 10 — Parity sign-off & micro-benchmarks (example data)**
- Run the full parity + e2e suite on x86_64 and record results in `PARITY.md`, with one row per
  R function: status (exact / within tolerance / deviation #n).
- Run §9.1 micro-benchmarks and write `benchmarks/micro/REPORT.md`. Profile and optimise any stage
  > 1.5× R, without changing results.
- ✅ Every exact-category output identical on the example data; all deviations documented
  (§8.3); benchmark report written. **Then stop and hand over to the human.** This is the end of
  the current pass.

**Phase 11 — Larger-dataset validation & benchmark — GATED**
- Start only on explicit human instruction. First interview the human about the dataset(s),
  reference database version, parameters, hardware/threads and acceptance thresholds. Then
  implement §9.2.
- ✅ Agreed acceptance criteria met; `benchmarks/dataset/REPORT.md` written; any new deviations
  are investigated and documented.

**Phase 12 — Remove the R reference (cleanup) — GATED**
- **Never start this phase without an explicit instruction from the human** (e.g. "run phase 12").
  Before starting, confirm that phases 10 (and 11, if run) are signed off.
- Steps:
  1. `git rm -r reference/r-dada2/` (pyDADA2's early history keeps it).
  2. Make `tests/reference/generate_reference.R` and the R benchmark scripts install/use
     Bioconductor `dada2==1.36` (pinned) instead of the local copy. Document the regeneration
     procedure in `tests/reference/README.md`.
  3. Remove R/Bioconductor from the `Dockerfile` and `devcontainer-post-create.sh`. Optionally
     keep a separate `Dockerfile.reference` / devcontainer profile for regenerating fixtures.
  4. Grep the repo for `reference/r-dada2`, `Rscript`, `R_LIBS_USER` and fix every remaining
     hit. Provenance headers already use the upstream SHA (§3).
  5. Remove the upstream R-specific files at the root that no longer apply (`.Rbuildignore`,
     `.travis.yml` if they remain).
  6. Re-run the full gate suite (fixtures are committed, so parity tests still run without R).
- ✅ No R files or R tooling are required to build, test or develop; all tests green; the human
  reviews and commits.

### 10.1 Quality gates (every phase)

`uv run ruff format --check . && uv run ruff check . && uv run mypy src && uv run pytest -m "not slow"`
Add these as a `pre-commit` config and a GitHub Actions workflow (`.github/workflows/ci.yml`:
Ubuntu + macOS (arm64 + x86_64) + Windows, Python 3.11–3.13, parity suite on every OS, `uv sync`, tests; build wheels with `cibuildwheel` on tags).
CI does not install R; it uses committed fixtures. The workflow is written but only runs once the
human pushes to raeslab/pyDADA2.

---

## 11. Licensing & attribution

Upstream DADA2 is LGPL-2 (see `LICENSE`, `DESCRIPTION`). The Python package includes derived C++
code (and C++ ported from R code), so it stays under the same license.

- Keep upstream copyright/author notices in every ported file.
- Cite the DADA2 paper in README and `CITATION.cff`.
- State that this is an unofficial port (unless/until upstream adopts it).
- Never copy code from R itself (GPL-2+); Boost.Math (BSL-1.0) is LGPL-compatible (§5.4).
- The R-to-C++ kernels (§6) are ports of DADA2's own LGPL R code, so they stay under LGPL with
  upstream notices.
- Flag any uncertainty about the exact SPDX identifier (`LGPL-2.0-only` vs `-or-later`) to the
  human instead of guessing.

---

## 12. Decisions & open questions

| # | Topic | Status |
|---|---|---|
| Q1 | Distribution/import name | Decided `pydada2` / `pydada2`, **but `pydada2` is already taken on PyPI** (omicverse/py-dada2, a pure-Python port). Placeholder kept until the human picks a distribution name (import name can stay or change) |
| Q2 | Commits | **Decided:** Claude never commits or pushes; the human does (§1.1) |
| Q3 | Repo history | **Decided:** fresh history (orphan branch → raeslab/pyDADA2), upstream credited by SHA |
| Q4 | ppois / RNG | **Decided:** Boost.Math + `std::mt19937_64`; no R code copied; taxonomy bootstrap parity statistical |
| Q5 | Kernels in scope this pass | **Decided:** all K1–K9 |
| Q6 | Platforms that must pass parity | **Decided:** Linux, macOS and Windows (CI runs the parity suite on all three) |
| Q7 | Minimum Python | **Decided:** 3.11 |
| Q8 | Larger dataset(s), reference DB, thresholds | Open: interview at the start of phase 11 |
| Q9 | License | **Decided:** LGPL (continues DADA2 code). Exact version open: upstream `DESCRIPTION` says `LGPL-2` but upstream `LICENSE` contains LGPL-**3** text. Placeholder `LGPL-2.0-or-later` (the bundled LGPL-3 text is then a valid choice); human to confirm before first release |

## 13. Working rules for Claude Code

- Before each phase, re-read the relevant section of this spec and the matching R source in
  `reference/r-dada2/R/`. The R code is the behavioural spec (defaults, edge cases, messages).
- Prefer small, verifiable steps: write/port a parity test first, then make it pass.
- Never weaken a tolerance or skip a parity test to get green; investigate and document.
- Don't add runtime dependencies beyond §4 without noting why in CLAUDE.md.
- When the spec and the R behaviour disagree, follow R behaviour and note it.
- If something is genuinely ambiguous (API shape, licensing), stop and ask the human.
- Git: see §1.1. No commits, no pushes, no PRs, nothing sent to `benjjneb/dada2`.
- Never start a GATED phase (11, 12) without an explicit instruction.
