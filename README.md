# pyDADA2

> [!WARNING]
> **Experimental: not yet ready for production use.** pyDADA2 reproduces R dada2 1.36 exactly on
> every dataset it has been validated on (Linux x86_64, see below), but it is new software that
> has not yet been tested by other users, on other platforms (macOS, Windows, ARM) or on real
> datasets of other amplicon types (ITS, PacBio full-length 16S, 18S). Check its results against
> R dada2 before relying on them, and please [report](https://github.com/raeslab/pyDADA2/issues)
> any difference you find.

Exact amplicon sequence variant (ASV) inference for Python: an **unofficial** port of the
[DADA2](https://github.com/benjjneb/dada2) R package (version 1.36). It reuses DADA2's C++
core and moves the performance-critical R-level code to C++ too. All of this is exposed
through a Pythonic, scikit-learn-style API.

* Home: <https://github.com/raeslab/pyDADA2>
* Upstream: <https://github.com/benjjneb/dada2> (the algorithms and their C++ implementation),
  ported at commit [`b61fe2e`](https://github.com/benjjneb/dada2/tree/b61fe2e90115be08e78a01a16af89fd21402ca42)
* License: LGPL-3.0-or-later, as the upstream code it derives from (see [License](#license))
* **All credit for the method and its implementation goes to the DADA2 authors.** If you use
  pyDADA2, please [cite the original DADA2 papers](#citation), not this port.

Results are validated against R dada2 1.36 on three datasets (see
[Benchmark](#benchmark-pydada2-vs-r-dada2)):
- the bundled example data;
- the MiSeq SOP dataset;
- the full 362-sample Stability dataset that MiSeq SOP comes from.

Filtering, dereplication, error learning, denoising, merging, sequence tables and chimera
removal give identical results on all three. Known differences are listed in
[`docs/porting-notes.md`](docs/porting-notes.md#3-known-differences-from-r-dada2-136); for example, taxonomy bootstrap values differ by random-number
noise only. pydada2 is **6–12x faster than R** on MiSeq SOP, and 4.5x on the larger dataset.

## Installation

pydada2 is not on PyPI yet: the name `pydada2` belongs to an unrelated project there (see
"Status"), so `pip install pydada2` does **not** install this package. Install it from source
or from the wheels built by CI. The plotting functions need matplotlib (the `plot` extra).
On Intel macOS, use scikit-misc 0.5.2 (`pip install "scikit-misc==0.5.2"`) unless you have a
Fortran compiler: later releases ship no wheels for that platform.

From source (needs a C++17 compiler and CMake; uses [uv](https://docs.astral.sh/uv/)):

```bash
git clone https://github.com/raeslab/pyDADA2 && cd pyDADA2
uv sync --all-groups
uv run pytest
```

Python ≥ 3.11 on Linux, macOS and Windows.

## Quickstart (paired-end Illumina)

```python
import pydada2 as dada2
from pydada2.datasets import example_path

fwd = [example_path("sam1F.fastq.gz"), example_path("sam2F.fastq.gz")]
rev = [example_path("sam1R.fastq.gz"), example_path("sam2R.fastq.gz")]
filt_f = ["filtered/sam1F.fq.gz", "filtered/sam2F.fq.gz"]
filt_r = ["filtered/sam1R.fq.gz", "filtered/sam2R.fq.gz"]

stats = dada2.filter_and_trim(fwd, filt_f, rev, filt_r, trunc_len=(240, 160), max_ee=(2, 2),
                              trunc_q=2, rm_phix=True, n_jobs=-1)

err_f = dada2.learn_errors(filt_f, n_jobs=-1)          # ErrorLearner; matrix in err_f.err_
err_r = dada2.learn_errors(filt_r, n_jobs=-1)

names = ["sam1", "sam2"]
drp_f = dict(zip(names, dada2.derep_fastq(filt_f)))
drp_r = dict(zip(names, dada2.derep_fastq(filt_r)))
dd_f = dada2.DADA(err=err_f.err_, n_jobs=-1).fit(drp_f)
dd_r = dada2.DADA(err=err_r.err_, n_jobs=-1).fit(drp_r)

mergers = dada2.merge_pairs(dd_f, list(drp_f.values()), dd_r, list(drp_r.values()))
seqtab = dada2.make_sequence_table(dict(zip(names, mergers)))
seqtab_nochim = dada2.BimeraRemover(method="consensus", n_jobs=-1).fit_transform(seqtab)

clf = dada2.TaxonomyClassifier(min_boot=50, random_state=0, n_jobs=-1)
clf.fit(example_path("example_train_set.fa.gz"))      # e.g. a DADA2-formatted SILVA file
taxa, boot = clf.predict(seqtab_nochim.columns, return_bootstrap=True)
```

Sequence tables are `pandas.DataFrame`s (samples × ASVs, int64). Indices are 0-based and `-1`
means NA in integer arrays. Estimators follow scikit-learn conventions: parameters go in
`__init__`, `fit` returns `self`, learned state ends in `_`, and they work with `sklearn.base.clone`.

## R → Python

| R | Python |
|---|---|
| `filterAndTrim`, `fastqFilter`, `fastqPairedFilter` | `filter_and_trim`, `fastq_filter`, `fastq_paired_filter` |
| `removePrimers` | `remove_primers` |
| `derepFastq` | `derep_fastq` → `Dereplication` |
| `learnErrors` | `learn_errors` / `ErrorLearner` |
| `loessErrfun`, `noqualErrfun`, `makeBinnedQualErrfun`, `PacBioErrfun` | `LoessErrorModel`, `NoQualErrorModel`, `BinnedQualErrorModel`, `PacBioErrorModel` |
| `getErrors`, `inflateErr` | `get_errors`, `inflate_err` |
| `dada`, `setDadaOpt`/`getDadaOpt` | `DADA(...).fit(...)` → `DadaResult`s; options in `DadaParams` |
| `mergePairs` | `merge_pairs` |
| `makeSequenceTable`, `mergeSequenceTables`, `collapseNoMismatch` | `make_sequence_table`, `merge_sequence_tables`, `collapse_no_mismatch` |
| `isBimera`, `isBimeraDenovo`, `isBimeraDenovoTable`, `removeBimeraDenovo` | `is_bimera`, `is_bimera_denovo`, `is_bimera_denovo_table`, `BimeraRemover` / `remove_bimera_denovo` |
| `isShiftDenovo` | `is_shift_denovo` |
| `assignTaxonomy`, `assignSpecies`, `addSpecies` | `TaxonomyClassifier` / `assign_taxonomy`, `SpeciesAssigner` / `assign_species`, `add_species` |
| `nwalign`, `nwhamming` | `pydada2.align.nwalign`, `pydada2.align.nwhamming` |
| `rc`, `seqComplexity`, `isPhiX` | `reverse_complement`, `seq_complexity`, `is_phix` |
| `getSequences`, `getUniques`, `uniquesToFasta` | `get_sequences`, `get_uniques`, `uniques_to_fasta` |
| `plotQualityProfile`, `plotErrors`, `plotComplexity` | `pydada2.plotting.plot_quality_profile`, `plot_errors`, `plot_complexity` |

## Benchmark: pydada2 vs R dada2

Both tools ran the full [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html)
pipeline with identical parameters, with SILVA 138.2 for taxonomy. Setup:

- **Machine:** Intel Core Ultra 7 265H (Linux), each run in a fresh process.
- **Versions:**
  - R 4.6.1 / dada2 1.36.0;
  - Python 3.13 / pydada2 0.1.0 at commit `5e6453a` ("now");
  - for comparison, the previous pydada2 at `dff04da` ("before"), before the fixes and
    optimisations that followed an independent review
    ([`docs/porting-notes.md`](docs/porting-notes.md#4-independent-review)).
- **Scripts:** every version ran the same pipeline scripts.
- **Multithreading:** R's `mergePairs` and `addSpecies` have no multithreading; pydada2 uses
  `n_jobs` there too.

### MiSeq SOP (20 samples, 16S V4, 2x250, 152,360 read pairs)

**Results are identical.** These outputs match R exactly:
- the read-tracking table;
- the denoised ASVs and abundances per sample;
- the merged pairs;
- the sequence table (232 non-chimeric ASVs, same order and counts);
- the chimera calls.

Learned error rates agree with R to a relative 1.4e-14. Every output, including taxonomy and
bootstraps at the same seed, is bit-identical between the previous and the current pydada2.

**Taxonomy.** pydada2 uses a different random-number generator for the bootstrap (by design,
see [`docs/porting-notes.md`](docs/porting-notes.md#3-known-differences-from-r-dada2-136)). 229 of 232 ASVs get the same assignment at every level.
The other 3 differ at a level where R's bootstrap was 45–48, right at the `min_boot=50` cutoff.

**pydada2 is 6.0x faster than R at 1 thread and 12.1x at 8 threads**, and uses 40–50 % less
memory (median of 3 runs; full report:
[`benchmarks/dataset/REPORT.md`](benchmarks/dataset/REPORT.md)):

| stage | R, 1 thr | before, 1 thr | **now, 1 thr** | R, 8 thr | before, 8 thr | **now, 8 thr** |
|---|---|---|---|---|---|---|
| filter_and_trim | 17.7 s | 6.7 s | **3.1 s** | 5.3 s | 1.5 s | **1.2 s** |
| learn_errors (F + R) | 104.7 s | 66.3 s | **16.6 s** | 30.2 s | 20.8 s | **3.6 s** |
| derep_fastq (F + R) | 3.3 s | 0.82 s | **0.77 s** | 3.5 s | 0.80 s | **0.77 s** |
| dada (F + R) | 17.5 s | 11.6 s | **9.6 s** | 4.2 s | 3.5 s | **1.6 s** |
| merge_pairs | 0.98 s | 0.65 s | **0.30 s** | 1.0 s | 0.18 s | **0.11 s** |
| make_sequence_table + remove_bimera_denovo | 0.28 s | 0.19 s | **0.15 s** | 0.05 s | 0.03 s | **0.03 s** |
| assign_taxonomy (SILVA) | 83.3 s | 57.9 s | **6.5 s** | 73.8 s | 47.3 s | **2.6 s** |
| add_species (SILVA) | 21.1 s | 4.8 s | **4.3 s** | 21.0 s | 2.1 s | **1.6 s** |
| **total** | **4.1 min** | **2.5 min** | **41.5 s** | **2.3 min** | **76.0 s** | **11.6 s** |
| peak memory (RSS) | 4.4 GiB | 2.9 GiB | **2.3 GiB** | 4.4 GiB | 3.0 GiB | **2.7 GiB** |

**Against the previous pydada2**, the current version is 3.6x faster at 1 thread and 6.6x at 8
threads, and no stage is slower. The main changes:

- **Taxonomy (8.9–18x):** a sparse, parallel build of the SILVA reference model with a k-mer-major
  table.
- **Error learning (4–6x):** alignments are reused across the self-consistency rounds,
  samples are denoised concurrently, and threads come from a persistent pool.
- **Filtering:** a faster phiX screen.
- **Merging:** pairs are batched and de-duplicated.

To reproduce (needs R with dada2 1.36; downloads about 250 MB into `benchmarks/dataset/data/`):

```bash
benchmarks/dataset/download.sh                       # MiSeq SOP + SILVA 138.2, SHA-256 checked
uv run python benchmarks/dataset/run.py 3 1 8        # 3 repeats x {1, 8} threads x {R, pydada2}
uv run python benchmarks/dataset/report.py           # compare outputs, write REPORT.md
# optional: an earlier pydada2 in another venv, then add its results dir as the 4th argument
uv run python benchmarks/dataset/run.py --tool py --python <venv>/bin/python --results <dir> 3 1 8
```

### Larger dataset: the full Stability data set (362 samples, 3.65M read pairs)

MiSeq SOP is a 20-sample subset of the mouse "Stability" study (Kozich et al. 2013). The full
data set is 24x larger, with the same chemistry. The same pipeline was used at 8 threads, with
2 runs per tool. Full report:
[`benchmarks/dataset/REPORT_stability.md`](benchmarks/dataset/REPORT_stability.md).

**Results are identical.** These outputs match R exactly:
- the read-tracking table, from 3.65M input pairs down to 2.97M non-chimeric reads;
- the dereplicated uniques;
- the denoised ASVs per sample, including their order;
- the merged pairs;
- the sequence table (752 ASVs);
- the chimera calls.

Error rates agree to a relative 9.2e-15.

**What changed since the previous version.** It created the clusters of 8 of the 724 denoising
runs in a different order than R (same ASVs and abundances). The cause was a Poisson tail that
underflowed to 0 where R returns a tiny subnormal number (deviation #1 in
[`docs/porting-notes.md`](docs/porting-notes.md#3-known-differences-from-r-dada2-136)). That
is fixed, so these runs now match R too.

**Taxonomy:** 412 of 420 ASVs are identical at every level. The other 8 differ at genus only,
where R's bootstrap was 41–52.

**pydada2 is 4.5x faster than R** (median of 2 runs, 8 threads):

| stage | R | before | **now** | now vs R |
|---|---|---|---|---|
| filter_and_trim | 53.3 s | 28.2 s | **20.4 s** | 2.6x |
| learn_errors (F + R) | 129.4 s | 92.6 s | **16.1 s** | 8.0x |
| derep_fastq (F + R) | 57.4 s | 18.8 s | **19.3 s** | 3.0x |
| dada (F + R) | 107.3 s | 84.2 s | **41.1 s** | 2.6x |
| merge_pairs | 21.4 s | 3.6 s | **1.5 s** | 13.9x |
| make_sequence_table + remove_bimera_denovo | 0.26 s | 0.22 s | **0.13 s** | 2.0x |
| assign_taxonomy (SILVA) | 75.7 s | 50.2 s | **3.0 s** | 25.2x |
| add_species (SILVA) | 24.3 s | 2.2 s | **1.8 s** | 13.8x |
| **total** | **7.8 min** | **4.7 min** | **103 s** | **4.5x** |
| peak memory (RSS) | 6.9 GiB | 5.9 GiB | **6.5 GiB** | 1.1x |

**Against the previous version**, this is 2.7x faster overall.

- **Dereplication** is unchanged code and runs at the same speed: the per-run ranges overlap.
- **Peak memory** is 0.6 GiB higher. That is the price of the alignment cache (bounded at
  1 GiB per fit) and of denoising samples concurrently.

To reproduce (downloads another 845 MB):

```bash
benchmarks/dataset/download.sh --stability
uv run python benchmarks/dataset/run.py --fastq-dir Stability --results benchmarks/dataset/results/stability 2 8
uv run python benchmarks/dataset/report.py benchmarks/dataset/results/stability \
    benchmarks/dataset/REPORT_stability.md Stability
```

### Micro-benchmarks

Per-function micro-benchmarks on the small bundled example data are in
[`benchmarks/micro/REPORT.md`](benchmarks/micro/REPORT.md).

- **Against R:** no stage is slower than R by more than 1.5x. `make_sequence_table` is level with
  R (~0.1 ms per call); every other stage is faster.
- **Against the previous pydada2:** stages are up to 10x faster. Two are marginally slower:
  - `seq_complexity` (~8 %), because it now sums exactly like R's long-double `sum()`;
  - `assign_species` (+0.05 ms per call), from a format check on the 14-sequence example
    reference. On real SILVA data that stage is faster.
  All other stages are faster or unchanged.
- **Fair comparison:** both pydada2 versions ran from the same kind of storage. In the dev
  container, the `/workspaces` bind mount adds ~1 ms per file operation.

## Status

Experimental (see the warning at the top). How the port was built, how it was validated against
R, every known difference from R and the outcome of an independent review are summarised in
[`docs/porting-notes.md`](docs/porting-notes.md).

Open items:
- CI (Linux, macOS, Windows; wheels) has not run yet, so parity on platforms other than Linux
  x86_64 is unverified;
- the distribution name `pydada2` is already used on PyPI by another project, so the package
  will be renamed before it is published there.

## Development

```bash
uv sync --all-groups                         # env + editable build (C++ rebuilds on import)
uv run pytest -n auto --cov-fail-under=90    # full suite incl. R parity (no R needed)
uv run ruff format . && uv run ruff check . && uv run mypy src
uv run pre-commit install                    # optional: run the same checks on commit/push
```

The dev container (`.devcontainer/`) provides the C++ toolchain and Python. R is only needed to
regenerate the parity fixtures or to benchmark against R; see
[`tests/reference/README.md`](tests/reference/README.md).

## Citation

**Please cite the original work, not this port.** The method, its C++ implementation and the
pipeline that pyDADA2 reproduces were developed by Benjamin Callahan, Paul McMurdie, Susan Holmes
and the DADA2 contributors ([benjjneb/dada2](https://github.com/benjjneb/dada2)). pyDADA2 only
translates their R package to Python, to fit Python-based workflows. The port was made with
extensive AI assistance (Claude Code) and validated against the R package (see
[`docs/porting-notes.md`](docs/porting-notes.md)).

If you use pyDADA2, cite the DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> DADA2: High-resolution sample inference from Illumina amplicon data.
> *Nature Methods*, 13, 581–583. doi:[10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)

Depending on what you use, also cite:

- **PacBio full-length amplicons** (`PacBioErrorModel`, PacBio parameters):
  Callahan BJ, Wong J, Heiner C, Oh S, Theriot CM, Gulati AS, McGill SK, Dougherty MK (2019).
  High-throughput amplicon sequencing of the full-length 16S rRNA gene with single-nucleotide
  resolution. *Nucleic Acids Research*, 47(18), e103.
  doi:[10.1093/nar/gkz569](https://doi.org/10.1093/nar/gkz569)
- **Taxonomy** (`TaxonomyClassifier` / `assign_taxonomy`, DADA2's implementation of the RDP
  naive Bayesian classifier): Wang Q, Garrity GM, Tiedje JM, Cole JR (2007). Naive Bayesian
  classifier for rapid assignment of rRNA sequences into the new bacterial taxonomy. *Applied
  and Environmental Microbiology*, 73(16), 5261–5267.
  doi:[10.1128/AEM.00062-07](https://doi.org/10.1128/AEM.00062-07)
- **The reference database** you classify against (e.g. SILVA, UNITE, GTDB), as its
  maintainers request.

Machine-readable citation metadata is in [`CITATION.cff`](CITATION.cff); GitHub's "Cite this
repository" button gives the DADA2 paper.

## License

pyDADA2 continues DADA2's code, which is licensed under the LGPL. It is distributed under the GNU
Lesser General Public License, version 3 or (at your option) any later version: see
[`LICENSE`](LICENSE) and the GNU GPL v3 it builds on, [`COPYING`](COPYING). Upstream DADA2
copyright notices are kept in every file derived from it. The compiled extension also contains
nanobind, robin-map and Boost.Math, under their own permissive licenses
([`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md)).
