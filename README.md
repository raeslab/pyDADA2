# pyDADA2

Exact amplicon sequence variant (ASV) inference for Python: an **unofficial** port of the
[DADA2](https://github.com/benjjneb/dada2) R package (version 1.36). It reuses DADA2's C++
core and moves the performance-critical R-level code to C++ too. All of this is exposed
through a Pythonic, scikit-learn-style API.

* Home: <https://github.com/raeslab/pyDADA2>
* Upstream: <https://github.com/benjjneb/dada2> (the algorithms and their C++ implementation)
* License: LGPL (see `LICENSE`), as the upstream code it derives from

Results are validated against R dada2 1.36 on the bundled example data, on the MiSeq SOP dataset
and on the full 362-sample Stability dataset it comes from (see
[Benchmark](#benchmark-pydada2-vs-r-dada2)). Filtering, dereplication, denoising, merging,
sequence tables and chimera removal give identical results on MiSeq SOP. On the larger dataset,
every table is identical, but in 8 of 724 denoising runs the clusters come out in a different
order (same ASVs and abundances; see [`REVIEW.md`](REVIEW.md) §1.1). Known differences are listed in
[`DEVIATIONS.md`](DEVIATIONS.md); for example, taxonomy bootstrap values differ by
random-number noise only.

## Installation

```bash
pip install pydada2            # placeholder name, see "Status" below
pip install "pydada2[plot]"    # with matplotlib for the plotting functions
```

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
pipeline with identical parameters on the MiSeq SOP dataset (20 samples, 16S V4, 2x250,
152,360 read pairs), with SILVA 138.2 for taxonomy. Each configuration ran 3 times in a fresh
process on an Intel Core Ultra 7 265H (Linux), with R 4.6.1 / dada2 1.36.0 and Python 3.13 /
pydada2 0.1.0. Full report: [`benchmarks/dataset/REPORT.md`](benchmarks/dataset/REPORT.md).

**Results are identical.** These outputs match R exactly: the read-tracking table (input,
filtered, denoised F/R, merged, non-chimeric reads per sample), the denoised ASVs and their
abundances, the merged pairs, the sequence table (232 non-chimeric ASVs, same order and
counts) and the chimera calls. Learned error rates agree to a relative 1.4e-14. Every run of
either tool, at 1 or 8 threads, gives the same output.

Taxonomy uses a different random-number generator for the bootstrap (by design, see
[`DEVIATIONS.md`](DEVIATIONS.md)). 229 of 232 ASVs get the same assignment at every level,
including species. The other 3 differ at genus or family level, where R's bootstrap was 45–48,
right at the `min_boot=50` cutoff. Bootstrap confidences differ by a median of 0 and a
maximum of 14, all within binomial noise.

**pydada2 is 1.7–1.8x faster and uses 35 % less memory** (median wall time per stage):

| stage | R, 1 thread | pydada2, 1 thread | R, 8 threads | pydada2, 8 threads |
|---|---|---|---|---|
| filter_and_trim | 17.8 s | 7.0 s | 5.6 s | 1.6 s |
| learn_errors (F + R) | 105.0 s | 68.3 s | 31.0 s | 21.0 s |
| derep_fastq (F + R) | 3.3 s | 0.8 s | 3.4 s | 0.8 s |
| dada (F + R) | 17.5 s | 11.8 s | 4.3 s | 3.5 s |
| merge_pairs | 1.0 s | 0.7 s | 1.0 s | 0.7 s |
| remove_bimera_denovo | 0.27 s | 0.18 s | 0.04 s | 0.03 s |
| assign_taxonomy (SILVA) | 84.1 s | 58.4 s | 75.8 s | 48.6 s |
| add_species (SILVA) | 21.1 s | 4.9 s | 20.9 s | 2.2 s |
| **total** | **4.2 min** | **2.5 min** | **2.4 min** | **78 s** |
| peak memory (RSS) | 4.4 GiB | 2.9 GiB | 4.4 GiB | 3.0 GiB |

To reproduce (needs R with dada2 1.36; downloads about 250 MB into `benchmarks/dataset/data/`):

```bash
benchmarks/dataset/download.sh                       # MiSeq SOP + SILVA 138.2, SHA-256 checked
uv run python benchmarks/dataset/run.py 3 1 8        # 3 repeats x {1, 8} threads x {R, pydada2}
uv run python benchmarks/dataset/report.py           # compare outputs, write REPORT.md
```

### Larger dataset: the full Stability data set (362 samples)

MiSeq SOP is a 20-sample subset of the mouse "Stability" study (Kozich et al. 2013). The full
data set has 362 samples and 3.65 million read pairs (24x more), with the same chemistry
(16S V4, 2x250). The same pipeline and parameters were used, at 8 threads, with 2 runs per tool.
Full report: [`benchmarks/dataset/REPORT_stability.md`](benchmarks/dataset/REPORT_stability.md).

**Results:**

- **Identical to R:**
  - the read-tracking table, from 3.65M input pairs down to 2.97M non-chimeric reads;
  - the dereplicated uniques;
  - the sequence table (752 ASVs: same columns, order and counts);
  - all chimera calls.
- **Error rates** agree to a relative 9.2e-15.
- **Denoising:** every sample has exactly the same ASVs and abundances as in R. In 8 of 724
  denoising runs (362 samples x forward/reverse), the clusters come out in a different
  order. As a result, the merged pairs of those 8 samples carry different forward/reverse
  cluster indices. Their sequences and abundances are unchanged, so the sequence table is not
  affected. The cause is a Poisson tail that underflows to 0 where R returns a tiny subnormal
  number; [`REVIEW.md`](REVIEW.md) §1.1 has the analysis and a fix.
- **Taxonomy:** 412 of 420 ASVs are identical at every level. The other 8 differ at genus only,
  where R's bootstrap was 41–52, at the `min_boot=50` cutoff. Bootstraps differ by a median of
  0 (maximum 23).

**Speed: pydada2 is 1.6x faster and uses 15 % less memory** (median of 2 runs, 8 threads):

| stage | R | pydada2 | speed-up |
|---|---|---|---|
| filter_and_trim | 56.0 s | 28.5 s | 2.0x |
| learn_errors (F + R) | 133.8 s | 93.3 s | 1.4x |
| derep_fastq (F + R) | 58.3 s | 19.1 s | 3.1x |
| dada (F + R) | 110.2 s | 85.6 s | 1.3x |
| merge_pairs | 21.6 s | 13.9 s | 1.5x |
| make_sequence_table + remove_bimera_denovo | 0.27 s | 0.23 s | 1.2x |
| assign_taxonomy (SILVA) | 76.8 s | 51.7 s | 1.5x |
| add_species (SILVA) | 25.0 s | 2.2 s | 11.4x |
| **total** | **8.0 min** | **4.9 min** | **1.6x** |
| peak memory (RSS) | 6.9 GiB | 5.9 GiB | 1.2x |

What scales with the data set:

- **Speed-up:** the speed-up shrinks from 1.8x (MiSeq SOP) to 1.6x.
- **Denoising and error learning:** these dominate on larger data, and both tools run the same
  per-sample, alignment-bound C++ algorithm, so pydada2 gains only 1.3–1.4x there.
- **Taxonomy:** the time is nearly constant (~1 min in R, ~50 s in pydada2), because building
  the SILVA reference model is a fixed cost.

[`REVIEW.md`](REVIEW.md) §5 identifies output-preserving changes to all of these, measured at
11x faster than R on MiSeq SOP.

To reproduce (downloads another 845 MB):

```bash
benchmarks/dataset/download.sh --stability
uv run python benchmarks/dataset/run.py --fastq-dir Stability --results benchmarks/dataset/results/stability 2 8
uv run python benchmarks/dataset/report.py benchmarks/dataset/results/stability \
    benchmarks/dataset/REPORT_stability.md Stability
```

Per-function micro-benchmarks on the small bundled example data are in
[`benchmarks/micro/REPORT.md`](benchmarks/micro/REPORT.md).

## Status

Work in progress, developed following [`spec.md`](spec.md). The distribution name `pydada2` is
already used on PyPI by another project, so the final name is still to be decided.

## Citation

If you use this software, please cite the DADA2 paper:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP (2016).
> DADA2: High-resolution sample inference from Illumina amplicon data.
> *Nature Methods*, 13, 581–583. doi:10.1038/nmeth.3869

and, for PacBio full-length amplicons, Callahan et al. (2019), *Nucleic Acids Research*, 47(18),
e103. doi:10.1093/nar/gkz569. Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).
