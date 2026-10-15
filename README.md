# pyDADA2

Exact amplicon sequence variant (ASV) inference for Python: an **unofficial** port of the
[DADA2](https://github.com/benjjneb/dada2) R package (version 1.36). It reuses DADA2's C++
core and moves the performance-critical R-level code to C++ too. All of this is exposed
through a Pythonic, scikit-learn-style API.

* Home: <https://github.com/raeslab/pyDADA2>
* Upstream: <https://github.com/benjjneb/dada2> (the algorithms and their C++ implementation)
* License: LGPL (see `LICENSE`), as the upstream code it derives from

Results are validated against R dada2 1.36 on the bundled example data. Filtering, dereplication,
denoising, merging, sequence tables and chimera removal give identical results. Known
differences are listed in [`DEVIATIONS.md`](DEVIATIONS.md); for example, taxonomy bootstrap
values differ by random-number noise only.

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
