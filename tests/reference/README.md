# R reference fixtures

`data/` holds fixtures produced by R dada2 1.36. The parity tests compare pydada2 against
them and never need R at test time.

## Regenerating

Needs R with Rcpp, RcppParallel, Biostrings, ShortRead, jsonlite, ggplot2 and reshape2 (the dev
container provides them when built with `INSTALL_R=true`).

```bash
scripts/install_r_dada2.sh                     # dada2 at upstream b61fe2e (1.36.0) into $R_LIBS_USER
Rscript tests/reference/generate_reference.R   # ~3.5 min
```

A Bioconductor install of dada2 1.36.x works too (the script checks the version).

`MANIFEST.json` records the dada2/R/package versions used. Regenerating should only change
`MANIFEST.json` (timestamp) and `taxonomy.json.gz` (upstream tie-breaking uses
`std::random_device`; see `docs/porting-notes.md`). Any other change means the R side changed
and must be investigated: `git checkout -- tests/reference/data/taxonomy.json.gz
tests/reference/data/MANIFEST.json` and review the rest of the diff.

## Format

- JSON written by `jsonlite` with `digits=I(17)` (exact double round-trip) and
  `auto_unbox=FALSE`: **every R value is a JSON array, scalars too**.
- `NA` → `null`, also for NaN/±Inf in doubles (compare with `rparity.assert_close_r`).
- Matrices are row-major `{rownames, colnames, values}`; data.frames are `{column: [...]}`.
- Files over 100 kB are `.json.gz`.
- R indices are stored 1-based as-is; the Python tests convert.
- Dereplication qualities are stored as integer sums (`qual_sums`): mean = sum / abundance.
- FASTQ outputs are stored as `ids`, `widths` and the md5 of the canonical text
  `"id\nseq\nqual\n"` per record (independent of file format and compression).
- `err/*.npy` are the bundled `errBalancedF/R` and `tperr1` (16 × 41, float64, C order).
- `inputs/` holds small synthetic FASTQ/FASTA inputs.
- `run_filter` records R's warnings; the parity tests require exactly those warnings.
- Loaders and comparison helpers: `tests/rparity.py`.
