# R reference fixtures

`data/` holds fixtures produced by R dada2 1.36. The parity tests compare pydada2 against
them and never need R at test time.

## Regenerating

```bash
# Option A (until spec phase 12): the bundled reference copy
R CMD INSTALL --no-docs --library="$R_LIBS_USER" reference/r-dada2
# Option B (always works): Bioconductor
#   install the Bioconductor release that ships dada2 1.36.x (the script checks the version)
Rscript tests/reference/generate_reference.R
```

`MANIFEST.json` records the dada2/R/package versions used. Regenerating should only change
`MANIFEST.json` (timestamp) and `taxonomy.json.gz` (upstream tie-breaking uses
`std::random_device`; see `DEVIATIONS.md`). Any other change means the R side changed and
must be investigated.

The format conventions are documented in `CLAUDE.md` ("Fixtures").
