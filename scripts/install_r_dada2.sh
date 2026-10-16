#!/usr/bin/env bash
# Install R dada2 at the upstream commit pyDADA2 was ported from (dada2 1.36.0) into $R_LIBS_USER.
# Only needed to regenerate the parity fixtures (tests/reference/) or to benchmark against R.
# Requires R with Rcpp, RcppParallel, Biostrings, ShortRead, ggplot2 and reshape2 installed.
set -euo pipefail
SHA=b61fe2e90115be08e78a01a16af89fd21402ca42
LIB="${R_LIBS_USER:?set R_LIBS_USER}"
mkdir -p "$LIB"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
curl -fsSL "https://github.com/benjjneb/dada2/archive/$SHA.tar.gz" -o "$WORK/dada2.tar.gz"
tar -xzf "$WORK/dada2.tar.gz" -C "$WORK"
R CMD INSTALL --no-docs --no-multiarch --library="$LIB" "$WORK/dada2-$SHA"
Rscript -e 'stopifnot(packageVersion("dada2") == "1.36.0")'
