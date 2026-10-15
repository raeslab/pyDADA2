#!/usr/bin/env Rscript
# Generate parity fixtures for pydada2 from the reference R package (dada2 1.36).
#
# Usage:  Rscript tests/reference/generate_reference.R
#
# Works against any installed dada2 1.36.x: the local copy in reference/r-dada2
# (installed by scripts/devcontainer-post-create.sh) or a Bioconductor install.
# Output: tests/reference/data/*.json[.gz] + inputs/ + err/*.npy + MANIFEST.json.
#
# Conventions (mirrored by tests/reference/fixtures.py):
#   * JSON written with digits=I(17) (exact double round-trip), auto_unbox=FALSE
#     (every R vector is a JSON array, scalars included), NA -> null,
#     matrices row-major, data.frames as {column: [values]}.
#   * Files larger than ~100 kB are gzipped (.json.gz).
#   * R indices (1-based) are stored as-is; the Python tests convert.
#   * FASTQ outputs are summarised as ids, widths and the md5 of the canonical
#     text "id\nseq\nqual\n" per record (independent of file formatting/compression).

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
  library(Biostrings)
  library(jsonlite)
})

stopifnot(packageVersion("dada2") >= "1.36.0", packageVersion("dada2") < "1.37.0")

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
OUT <- file.path(script_dir, "data")
INP <- file.path(OUT, "inputs")
ERR <- file.path(OUT, "err")
for (d in c(OUT, INP, ERR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
TMP <- tempfile("dada2fx")
dir.create(TMP)
on.exit(unlink(TMP, recursive = TRUE), add = TRUE)

ext <- function(f) system.file("extdata", f, package = "dada2")
log_msg <- function(...) message(">> ", ...)

# ---------------------------------------------------------------- writers ---
wj <- function(obj, name) {
  txt <- toJSON(obj, digits = I(17), auto_unbox = FALSE, na = "null", null = "null",
                dataframe = "columns", matrix = "rowmajor")
  if (nchar(txt, type = "bytes") > 100000) {
    path <- file.path(OUT, paste0(name, ".json.gz"))
    con <- gzfile(path, "w")
  } else {
    path <- file.path(OUT, paste0(name, ".json"))
    con <- file(path, "w")
  }
  writeLines(txt, con)
  close(con)
  log_msg("wrote ", basename(path))
}

# Minimal .npy (v1.0) writer for float64 matrices, C order.
write_npy <- function(m, path) {
  m <- unname(as.matrix(m))
  storage.mode(m) <- "double"
  hdr <- sprintf("{'descr': '<f8', 'fortran_order': False, 'shape': (%d, %d), }", nrow(m), ncol(m))
  total <- 10 + nchar(hdr) + 1
  pad <- (64 - total %% 64) %% 64
  hdr <- paste0(hdr, strrep(" ", pad), "\n")
  con <- file(path, "wb")
  writeBin(as.raw(c(0x93)), con)
  writeChar("NUMPY", con, eos = NULL)
  writeBin(as.raw(c(1, 0)), con)
  writeBin(as.integer(nchar(hdr)), con, size = 2, endian = "little")
  writeChar(hdr, con, eos = NULL)
  writeBin(as.vector(t(m)), con, size = 8, endian = "little")
  close(con)
  log_msg("wrote ", basename(path))
}

# ------------------------------------------------------------- converters ---
unq_list <- function(u) list(sequence = names(u), abundance = as.integer(unname(u)))
mat_list <- function(m) {
  if (is.null(m)) return(NULL)
  list(rownames = rownames(m), colnames = colnames(m), values = unname(m))
}
df_list <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (cn in names(df)) if (is.factor(df[[cn]])) df[[cn]] <- as.character(df[[cn]])
  rownames(df) <- NULL
  df
}
dada_list <- function(d) {
  err_in <- d$err_in
  err_in <- if (is.list(err_in) && !is.matrix(err_in)) lapply(err_in, mat_list) else mat_list(err_in)
  list(
    denoised = unq_list(d$denoised),
    clustering = df_list(d$clustering),
    map = d$map,
    pval = d$pval,
    quality = mat_list(d$quality),
    birth_subs = df_list(d$birth_subs),
    trans = mat_list(d$trans),
    err_in = err_in,
    err_out = mat_list(d$err_out),
    opts = d$opts
  )
}
dada_multi <- function(dd) {
  if (is(dd, "dada")) dd <- list(dd)
  lapply(dd, dada_list)
}
derep_list <- function(dr) {
  # Store cumulative (integer) quality sums; Python divides by abundance exactly like R.
  qsum <- round(dr$quals * as.numeric(dr$uniques))
  storage.mode(qsum) <- "integer"
  list(uniques = unq_list(dr$uniques), qual_sums = unname(qsum), quals_ncol = ncol(dr$quals),
       map = dr$map)
}
fq_digest <- function(path) {
  if (!file.exists(path)) return(list(exists = FALSE))
  fq <- readFastq(path)
  ids <- as.character(id(fq))
  sq <- as.character(sread(fq))
  ql <- as.character(quality(quality(fq)))
  txt <- paste0(ids, "\n", sq, "\n", ql, "\n", collapse = "")
  f <- tempfile(tmpdir = TMP)
  con <- file(f, "wb")
  writeChar(txt, con, eos = NULL, useBytes = TRUE)
  close(con)
  list(exists = TRUE, n = length(fq), ids = ids, widths = width(fq),
       md5 = unname(tools::md5sum(f)))
}
write_fq <- function(fq, path) {
  if (file.exists(path)) file.remove(path)
  writeFastq(fq, path, mode = "w", compress = TRUE)
}

# ================================================================ MANIFEST ===
log_msg("manifest")
pkgs <- c("dada2", "Rcpp", "RcppParallel", "ShortRead", "Biostrings", "jsonlite", "BiocGenerics", "IRanges", "XVector")
manifest <- list(
  dada2_version = as.character(packageVersion("dada2")),
  upstream_repo = "https://github.com/benjjneb/dada2",
  upstream_commit = "b61fe2e",
  r_version = R.version.string,
  platform = R.version$platform,
  packages = setNames(lapply(pkgs, function(p) as.character(packageVersion(p))), pkgs),
  generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  seed = 100L
)
writeLines(toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), file.path(OUT, "MANIFEST.json"))
exports <- sort(getNamespaceExports("dada2"))
exports <- exports[!startsWith(exports, ".")]  # drop S4 method tables
writeLines(exports, file.path(OUT, "NAMESPACE_exports.txt"))

# ================================================== bundled error matrices ===
log_msg("bundled data")
data("errBalancedF", "errBalancedR", "tperr1", package = "dada2", envir = environment())
write_npy(errBalancedF, file.path(ERR, "errBalancedF.npy"))
write_npy(errBalancedR, file.path(ERR, "errBalancedR.npy"))
write_npy(tperr1, file.path(ERR, "tperr1.npy"))
wj(list(errBalancedF = mat_list(errBalancedF), errBalancedR = mat_list(errBalancedR),
        tperr1 = mat_list(tperr1)), "bundled_err")

# ================================================= synthetic input files ===
log_msg("synthetic inputs")
set.seed(100)
fq1F <- readFastq(ext("sam1F.fastq.gz"))
fq1R <- readFastq(ext("sam1R.fastq.gz"))
rc_fq <- function(fq) ShortReadQ(sread = reverseComplement(sread(fq)),
                                  quality = FastqQuality(reverse(quality(quality(fq)))),
                                  id = id(fq))
# Single-end orientation: reads 101..200 reverse-complemented.
mixed <- append(fq1F[1:100], rc_fq(fq1F[101:200]))
write_fq(mixed, file.path(INP, "orient_mixed.fastq.gz"))
# Paired orientation: F/R swapped for pairs 101..150.
pF <- append(append(fq1F[1:100], fq1R[101:150]), fq1F[151:200])
pR <- append(append(fq1R[1:100], fq1F[101:150]), fq1R[151:200])
write_fq(pF, file.path(INP, "orient_pair_F.fastq.gz"))
write_fq(pR, file.path(INP, "orient_pair_R.fastq.gz"))
# matchIDs: differential membership in F and R (same relative order).
dropF <- seq(3, 300, by = 11)
dropR <- seq(5, 300, by = 7)
write_fq(fq1F[setdiff(1:300, dropF)], file.path(INP, "matchids_F.fastq.gz"))
write_fq(fq1R[setdiff(1:300, dropR)], file.path(INP, "matchids_R.fastq.gz"))
inp <- function(f) file.path(INP, f)

# ============================================================= filtering ===
log_msg("filterAndTrim")
run_filter <- function(.name, .fwd, .rev = NULL, ...) {
  # Dot-prefixed formals: filterAndTrim's `n` must not partially match them.
  log_msg("  filter ", .name)
  odir <- file.path(TMP, paste0("filt_", .name))
  dir.create(odir)
  filt <- file.path(odir, sub("\\.gz$", "", basename(.fwd)))
  filt.rev <- if (is.null(.rev)) NULL else file.path(odir, paste0("R_", sub("\\.gz$", "", basename(.rev))))
  out <- suppressMessages(filterAndTrim(.fwd, filt, rev = .rev, filt.rev = filt.rev,
                                        compress = FALSE, multithread = FALSE, ...))
  res <- list(
    inputs = basename(.fwd), inputs_rev = if (is.null(.rev)) NULL else basename(.rev),
    params = list(...),
    reads_in = unname(out[, 1]), reads_out = unname(out[, 2]),
    outputs = lapply(filt, fq_digest),
    outputs_rev = if (is.null(.rev)) NULL else lapply(filt.rev, fq_digest)
  )
  list(res = res, filt = filt, filt.rev = filt.rev)
}
samF <- c(ext("sam1F.fastq.gz"), ext("sam2F.fastq.gz"))
samR <- c(ext("sam1R.fastq.gz"), ext("sam2R.fastq.gz"))

filters <- list()
f1 <- run_filter("F1", samF, truncLen = 240, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE)
filters$F1 <- f1$res
filters$F2 <- run_filter("F2", samF, trimLeft = 10, trimRight = 20, truncQ = 11, minLen = 150,
                         maxEE = 1, rm.phix = FALSE)$res
filters$F3 <- run_filter("F3", samF, maxLen = 250, minQ = 10, truncQ = 2, rm.lowcomplex = 10)$res
filters$F4 <- run_filter("F4", samF, truncQ = 0, maxN = 2, truncLen = 200, rm.phix = FALSE)$res
filters$F5_chunked <- run_filter("F5", samF, truncLen = 240, maxN = 0, maxEE = 2, truncQ = 2,
                                 rm.phix = TRUE, n = 300)$res
filters$F6_orient <- run_filter("F6", inp("orient_mixed.fastq.gz"), orient.fwd = "TACG",
                                truncQ = 2, rm.phix = FALSE)$res
filters$F7_orient_chunked <- run_filter("F7", inp("orient_mixed.fastq.gz"), orient.fwd = "TACG",
                                        truncQ = 2, rm.phix = FALSE, n = 64)$res
filters$F8_none_pass <- run_filter("F8", samF, truncLen = 400)$res
fpb <- run_filter("FPB", ext("samPB.fastq.gz"), minQ = 3, minLen = 1000, maxLen = 1600, maxN = 0,
                  rm.phix = FALSE, maxEE = 2)
filters$FPB <- fpb$res
p1 <- run_filter("P1", samF, samR, truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2), truncQ = 2,
                 rm.phix = TRUE)
filters$P1 <- p1$res
filters$P2_matchids <- run_filter("P2", inp("matchids_F.fastq.gz"), inp("matchids_R.fastq.gz"),
                                  matchIDs = TRUE, truncLen = c(240, 160), rm.phix = FALSE)$res
filters$P3_matchids_chunked <- run_filter("P3", inp("matchids_F.fastq.gz"), inp("matchids_R.fastq.gz"),
                                          matchIDs = TRUE, truncLen = c(240, 160), rm.phix = FALSE,
                                          n = 50)$res
filters$P4_mixed <- run_filter("P4", samF, samR, trimLeft = c(10, 0), maxEE = c(1, 3),
                               minLen = c(50, 50), rm.lowcomplex = c(0, 8), truncQ = c(2, 11),
                               rm.phix = c(TRUE, FALSE))$res
filters$P5_orient <- run_filter("P5", inp("orient_pair_F.fastq.gz"), inp("orient_pair_R.fastq.gz"),
                                orient.fwd = "TACG", rm.phix = FALSE)$res
wj(filters, "filter_and_trim")
filtF <- f1$filt
filtF_P <- p1$filt
filtR_P <- p1$filt.rev

# ========================================================== removePrimers ===
log_msg("removePrimers")
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
run_primers <- function(.name, ...) {
  fout <- file.path(TMP, paste0("prim_", .name, ".fastq"))
  out <- suppressMessages(removePrimers(ext("samPBprimers.fastq.gz"), fout, compress = FALSE, ...))
  list(params = list(...), reads_in = unname(out[, 1]), reads_out = unname(out[, 2]),
       output = fq_digest(fout))
}
primers <- list(
  F27 = F27, R1492_rc = rc(R1492),
  default = run_primers("default", primer.fwd = F27, primer.rev = rc(R1492), orient = TRUE),
  indels = run_primers("indels", primer.fwd = F27, primer.rev = rc(R1492), orient = TRUE,
                       allow.indels = TRUE),
  fwd_only = run_primers("fwdonly", primer.fwd = F27, orient = TRUE),
  no_trim = run_primers("notrim", primer.fwd = F27, primer.rev = rc(R1492), trim.fwd = FALSE,
                        trim.rev = FALSE, orient = FALSE),
  mm0 = run_primers("mm0", primer.fwd = F27, primer.rev = rc(R1492), max.mismatch = 0)
)
wj(primers, "remove_primers")

# ====================================================== Biostrings matching ===
log_msg("matchPattern panel")
set.seed(100)
rand_dna <- function(n, alphabet = c("A", "C", "G", "T")) paste(sample(alphabet, n, replace = TRUE), collapse = "")
mutate <- function(s, nsub, nins, ndel) {
  x <- strsplit(s, "")[[1]]
  for (k in seq_len(nsub)) { i <- sample(length(x), 1); x[i] <- sample(setdiff(c("A","C","G","T"), x[i]), 1) }
  for (k in seq_len(nins)) { i <- sample(length(x) + 1, 1) - 1; x <- append(x, sample(c("A","C","G","T"), 1), after = i) }
  for (k in seq_len(ndel)) { if (length(x) > 3) x <- x[-sample(length(x), 1)] }
  paste(x, collapse = "")
}
mp_cases <- list()
for (case in 1:300) {
  plen <- sample(6:14, 1)
  pat <- rand_dna(plen)
  iupac <- case %% 3 == 0
  if (iupac) {  # put 1-3 ambiguity codes in the pattern
    x <- strsplit(pat, "")[[1]]
    x[sample(plen, sample(1:3, 1))] <- sample(c("R","Y","S","W","K","M","B","D","H","V","N"), 1)
    pat <- paste(x, collapse = "")
  }
  core <- mutate(gsub("[^ACGT]", "A", pat), sample(0:2, 1), sample(0:1, 1), sample(0:1, 1))
  subj <- paste0(rand_dna(sample(0:12, 1)), core, rand_dna(sample(0:12, 1)),
                 if (case %% 5 == 0) "NN" else "", if (case %% 7 == 0) core else "")
  if (case %% 11 == 0) subj <- paste0(substr(core, 3, nchar(core)), rand_dna(10))  # partial at left edge
  k <- sample(0:3, 1)
  fixed <- !iupac
  res <- list(subject = subj, pattern = pat, max_mismatch = k, fixed = fixed)
  for (ind in c(FALSE, TRUE)) {
    m <- tryCatch(matchPattern(pat, DNAString(subj), max.mismatch = k, with.indels = ind, fixed = fixed),
                  error = function(e) NULL)
    res[[if (ind) "indels" else "mismatch"]] <- if (is.null(m)) NULL else list(start = start(m), end = end(m))
  }
  mp_cases[[case]] <- res
}
wj(mp_cases, "match_pattern")

# =========================================================== dereplication ===
log_msg("derepFastq")
dF <- derepFastq(filtF)          # list of 2
dR <- derepFastq(filtR_P)
dF_P <- derepFastq(filtF_P)
names(dF) <- names(dR) <- names(dF_P) <- c("sam1", "sam2")
dPB <- derepFastq(fpb$filt)
derep <- list(
  filtF_sam1 = derep_list(dF$sam1), filtF_sam2 = derep_list(dF$sam2),
  filtF_sam1_n200 = derep_list(derepFastq(filtF[[1]], n = 200)),
  filtR_P_sam1 = derep_list(dR$sam1), filtR_P_sam2 = derep_list(dR$sam2),
  filtF_P_sam1 = derep_list(dF_P$sam1), filtF_P_sam2 = derep_list(dF_P$sam2),
  raw_sam1F = derep_list(derepFastq(ext("sam1F.fastq.gz"))),
  raw_samPB = derep_list(derepFastq(ext("samPB.fastq.gz"))),
  raw_samPB_n100 = derep_list(derepFastq(ext("samPB.fastq.gz"), n = 100)),
  filt_samPB = derep_list(dPB)
)
wj(derep, "derep")

# ========================================================= low-level panel ===
log_msg("alignment panel")
d0 <- suppressMessages(dada(dF$sam1, err = errBalancedF, multithread = FALSE, verbose = FALSE))
asv <- getSequences(d0)[1:8]
s1 <- asv[[1]]
mut <- function(s, pos, ch) paste0(substr(s, 1, pos - 1), ch, substr(s, pos + 1, nchar(s)))
panel <- c(asv,
           del = paste0(substr(s1, 1, 99), substr(s1, 101, 240)),
           ins = paste0(substr(s1, 1, 120), "G", substr(s1, 121, 240)),
           mm2 = mut(mut(s1, 50, ifelse(substr(s1, 50, 50) == "A", "C", "A")), 180,
                     ifelse(substr(s1, 180, 180) == "G", "T", "G")),
           shift = paste0(substr(s1, 6, 240), "ACGTA"),
           short = substr(s1, 30, 160),
           homo = sub("AAAA", "AAAAA", s1),
           short2 = "ACGTACGTTTGACCA",
           s_ex1 = "CTAATACATGCAAGTCGAGCGAGTCTGCCTTGAAGATCGGAGTGCTTGCACTCTGTGAAACAAGATA",
           s_ex2 = "TTAACACATGCAAGTCGAACGGAAAGGCCAGTGCTTGCACTGGTACTCGAGTGGCGAACGGGTGAGT")
panel <- unname(panel)
np <- length(panel)
pairs <- which(upper.tri(diag(np)), arr.ind = TRUE)
pairs <- rbind(pairs, cbind(pairs[1:10, 2], pairs[1:10, 1]))  # some reversed pairs
nw_variants <- list(
  default = list(),
  band16 = list(band = 16),
  noends = list(endsfree = FALSE),
  homo = list(homo_gap = -1),
  vec = list(vec = TRUE),
  vec_band16 = list(vec = TRUE, band = 16),
  merge0 = list(match = 1, mismatch = -64, gap = -64),
  merge8 = list(match = 1, mismatch = -8, gap = -8)
)
align <- list(panel = panel, pairs = pairs, variants = nw_variants, results = list())
for (vn in names(nw_variants)) {
  als <- lapply(seq_len(nrow(pairs)), function(k) {
    do.call(nwalign, c(list(panel[pairs[k, 1]], panel[pairs[k, 2]]), nw_variants[[vn]]))
  })
  ev <- t(sapply(als, function(a) dada2:::C_eval_pair(a[1], a[2])))
  align$results[[vn]] <- list(al1 = sapply(als, `[`, 1), al2 = sapply(als, `[`, 2),
                              match = ev[, 1], mismatch = ev[, 2], indel = ev[, 3])
}
align$nwhamming <- nwhamming(panel[pairs[, 1]], panel[pairs[, 2]])
align$nwhamming_band16 <- nwhamming(panel[pairs[, 1]], panel[pairs[, 2]], band = 16)
align$kmer <- list(
  kmer_dist5 = mapply(function(a, b) dada2:::kmer_dist(a, b, 5), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE),
  kord_dist5_sse0 = mapply(function(a, b) dada2:::kord_dist(a, b, 5, 0), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE),
  kord_dist5_sse1 = mapply(function(a, b) dada2:::kord_dist(a, b, 5, 1), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE),
  kord_dist5_sse2 = mapply(function(a, b) dada2:::kord_dist(a, b, 5, 2), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE),
  kmer_matches5 = mapply(function(a, b) dada2:::kmer_matches(a, b, 5), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE),
  kdist_matches5 = mapply(function(a, b) dada2:::kdist_matches(a, b, 5), panel[pairs[, 1]], panel[pairs[, 2]], USE.NAMES = FALSE)
)
# pair_consensus on merge-style alignments (F vs rc(R) of sample 1 ASVs)
dR0 <- suppressMessages(dada(dR$sam1, err = errBalancedR, multithread = FALSE, verbose = FALSE))
dF0p <- suppressMessages(dada(dF_P$sam1, err = errBalancedF, multithread = FALSE, verbose = FALSE))
fseq <- getSequences(dF0p)[1:5]
rseq <- rc(getSequences(dR0)[1:5])
pc <- list(fwd = fseq, rev_rc = rseq, results = list())
for (i in seq_along(fseq)) for (j in seq_along(rseq)) {
  a <- nwalign(fseq[[i]], rseq[[j]], match = 1, mismatch = -8, gap = -8, band = -1)
  pc$results[[length(pc$results) + 1]] <- list(
    i = i, j = j, al1 = a[1], al2 = a[2],
    eval = unname(dada2:::C_eval_pair(a[1], a[2])),
    p1 = dada2:::C_pair_consensus(a[1], a[2], 1L, FALSE),
    p2 = dada2:::C_pair_consensus(a[1], a[2], 2L, FALSE),
    p1_trim = dada2:::C_pair_consensus(a[1], a[2], 1L, TRUE),
    p2_trim = dada2:::C_pair_consensus(a[1], a[2], 2L, TRUE))
}
align$pair_consensus <- pc
align$is_acgt <- list(seqs = c(panel[1:3], "ACGN", "acgt", "", "ACGU", "RYACGT"),
                      result = dada2:::C_isACGT(c(panel[1:3], "ACGN", "acgt", "", "ACGU", "RYACGT")))
wj(align, "align")

# matrixEE: first 300 reads of sam1F, first 100 of samPB (variable length -> NA padding)
qmF <- as(quality(fq1F[1:300]), "matrix")
fqPB <- readFastq(ext("samPB.fastq.gz"))
qmPB <- as(quality(fqPB[1:100]), "matrix")
wj(list(sam1F_300 = dada2:::C_matrixEE(qmF), samPB_100 = dada2:::C_matrixEE(qmPB),
        samPB_100_ncol = ncol(qmPB)), "matrix_ee")

# ppois upper tail grid (used by the C++ abundance p-value)
log_msg("ppois")
lams <- c(1e-12, 1e-8, 1e-6, 1e-4, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 100, 500, 1000, 1e4, 1e5, 1e6)
ns <- c(0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 1e4, 5e4, 1e5, 1e6, 2e6)
grid <- expand.grid(n = ns, lambda = lams)
grid$p_upper <- ppois(grid$n, grid$lambda, lower.tail = FALSE)
wj(df_list(grid), "ppois")

# ============================================== sequences / complexity / phiX ===
log_msg("sequences")
phix <- getSequences(ext("phix_genome.fa"))[[1]]
phix_frag <- c(substr(phix, 1001, 1250), substr(phix, 3001, 3250),
               rc(substr(phix, 2001, 2250)), mut(substr(phix, 4001, 4250), 100, "A"),
               substr(phix, 5001, 5040))
raw200 <- as.character(sread(fq1F[1:200]))
cx_seqs <- c(raw200, strrep("A", 100), strrep("AC", 60), strrep("ACGT", 30),
             paste0(strrep("A", 80), raw200[[1]]), "ACGTNNNNACGT", "AC")
sequences <- list(
  phix_frag = phix_frag,
  is_phix = list(seqs = c(raw200[1:50], phix_frag),
                 default = isPhiX(c(raw200[1:50], phix_frag)),
                 ws12_mm1 = isPhiX(c(raw200[1:50], phix_frag), wordSize = 12, minMatches = 1),
                 match_ref = dada2:::C_matchRef(c(raw200[1:50], phix_frag), phix, 16, TRUE),
                 match_ref_overlap = dada2:::C_matchRef(c(raw200[1:50], phix_frag), phix, 16, FALSE)),
  complexity = list(seqs = cx_seqs,
                    k2 = suppressWarnings(seqComplexity(cx_seqs)),
                    k3 = suppressWarnings(seqComplexity(cx_seqs, kmerSize = 3)),
                    win50 = suppressWarnings(seqComplexity(cx_seqs, window = 50, by = 5))),
  rc = list(seqs = c(panel[1:3], "ACGTN", "RYKMSWBDHVN", "acgt", ""),
            result = rc(c(panel[1:3], "ACGTN", "RYKMSWBDHVN", "acgt", "")))
)
wj(sequences, "sequences")

# quality profile counts (ShortRead::qa, as used by plotQualityProfile)
log_msg("quality profile")
qprof <- list()
for (f in c("sam1F.fastq.gz", "samPB.fastq.gz")) {
  q <- qa(ext(f), n = 5e5)
  st <- plotQualityProfile(ext(f))$layers[[2]]$data  # per-cycle summary drawn by plotQualityProfile
  qprof[[f]] <- list(per_cycle = df_list(q[["perCycle"]]$quality[, c("Cycle", "Score", "Count")]),
                     read_count = unname(q[["readCounts"]]$read),
                     stats = df_list(st[, c("Cycle", "Mean", "Q25", "Q50", "Q75", "Cum")]))
}
wj(qprof, "quality_profile")

# =========================================================== error models ===
log_msg("error models")
trans <- d0$trans
binnedQ <- sort(unique(as.integer(colnames(trans)[colSums(trans) > 0])))
binnedQ <- unique(c(min(binnedQ), binnedQ[seq(1, length(binnedQ), by = 3)], max(binnedQ)))
errm <- list(
  trans = mat_list(trans),
  loess = mat_list(loessErrfun(trans)),
  noqual = mat_list(noqualErrfun(trans)),
  noqual_pc5 = mat_list(noqualErrfun(trans, pseudocount = 5)),
  binned_q = binnedQ,
  binned = mat_list(suppressWarnings(makeBinnedQualErrfun(binnedQ)(trans))),
  inflate3 = mat_list(inflateErr(errBalancedF, 3)),
  inflate3_self = mat_list(inflateErr(errBalancedF, 3, inflateSelfTransitions = TRUE)),
  get_errors = mat_list(getErrors(d0)),
  get_errors_detailed = list(err_out = mat_list(getErrors(d0, detailed = TRUE)$err_out),
                             trans = mat_list(getErrors(d0, detailed = TRUE)$trans))
)
wj(errm, "error_models")

log_msg("learnErrors")
set.seed(100)
leF <- learnErrors(filtF_P, multithread = FALSE, verbose = 0)
leR <- learnErrors(filtR_P, multithread = FALSE, verbose = 0)
le_noqual <- learnErrors(filtF_P, errorEstimationFunction = noqualErrfun, multithread = FALSE, verbose = 0)
le_small <- learnErrors(filtF_P, nbases = 1e5, multithread = FALSE, verbose = 0)
le_pb <- learnErrors(fpb$filt, errorEstimationFunction = PacBioErrfun, BAND_SIZE = 32,
                     multithread = FALSE, verbose = 0)
le_out <- function(le) list(err_out = mat_list(le$err_out),
                            err_in = if (is.list(le$err_in) && !is.matrix(le$err_in)) lapply(le$err_in, mat_list) else mat_list(le$err_in),
                            trans = mat_list(le$trans))
wj(list(fwd = le_out(leF), rev = le_out(leR), noqual = le_out(le_noqual),
        nbases_1e5 = le_out(le_small), pacbio = le_out(le_pb)), "learn_errors")

# ================================================================== dada ===
log_msg("dada")
run_dada <- function(drp, err, ...) suppressMessages(dada(drp, err = err, multithread = FALSE, verbose = FALSE, ...))
prior_seqs <- getSequences(d0)[c(3, 5)]
dd <- list()
ddF <- run_dada(dF_P, errBalancedF)
ddR <- run_dada(dR, errBalancedR)
dd$default_F <- dada_multi(ddF)
dd$default_R <- dada_multi(ddR)
dd$single_filtF <- dada_multi(run_dada(dF$sam1, errBalancedF))
dd$pool_true <- dada_multi(run_dada(dF_P, errBalancedF, pool = TRUE))
dd$pool_pseudo <- dada_multi(run_dada(dF_P, errBalancedF, pool = "pseudo"))
dd$no_quals <- dada_multi(run_dada(dF_P, errBalancedF, USE_QUALS = FALSE))
dd$not_gapless <- dada_multi(run_dada(dF_P, errBalancedF, GAPLESS = FALSE))
dd$band32 <- dada_multi(run_dada(dF_P, errBalancedF, BAND_SIZE = 32))
dd$priors <- dada_multi(run_dada(dF_P, errBalancedF, priors = prior_seqs))
dd$detect_singletons <- dada_multi(run_dada(dF_P, errBalancedF, DETECT_SINGLETONS = TRUE))
dd$omega_a_1e10 <- dada_multi(run_dada(dF_P, errBalancedF, OMEGA_A = 1e-10))
dd$no_kmers <- dada_multi(run_dada(dF_P, errBalancedF, USE_KMERS = FALSE))
dd$not_vectorized <- dada_multi(run_dada(dF_P, errBalancedF, VECTORIZED_ALIGNMENT = FALSE))
dd$homo_gap <- dada_multi(run_dada(dF_P, errBalancedF, HOMOPOLYMER_GAP_PENALTY = -1, VECTORIZED_ALIGNMENT = FALSE))
dd$sse0 <- dada_multi(run_dada(dF_P, errBalancedF, SSE = 0))
dd$sse1 <- dada_multi(run_dada(dF_P, errBalancedF, SSE = 1))
dd$not_greedy <- dada_multi(run_dada(dF_P, errBalancedF, GREEDY = FALSE))
dd$min_fold_hamming <- dada_multi(run_dada(dF_P, errBalancedF, MIN_FOLD = 2, MIN_HAMMING = 2, MIN_ABUNDANCE = 2))
dd$max_clust5 <- dada_multi(run_dada(dF_P, errBalancedF, MAX_CLUST = 5))
dd$self_consist <- dada_multi(run_dada(dF_P, errBalancedF, selfConsist = TRUE))
dd$pacbio <- dada_multi(run_dada(dPB, le_pb$err_out, BAND_SIZE = 32))
dd$prior_seqs <- prior_seqs
wj(dd, "dada")

# ============================================================ mergePairs ===
log_msg("mergePairs")
mm <- function(...) suppressMessages(mergePairs(ddF, dF_P, ddR, dR, ...))
mergers <- mm()
mg <- list(
  default = lapply(mergers, df_list),
  rejects = lapply(mm(returnRejects = TRUE), df_list),
  mismatch1 = lapply(mm(maxMismatch = 1, returnRejects = TRUE), df_list),
  min_overlap20 = lapply(mm(minOverlap = 20), df_list),
  concat = lapply(mm(justConcatenate = TRUE), df_list),
  trim_overhang = lapply(mm(trimOverhang = TRUE), df_list),
  propagate = lapply(mm(propagateCol = c("n0", "birth_fold")), df_list)
)
wj(mg, "merge_pairs")

# ===================================================== sequence tables ===
log_msg("sequence tables")
seqtab <- makeSequenceTable(mergers)
seqtabF <- makeSequenceTable(ddF)
# Variable-length table for collapseNoMismatch: add prefix / shifted / duplicate columns.
top <- colnames(seqtabF)[1:6]
extra <- c(substr(top[1], 1, 230), substr(top[2], 5, 240), paste0(top[3], "ACG"),
           substr(top[4], 1, 200), paste0("TT", top[5]))
ex_mat <- matrix(c(7L, 3L, 0L, 12L, 5L, 2L, 9L, 0L, 1L, 4L), nrow = 2, dimnames = list(rownames(seqtabF), extra))
st_var <- cbind(seqtabF, ex_mat)
st_dup <- cbind(seqtabF[, 1:5], seqtabF[, 1:2])
colnames(st_dup) <- c(colnames(seqtabF)[1:5], colnames(seqtabF)[1:2])
stB <- seqtab[, 1:min(6, ncol(seqtab))]
rownames(stB) <- c("sam3", "sam4")
stC <- seqtab[1, 3:ncol(seqtab), drop = FALSE]
stRC <- seqtab[, 1:min(6, ncol(seqtab))]
colnames(stRC)[1:3] <- rc(colnames(stRC)[1:3])
rownames(stRC) <- c("sam5", "sam6")
st <- list(
  seqtab = mat_list(seqtab),
  seqtabF = mat_list(seqtabF),
  seqtab_nsamples = mat_list(makeSequenceTable(mergers, orderBy = "nsamples")),
  seqtab_unordered = mat_list(makeSequenceTable(mergers, orderBy = NULL)),
  st_var = mat_list(st_var),
  collapse = mat_list(collapseNoMismatch(st_var)),
  collapse_ov100 = mat_list(collapseNoMismatch(st_var, minOverlap = 100)),
  collapse_novec = mat_list(collapseNoMismatch(st_var, vec = FALSE)),
  st_dup = mat_list(st_dup),
  collapse_identical = mat_list(suppressMessages(collapseNoMismatch(st_dup, identicalOnly = TRUE))),
  stB = mat_list(stB), stC = mat_list(stC), stRC = mat_list(stRC),
  merge_AB = mat_list(mergeSequenceTables(seqtab, stB)),
  merge_AC_sum = mat_list(mergeSequenceTables(seqtab, stC, repeats = "sum")),
  merge_AB_unordered = mat_list(mergeSequenceTables(seqtab, stB, orderBy = NULL)),
  merge_ARC_tryrc = mat_list(suppressMessages(mergeSequenceTables(seqtab, stRC, tryRC = TRUE)))
)
wj(st, "seqtab")

# ============================================================== chimeras ===
log_msg("chimeras")
tp <- colnames(seqtab)[1:6]
chims <- c(paste0(substr(tp[1], 1, 120), substr(tp[2], 121, nchar(tp[2]))),
           paste0(substr(tp[3], 1, 80), substr(tp[1], 81, nchar(tp[1]))),
           paste0(substr(tp[2], 1, 170), substr(tp[4], 171, nchar(tp[4]))),
           paste0(substr(tp[5], 1, 140), substr(tp[6], 141, nchar(tp[6]))))
chims <- unique(chims[!chims %in% colnames(seqtab)])
ch_mat <- matrix(c(3L, 2L, 4L, 0L, 2L, 2L, 5L, 1L)[seq_len(2 * length(chims))], nrow = 2,
                 dimnames = list(rownames(seqtab), chims))
st_chim <- cbind(seqtab, ch_mat)
# one-off chimera: chimera of tp1/tp2 with a single mismatch
oo <- chims[[1]]
oo <- mut(oo, 200, ifelse(substr(oo, 200, 200) == "A", "C", "A"))
unq_chim <- getUniques(setNames(as.integer(colSums(st_chim)), colnames(st_chim)))
bimtab <- dada2:::C_table_bimera2(st_chim, colnames(st_chim), 1.5, 2, FALSE, 4, 5L, -4L, -8L, 16)
chim <- list(
  chims = chims, st_chim = mat_list(st_chim),
  table_bimera = df_list(bimtab),
  is_bimera_denovo_table = unname(isBimeraDenovoTable(st_chim)),
  is_bimera_denovo_table_oneoff = unname(isBimeraDenovoTable(st_chim, allowOneOff = TRUE, minSampleFraction = 0.5)),
  remove_consensus = colnames(removeBimeraDenovo(st_chim, method = "consensus")),
  remove_pooled = colnames(removeBimeraDenovo(st_chim, method = "pooled")),
  remove_per_sample = mat_list(removeBimeraDenovo(st_chim, method = "per-sample")),
  is_bimera_denovo = unname(isBimeraDenovo(unq_chim)),
  is_bimera_denovo_oneoff = unname(isBimeraDenovo(unq_chim, allowOneOff = TRUE, minFoldParentOverAbundance = 1, minParentAbundance = 2)),
  is_bimera_denovo_dada = unname(isBimeraDenovo(ddF$sam1)),
  is_bimera = list(
    seqs = c(chims, oo, tp[1]), parents = tp,
    default = sapply(c(chims, oo, tp[1]), function(s) isBimera(s, tp[tp != s]), USE.NAMES = FALSE),
    one_off = sapply(c(chims, oo, tp[1]), function(s) isBimera(s, tp[tp != s], allowOneOff = TRUE), USE.NAMES = FALSE),
    shift2 = sapply(c(chims, oo, tp[1]), function(s) isBimera(s, tp[tp != s], maxShift = 2), USE.NAMES = FALSE)
  )
)
# isShiftDenovo: shifted copies of abundant sequences
sh_unq <- getUniques(ddF$sam1)[1:8]
shifted <- c(paste0(substr(names(sh_unq)[1], 6, 240), "ACGTA"), paste0("GGT", substr(names(sh_unq)[2], 1, 237)),
             substr(names(sh_unq)[3], 1, 200))
sh_all <- c(sh_unq, setNames(c(3L, 2L, 2L), shifted))
chim$shift <- list(uniques = unq_list(sh_all),
                   default = unname(isShiftDenovo(sh_all)),
                   ov100 = unname(isShiftDenovo(sh_all, minOverlap = 100)),
                   subseqs = unname(isShiftDenovo(sh_all, flagSubseqs = TRUE)))
wj(chim, "chimeras")

# ============================================================== taxonomy ===
log_msg("taxonomy")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
queries <- colnames(seqtab.nochim)
ten <- getSequences(ext("ten_16s.100.fa.gz"))
q_full <- unname(ten[seq(1, length(ten), by = 40)])  # 100 full-length queries
q_spec <- q_full[dada2:::C_isACGT(q_full)]          # assignSpecies requires ACGT-only
train <- ext("example_train_set.fa.gz")
spec <- ext("example_species_assignment.fa.gz")
tx <- function(seqs, seed = 100, ...) { set.seed(seed); assignTaxonomy(seqs, train, multithread = FALSE, ...) }
boot0 <- tx(queries, minBoot = 0, outputBootstraps = TRUE)
boot0_full <- tx(q_full, minBoot = 0, outputBootstraps = TRUE)
# Ties are broken with std::random_device upstream (not reproducible in R): record the
# set of genus-level assignments observed over repeated runs.
tie_sets <- function(seqs) {
  runs <- replicate(15, unname(tx(seqs, seed = sample.int(1e6, 1), minBoot = 0)[, "Genus"]))
  apply(runs, 1, function(r) sort(unique(ifelse(is.na(r), "NA", r))))
}
set.seed(7)
tq <- tie_sets(queries)
tf <- tie_sets(q_full)
rc_queries <- queries
rc_queries[seq(1, length(rc_queries), by = 2)] <- rc(rc_queries[seq(1, length(rc_queries), by = 2)])
taxa50 <- tx(queries)
taxo <- list(
  queries = queries, queries_full = q_full,
  boot0 = list(tax = df_list(as.data.frame(unname(boot0$tax))), boot = unname(boot0$boot)),
  boot0_full = list(tax = df_list(as.data.frame(unname(boot0_full$tax))), boot = unname(boot0_full$boot)),
  minboot50 = df_list(as.data.frame(unname(taxa50))),
  minboot50_full = df_list(as.data.frame(unname(tx(q_full)))),
  tryrc = list(queries = rc_queries, tax = df_list(as.data.frame(unname(tx(rc_queries, tryRC = TRUE, minBoot = 0))))),
  tie_sets_queries = lapply(tq, I), tie_sets_full = lapply(tf, I),
  queries_species = q_spec,
  species = list(
    default = df_list(as.data.frame(unname(assignSpecies(q_spec, spec)))),
    multiple = df_list(as.data.frame(unname(assignSpecies(q_spec, spec, allowMultiple = TRUE)))),
    multiple3 = df_list(as.data.frame(unname(assignSpecies(q_spec, spec, allowMultiple = 3)))),
    tryrc = df_list(as.data.frame(unname(assignSpecies(c(q_spec[1:20], rc(q_spec[21:40])), spec, tryRC = TRUE)))),
    asv = df_list(as.data.frame(unname(assignSpecies(queries, spec))))
  ),
  add_species = df_list(as.data.frame(unname(addSpecies(tx(q_spec), spec)))),
  add_species_multiple = df_list(as.data.frame(unname(addSpecies(tx(q_spec), spec, allowMultiple = TRUE))))
)
wj(taxo, "taxonomy")

# ===================================================== misc / helpers ===
log_msg("misc")
ufa <- file.path(TMP, "uniques.fasta")
uniquesToFasta(getUniques(ddF$sam1), ufa)
wj(list(uniques_to_fasta = readLines(ufa),
        get_uniques_seqtab_row = unq_list(getUniques(seqtab[1, ])),
        get_sequences_dada = getSequences(ddF$sam1)), "misc")

# ===================================================== e2e (vignette) ===
log_msg("e2e")
getN <- function(x) sum(getUniques(x))
track <- cbind(p1$res$reads_in, p1$res$reads_out, sapply(ddF, getN), sapply(ddR, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- c("sam1", "sam2")
wj(list(track = mat_list(track), seqtab_nochim = mat_list(seqtab.nochim),
        taxa = df_list(as.data.frame(unname(taxa50)))), "e2e")

log_msg("done")
