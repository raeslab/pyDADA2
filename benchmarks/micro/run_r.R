#!/usr/bin/env Rscript
# Micro-benchmarks of R dada2 1.36 on the bundled example data.
#
# Usage: Rscript benchmarks/micro/run_r.R [repeats] [outfile]
#
# Method (identical in benchmarks/micro/run_py.py):
#   * every stage runs once as warm-up; stages faster than MIN_SAMPLE seconds are then timed in
#     loops of k calls (k chosen from the warm-up), so the timer resolution does not distort
#     sub-millisecond stages; the reported time is the median over `repeats` samples of
#     (loop time / k);
#   * a full gc() runs before every timed sample (outside the timed region); garbage collections
#     that R's own functions trigger internally (assignSpecies, isBimeraDenovo) are part of their
#     time and are reported as caveats by report.py, using the gc() cost measured here;
#   * thread counts are equal on both sides: filter_and_trim is timed at 1 and 2 threads (R
#     parallelises over the 2 files only), the multithreaded R stages at 1 and 4;
#   * inputs are prepared outside the timed regions (e.g. seq_complexity gets the sequences,
#     not the FASTQ file).
suppressPackageStartupMessages({ library(dada2); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
REPEATS <- if (length(args) >= 1) as.integer(args[[1]]) else 7
OUT <- if (length(args) >= 2) args[[2]] else "benchmarks/micro/results_r.json"
MIN_SAMPLE <- 0.2  # seconds

ext <- function(f) system.file("extdata", f, package = "dada2")
tmp <- tempfile("bench"); dir.create(tmp)
fwd <- c(ext("sam1F.fastq.gz"), ext("sam2F.fastq.gz"))
rev <- c(ext("sam1R.fastq.gz"), ext("sam2R.fastq.gz"))
filtF <- file.path(tmp, c("F1.fq.gz", "F2.fq.gz"))
filtR <- file.path(tmp, c("R1.fq.gz", "R2.fq.gz"))
data(errBalancedF, errBalancedR)
mt <- function(n) if (n == 1) FALSE else n
now <- function() as.numeric(Sys.time())  # microsecond resolution

LOOPS <- list()
bench <- function(key, f) {
  q <- function() suppressWarnings(suppressMessages(f()))
  gc(verbose = FALSE)
  t0 <- now(); invisible(q()); warm <- now() - t0
  k <- as.integer(max(1, min(1000, ceiling(MIN_SAMPLE / max(warm, 1e-6)))))
  times <- vapply(seq_len(REPEATS), function(r) {
    gc(verbose = FALSE)
    t0 <- now()
    for (i in seq_len(k)) q()
    (now() - t0) / k
  }, numeric(1))
  LOOPS[[key]] <<- k
  median(times)
}

# Shared inputs (not timed)
suppressMessages(filterAndTrim(fwd, filtF, rev, filtR, truncLen = c(240, 160), maxN = 0,
                               maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE))
drpF <- derepFastq(filtF); drpR <- derepFastq(filtR)
ddF <- dada(drpF, err = errBalancedF, verbose = FALSE)
ddR <- dada(drpR, err = errBalancedR, verbose = FALSE)
mergers <- mergePairs(ddF, drpF, ddR, drpR)
st <- makeSequenceTable(mergers)
# Chimera table: st plus chimeras of its 6 most abundant sequences (breakpoint 120).
top <- colnames(st)[1:6]
chims <- unlist(lapply(seq_along(top), function(i) lapply(seq_along(top), function(j)
  if (i != j) paste0(substr(top[i], 1, 120), substr(top[j], 121, nchar(top[j]))))))
chims <- unique(chims[!chims %in% colnames(st)])
st_chim <- cbind(st, matrix(rep(c(5L, 3L), length(chims)), nrow = 2,
                            dimnames = list(rownames(st), chims)))
ten <- getSequences(ext("ten_16s.100.fa.gz"))
qfull <- unname(ten[seq(1, length(ten), by = 10)])                 # 400 full-length 16S
# assignSpecies queries that hit: 250-nt windows (every 50 nt) of the species reference.
sref <- unname(getSequences(ext("example_species_assignment.fa.gz")))
qspec <- unlist(lapply(sref, function(s) {
  starts <- seq(1, nchar(s) - 249, by = 50); substring(s, starts, starts + 249) }))
qspec <- qspec[dada2:::C_isACGT(qspec)]
F27 <- "AGRGTTYGATYMTGGCTCAG"; R1492rc <- rc("RGYTACCTTGTTACGACTT")
sam1_seqs <- getSequences(ext("sam1F.fastq.gz"))
spec_hits <- assignSpecies(qspec, ext("example_species_assignment.fa.gz"))

gc_cost <- median(replicate(5, { t0 <- now(); gc(verbose = FALSE); now() - t0 }))
res <- list(meta = list(r_version = R.version.string, dada2 = as.character(packageVersion("dada2")),
                        repeats = REPEATS, cores = parallel::detectCores(), min_sample = MIN_SAMPLE,
                        gc_seconds = gc_cost,
                        # gc() calls made inside the R function, per call of the stage
                        internal_gc = list(assign_species = length(unique(nchar(qspec))) + 1,
                                           remove_bimera_pooled = 1,
                                           remove_bimera_per_sample = nrow(st_chim)),
                        n_chimeras = length(chims), n_species_queries = length(qspec),
                        species_hits = sum(!is.na(spec_hits[, "Species"]))))
for (nt in c(1, 2)) {
  key <- paste0("filter_and_trim@", nt)
  res[[key]] <- bench(key, function() filterAndTrim(fwd, filtF, rev, filtR,
    truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, multithread = mt(nt)))
}
for (nt in c(1, 4)) {
  k <- function(name) paste0(name, "@", nt)
  res[[k("learn_errors")]] <- bench(k("learn_errors"), function() learnErrors(filtF, multithread = mt(nt)))
  res[[k("dada")]] <- bench(k("dada"), function() dada(drpF, err = errBalancedF, multithread = mt(nt), verbose = FALSE))
  res[[k("dada_pool")]] <- bench(k("dada_pool"), function() dada(drpF, err = errBalancedF, pool = TRUE, multithread = mt(nt), verbose = FALSE))
  res[[k("dada_pseudo")]] <- bench(k("dada_pseudo"), function() dada(drpF, err = errBalancedF, pool = "pseudo", multithread = mt(nt), verbose = FALSE))
  res[[k("remove_bimera_consensus")]] <- bench(k("remove_bimera_consensus"), function() removeBimeraDenovo(st_chim, method = "consensus", multithread = mt(nt)))
  res[[k("remove_bimera_pooled")]] <- bench(k("remove_bimera_pooled"), function() removeBimeraDenovo(st_chim, method = "pooled", multithread = mt(nt)))
  res[[k("remove_bimera_per_sample")]] <- bench(k("remove_bimera_per_sample"), function() removeBimeraDenovo(st_chim, method = "per-sample", multithread = mt(nt)))
  res[[k("assign_taxonomy")]] <- bench(k("assign_taxonomy"), function() assignTaxonomy(qfull, ext("example_train_set.fa.gz"), multithread = mt(nt)))
}
# single-threaded only in R
res[["derep_fastq@1"]] <- bench("derep_fastq@1", function() derepFastq(filtF))
res[["merge_pairs@1"]] <- bench("merge_pairs@1", function() mergePairs(ddF, drpF, ddR, drpR))
res[["make_sequence_table@1"]] <- bench("make_sequence_table@1", function() makeSequenceTable(mergers))
res[["collapse_no_mismatch@1"]] <- bench("collapse_no_mismatch@1", function() collapseNoMismatch(st_chim))
res[["assign_species@1"]] <- bench("assign_species@1", function() assignSpecies(qspec, ext("example_species_assignment.fa.gz")))
res[["remove_primers@1"]] <- bench("remove_primers@1", function() removePrimers(ext("samPBprimers.fastq.gz"), file.path(tmp, "p.fq.gz"), F27, R1492rc))
res[["remove_primers_indels@1"]] <- bench("remove_primers_indels@1", function() removePrimers(ext("samPBprimers.fastq.gz"), file.path(tmp, "p2.fq.gz"), F27, R1492rc, allow.indels = TRUE))
res[["seq_complexity@1"]] <- bench("seq_complexity@1", function() seqComplexity(sam1_seqs))
res$meta$loops <- LOOPS
writeLines(toJSON(res, auto_unbox = TRUE, digits = NA, pretty = TRUE), OUT)
message("wrote ", OUT)
