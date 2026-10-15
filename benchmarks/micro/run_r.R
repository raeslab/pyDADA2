#!/usr/bin/env Rscript
# Micro-benchmarks of R dada2 1.36 on the bundled example data (spec §9.1).
#
# Usage: Rscript benchmarks/micro/run_r.R [repeats] [outfile]
# Each stage is run once as warm-up, then `repeats` times; the median elapsed time is reported.
# Stage names and inputs mirror benchmarks/micro/run_py.py.
suppressPackageStartupMessages({ library(dada2); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
REPEATS <- if (length(args) >= 1) as.integer(args[[1]]) else 7
OUT <- if (length(args) >= 2) args[[2]] else "benchmarks/micro/results_r.json"

ext <- function(f) system.file("extdata", f, package = "dada2")
tmp <- tempfile("bench"); dir.create(tmp)
fwd <- c(ext("sam1F.fastq.gz"), ext("sam2F.fastq.gz"))
rev <- c(ext("sam1R.fastq.gz"), ext("sam2R.fastq.gz"))
filtF <- file.path(tmp, c("F1.fq.gz", "F2.fq.gz"))
filtR <- file.path(tmp, c("R1.fq.gz", "R2.fq.gz"))
data(errBalancedF, errBalancedR)
mt <- function(n) if (n == 1) FALSE else n

bench <- function(f) {
  q <- function() suppressWarnings(suppressMessages(f()))
  invisible(q())
  median(replicate(REPEATS, system.time(q())[["elapsed"]]))
}

# Shared inputs
suppressMessages(filterAndTrim(fwd, filtF, rev, filtR, truncLen = c(240, 160), maxN = 0,
                               maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE))
drpF <- derepFastq(filtF); drpR <- derepFastq(filtR)
ddF <- dada(drpF, err = errBalancedF, verbose = FALSE)
ddR <- dada(drpR, err = errBalancedR, verbose = FALSE)
mergers <- mergePairs(ddF, drpF, ddR, drpR)
st <- makeSequenceTable(mergers)
ten <- getSequences(ext("ten_16s.100.fa.gz"))
qfull <- unname(ten[seq(1, length(ten), by = 10)])                 # 400 full-length 16S
qspec <- qfull[dada2:::C_isACGT(qfull)]
F27 <- "AGRGTTYGATYMTGGCTCAG"; R1492rc <- rc("RGYTACCTTGTTACGACTT")

res <- list(meta = list(r_version = R.version.string, dada2 = as.character(packageVersion("dada2")),
                        repeats = REPEATS, cores = parallel::detectCores()))
for (nt in c(1, 4)) {
  k <- function(name) paste0(name, "@", nt)
  res[[k("filter_and_trim")]] <- bench(function() filterAndTrim(fwd, filtF, rev, filtR,
    truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, multithread = mt(nt)))
  res[[k("learn_errors")]] <- bench(function() learnErrors(filtF, multithread = mt(nt)))
  res[[k("dada")]] <- bench(function() dada(drpF, err = errBalancedF, multithread = mt(nt), verbose = FALSE))
  res[[k("dada_pool")]] <- bench(function() dada(drpF, err = errBalancedF, pool = TRUE, multithread = mt(nt), verbose = FALSE))
  res[[k("dada_pseudo")]] <- bench(function() dada(drpF, err = errBalancedF, pool = "pseudo", multithread = mt(nt), verbose = FALSE))
  res[[k("remove_bimera_consensus")]] <- bench(function() removeBimeraDenovo(st, method = "consensus", multithread = mt(nt)))
  res[[k("remove_bimera_pooled")]] <- bench(function() removeBimeraDenovo(st, method = "pooled", multithread = mt(nt)))
  res[[k("assign_taxonomy")]] <- bench(function() assignTaxonomy(qfull, ext("example_train_set.fa.gz"), multithread = mt(nt)))
}
# single-threaded only in R
res[["derep_fastq@1"]] <- bench(function() derepFastq(filtF))
res[["merge_pairs@1"]] <- bench(function() mergePairs(ddF, drpF, ddR, drpR))
res[["make_sequence_table@1"]] <- bench(function() makeSequenceTable(mergers))
res[["collapse_no_mismatch@1"]] <- bench(function() collapseNoMismatch(st))
res[["assign_species@1"]] <- bench(function() assignSpecies(qspec, ext("example_species_assignment.fa.gz")))
res[["remove_primers@1"]] <- bench(function() removePrimers(ext("samPBprimers.fastq.gz"), file.path(tmp, "p.fq.gz"), F27, R1492rc))
res[["remove_primers_indels@1"]] <- bench(function() removePrimers(ext("samPBprimers.fastq.gz"), file.path(tmp, "p2.fq.gz"), F27, R1492rc, allow.indels = TRUE))
res[["seq_complexity@1"]] <- bench(function() seqComplexity(getSequences(ext("sam1F.fastq.gz"))))
writeLines(toJSON(res, auto_unbox = TRUE, digits = 6, pretty = TRUE), OUT)
message("wrote ", OUT)
