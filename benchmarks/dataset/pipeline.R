#!/usr/bin/env Rscript
# DADA2 tutorial pipeline on the MiSeq SOP data with R dada2 1.36 (spec §9.2, phase 11).
#
# Usage: Rscript benchmarks/dataset/pipeline.R <threads> <outdir> [datadir] [fastq_subdir]
# fastq_subdir defaults to MiSeq_SOP; files ending in _R1_001.fastq(.gz) / _R2_001.fastq(.gz) are
# used. Samples with no reads left after filtering are dropped (neither tool writes their files).
# Mirrors benchmarks/dataset/pipeline.py stage by stage, with the same parameters. Writes
# <outdir>/results.json: per-stage wall times plus every intermediate result that is compared.
suppressPackageStartupMessages({ library(dada2); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
THREADS <- as.integer(args[[1]])
OUT <- args[[2]]
DATA <- if (length(args) >= 3) args[[3]] else "benchmarks/dataset/data"
mt <- if (THREADS == 1) FALSE else THREADS
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

times <- list()
stage <- function(name, expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- suppressMessages(expr)
  times[[name]] <<- proc.time()[["elapsed"]] - t0
  val
}

sop <- file.path(DATA, if (length(args) >= 4) args[[4]] else "MiSeq_SOP")
fnFs <- sort(list.files(sop, pattern = "_R1_001[.]fastq([.]gz)?$", full.names = TRUE), method = "radix")
fnRs <- sort(list.files(sop, pattern = "_R2_001[.]fastq([.]gz)?$", full.names = TRUE), method = "radix")
samples <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(OUT, "filtered", paste0(samples, "_F_filt.fastq.gz"))
filtRs <- file.path(OUT, "filtered", paste0(samples, "_R_filt.fastq.gz"))

t_total <- proc.time()[["elapsed"]]
out <- stage("filter_and_trim", filterAndTrim(
  fnFs, filtFs, fnRs, filtRs, truncLen = c(240, 160), maxN = 0, maxEE = c(2, 2), truncQ = 2,
  rm.phix = TRUE, compress = TRUE, multithread = mt))
# Neither tool writes the files of samples without reads left: drop those samples.
keep <- out[, 2] > 0
samples <- samples[keep]; filtFs <- filtFs[keep]; filtRs <- filtRs[keep]
errF <- stage("learn_errors_F", learnErrors(filtFs, multithread = mt))
errR <- stage("learn_errors_R", learnErrors(filtRs, multithread = mt))
drpF <- stage("derep_F", derepFastq(filtFs)); names(drpF) <- samples
drpR <- stage("derep_R", derepFastq(filtRs)); names(drpR) <- samples
ddF <- stage("dada_F", dada(drpF, err = errF, multithread = mt, verbose = FALSE))
ddR <- stage("dada_R", dada(drpR, err = errR, multithread = mt, verbose = FALSE))
mergers <- stage("merge_pairs", mergePairs(ddF, drpF, ddR, drpR))
seqtab <- stage("make_sequence_table", makeSequenceTable(mergers))
bim <- stage("remove_bimera_denovo",
             isBimeraDenovoTable(seqtab, multithread = mt))   # == removeBimeraDenovo("consensus")
nochim <- seqtab[, !bim, drop = FALSE]
set.seed(100)
taxa <- stage("assign_taxonomy", assignTaxonomy(
  nochim, file.path(DATA, "silva_nr99_v138.2_toGenus_trainset.fa.gz"), multithread = mt,
  outputBootstraps = TRUE))
species <- stage("add_species", addSpecies(
  taxa$tax, file.path(DATA, "silva_v138.2_assignSpecies.fa.gz")))
times$total <- proc.time()[["elapsed"]] - t_total

getN <- function(x) sum(getUniques(x))
uniq <- function(dd) lapply(dd, function(d) list(seqs = names(getUniques(d)),
                                                 abund = unname(getUniques(d))))
rowmajor <- function(m) as.vector(t(m))
res <- list(
  meta = list(tool = "R", version = paste(R.version.string, "/ dada2",
                                          as.character(packageVersion("dada2"))),
              threads = THREADS),
  times = times,
  samples = samples,
  track = list(input = unname(out[, 1]), filtered = unname(out[, 2]),
               denoisedF = unname(sapply(ddF, getN)), denoisedR = unname(sapply(ddR, getN)),
               merged = unname(sapply(mergers, getN)), nonchim = unname(rowSums(nochim))),
  err_F = rowmajor(getErrors(errF)), err_R = rowmajor(getErrors(errR)),
  derep_F = unname(sapply(drpF, function(d) length(d$uniques))),
  derep_R = unname(sapply(drpR, function(d) length(d$uniques))),
  dada_F = uniq(ddF), dada_R = uniq(ddR),
  mergers = lapply(mergers, function(m) list(sequence = m$sequence, abundance = m$abundance,
                                             forward = m$forward, reverse = m$reverse)),
  seqtab = list(seqs = colnames(seqtab), counts = rowmajor(unname(seqtab))),
  bimera = unname(bim),
  taxa = list(seqs = rownames(species), levels = colnames(species),
              values = rowmajor(unname(species))),
  boot = list(levels = colnames(taxa$boot), values = rowmajor(unname(taxa$boot)))
)
write_json(res, file.path(OUT, "results.json"), digits = I(17), auto_unbox = TRUE,
           na = "null")
