#!/usr/bin/env Rscript
# Extract karyoploteR's bundled cytoband and chromosome-length data to TSV
# files that pykaryoplot can ship without needing R at runtime.
#
# Usage:  Rscript scripts/extract_sysdata.R /path/to/karyoploteR/R/sysdata.rda
suppressPackageStartupMessages({
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
rda <- if (length(args) >= 1) args[1] else "../karyoploteR/R/sysdata.rda"
out_genomes <- "src/pykaryoplot/data/genomes"
out_cyto    <- "src/pykaryoplot/data/cytobands"
dir.create(out_genomes, showWarnings = FALSE, recursive = TRUE)
dir.create(out_cyto,    showWarnings = FALSE, recursive = TRUE)

env <- new.env()
load(rda, envir = env)
dc <- env$data.cache
if (is.null(dc)) stop("data.cache not found in ", rda)

# genomes: each entry is a GRanges with one range per chromosome
for (g in names(dc$genomes)) {
  gr <- dc$genomes[[g]]
  df <- data.frame(chrom = as.character(seqnames(gr)),
                   start = start(gr),
                   end   = end(gr),
                   stringsAsFactors = FALSE)
  write.table(df, file = file.path(out_genomes, paste0(g, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# cytobands: each entry is a GRanges with name + gieStain mcols
for (g in names(dc$cytobands)) {
  gr <- dc$cytobands[[g]]
  if (length(gr) == 0) next
  m <- mcols(gr)
  name      <- if ("name"     %in% names(m)) as.character(m$name)     else rep("", length(gr))
  gie_stain <- if ("gieStain" %in% names(m)) as.character(m$gieStain) else rep("gpos50", length(gr))
  df <- data.frame(chrom = as.character(seqnames(gr)),
                   start = start(gr),
                   end   = end(gr),
                   name  = name,
                   gie_stain = gie_stain,
                   stringsAsFactors = FALSE)
  write.table(df, file = file.path(out_cyto, paste0(g, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Extracted", length(dc$genomes), "genomes and", length(dc$cytobands), "cytoband sets\n")
