#!/bin/env Rscript
#
# Extract data processing information from raw GEO datasets
# KH (April 2020)
#
infiles <- Sys.glob("data_processing/GSE*_data_processing.txt")

for (infile in infiles) {
  x <- read.table(infile)
  cat(basename(infile))
  cat("\n")
  cat(levels(x[1, 2]))
  cat("\n")
  cat("\n")
}

