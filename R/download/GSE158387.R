#!/bin/env/Rscript
###############################################################################
#
# Download GSE158387 
#
###############################################################################
library(GEOquery)
library(tidyverse)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expression data is missing from getGEO() query result and must be downloaded
# separately
supp_file1 <- file.path(cache_dir, 'GSE158387_RawCounts.tsv.gz')

if (!file.exists(supp_file1)) {
  getGEOSuppFiles(acc, baseDir = cache_dir)
}

pdata <- pData(eset)

# load raw counts
#   Symbol              Biotype IGF106703
# A1BG-AS1            antisense        53
#     A1BG processed_transcript        13
#     A1BG       protein_coding         0
expr_dat <- read.delim(gzfile(supp_file1))

# save gene metadata (just symbol & biotype in this case)
fdata <- expr_dat[, 1:2]

# drop biotype column
expr_dat <- expr_dat[, -2]

# sanity check
# all(colnames(expr_dat)[-1] == pdata$title)
# TRUE

colnames(expr_dat) <- c("symbol", pdata$geo_accession)

write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
