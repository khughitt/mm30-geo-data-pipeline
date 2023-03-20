#!/bin/env/Rscript
###############################################################################
#
# Download GSE106218
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
supp_file <- file.path(cache_dir, "GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz")

if (!file.exists(supp_file)) {
  getGEOSuppFiles(acc, baseDir = cache_dir, makeDirectory = FALSE)
}

# load TPM counts
#         MM02_38 MM02_48 MM02_50
# 5S_rRNA       0       0       0
# 7SK           0       0       0
# A1BG          0       0       0
expr_dat <- read.delim(gzfile(supp_file), row.names = 1) %>%
  rownames_to_column("feature")

pdata <- pData(eset)

# save a placeholder "row metadata" dataframe with just the gene identifiers (no
# metadata is available via fData)
fdata <- data.frame("feature" = expr_dat$feature)

write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
