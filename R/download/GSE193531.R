#!/bin/env/Rscript
###############################################################################
#
# GSE193531
#
###############################################################################
library(GEOquery)
library(tidyverse)
library(arrow)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expression data is missing from getGEO() query result and must be downloaded
# separately;
# this will also download cell level metadata used in a later step
supp_file <- file.path(cache_dir, "GSE193531_umi-count-matrix.csv.gz")

if (!file.exists(supp_file)) {
  getGEOSuppFiles(acc, baseDir = cache_dir, makeDirectory = FALSE)
}

# load UMI counts
#       feature AAACCTGAGGTAAACT.1.MM.5.138P AAACCTGGTCCCTACT.1.MM.5.138P
# 1 RP11-34P13.7                            0                            0
# 2 RP11-34P13.8                            0                            0
# 3   FO538757.3                            0                            0
expr_dat <- read.csv(gzfile(supp_file), header = TRUE, row.names = 1) %>%
  rownames_to_column("feature")

# [1] 22273 29388
# dim(expr_dat)

# table(ids %in% grch38$symbol)
# FALSE  TRUE
# 4747 17526

#> table(ids %in% grch37$symbol)
#FALSE  TRUE
# 1439 20834

pdata <- pData(eset)

# save a placeholder "row metadata" dataframe with just the gene identifiers (no
# metadata is available via fData)
fdata <- data.frame("feature" = expr_dat$feature)

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])

