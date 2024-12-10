#!/bin/env/Rscript
###############################################################################
#
# Download GSE128251
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
expr_dat <- exprs(eset)

# shift the normalized expression measurements so that the minimum value is "0"
expr_dat <- expr_dat + abs(min(expr_dat))

expr_dat <- expr_dat %>%
    as.data.frame() %>%
    rownames_to_column("feature")

pdata <- pData(eset)
fdata <- fData(eset)

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
