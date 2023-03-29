#!/bin/env/Rscript
###############################################################################
#
# GSE144249
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

# GSE144249 includes two esets: the first one for the human samples, and second one for
# mouse samples
eset <- getGEO(acc, destdir = cache_dir)[[1]]

# expr data is provided as a supplementary file on GEO
# https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144249&format=file
# for convenience, the per-sample counts have been assembled into a single table, using
# the first column from each sample and dropping the first four QA-related rows.
expr_dat <- read_tsv("supp/clean/GSE144249_raw_counts.tsv", col_types = cols()) %>%
  rename(feature = symbol)

pdata <- pData(eset)
fdata <- data.frame("feature" = expr_dat$feature)

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
