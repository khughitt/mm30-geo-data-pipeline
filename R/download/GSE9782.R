#!/bin/env/Rscript
###############################################################################
#
# Download GSE9782
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

# GSE7039 includes two esets/arrays for each patient, Affy HG-U133A & HG-U133B, which
# overlap in a very small number of probes (n = 168).

# Below, the two arrays are processed separately following the usual approach, and then combined
esets <- getGEO(acc, destdir = cache_dir)

e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

e1 <- e1 %>%
  as.data.frame() %>%
  rownames_to_column("feature")

e2 <- e2 %>%
  as.data.frame() %>%
  rownames_to_column("feature")

colnames(e2) <- colnames(e1)

# join expression datasets
expr_dat <- rbind(e1, e2)

# get sample & feature metadta
pdata <- pData(esets[[1]])

# each patient is associated with two sample ids: one for each microarray
pdata$geo_accession2 <- pData(esets[[2]])$geo_accession

fdata <- rbind(fData(esets[[1]]), fData(esets[[2]]))

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
