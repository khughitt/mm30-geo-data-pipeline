#!/bin/env/Rscript
###############################################################################
#
# Download GSE134598
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
supp_file <- file.path(cache_dir, "GSE134598_All_Processed_Data.txt.gz")

if (!file.exists(supp_file)) {
  getGEOSuppFiles(acc, baseDir = cache_dir, makeDirectory = FALSE)
}

# get expression data;
expr_dat <- read_tsv(supp_file, col_types = cols()) %>%
  select(-Chromosome, -Start, -End, -Length, -GeneBiotype, -GeneName) %>%
  rename(feature = GeneId)

pdata <- pData(eset)

# save a placeholder "row metadata" dataframe with just the gene identifiers (no
# metadata is available via fData)
fdata <- data.frame("feature" = expr_dat$feature)

write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
