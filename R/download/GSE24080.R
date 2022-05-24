#!/bin/env/Rscript
###############################################################################
#
# GSE24080
#
# Note: this dataset includes some of the same samples from GSE2658, but processed in a
# different manner, and with different metadata.
# This verion of the dataset (MAQC-II) includes ~2x samples before filtering and also
# includes survival-related metadata.
#
###############################################################################
library(GEOquery)
library(tidyverse)

acc <- snakemake@wildcards[["acc"]]
cache_dir <- file.path("/data/raw/geo", acc)

# create cache dir
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE, mode = "0755")
}

eset <- getGEO(acc, destdir = cache_dir)[[1]]

expr_dat <- exprs(eset) %>%
    as.data.frame() %>%
    rownames_to_column("feature")

# load supplemental clinical metadata;
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24080/suppl/GSE24080%5FMM%5FUAMS565%5FClinInfo%5F27Jun2008%5FLS%5Fclean%2Exls%2Egz
clinical_metadata <- read_tsv('supp/clean/GSE24080_MM_UAMS565_ClinInfo_27Jun2008_LS_clean.tsv', 
                              col_types = cols())

clinical_metadata <- clinical_metadata %>%
  mutate(title = sub('.CEL', '', cel_filename)) %>%
  select(-cel_filename)

# add to geo sample metadata
pdata <- pData(eset) %>%
  inner_join(clinical_metadata, by = 'title')

fdata <- fData(eset)

write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(pdata, snakemake@output[[3]])
