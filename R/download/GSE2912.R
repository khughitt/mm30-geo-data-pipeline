#!/bin/env/Rscript
###############################################################################
#
# GSE2912
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

expr_dat <- exprs(eset) %>%
    as.data.frame() %>%
    rownames_to_column("feature")

pdata <- pData(eset)

# load additional metadata from Agnelli et al. (2005), Appendix A
# https://pubmed.ncbi.nlm.nih.gov/16129847/
supp_mdat <- read_csv("supp/clean/agnelli2005.csv", show_col_types = FALSE) %>%
  rename(patient_id = Patient)

# reoder patient metadata
pdata$patient_id <- str_extract(pdata$title, "MM-[0-9]+")
supp_mdat <- supp_mdat[match(supp_mdat$patient_id, pdata$patient_id), ]

pdata <- pdata %>%
  inner_join(supp_mdat, by = "patient_id")

fdata <- fData(eset)

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
