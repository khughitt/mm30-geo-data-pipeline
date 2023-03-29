#!/bin/env/Rscript
###############################################################################
#
# GSE26760
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

# load additional metadata from
# http://portals.broadinstitute.org/mmgp/data/browseData?conversationPropagation=begin
mmrc_sample_info <- read_csv("supp/clean/mmrc.sample.information.csv", col_types = cols()) %>%
  rename(patient_id = Array)

pdata <- pData(eset)

# parse out patient ids
pdata$patient_id <- str_split(pdata$source_name_ch1, " ", simplify = TRUE)[, 4] 

pdata <- pdata %>%
  inner_join(mmrc_sample_info, by = "patient_id")

fdata <- fData(eset)

write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(pdata, snakemake@output[[3]])
