#!/bin/env/Rscript
###############################################################################
#
# GSE117847
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, disease_stage = `patient diagnosis:ch1`) %>%
  mutate(disease_stage = recode(disease_stage, `active MM` = "MM")) %>%
  mutate(disease_stage = recode(disease_stage, `progressed SMM` = "SMM",
                                `non-progressed SMM` = "SMM"))

sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

# map from ensgenes to gene symbols ("SPOT_ID" column for this dataset contains ensembl
# gene identifiers..)
fdata$symbol <- grch38$symbol[match(fdata$SPOT_ID, grch38$ensgene)]

# get expression data and swap ensgenes for gene symbols
expr_dat <- dat %>%
  add_column(symbol = fdata$symbol, .before = 1) %>%
  filter(symbol != "") %>%
  as_tibble()

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
