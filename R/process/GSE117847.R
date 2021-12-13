#!/bin/env/Rscript
###############################################################################
#
# GEO Dataset processing - GSE117847
#
# A retained transcriptomic profile characterizes the CD138+ cells in the progression
# from smoldering to active multiple myeloma
#
# Storti et al. (2019)
#
###############################################################################
library(annotables)
library(tidyverse)
library(iodag)

# GEO accession
accession <- 'GSE117847'

# directory to store raw and processed data
data_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat <- read_csv(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

fdata <- read_csv(snakemake@input[[2]])
pdata <- read_csv(snakemake@input[[3]])

# size factor normalization
dat <- sweep(dat, 2, colSums(dat), '/') * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, mm_stage = `patient diagnosis:ch1`) %>%
  mutate(mm_stage = recode(mm_stage, `active MM` = 'MM')) %>%
  mutate(disease_stage = recode(mm_stage, `progressed SMM` = 'SMM', `non-progressed SMM` = 'SMM'))

sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# map from ensgenes to gene symbols
fdata$symbol <- grch38$symbol[match(fdata$SPOT_ID, grch38$ensgene)]

# get expression data and swap ensgenes for gene symbols
expr_dat <- dat %>%
  add_column(symbol = fdata$symbol, .before = 1) %>%
  filter(symbol != '') %>%
  as_tibble()

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

pkg_dir <- dirname(snakemake@output[[1]])
setwd(pkg_dir)

# create a new data package, based off the old one
pkgr <- Packager$new()

updates <- list(
  data=list(
    processing="reprocessed"
  ),
  rows="symbol"
)

pkg <- pkgr$update_package(snakemake@input[[4]], updates, "reprocess data",
                           expr_dat, fdata, sample_metadata)

pkg %>%
  write_package(pkg_dir)
