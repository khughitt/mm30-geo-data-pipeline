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
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE) %>%
  column_to_rownames("feature")

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization
dat <- sweep(dat, 2, colSums(dat), '/') * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, mm_stage = `patient diagnosis:ch1`) %>%
  mutate(mm_stage = recode(mm_stage, `active MM` = 'MM')) %>%
  mutate(disease_stage = recode(mm_stage, `progressed SMM` = 'SMM', `non-progressed SMM` = 'SMM'))

sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# map from ensgenes to gene symbols ("SPOT_ID" column for this dataset contains ensembl
# gene identifiers..)
fdata$symbol <- grch38$symbol[match(fdata$SPOT_ID, grch38$ensgene)]

# get expression data and swap ensgenes for gene symbols
expr_dat <- dat %>%
  add_column(symbol = fdata$symbol, .before = 1) %>%
  filter(symbol != '') %>%
  as_tibble()

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a new data package, based off the old one
pkgr <- Packager$new()

# resource list (consider converting values to lists and including "data_type" field for
# each resource? i.e. "data"/"row metadata"/"column metadata")
resources <- list(
  "data" = expr_dat,
  "row-metadata" = fdata,
  "column-metadata" = sample_metadata
)

# changes to make to metadata
mdata <- list(
  data=list(
    processing="reprocessed"
  ),
  rows="symbol"
)

# annotations
annot <- list("data-prep" = read_file("annot/prepare-data/GSE117847.md"))

pkg <- pkgr$update_package(snakemake@input[[4]], mdata, "Reprocess data", 
                           resources, annotations=annot)

pkg %>%
  write_package(dirname(snakemake@output[[1]]))
