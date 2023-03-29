#!/bin/env/Rscript
###############################################################################
#
# GSE128251
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]])

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$Symbol, .before = 1) %>%
  as_tibble()

# drop genes with missing symbols
mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]

# exclude erroneous "1-Mar", etc. excel gene names..
mask <- !grepl("^[0-9]+", expr_dat$symbol)
expr_dat <- expr_dat[mask, ]

# drop rows with missing values
expr_dat <- expr_dat[complete.cases(expr_dat), ]

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         treatment = `treatment:ch1`, replicate = `replicate:ch1`,
         cell_line = source_name_ch1)

sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Cell Line"

sample_metadata$disease_stage <- NA

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
