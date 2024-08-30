#!/bin/env/Rscript
###############################################################################
#
# GSE158387
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]])
pdata <- read_feather(snakemake@input[[3]])

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), "/") * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, title, description, platform_id,
         cell_line = `cell line:ch1`, replicate = `replicate:ch1`,
         time_days = `time point:ch1`,
         treatment = `treatment:ch1`)

sample_metadata$replicate <- as.numeric(substr(sample_metadata$replicate, 10, 11))
sample_metadata$time_days <- as.numeric(str_split(sample_metadata$time_days, " ", simplify = TRUE)[, 2])

sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Cell Line"
sample_metadata$disease_stage <- NA

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(dat$symbol, grch38$symbol), ]

# store results
write_feather(dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
