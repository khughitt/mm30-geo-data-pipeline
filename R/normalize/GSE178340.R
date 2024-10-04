#!/bin/env/Rscript
###############################################################################
#
# GSE178340
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]])
pdata <- read_feather(snakemake@input[[3]])

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, description, 
         cell_line = `cell line:ch1`, tissue = `tissue:ch1`, treatment = `treatment:ch1`)

sample_metadata$replicate <- as.numeric(str_split(pdata$description, " ", simplify = TRUE)[, 2])

sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Cell Line"
sample_metadata$disease_stage <- NA

# drop rows with missing gene symbols
dat <- dat[!is.na(dat$symbol), ]

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(dat$symbol, grch38$symbol), ]

# store results
write_feather(dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
