#!/bin/env/Rscript
###############################################################################
#
# GSE178340
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

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
write_csv(dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
