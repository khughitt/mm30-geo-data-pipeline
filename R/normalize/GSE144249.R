#!/bin/env/Rscript
###############################################################################
#
# GSE144249
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

# rename feature column
expr_dat <- dat %>%
  rename(symbol = feature)

# sample metadata columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, description,
         drug_resistance=`drug resistance:ch1`)

sample_metadata$replicate <- as.numeric(substr(sample_metadata$description, 12, 12))

# add platform, etc.
sample_metadata$platform_type <- 'Microarray'
sample_metadata$sample_type <- "Cell Line"
sample_metadata$cell_line <- "RPMI-8226"
sample_metadata$disease_stage <- NA

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch37[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
