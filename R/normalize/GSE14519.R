#!/bin/env/Rscript
###############################################################################
#
# GSE14519
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

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title)

# treatment status
sample_metadata$treatment <- 'Untreated'
sample_metadata$treatment[grepl('DAR 6', sample_metadata$title)] <- 'DAR (6 hrs)'
sample_metadata$treatment[grepl('DAR 24', sample_metadata$title)] <- 'DAR (24 hrs)'
sample_metadata$treatment[grepl('ATO 6', sample_metadata$title)] <- 'ATO (6 hrs)'
sample_metadata$treatment[grepl('ATO 24', sample_metadata$title)] <- 'ATO (24 hrs)'
sample_metadata$treatment[grepl('ATO 48', sample_metadata$title)] <- 'ATO (48 hrs)'

# cell line
sample_metadata$cell_type <- 'U266'
sample_metadata$cell_type[grepl('MM.1s', sample_metadata$title)] <- 'MM.1s'
sample_metadata$cell_type[grepl('KMS11', sample_metadata$title)] <- 'KMS11'
sample_metadata$cell_type[grepl('8226', sample_metadata$title)] <- '8226'

# disease stage
sample_metadata$disease_stage <- 'MM'

# add platform
sample_metadata$platform_type <- 'Microarray'

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
