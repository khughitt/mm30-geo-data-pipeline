#!/bin/env/Rscript
###############################################################################
#
# GSE13591
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
  select(geo_accession, platform_id, title,
         sample_type = `sample type:ch1`)

sample_metadata$mm_stage <- sample_metadata$sample_type
sample_metadata$mm_stage[startsWith(sample_metadata$mm_stage, 'TC')] <- 'MM'
sample_metadata$mm_stage[startsWith(sample_metadata$mm_stage, 'N')] <- 'Healthy'

sample_metadata$disease_stage <- sample_metadata$mm_stage

#table(sample_metadata$sample_type)
#
# MGUS    N  PCL  TC1  TC2  TC3  TC4  TC5
#   11    5    9   29   25   49   24    6

#table(sample_metadata$disease_stage)
#
# Healthy    MGUS      MM     PCL
#       5      11     133       9

# drop PCL samples
mask <- sample_metadata$mm_stage  != 'PCL'

sample_metadata <- sample_metadata[mask, ]
expr_dat <- expr_dat[, c(TRUE, mask)]

# add platform
sample_metadata$platform_type <- 'Microarray'

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
