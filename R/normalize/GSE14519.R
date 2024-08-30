#!/bin/env/Rscript
###############################################################################
#
# GSE14519
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
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# drop rows with missing values
expr_dat <- expr_dat[complete.cases(expr_dat), ]

expr_dat <- expr_dat[expr_dat$symbol != "", ]

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title)

# treatment status
sample_metadata$treatment <- "Untreated"
sample_metadata$treatment[grepl("DAR 6", sample_metadata$title)] <- "DAR (6 hrs)"
sample_metadata$treatment[grepl("DAR 24", sample_metadata$title)] <- "DAR (24 hrs)"
sample_metadata$treatment[grepl("ATO 6", sample_metadata$title)] <- "ATO (6 hrs)"
sample_metadata$treatment[grepl("ATO 24", sample_metadata$title)] <- "ATO (24 hrs)"
sample_metadata$treatment[grepl("ATO 48", sample_metadata$title)] <- "ATO (48 hrs)"

# cell line
sample_metadata$cell_line <- "U266"
sample_metadata$cell_line[grepl("MM.1s", sample_metadata$title)] <- "MM.1s"
sample_metadata$cell_line[grepl("KMS11", sample_metadata$title)] <- "KMS11"
sample_metadata$cell_line[grepl("8226", sample_metadata$title)] <- "8226"

# disease stage
sample_metadata$disease_stage <- "MM"

# add platform
sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Cell Line"

sample_metadata$disease_stage <- NA

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
