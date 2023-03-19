#!/bin/env/Rscript
###############################################################################
#
# GSE47552
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$`Gene symbol`, .before = 1) %>%
  as_tibble()

mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# drop rows with missing values
expr_dat <- expr_dat[complete.cases(expr_dat), ]

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, 
         disease_stage_raw = `cell type:ch1`)

# add disease stage
sample_metadata$disease_stage <- rep("MM", length(sample_metadata$disease_stage_raw))

sample_metadata$disease_stage[grepl("Normal", sample_metadata$disease_stage_raw)] <- "Healthy"
sample_metadata$disease_stage[grepl("MGUS", sample_metadata$disease_stage_raw)] <- "MGUS"
sample_metadata$disease_stage[grepl("SMM", sample_metadata$disease_stage_raw)] <- "SMM"

sample_metadata <- sample_metadata %>%
  select(-disease_stage_raw)

# add disease stage
sample_metadata$disease_stage <- sample_metadata$disease_stage

# add platform
sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
