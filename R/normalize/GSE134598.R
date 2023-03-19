#!/bin/env/Rscript
###############################################################################
#
# GSE134598
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
expr_dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# map from ensgenes -> gene symbols
ind <- match(expr_dat$feature, grch38$ensgene)
gene_symbols <- grch38$symbol[ind]

expr_dat <- expr_dat %>%
  select(-feature) %>%
  add_column(symbol = gene_symbols, .before = 1)

# drop genes with missing symbols
mask <- !is.na(gene_symbols)
expr_dat <- expr_dat[mask, ]

# drop empty rows
expr_dat <- expr_dat[rowSums(expr_dat[, -1]) > 0, ]

# get relevant sample metadata
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         cell_line = `cell line:ch1`,
         sample_type = `sample type:ch1`,
         treatment = `treatment:ch1`, dose = `dose:ch1`)

sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Mixed"

sample_metadata$disease_stage <- NA

if (!all(colnames(expr_dat)[-1] == sample_metadata$title)) {
  stop("Sample ID mismatch!")
}

# for consistency, use GEO sample accessions as ids
colnames(expr_dat) <- c("symbol", sample_metadata$geo_accession)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
