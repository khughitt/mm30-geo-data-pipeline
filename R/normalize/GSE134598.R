#!/bin/env/Rscript
###############################################################################
#
# GSE134598
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore feature column)
expr_dat <- dat

expr_dat[, -1] <- sweep(expr_dat[, -1], 2, colSums(expr_dat[, -1]), '/') * 1E6

# map from ensgenes -> gene symbols
ind <- match(expr_dat$feature, grch38$ensgene)
gene_symbols <- grch38$symbol[ind]

expr_dat <- expr_dat %>%
  select(-feature) %>%
  add_column(symbol = gene_symbols, .before = 1)

fdata$symbol <- gene_symbols

# drop genes with missing symbols
mask <- !is.na(gene_symbols)

expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# drop empty rows
expr_dat <- expr_dat[rowSums(expr_dat[, -1]) > 0, ]

# get relevant sample metadata
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         cell_type = `cell line:ch1`,
         sample_type = `sample type:ch1`,
         treatment = `treatment:ch1`, dose = `dose:ch1`)

sample_metadata$disease_stage <- "MM"
sample_metadata$platform_type <- 'RNA-Seq'

if (!all(colnames(expr_dat)[-1] == sample_metadata$title)) {
  stop("Sample ID mismatch!")
}

# for consistency, use GEO sample accessions as ids
colnames(expr_dat) <- c('symbol', sample_metadata$geo_accession)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
