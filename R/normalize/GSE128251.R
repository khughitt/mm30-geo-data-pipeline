#!/bin/env/Rscript
###############################################################################
#
# GSE128251
#
###############################################################################
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
  add_column(symbol = fdata$Symbol, .before = 1) %>%
  as_tibble()

# drop genes with missing symbols
mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# exclude erroneous "1-Mar", etc. excel gene names..
mask <- !grepl("^[0-9]+", expr_dat$symbol)
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         treatment = `treatment:ch1`, replicate = `replicate:ch1`,
         cell_line = source_name_ch1)

sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'H929'
sample_metadata$platform_type <- 'Microarray'

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
