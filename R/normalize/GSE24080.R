#!/bin/env/Rscript
###############################################################################
#
# GSE24080
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
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         maqc_status = `maqc_distribution_status:ch1`) 

# add platform & disease stage
sample_metadata$platform_type <- 'Microarray'
sample_metadata$disease_stage <- 'MM'

# remove five samples with MAQC_Remove flag
# MAQC_Remove    Training  Validation
#          5         340         214
bad_samples <- sample_metadata %>%
  filter(maqc_status == 'MAQC_Remove') %>%
  pull(geo_accession)

sample_metadata <- sample_metadata %>%
  filter(maqc_status != 'MAQC_Remove')

expr_dat <- expr_dat[, !colnames(expr_dat) %in% bad_samples]

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
