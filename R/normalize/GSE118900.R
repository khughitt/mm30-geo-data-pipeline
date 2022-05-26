#!/bin/env/Rscript
###############################################################################
#
# GSE118900
#
###############################################################################
library(tidyverse)

# load data & metadata
expr_dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore gene symbol column)
expr_dat[, -1] <- sweep(expr_dat[, -1], 2, colSums(expr_dat[, -1]), '/') * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, title,  platform_id,
         sample_name = description,
         disease_stage = `tumor stage:ch1`,
         patient = `patient:ch1`,
         disease_stage = `tumor stage:ch1`) %>%
  mutate(disease_stage = recode(disease_stage, `IgM-MGUS` = 'MGUS', `NDMM` = 'MM'))

# add cell and platform type
sample_metadata$platform_type <- 'RNA-Seq'
sample_metadata$sample_type <- "Patient"

# normalize sample ids, ex: "IgM-MGUS1_C37" -> "IgM.MGUS1_C37"
sample_metadata$sample_name <- gsub('-', '.', sample_metadata$sample_name)

# drop single outlier patient NDMM7 whose samples had a low correlation with all other
# patient samples (median pairwise correlation ~0.06 vs. 0.6 for other patients)
exclude_patient <- 'NDMM7'

exclude_samples <- sample_metadata %>%
  filter(patient == exclude_patient) %>%
  pull(geo_accession)

expr_dat <- expr_dat[, !startsWith(colnames(expr_dat), exclude_patient)]

sample_metadata <- sample_metadata %>%
  filter(patient != exclude_patient)

# exclude any zero variance genes present
mask <- apply(expr_dat[, -1], 1, var) > 0
expr_dat <- expr_dat[mask, ]

# table(mask)
# mask
# FALSE  TRUE
#  3925 19473

# match sample order to metadata
expr_dat <- expr_dat[, c('symbol', sample_metadata$geo_accession)]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
