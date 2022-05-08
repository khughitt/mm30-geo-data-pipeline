#!/bin/env/Rscript
###############################################################################
#
# GSE178340 (Wang et al., 2022)
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

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, description, 
         cell_line = `cell line:ch1`, tissue = `tissue:ch1`, treatment = `treatment:ch1`)

sample_metadata$replicate <- as.numeric(str_split(pdata$description, " ", simplify = TRUE)[, 2])

  # mutate(mm_stage = recode(mm_stage, `active MM` = 'MM')) %>%
  # mutate(disease_stage = recode(mm_stage, `progressed SMM` = 'SMM', `non-progressed SMM` = 'SMM'))

sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# drop rows with missing gene symbols
dat <- dat[!is.na(dat$symbol), ]

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])