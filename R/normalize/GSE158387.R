#!/bin/env/Rscript
###############################################################################
#
# GSE158387
#
###############################################################################
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, title, description, platform_id, 
         cell_line = `cell line:ch1`, replicate = `replicate:ch1`,
         time_days = `time point:ch1`,
         treatment = `treatment:ch1`)

sample_metadata$replicate <- as.numeric(substr(sample_metadata$replicate, 10, 11))

sample_metadata$time_days <- as.numeric(str_split(sample_metadata$time_days, " ", simplify = TRUE)[, 2])

  # mutate(mm_stage = recode(mm_stage, `active MM` = 'MM')) %>%
  # mutate(disease_stage = recode(mm_stage, `progressed SMM` = 'SMM', `non-progressed SMM` = 'SMM'))

sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'RNA-Seq'

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
