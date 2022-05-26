#!/bin/env/Rscript
###############################################################################
#
# GSE162205
#
###############################################################################
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

# replace combined ensgene/symbol id column with symbols alone
dat$symbol <- fdata$symbol

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title, 
         cell_line = `cell type:ch1`, 
         time_hours =`time point:ch1`,
         treatment = `treatment:ch1`)

sample_metadata$time_hours <- as.numeric(sub("h", "", sample_metadata$time_hours))
sample_metadata$replicate <- as.numeric(endsWith(sample_metadata$title, "B")) + 1

sample_metadata$platform_type <- 'RNA-Seq'
sample_metadata$sample_type <- "Cell Line"

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
