#!/bin/env/Rscript
###############################################################################
#
# GSE162205
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]])

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# size factor normalization (ignore gene symbol column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), "/") * 1E6

# replace combined ensgene/symbol id column with symbols alone
dat$symbol <- fdata$symbol

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, platform_id, title,
         cell_line = `cell type:ch1`,
         time_hours = `time point:ch1`,
         treatment = `treatment:ch1`)

sample_metadata$time_hours <- as.numeric(sub("h", "", sample_metadata$time_hours))

# GSE162205 includes three replicates and two technical replicates for each condition
sample_metadata$replicate <- letters[as.numeric(factor(substr(sample_metadata$title, 7, 7)))]
sample_metadata$technical_replicate <- factor(substr(sample_metadata$title,
                                                     nchar(sample_metadata$title),
                                                     nchar(sample_metadata$title)))

# replace NA values for control samples; for "time", the earlier time point of 4 hours
# is arbitrarily used
is_control <- is.na(sample_metadata$treatment)

sample_metadata$treatment[is_control] <- "Control"
sample_metadata$time_hours[is_control] <- 4

sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Cell Line"
sample_metadata$disease_stage <- NA

if (!all(colnames(dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(dat$symbol, grch38$symbol), ]

# store results
write_feather(dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
