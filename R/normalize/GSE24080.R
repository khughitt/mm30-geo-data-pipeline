#!/bin/env/Rscript
###############################################################################
#
# GSE24080
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
dat <- read_feather(snakemake@input[[1]])

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
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
         age, race, sex, isotype, os_time,
         os_censor,
         pfs_censor = efs_censor,
         pfs_time = efs_time,
         maqc_status = `maqc_distribution_status:ch1`)

sample_metadata <- sample_metadata %>%
  mutate(os_time = as.numeric(os_time),
         os_censor = as.logical(as.numeric(os_censor)),
         pfs_time = as.numeric(pfs_time),
         pfs_censor = as.logical(as.numeric(pfs_censor)),
         race = as.factor(race),
         sex = as.factor(sex),
         isotype = as.factor(isotype))


# add platform & disease stage
sample_metadata$platform_type <- "Microarray"
sample_metadata$disease_stage <- "MM"
sample_metadata$sample_type <- "Patient"

# remove five samples with MAQC_Remove flag
# MAQC_Remove    Training  Validation
#          5         340         214
bad_samples <- sample_metadata %>%
  filter(maqc_status == "MAQC_Remove") %>%
  pull(geo_accession)

sample_metadata <- sample_metadata %>%
  filter(maqc_status != "MAQC_Remove")

expr_dat <- expr_dat[, !colnames(expr_dat) %in% bad_samples]

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
