#!/bin/env/Rscript
###############################################################################
#
# GSE31161
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# output directory to store data packages to
out_dir <- dirname(snakemake@output[[1]])

# load data & metadata
expr_dat <- read_feather(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# there are three samples in GSE31161 with all missing data, which are excluded below..
num_missing <- apply(expr_dat, 2, function(x) {
  sum(is.na(x))
})

# table(num_missing)
# num_missing
#     0 54675
#  1035     3

expr_dat <- expr_dat[, num_missing == 0]

# exclude outlier samples;
# these samples were found to have very median pairwise correlations (0.38 - 0.68)
# compared to all other samples (>0.9 for most)
exclude_samples <- c("GSM771497", "GSM772341", "GSM772335")

expr_dat <- expr_dat[, !colnames(expr_dat) %in% exclude_samples]

# drop genes with missing symbols
mask <- !is.na(fdata$`Gene Symbol`)
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# add gene symbol column
expr_dat <- expr_dat %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

expr_dat <- expr_dat[expr_dat$symbol != "", ]

# columns to include
sample_metadata <- pdata %>%
  filter(geo_accession %in% colnames(expr_dat)) %>%
  select(geo_accession, platform_id,
         treatment = `treatment:ch1`, time_of_testing = `time of testing:ch1`) %>%
  mutate(relapsed = time_of_testing == "relapse")

# add disease stage
sample_metadata$disease_stage <- ifelse(sample_metadata$relapsed, "RRMM", "MM")

# add platform and cell type (same for all samples)
sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

platform <- pdata$platform_id[1]

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
