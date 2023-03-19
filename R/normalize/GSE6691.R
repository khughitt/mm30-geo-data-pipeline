#!/bin/env/Rscript
###############################################################################
#
# GSE6691
#
###############################################################################
library(annotables)
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

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
  select(geo_accession, platform_id, title)

patient_desc <- pdata$characteristics_ch1

disease <- rep("MM", length(patient_desc))
disease[grepl("healthy", patient_desc)] <- "Healthy"
disease[grepl("WM", patient_desc)] <- "Waldenström's Macroglobulinemia"
disease[grepl("CLL", patient_desc)] <- "Chronic Lymphocytic Leukemia"

sample_metadata$disease_stage <- disease
sample_metadata$disease_stage <- disease

sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

# keep only healthy / myeloma samples
mask <- sample_metadata$disease_stage %in% c("Healthy", "MM")

#table(mask)
# mask
# FALSE  TRUE
#    31    25

expr_dat <- expr_dat[, c(TRUE, mask)]
sample_metadata <- sample_metadata[mask, ]

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
