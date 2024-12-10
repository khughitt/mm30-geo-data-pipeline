#!/bin/env/Rscript
###############################################################################
#
# GSE118900
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# load data & metadata
expr_dat <- read_feather(snakemake@input[[1]])
pdata <- read_feather(snakemake@input[[3]])

# columns to include
sample_metadata <- pdata %>%
  select(geo_accession, title,  platform_id,
         sample_name = description,
         disease_stage = `tumor stage:ch1`,
         patient = `patient:ch1`,
         disease_stage = `tumor stage:ch1`) %>%
  mutate(disease_stage = recode(disease_stage, `IgM-MGUS` = "MGUS", `NDMM` = "MM"))

# add cell and platform type
sample_metadata$platform_type <- "RNA-Seq"
sample_metadata$sample_type <- "Patient"

# normalize sample ids, ex: "IgM-MGUS1_C37" -> "IgM.MGUS1_C37"
sample_metadata$patient <- gsub("-", ".", sample_metadata$patient)
sample_metadata$sample_name <- gsub("-", ".", sample_metadata$sample_name)

# drop single outlier patient NDMM7 whose samples had a low correlation with all other
# patient samples (median pairwise correlation ~0.06 vs. 0.6 for other patients)
exclude_patient <- "NDMM7"

exclude_samples <- sample_metadata %>%
  filter(patient == exclude_patient) %>%
  pull(geo_accession)

expr_dat <- expr_dat[, !startsWith(colnames(expr_dat), exclude_patient)]

sample_metadata <- sample_metadata %>%
  filter(patient != exclude_patient)

# combine samples within each patient
patient_ids <- unique(sample_metadata$patient)

expr_list <- list(
  "symbol" = expr_dat$symbol
)

for (patient_id in patient_ids) {
  patient_samples <- sample_metadata %>%
    filter(patient == patient_id) %>%
    pull(geo_accession)

  expr_list[[patient_id]] <- rowSums(expr_dat[, patient_samples])
}
expr_dat <- data.frame(expr_list)

# exclude any zero variance genes present
mask <- apply(expr_dat[, -1], 1, var) > 0
expr_dat <- expr_dat[mask, ]

# adjust sample metadata to be at the patient level
sample_metadata <- sample_metadata %>%
  group_by(patient) %>%
  slice(1) %>%
  select(-geo_accession) %>%
  ungroup()

# normalize expr data / metadata order
expr_dat <- expr_dat[, c("symbol", sample_metadata$patient)]

if (!all(sample_metadata$patient == colnames(expr_dat[, -1]))) {
  stop("patient id mismatch!")
}

sample_metadata <- sample_metadata %>%
  select(patient, everything())

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
