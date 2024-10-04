#!/bin/env/Rscript
###############################################################################
#
# GSE19784
#
###############################################################################
library(annotables)
library(tidyverse)
library(arrow)

# output directory to store data packages to
out_dir <- dirname(snakemake@output[[1]])

# load data & metadata
dat <- read_feather(snakemake@input[[1]]) %>%
  column_to_rownames("feature")

fdata <- read_feather(snakemake@input[[2]])
pdata <- read_feather(snakemake@input[[3]])

# add gene symbol column
expr_dat <- dat %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

# drop rows with missing values
expr_dat <- expr_dat[complete.cases(expr_dat), ]

# get relevant sample metadata
sample_metadata <- pdata %>%
  select(geo_accession,
         platform_id,
         iss_stage = `iss:ch1`,
         patient_subgroup = `cluster:ch1`)

# add platform, cell type and disease (same for all samples)
sample_metadata$disease_stage <- "MM"
sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

# Note; there is not sufficient information provided to link patients in table S11 to
# GSM sample identifiers; skipping.
# table_s11 <- read_csv(file.path(base_dir, "metadata", "broyl2010_supp_table_s11.csv"))

# load survival metadata from Kuiper et al. (2012)
survival_mdata <- read_csv("supp/clean/kuiper2012_supp_patient_survival.csv", col_types = cols())

survival_mdata <- survival_mdata %>%
  rename(geo_accession = Patient) %>%
  filter(geo_accession %in% colnames(dat))

colnames(survival_mdata) <- c("geo_accession", "os_time", "os_censor", "pfs_time", "pfs_censor")

# exclude samples without metadata
mask <- sample_metadata$geo_accession %in% survival_mdata$geo_accession

# table(mask)
# mask
# FALSE  TRUE
#    46   282

expr_dat <- expr_dat[, c(TRUE, mask)]
sample_metadata <- sample_metadata[mask, ]

#all(colnames(dat) == sample_metadata$geo_accession)
# [1] TRUE

# combine metadata
sample_metadata <- sample_metadata %>%
  inner_join(survival_mdata, by = "geo_accession")

sample_metadata <- sample_metadata %>%
  mutate(os_censor = as.logical(os_censor),
         os_time = as.numeric(os_time),
         pfs_censor = as.logical(pfs_censor),
         pfs_time = as.numeric(pfs_time),
         iss_stage = factor(iss_stage, levels = c("I", "II", "III"), ordered = TRUE),
         patient_subgroup = as.factor(patient_subgroup))

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_feather(expr_dat, snakemake@output[[1]])
write_feather(fdata, snakemake@output[[2]])
write_feather(sample_metadata, snakemake@output[[3]])
