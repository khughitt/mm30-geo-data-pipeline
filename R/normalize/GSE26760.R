#!/bin/env/Rscript
###############################################################################
#
# GSE26760
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
  select(geo_accession, platform_id, title, patient_id,
         age = `Age at Diagnosis`, sex = Gender,
         race = Race, disease_stage = Diagnosis)

# add platform & disease stage
sample_metadata$platform_type <- "Microarray"
sample_metadata$sample_type <- "Patient"

# drop samples with unknown stage and normalize stage names
sample_metadata <- sample_metadata %>%
  filter(disease_stage != "Unknown") %>%
  mutate(disease_stage = recode(disease_stage, 
                           `Multiple Myeloma` = "MM", 
                           `Primary Plasma Cell Leukemia` = "PCL",
                           `Smoldering Myeloma` = "SMM"))

# update expression data to match samples & order in metadata
expr_dat <- expr_dat[, c("symbol", sample_metadata$geo_accession)] 

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# update feature annotations
fdata <- grch38[match(expr_dat$symbol, grch38$symbol), ]

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
