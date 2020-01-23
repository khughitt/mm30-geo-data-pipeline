#!/bin/env/Rscript
#
# Identification of PIKfyve kinase as a target in multiple myeloma
#
# Bonolo et al. (2019)
#
library(GEOquery)
library(tidyverse)

# GEO accession
accession <- 'GSE134598'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
clean_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, clean_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
# result is a list with a single entry containing an ExpressionSet instance
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# eset is missing expression data; download separately..
supp_file <- file.path(raw_data_dir, accession, 'GSE134598_All_Processed_Data.txt.gz')

if (!file.exists(supp_file)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir)[[1]]
}

expr <- read_tsv(supp_file) %>%
  select(-Chromosome, -Start, -End, -Length, -GeneBiotype)

# in order to normalize downstream comparisons across datasets, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
expr[, -(1:2)] <- sweep(expr[, -(1:2)], 2, colSums(expr[, -(1:2)]), '/') * 1E6

# drop empty rows
expr <- expr[rowSums(expr[, -(1:2)]) > 0, ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         cell_line = `cell line:ch1`,
         sample_type = `sample type:ch1`,
         treatment = `treatment:ch1`, dose = `dose:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# store cleaned expression data and metadata
expr_outfile <- file.path(clean_data_dir, sprintf('%s_gene_expr.csv', accession))
mdat_outfile <- file.path(clean_data_dir, sprintf('%s_sample_metadata.csv', accession))

write_csv(expr, expr_outfile)
write_csv(sample_metadata, mdat_outfile)

sessionInfo()
