#!/bin/env/Rscript
#
# Identification of PIKfyve kinase as a target in multiple myeloma
#
# Bonolo et al. (2019)
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE134598'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
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

expr_dat <- read_tsv(supp_file) %>%
  select(-Chromosome, -Start, -End, -Length, -GeneBiotype) %>%
  rename(ensgene = GeneId, symbol = GeneName)

# in order to normalize downstream comparisons across datasets, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
expr_dat[, -(1:2)] <- sweep(expr_dat[, -(1:2)], 2, colSums(expr_dat[, -(1:2)]), '/') * 1E6

# drop empty rows
expr_dat <- expr_dat[rowSums(expr_dat[, -(1:2)]) > 0, ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         cell_line = `cell line:ch1`,
         sample_type = `sample type:ch1`,
         treatment = `treatment:ch1`, dose = `dose:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  select(-ensgene) %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_nr_outfile <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_nr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))

sessionInfo()
