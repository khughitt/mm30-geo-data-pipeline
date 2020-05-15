#!/bin/env/Rscript
#
#	Gene Expression Profiles of Multiple Myeloma (N=414) Before Treatment
#
# NOTE: samples in this series are included as part of GSE24080 (MAQC-II) and as such,
# this script is no longer needed
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE4581'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/1.1', accession)

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
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# GSE4581 includes a small number (131) of missing values all relating to a single
# probe: 1552256_a_at / SCARB1
eset <- eset[complete.cases(exprs(eset)), ]

# exclude control sequences
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# sample size range
#range(colSums(exprs(eset)))
# [1] 32986777 54993524

# extract sample metadata including survival time and censor status from metadata
# $ title                     <chr> "U133Plus-P002 (TT2 pre-treatment)"
# $ characteristics_ch1       <chr> "[SURIND=0 (Indicator of disease-related death; integer, 0=alive or death by other cause, 1=disease related death, na=death cause undetermined)]"
# $ characteristics_ch1.2     <chr> "[SURTIM=69.24 (Follow-up time in months from Pre-Treatment baseline; integer)]"
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, characteristics_ch1, characteristics_ch1.2, title) %>%
  mutate(condition = gsub('\\(|\\)', '', str_match(title, '\\(.*\\)')))

sample_metadata$os_censor <- as.numeric(sub('=', '', str_match(sample_metadata$characteristics_ch1, '=[01]')))
sample_metadata$os_time <- as.numeric(sub('=', '', str_match(sample_metadata$characteristics_ch1.2, '=[0-9\\.]+')))

sample_metadata <- sample_metadata %>%
  select(-characteristics_ch1, -characteristics_ch1.2, -title)

#
# conditions present:
#
# 1. CD-138-selected plasma cells from bone marrow of patients with newly diagnosed
#    multiple myeloma subsequently treated with high dose therapy and stem cell
#    transplants termed Total Therapy 2 (pre-treatment TT2)
# 2. CD-138-selected plasma cells from bone marrow of patients with newly diagnosed
#    multiple myeloma subsequently treated with Total Therapy 3 (pre treatment TT3).
#

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# load expression data and add gene symbol column
expr_dat <- as.data.frame(exprs(eset))

expr_dat <- expr_dat %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = fData(eset)$`Gene symbol`, .after = 1)

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  select(-probe_id) %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_nr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, path = expr_outfile)
write_feather(expr_dat_nr, path = expr_nr_outfile)
write_tsv(sample_metadata, path = mdat_outfile)

sessionInfo()
