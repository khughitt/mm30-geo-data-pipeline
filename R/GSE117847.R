#!/bin/env/Rscript
#
# A retained transcriptomic profile characterizes the CD138+ cells in the progression
# from smoldering to active multiple myeloma
#
# Storti et al. (2019)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE117847'

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

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# columns to include (GSE31161)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, mm_stage = `patient diagnosis:ch1`) %>%
  mutate(mm_stage = recode(mm_stage, `active MM` = 'MM'))

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'CD138+'

# map from ensgenes to gene symbols
ensgenes <- fData(eset)$SPOT_ID
gene_symbols <- grch38$symbol[match(ensgenes, grch38$ensgene)]

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = gene_symbols, .after = 1)

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  select(-probe_id) %>%
  group_by(symbol) %>%
  summarize_all(median)

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_outfile_nr <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_outfile_nr)
write_tsv(sample_metadata, mdat_outfile)

sessionInfo()
