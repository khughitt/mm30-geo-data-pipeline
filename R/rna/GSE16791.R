#!/bin/env/Rscript
#
# Expression data from CD138+ cells obtained from MM patients at diagnosis
#
# Terragna et al. (2009)
#
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE16791'

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

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, mm_stage = `stage:ch1`, age = `age:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# get expression data and add gene symbol column
gene_symbols <- gsub('///', ' // ', fData(eset)$`Gene symbol`)

expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = gene_symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
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
