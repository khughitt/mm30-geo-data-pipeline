#!/bin/env/Rscript
#
# Potent Anti-Myeloma Activity of the TOPK inhibitor OTS514 in Pre-Clinical Models
#
# Stefka et al. (2020)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("../util/eset.R")

# GEO accession
accession <- 'GSE128251'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.1', accession)

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

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         treatment = `treatment:ch1`, replicate = `replicate:ch1`,
         cell_line = source_name_ch1)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'H929'

sample_metadata$platform_type <- 'Microarray'

# extract gene expression data
expr_dat <- process_eset(eset)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create non-redundant version
expr_dat_nr <- expr_dat %>%
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
