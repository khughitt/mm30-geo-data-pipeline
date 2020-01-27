#!/bin/env/Rscript
#
# Gene Expression of Bone Marrow Plasma Cells from Healthy Donors (N=22), MGUS (N=44),
# and Smoldering Myeloma (N=12)
#
# Zhan et al. (2006)
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE5900'

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
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
#eset <- eset[apply(exprs(eset), 1, var, na.rm = TRUE) > 0, ]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

group <- pData(eset)$source_name_ch1

sample_metadata$mm_stage <- rep('Healthy', nrow(sample_metadata))
sample_metadata$mm_stage[grepl('MGUS', group)] <- 'MGUS'
sample_metadata$mm_stage[grepl('smoldering', group)] <- 'SMM'

sample_metadata$disease <- rep('Healthy', nrow(sample_metadata))
sample_metadata$disease[!grepl('healthy', group)] <- 'Multiple Myeloma'


# add cell type and disease (same for all samples)
sample_metadata$cell_type <- 'CD138+'

table(sample_metadata$mm_stage)
#
# Healthy    MGUS     SMM
#      22      44      12
#

# get gene symbols
symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
  add_column(symbol = symbols, .before = 1)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  separate_rows(symbol, sep = " ?//+ ?") %>%
  group_by(symbol) %>%
  summarize_all(median)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_outfile_nr <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_outfile_nr))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
