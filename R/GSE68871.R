#!/bin/env/Rscript
#
# Expression data from BM-CD138+, obtained from newly diagnosed Multiple Myeloma
# patients [response to VTD therapy]
#
# Terragna et al. (2016)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE68871'

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
eset <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)[[1]]

# perform size-factor adjustment
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         treatment_response = `response to vtd therapy:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  add_column(symbol = fData(eset)$`Gene symbol`, .before = 1)

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

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_nr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_nr_outfile)
write_tsv(sample_metadata, mdat_outfile)
