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
library(arrow)
source("R/util/eset.R")

# GEO accession
accession <- 'GSE68871'

# directory to store raw and processed data
raw_data_dir <- file.path('/data/raw', accession)
processed_data_dir <- sub('raw', 'clean', raw_data_dir)

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
         treatment_response = `response to vtd therapy:ch1`)

# add cell type and disease stage (same for all samples)
sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'CD138+'

sample_metadata$platform_type <- 'Microarray'

# extract gene expression data
expr_dat <- process_eset(eset)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  group_by(symbol) %>%
  summarize_all(median)

# store cleaned expression data and metadata
expr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr.feather', accession))
expr_nr_outfile <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

print("Final dimensions:")
print(paste0("- Num rows: ", nrow(expr_dat_nr)))
print(paste0("- Num cols: ", ncol(expr_dat_nr)))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_nr_outfile)
write_tsv(sample_metadata, mdat_outfile)
