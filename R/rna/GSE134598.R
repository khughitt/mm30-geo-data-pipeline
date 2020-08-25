#!/bin/env/Rscript
#
# Identification of PIKfyve kinase as a target in multiple myeloma
#
# Bonolo et al. (2019)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE134598'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.0', accession)

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

# get expression data;
expr_dat <- read_tsv(supp_file, col_types = cols()) %>%
  select(-Chromosome, -Start, -End, -Length, -GeneBiotype, -GeneName) 

# mistake in supplemental file when attempting to use the "GeneName" field (excel gene symbol issue...):
# > expr_dat_nr[1:3, 1:3]
# A tibble: 3 x 3                            
#   symbol KMS26_Control KMS26_APY0201_100nM
#   <chr>          <dbl>               <dbl>
# 1 1-Dec        0                   0
# 2 1-Mar        0.00229             0.00255
# 3 1-Sep        1.57                1.61
# (reported upstream jan 26, 2020)

# map from ensgene -> grch38
ind <- match(expr_dat$GeneId, grch38$ensgene)
gene_symbols <- grch38$symbol[ind]

expr_dat <- expr_dat %>%
  select(-GeneId) %>%
  add_column(symbol = gene_symbols, .before = 1)

# size factor normalization
expr_dat[, -1] <- sweep(expr_dat[, -1], 2, colSums(expr_dat[, -1]), '/') * 1E6

# drop entries that could not be mapped to gene symbols
expr_dat <- expr_dat[!is.na(expr_dat$symbol), ]

# drop empty rows
expr_dat <- expr_dat[rowSums(expr_dat[, -1]) > 0, ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         sample_name = title,
         cell_type = `cell line:ch1`,
         sample_type = `sample type:ch1`,
         treatment = `treatment:ch1`, dose = `dose:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$cell_type[is.na(sample_metadata$cell_type)] <- 'CD138+'
sample_metadata$disease_stage <- "MM"

if (!all(colnames(expr_dat)[-1] == sample_metadata$sample_name)) {
  stop("Sample ID mismatch!")
}

# for consistency, use GEO sample accessions as ids
colnames(expr_dat) <- c('symbol', sample_metadata$geo_accession)

# only entries which could be mapped to a known gene symbol
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
