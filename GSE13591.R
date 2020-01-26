#!/bin/env/Rscript
#
# Integrated genomics approach to detect allelic imbalances in multiple myeloma#
#
# Agnelli et al. (2009)
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE13591'

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

# in order to normalize downstream comparisons across datasets, we will
# aply a size-factor normalization so that the sample sizes all sum to exactly the
# same amount..
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# get relevant sample metadata
# "treatment_protocol_ch1" is another alias for death;
# "grow_protocol_ch1" is an alias for relapse;
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, sample_type = `sample type:ch1`)

sample_metadata$mm_stage <- sample_metadata$sample_type
sample_metadata$mm_stage[startsWith(sample_metadata$mm_stage, 'TC')] <- 'MM'
sample_metadata$mm_stage[startsWith(sample_metadata$mm_stage, 'N')] <- 'Healthy'

#table(sample_metadata$sample_type)
# 
# MGUS    N  PCL  TC1  TC2  TC3  TC4  TC5 
#   11    5    9   29   25   49   24    6 

# table(sample_metadata$mm_stage)
# 
# MGUS   MM    N  PCL 
#   11  133    5    9 

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$disease[sample_metadata$mm_stage == 'Healthy'] <- 'Healthy'

sample_metadata$cell_type <- 'BM-CD138+'

# normalize divider used for multi-mapped probe entries
gene_symbols <- fData(eset)$`Gene symbol`
gene_symbols <- gsub('///', ' // ', gene_symbols)

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = gene_symbols, .after = 1)

# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  filter(symbol != '') %>%
  select(-probe_id) %>%
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
