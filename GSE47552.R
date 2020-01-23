#!/bin/env/Rscript
#
#	Transcriptome analysis reveals molecular profiles associated with evolving steps of
#	monoclonal gammopathies
#
# López-Corral et al. (2014)
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE47552'

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

# columns to include (GSE47552)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         mm_stage_raw = `cell type:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'CD138+'

sample_metadata$disease[grepl('NPC', sample_metadata$mm_stage_raw)] <- 'Healthy'

sample_metadata$mm_stage <- rep('MM', length(sample_metadata$mm_stage_raw))

sample_metadata$mm_stage[grepl('Normal', sample_metadata$mm_stage_raw)] <- 'Healthy'
sample_metadata$mm_stage[grepl('MGUS', sample_metadata$mm_stage_raw)] <- 'MGUS'
sample_metadata$mm_stage[grepl('SMM', sample_metadata$mm_stage_raw)] <- 'SMM'

sample_metadata <- sample_metadata %>%
  select(-mm_stage_raw)

# get gene symbols
symbols <- fData(eset)$`Gene symbol`

# get expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame %>%
  rownames_to_column('probe_id') %>%
  add_column(symbol = symbols, .after = 1)

# table(expr_dat$symbol == '')
#
# FALSE  TRUE
# 22195 11102

# drop unmapped genes
expr_dat <- expr_dat %>%
  filter(symbol != '')

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))

sessionInfo()
