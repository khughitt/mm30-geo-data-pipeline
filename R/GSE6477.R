#!/bin/env/Rscript
#
# Expression data from different stages of plasma cell neoplasm
#
# Chng et al. (2007)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE6477'

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

# perform size-factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# exclude control sequences present in some datasets (GSE6477)
eset <- eset[!startsWith(rownames(eset), 'AFFX-'), ]

# exclude any probes with zero variance (uninformative)
eset <- eset[apply(exprs(eset), 1, var) > 0, ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id, 
         title, ploidy = characteristics_ch1,
         ch13_status = characteristics_ch1.1) %>%
  mutate(mm_stage = ifelse(grepl('Normal', title), 'Normal',
                    ifelse(grepl('New', title), 'New',
                    ifelse(grepl('MGUS', title), 'MGUS',
                    ifelse(grepl('Relapsed', title), 'Relapsed',
                    ifelse(grepl('Smoldering', title), 'Smoldering', 'Normal')))))) %>%
  select(-title)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'BM-CD138+'

# load expression data and add gene symbol column
expr_dat <- exprs(eset) %>%
  as.data.frame() %>%
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

write_feather(expr_dat, path = expr_outfile)
write_feather(expr_dat_nr, path = expr_nr_outfile)
write_tsv(sample_metadata, path = mdat_outfile)
