#!/bin/env/Rscript
#
# MAQC-II Project: Multiple myeloma (MM) data set
#
# Popovici et al. (2010)
#
# Note: this dataset appears to use the same samples from GSE2658, but processed in
# a different manner, and with different metadata.
#
# This verion of the dataset (MAQC-II) includes ~2x samples before filtering and also
# includes survival-related metadata.
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("../util/eset.R")

options(stringsAsFactors = FALSE)

# GEO accession
accession <- 'GSE24080'

# directory to store raw and processed data
raw_data_dir <- file.path('/data/raw/geo/3.1', accession)
processed_data_dir <- sub('raw', 'clean', raw_data_dir)

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data;
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         sample_name = title,
         maqc_status = `maqc_distribution_status:ch1`) 

sample_metadata$sample_name <- as.character(sample_metadata$sample_name)

# add platform, cell type and disease stage (same for all samples)
sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# load supplemental clinical metadata;
# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24080/suppl/GSE24080%5FMM%5FUAMS565%5FClinInfo%5F27Jun2008%5FLS%5Fclean%2Exls%2Egz
clinical_metadata <- read_tsv('../../supp/clean/GSE24080_MM_UAMS565_ClinInfo_27Jun2008_LS_clean.tsv', 
                              col_types = cols())

clinical_metadata <- clinical_metadata %>%
  mutate(sample_name = sub('.CEL', '', cel_filename)) %>%
  select(-cel_filename)

# add to sample metadata table
sample_metadata <- sample_metadata %>%
  inner_join(clinical_metadata, by = 'sample_name')

# extract gene expression data
expr_dat <- process_eset(eset)

# remove five samples with MAQC_Remove flag
# MAQC_Remove    Training  Validation
#          5         340         214
bad_samples <- sample_metadata %>%
  filter(maqc_status == 'MAQC_Remove') %>%
  pull(geo_accession)

sample_metadata <- sample_metadata %>%
  filter(maqc_status != 'MAQC_Remove')

expr_dat <- expr_dat[, !colnames(expr_dat) %in% bad_samples]

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
  group_by(symbol) %>%
  summarize_all(median)

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
expr_nr_outfile <- sprintf('%s_gene_expr_nr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

print("Final dimensions:")
print(paste0("- Num rows: ", nrow(expr_dat_nr)))
print(paste0("- Num cols: ", ncol(expr_dat_nr)))

# store cleaned expression data and metadata
write_feather(expr_dat, file.path(processed_data_dir, expr_outfile))
write_feather(expr_dat_nr, file.path(processed_data_dir, expr_nr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))
