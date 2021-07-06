#!/bin/env/Rscript
#
# Molecular classification of multiple myeloma
#
# Agnelli et al. (2005)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("util/eset.R")

# GEO accession
accession <- 'GSE2912'

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

# metadata stored separately (source: Agnelli et al, 2005, Appendix A)
mdat <- read.csv('/data/raw/agnelli2005/Agnelli2005.csv')

patient_ids <- str_match(pData(eset)$title, 'MM-[0-9]+')

# reoder patient metadata
mdat <- mdat[match(mdat$Patient, patient_ids), ]

# get relevant sample metadata
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

GENDER_IND <- 3
AGE_IND <- 5
STAGE_IND <- 6

sample_metadata$gender <- mdat[, GENDER_IND]
sample_metadata$age <- mdat[, AGE_IND]
sample_metadata$mm_stage <- mdat[, STAGE_IND]

# add cell type and disease stage (same for all samples)
sample_metadata$disease_stage <- 'MM'
sample_metadata$cell_type <- 'CD138+'

sample_metadata$platform_type <- 'Microarray'

expr_dat <- process_eset(eset)

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
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
