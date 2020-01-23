#!/bin/env/Rscript
#
# Molecular Signatures of Multiple Myeloma Progression through Single Cell RNA-Seq
#
# Jang et al. (2019)
#
library(GEOquery)
library(tidyverse)
library(feather)

# GEO accession
accession <- 'GSE118900'

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

# download GEO data
#
# eset for this dataset only include metadata; expression data is empty and will be
# loaded separately..
#
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

supp_file <- file.path(raw_data_dir, accession, 'GSE118900_MM.scrna-seq.tpm.pass.txt.gz')

if (!file.exists(supp_file)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir, filter_regex = 'tpm.pass')
}

expr <- read.delim(gzfile(supp_file), row.names = 1)

# replace rownames with gene symbol

#head(rownames(expr))
# [1] "AADACL3|chr1|12776118" "AADACL4|chr1|12704566" "ABCA4|chr1|94458394"   "ABCB10|chr1|229652329"
# [5] "ABCD3|chr1|94883933"   "ABL2|chr1|179068462"

symbols <- str_split(rownames(expr), '\\|', simplify = TRUE)[, 1]

# remove small number of duplicated gene symbol entries
mask <- !duplicated(symbols)

expr <- expr[mask, ]
rownames(expr) <- symbols[mask]

# perform size-factor normalization
expr <- sweep(expr, 2, colSums(expr), '/') * 1E6

# exclude any zero variance genes present
row_vars <- apply(expr, 1, var)
expr <- expr[row_vars != 0, ]

# columns to include (GSE118900)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         patient = `patient:ch1`,
         mm_stage = `tumor stage:ch1`)

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- 'CD138+'

expr <- expr %>%
  rownames_to_column('symbol')

# determine filenames to use for outputs and save to disk
expr_outfile <- sprintf('%s_gene_expr.feather', accession)
mdat_outfile <- sprintf('%s_sample_metadata.tsv', accession)

# store cleaned expression data and metadata
write_feather(expr, file.path(processed_data_dir, expr_outfile))
write_tsv(sample_metadata, file.path(processed_data_dir, mdat_outfile))

sessionInfo()
