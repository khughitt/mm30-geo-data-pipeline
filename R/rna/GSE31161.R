#!/bin/env/Rscript
#
# Identification of multiple risk loci and regulatory mechanisms influencing
# susceptibility to multiple myeloma
#
# Went et al. (2018)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("../util/eset.R")

# GEO accession
accession <- 'GSE31161'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.1', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# there are three samples in GSE31161 with all missing data..
num_missing <- apply(exprs(eset), 2, function(x) {
  sum(is.na(x))
})

#table(num_missing)
# num_missing
#     0 54675
#  1035     3

eset <- eset[, num_missing == 0]

# exclude outlier samples;
# these samples were found to have very median pairwise correlations (0.38 - 0.68) 
# compared to all other samples (>0.9 for most)
exclude_samples <- c("GSM771497", "GSM772341", "GSM772335")

eset <- eset[, !colnames(eset) %in% exclude_samples]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# columns to include (GSE31161)
sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id,
         treatment = `treatment:ch1`, time_of_testing = `time of testing:ch1`) %>%
  mutate(relapsed = time_of_testing == 'relapse')

# add disease stage
sample_metadata$disease_stage <- ifelse(sample_metadata$relapsed, 'RRMM', 'MM')

# add platform and cell type (same for all samples)
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
expr_outfile_nr <- file.path(processed_data_dir, sprintf('%s_gene_expr_nr.feather', accession))
mdat_outfile <- file.path(processed_data_dir, sprintf('%s_sample_metadata.tsv', accession))

write_feather(expr_dat, expr_outfile)
write_feather(expr_dat_nr, expr_outfile_nr)
write_tsv(sample_metadata, mdat_outfile)
