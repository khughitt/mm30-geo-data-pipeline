#!/bin/env/Rscript
#
# Expression data from different stages of plasma cell neoplasm
#
# Chng et al. (2007)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("util/eset.R")

# GEO accession
accession <- 'GSE6477'

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
         title, ploidy = characteristics_ch1,
         ch13_status = characteristics_ch1.1) %>%
  mutate(mm_stage = ifelse(grepl('Normal', title), 'Normal',
                    ifelse(grepl('New', title), 'New',
                    ifelse(grepl('MGUS', title), 'MGUS',
                    ifelse(grepl('Relapsed', title), 'Relapsed',
                    ifelse(grepl('Smoldering', title), 'Smoldering', 'Normal')))))) %>%
  select(-title)

# add disease stage
sample_metadata <- sample_metadata %>%
  mutate(disease_stage = recode(mm_stage, Normal = 'Healthy', New = 'MM',
                                Relapsed = 'RRMM', Smoldering = 'SMM'))

# add cell and platform type
sample_metadata$cell_type <- 'CD138+'
sample_metadata$platform_type <- 'Microarray'

# exclude outlier samples with low median pairwise correlations (0.63-0.64 vs. ~0.9 for
# most other samples)
exclude_samples <- c("GSM149035", "GSM149037")

eset <- eset[, !colnames(eset) %in% exclude_samples]

sample_metadata <- sample_metadata %>%
  filter(!geo_accession %in% exclude_samples)

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
