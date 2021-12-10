#!/bin/env/Rscript
#
# MMRC expression reference collection
#
# Chapman et al. (2011)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("R/util/eset.R")

# GEO accession
accession <- 'GSE26760'

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
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

# size factor normalization
exprs(eset) <- sweep(exprs(eset), 2, colSums(exprs(eset)), '/') * 1E6

# "Sample from patient MMRC0091" -> "MMRC0091"
patient_ids <- str_split(pData(eset)$source_name_ch1, ' ', simplify = TRUE)[, 4]

sample_metadata <- pData(eset) %>%
  select(geo_accession, platform_id)

sample_metadata$patient_id <- patient_ids

# add cell type (same for all samples)
sample_metadata$cell_type <- 'CD138+'

sample_metadata$platform_type <- 'Microarray'

# add additional metadata from
# http://portals.broadinstitute.org/mmgp/data/browseData?conversationPropagation=begin
mdat <- read_csv('supp/clean/mmrc.sample.information.csv', col_types = cols()) %>%
  select(patient_id = Array, age = `Age at Diagnosis`, gender = Gender,
         race = Race, mm_stage = Diagnosis)

sample_metadata <- sample_metadata %>%
  inner_join(mdat, by = 'patient_id') %>%
  filter(mm_stage != 'Unknown')

sample_metadata$mm_stage[grepl('Multiple Myeloma', sample_metadata$mm_stage)] <- 'MM'
sample_metadata$mm_stage[grepl('Leukemia', sample_metadata$mm_stage)] <- 'PCL'
sample_metadata$mm_stage[sample_metadata$mm_stage == 'Smoldering Myeloma'] <- 'SMM'

#table(sample_metadata$mm_stage)
#
# MGUS   MM  PCL  SMM
#    2  224    3   10
#

# drop samples with no available metadata and normalize order
eset <- eset[, sample_metadata$geo_accession] 

if (!all(colnames(eset) == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# drop PCL samples
mask <- sample_metadata$mm_stage != 'PCL'

sample_metadata <- sample_metadata[mask, ]
eset <- eset[, mask]

sample_metadata$disease_stage <- sample_metadata$mm_stage

# extract gene expression data
expr_dat <- process_eset(eset)

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
