#!/bin/env/Rscript
#
# Gene expression profiling and correlation with outcome in clinical trials of the
# proteasome inhibitor bortezomib
#
# Mulligan et al. (2007)
#
library(GEOquery)
library(tidyverse)
library(arrow)
source("../util/eset.R")

# GEO accession
accession <- 'GSE9782'

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

# download GEO data;
# for GSE9782, data includes two separate esets with a small number of overlapping
# probes
esets <- getGEO(accession, destdir = raw_data_dir)

#
# Note: some of the metadata fields appear to be incorrectly encoded, e.g.:
#
# > pData(eset)[, 'characteristics_ch1.9']
# [1] PGx_Days_To_Progression = 87           PGx_Progression(0=No,1=Yes) = 1
# [3] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [5] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [7] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [9] PGx_CensorReason = No more data        PGx_Progression(0=No,1=Yes) = 1
#
# For the case of ch1.9, the column contains a mix of "pgx progression",
# "days to progress", and "censor" fields.
#

#
# ch1.7 - PGx_Response   (-> treatment_response)
# ch1.8 - PGx_Responder  (-> patient_subgroup)
#

# to get around this, we will manually detect and parse those fields..
sample_metadata <- pData(esets[[1]]) %>%
  mutate(
    patient_id         = title,
    study_code         = as.numeric(str_match(characteristics_ch1, '\\d+$')),
    treatment          = str_match(characteristics_ch1.1, '\\w+$'),
    gender             = str_match(characteristics_ch1.2, '\\w+$'),
    ethnicity          = str_match(characteristics_ch1.3, '\\w+$'),
    age                = as.numeric(str_match(characteristics_ch1.4, '\\d+$')),
    treatment_response = str_match(characteristics_ch1.7, '\\w+$'),
    patient_subgroup   = str_match(characteristics_ch1.8, '\\w+$')) %>%
  select(geo_accession, patient_id, platform_id, study_code, treatment, gender,
          ethnicity, age, treatment_response, patient_subgroup) %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1)

# sanity check (comparing patient ids for two 133A/B to make sure they match)
if (!all(pData(esets[[1]])$title == pData(esets[[2]])$title)) {
  stop("Sample mismatch!")
}

# convert metadata from factor to character columns for parsing
str_mdat <- pData(esets[[1]])

for (cname in colnames(str_mdat)) {
  str_mdat[, cname] <- as.character(str_mdat[, cname])
}

pfs_event <- c()
pfs_time <- c()
pfs_event_reason <- c()
patient_died <- c()
os_time <- c()

# iterate over samples and find relevant information
for (sample_num in 1:nrow(pData(esets[[1]]))) {
  # PGx Progression
  ind <- which(grepl('PGx_Prog', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    pfs_event <- c(pfs_event, 0)
  } else {
    pfs_event <- c(pfs_event, ifelse(endsWith(str_mdat[sample_num, ind], '1'), 1, 0))
  }

  # PGx days
  ind <- which(grepl('PGx_Days', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    pfs_time <- c(pfs_time, NA)
  } else {
    pfs_time <- c(pfs_time, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
  }

  # PGx Censor Reason
  ind <- which(grepl('PGx_CensorReason', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
      pfs_event_reason <- c(pfs_event_reason, NA)
  } else {
    pfs_event_reason <- c(pfs_event_reason, trimws(str_match(str_mdat[sample_num, ind], '[ \\w]+$')))
  }

  # Deceased
  ind <- which(grepl('Did_Patient_Die', str_mdat[sample_num, ]))

  if (length(ind) == 0) {
    patient_died <- c(patient_died, 0)
  } else {
    patient_died <- c(patient_died, 1)
  }

  # Days Survived
  ind <- which(grepl('Days_Survived', str_mdat[sample_num, ]))
  os_time <- c(os_time, as.numeric(str_match(str_mdat[sample_num, ind], '\\d+$')))
}

# add to sample metadata
sample_metadata$pfs_time <- pfs_time
sample_metadata$pfs_event <- pfs_event
sample_metadata$pfs_event_reason <- pfs_event_reason
sample_metadata$os_time <- os_time
sample_metadata$patient_died <- patient_died

# add disease stage
sample_metadata$disease_stage <- 'RRMM'

# add cell type
sample_metadata$cell_type <- "CD138+"

sample_metadata$platform_type <- 'Microarray'

# GSE7039 includes two esets/arrays for each patient, Affy HG-U133A & HG-U133B, which
# overlap in a very small number of probes (n = 168).
# Below, the two arrays are processed separately following the usual approach, and then
# bound together at the end before deduplication.

# size factor normalization
exprs(esets[[1]]) <- sweep(exprs(esets[[1]]), 2, colSums(exprs(esets[[1]])), '/') * 1E6
exprs(esets[[2]]) <- sweep(exprs(esets[[2]]), 2, colSums(exprs(esets[[2]])), '/') * 1E6

# extract and process gene expression data from each, separately
expr_dat1 <- process_eset(esets[[1]])
expr_dat2 <- process_eset(esets[[2]])

colnames(expr_dat1) <- c('symbol', sample_metadata$geo_accession)
colnames(expr_dat2) <- c('symbol', sample_metadata$geo_accession)

# join datasets
expr_dat <- rbind(expr_dat1, expr_dat2)

# create a version of gene expression data with a single entry per gene, including
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
