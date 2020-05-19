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

# GEO accession
accession <- 'GSE9782'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/2.0', accession)

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
esets <- getGEO(accession, destdir = raw_data_dir, AnnotGPL = TRUE)

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
    study_code         = as.numeric(str_match(characteristics_ch1, '\\d+$')),
    treatment          = str_match(characteristics_ch1.1, '\\w+$'),
    gender             = str_match(characteristics_ch1.2, '\\w+$'),
    ethnicity          = str_match(characteristics_ch1.3, '\\w+$'),
    age                = as.numeric(str_match(characteristics_ch1.4, '\\d+$')),
    treatment_response = str_match(characteristics_ch1.7, '\\w+$'),
    patient_subgroup   = str_match(characteristics_ch1.8, '\\w+$')) %>%
  select(geo_accession, platform_id, study_code, treatment, gender,
          ethnicity, age, treatment_response, patient_subgroup) %>%
  add_column(geo_accession2 = pData(esets[[2]])$geo_accession, .after = 1)

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

# add cell type and disease (same for all samples)
sample_metadata$disease <- 'Multiple Myeloma'
sample_metadata$cell_type <- NA # likely CD138+, but not 100% certain

# for GSE7039, data includes two separate esets with a small number of overlapping probes
e1 <- exprs(esets[[1]])
e2 <- exprs(esets[[2]])

symbols1 <- fData(esets[[1]])[, 'Gene symbol']
symbols2 <- fData(esets[[2]])[, 'Gene symbol']

mask1 <- !startsWith(rownames(e1), 'AFFX-')
mask2 <- !startsWith(rownames(e2), 'AFFX-')

e1 <- e1[mask1, ]
e2 <- e2[mask2, ]

symbols1 <- symbols1[mask1]
symbols2 <- symbols2[mask2]

# small number of shared probes exist...
#length(intersect(rownames(e1), rownames(e2)))
# [1] 168

#head(intersect(rownames(e1), rownames(e2)))
# [1] "200000_s_at" "200001_at"   "200002_at"   "200003_s_at" "200004_at"   "200005_at"

# use average values for those probes..
shared_probes <- intersect(rownames(e1), rownames(e2))

e1[shared_probes, ] <- (e1[shared_probes, ] + e2[shared_probes, ]) / 2

mask <- !rownames(e2) %in% shared_probes
e2 <- e2[mask, ]
symbols2 <- symbols2[mask]

expr_dat <- rbind(e1, e2)
symbols <- c(symbols1, symbols2)

# size factor normalization
expr_dat <- sweep(expr_dat, 2, colSums(expr_dat), '/') * 1E6

# get expression data and add gene symbol column
expr_dat <- expr_dat %>%
  as.data.frame() %>%
  add_column(symbol = symbols, .before = 1) %>%
  filter(symbol != '')

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# create a version of gene expression data with a single entry per gene, including
# only entries which could be mapped to a known gene symbol
expr_dat_nr <- expr_dat %>%
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
