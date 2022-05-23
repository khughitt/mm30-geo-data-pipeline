#!/bin/env/Rscript
###############################################################################
#
# GSE9782
#
###############################################################################
library(tidyverse)

# load data & metadata
dat <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

fdata <- read_csv(snakemake@input[[2]], show_col_types = FALSE)
pdata <- read_csv(snakemake@input[[3]], show_col_types = FALSE)

# size factor normalization (ignore feature column)
dat[, -1] <- sweep(dat[, -1], 2, colSums(dat[, -1]), '/') * 1E6

# add gene symbol column
expr_dat <- dat %>%
  select(-feature) %>%
  add_column(symbol = fdata$`Gene Symbol`, .before = 1) %>%
  as_tibble()

mask <- expr_dat$symbol != ""
expr_dat <- expr_dat[mask, ]
fdata <- fdata[mask, ]

# split multi-mapped symbols
expr_dat <- expr_dat %>%
  separate_rows(symbol, sep = " ?//+ ?")

#
# Note: some of the metadata fields appear to be incorrectly encoded, e.g.:
#
# > pdata[, 'characteristics_ch1.9']
# [1] PGx_Days_To_Progression = 87           PGx_Progression(0=No,1=Yes) = 1
# [3] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [5] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [7] PGx_Progression(0=No,1=Yes) = 1        PGx_Progression(0=No,1=Yes) = 1
# [9] PGx_CensorReason = No more data        PGx_Progression(0=No,1=Yes) = 1
#
# For the case of ch1.9, the column contains a mix of "pgx progression",
# "days to progress", and "censor" fields.
#

# ch1.7 - PGx_Response   (-> treatment_response)
# ch1.8 - PGx_Responder  (-> patient_subgroup)

# to get around this, we will manually detect and parse those fields..
sample_metadata <- pdata %>%
  mutate(
    patient_id         = title,
    study_code         = as.numeric(str_extract(characteristics_ch1, '\\d+$')),
    treatment          = str_extract(characteristics_ch1.1, '\\w+$'),
    gender             = str_extract(characteristics_ch1.2, '\\w+$'),
    ethnicity          = str_extract(characteristics_ch1.3, '\\w+$'),
    age                = as.numeric(str_extract(characteristics_ch1.4, '\\d+$')),
    treatment_response = str_extract(characteristics_ch1.7, '\\w+$'),
    patient_subgroup   = str_extract(characteristics_ch1.8, '\\w+$')) %>%
  select(geo_accession, patient_id, platform_id, study_code, title, treatment, gender,
          ethnicity, age, treatment_response, patient_subgroup)

pfs_event <- c()
pfs_time <- c()
pfs_event_reason <- c()
patient_died <- c()
os_time <- c()

pdata_df <- as.data.frame(pdata)

# iterate over samples and find relevant information
for (sample_num in 1:nrow(pdata_df)) {
  # PGx Progression
  ind <- which(grepl('PGx_Prog', pdata_df[sample_num, ]))

  if (length(ind) == 0) {
    pfs_event <- c(pfs_event, 0)
  } else {
    pfs_event <- c(pfs_event, ifelse(endsWith(pdata_df[sample_num, ind], '1'), 1, 0))
  }

  # PGx days
  ind <- which(grepl('PGx_Days', pdata_df[sample_num, ]))

  if (length(ind) == 0) {
    pfs_time <- c(pfs_time, NA)
  } else {
    pfs_time <- c(pfs_time, as.numeric(str_extract(pdata_df[sample_num, ind], '\\d+$')))
  }

  # PGx Censor Reason
  ind <- which(grepl('PGx_CensorReason', pdata_df[sample_num, ]))

  if (length(ind) == 0) {
    pfs_event_reason <- c(pfs_event_reason, NA)
  } else {
    pfs_event_reason <- c(pfs_event_reason, 
                          trimws(str_extract(pdata_df[sample_num, ind], '[ \\w]+$')))
  }

  # Deceased
  ind <- which(grepl('Did_Patient_Die', pdata_df[sample_num, ]))

  if (length(ind) == 0) {
    patient_died <- c(patient_died, 0)
  } else {
    patient_died <- c(patient_died, 1)
  }

  # Days Survived
  ind <- which(grepl('Days_Survived', pdata_df[sample_num, ]))
  os_time <- c(os_time, as.numeric(str_extract(pdata_df[sample_num, ind], '\\d+$')))
}

# add to sample metadata
sample_metadata$pfs_time <- pfs_time
sample_metadata$pfs_event <- pfs_event
sample_metadata$pfs_event_reason <- pfs_event_reason
sample_metadata$os_time <- os_time
sample_metadata$patient_died <- patient_died

sample_metadata$disease_stage <- 'RRMM'
sample_metadata$platform_type <- 'Microarray'

if (!all(colnames(expr_dat)[-1] == sample_metadata$geo_accession)) {
  stop("Sample ID mismatch!")
}

# store results
write_csv(expr_dat, snakemake@output[[1]])
write_csv(fdata, snakemake@output[[2]])
write_csv(sample_metadata, snakemake@output[[3]])
