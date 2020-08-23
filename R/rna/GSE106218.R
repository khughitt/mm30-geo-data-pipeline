#!/bin/env/Rscript
#
# GSE106218
#
# Single cell RNA sequencing of multiple myeloma I
#
# Ryu et al. (2019)
#
library(annotables)
library(GEOquery)
library(tidyverse)
library(arrow)

# GEO accession
accession <- 'GSE106218'

# directory to store raw and processed data
base_dir <- file.path('/data/human/geo/3.0', accession)

raw_data_dir <- file.path(base_dir, 'raw')
processed_data_dir <- file.path(base_dir, 'processed')

# create output directories if they don't already exist
for (dir_ in c(raw_data_dir, processed_data_dir)) {
  if (!dir.exists(dir_)) {
      dir.create(dir_, recursive = TRUE)
  }
}

# expression data is missing from getGEO() query result and must be downloaded
# separately; further, some useful clinical metadata is available as a separate file
# which we will also download using the getGEOSuppFiles() function

# download GEO data
eset <- getGEO(accession, destdir = raw_data_dir)[[1]]

supp_file1 <- file.path(raw_data_dir, accession, 'GSE106218_GEO_processed_MM_raw_TPM_matrix.txt.gz')
supp_file2 <- file.path(raw_data_dir, accession, 'GSE106218_GEO_clinical_information_MM.txt.gz')

if (!file.exists(supp_file1)) {
  getGEOSuppFiles(accession, baseDir = raw_data_dir)
}

# load TPM counts
tpm_counts <- read.delim(gzfile(supp_file1), row.names = 1)

# load sample metadata
clinical_metadata <- as.data.frame(t(read.delim(gzfile(supp_file2), row.names = 1))) %>%
  rownames_to_column('patient_id') %>%
  select(patient_id,
         iss_stage = `ISS stage`, os_event = `Death(alive=0; death=1)`,
         os_time = `Survival time(EM;month)`,
         heavy_chain = `Heavy chain`,
         light_chain = `Light chain`) %>%
  mutate(os_event = os_event == "1",
         disease_stage = 'MM')

# drop the month component of survival time
clinical_metadata$os_time <- as.numeric(str_match(clinical_metadata$os_time, '[0-9]+'))

sample_metadata <- pData(eset) %>%
  select(patient_id = `patient id:ch1`,
         geo_accession, platform_id,
         disease_stage,
         gender = `gender:ch1`,
         prep_site = `prep-site:ch1`)

sample_metadata <- sample_metadata %>%
  inner_join(clinical_metadata, by = 'patient_id')

# tpm_counts
#         MM02_38 MM02_48 MM02_50
# 5S_rRNA       0       0       0
# 7SK           0       0       0
# A1BG          0       0       0

# clinical_metadata
#                         MM02 MM16 MM17
# Sex                        M    F    F
# ISS stage                 II    I    I
# Death(alive=0; death=1)    1    1    1

#table(rowSums(tpm_counts) == 0)
#
# FALSE  TRUE
# 35582 20278

# remove empty rows
tpm_counts <- tpm_counts[rowSums(tpm_counts) > 0, ]

sample_metadata$disease <- "Multiple Myeloma"
sample_metadata$cell_type <- 'CD138+'

# collapse patient samples
expr_patient_ids <- str_match(colnames(tpm_counts), 'MM[0-9]+')

# sort(table(expr_patient_ids))
# expr_patient_ids
# MM33 MM16 MM28 MM30 MM17 MM36 MM38 MM25 MM02 MM34
#   13   24   28   31   38   47   49   52   77  129

expr_dat <- NULL

for (patient_id in unique(expr_patient_ids)) {
  dat <- rowMeans(tpm_counts[, expr_patient_ids == patient_id])
  expr_dat <- cbind(expr_dat, dat)
}
colnames(expr_dat) <- unique(expr_patient_ids)

expr_dat <- expr_dat %>%
  as.data.frame() %>%
  rownames_to_column('symbol')

# update GRCh37 symbols -> GRCh38
grch37_mask <- expr_dat$symbol %in% grch37$symbol

gene_symbols <- expr_dat$symbol[grch37_mask]
ensgenes <- grch37$ensgene[match(gene_symbols, grch37$symbol)]

mapped_mask <- ensgenes %in% grch38$ensgene
gene_symbols[mapped_mask] <- grch38$symbol[match(ensgenes[mapped_mask], grch38$ensgene)]

expr_dat$symbol[grch37_mask] <- gene_symbols

# table(expr_dat$symbol %in% grch38$symbol)
#
# FALSE  TRUE
#  1242 34340

# limit sample metadata to first sample for each patient
sample_metadata <- sample_metadata %>%
  group_by(patient_id) %>%
  slice(1)

if (!all(colnames(expr_dat)[-1] == sample_metadata$patient_id)) {
  stop("Sample ID mismatch!")
}

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

